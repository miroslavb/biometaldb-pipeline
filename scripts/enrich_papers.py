#!/usr/bin/env python3
"""
Enrich papers table with title and authors from CrossRef API.
Processes all DOIs in the measurements table.
"""
import sqlite3
import json
import time
import sys
import os
from urllib.request import Request, urlopen
from urllib.error import URLError, HTTPError
from urllib.parse import quote

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite")

def fetch_crossref(doi, retries=2):
    """Fetch paper metadata from CrossRef API."""
    url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
    headers = {
        "User-Agent": "MetalCytoToxDB/1.0 (mailto:serebryanskaya@example.com)",
        "Accept": "application/json"
    }
    
    for attempt in range(retries + 1):
        try:
            req = Request(url, headers=headers)
            with urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode())
                msg = data.get("message", {})
                
                # Extract title
                title_list = msg.get("title", [])
                title = title_list[0] if title_list else None
                
                # Extract authors
                authors_list = msg.get("author", [])
                if authors_list:
                    authors = "; ".join(
                        f"{a.get('family', '')}, {a.get('given', '')}" 
                        for a in authors_list
                        if a.get('family')
                    )
                else:
                    authors = None
                
                # Extract journal
                container = msg.get("container-title", [])
                journal = container[0] if container else None
                
                # Extract year
                year = None
                for date_field in ("published-print", "published-online", "created"):
                    date_parts = msg.get(date_field, {}).get("date-parts", [[]])
                    if date_parts and date_parts[0]:
                        year = date_parts[0][0]
                        break
                
                return {
                    "title": title,
                    "authors": authors,
                    "journal": journal,
                    "year": year
                }
        except HTTPError as e:
            if e.code == 404:
                return None  # DOI not found
            if attempt < retries:
                time.sleep(1)
                continue
            return None
        except (URLError, TimeoutError):
            if attempt < retries:
                time.sleep(2)
                continue
            return None
        except Exception:
            if attempt < retries:
                time.sleep(1)
                continue
            return None
    return None


def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Get all DOIs that need enrichment
    cur.execute("""
        SELECT DISTINCT doi FROM measurements 
        WHERE doi IS NOT NULL AND doi != ''
        ORDER BY doi
    """)
    all_dois = [row[0] for row in cur.fetchall()]
    print(f"Total unique DOIs: {len(all_dois)}")
    
    # Check which already have titles
    cur.execute("SELECT doi FROM papers WHERE title IS NOT NULL AND title != ''")
    already_done = set(row[0] for row in cur.fetchall())
    print(f"Already enriched: {len(already_done)}")
    
    todo = [d for d in all_dois if d not in already_done]
    print(f"To process: {len(todo)}")
    
    if not todo:
        print("All DOIs already enriched!")
        conn.close()
        return
    
    # Process in batches
    success = 0
    failed = 0
    not_found = 0
    batch_size = 50  # CrossRef polite rate: ~50 req/s
    
    for i, doi in enumerate(todo):
        if i % 10 == 0:
            print(f"\n[{i+1}/{len(todo)}] success={success} failed={failed} not_found={not_found}")
            sys.stdout.flush()
        
        result = fetch_crossref(doi)
        
        if result and (result["title"] or result["authors"]):
            # Update papers table
            cur.execute("""
                INSERT INTO papers (doi, year, title, authors)
                VALUES (?, ?, ?, ?)
                ON CONFLICT(doi) DO UPDATE SET
                    title = COALESCE(excluded.title, title),
                    authors = COALESCE(excluded.authors, authors),
                    year = COALESCE(excluded.year, year)
            """, (doi, result.get("year"), result.get("title"), result.get("authors")))
            
            # Also update year in measurements if we got a better one
            if result.get("year"):
                cur.execute("""
                    UPDATE measurements SET year = ? 
                    WHERE doi = ? AND (year IS NULL OR year = 0)
                """, (result["year"], doi))
            
            success += 1
        elif result is None:
            not_found += 1
        else:
            failed += 1
        
        # Commit every 50
        if (i + 1) % 50 == 0:
            conn.commit()
        
        # Rate limiting: ~10 req/s to be polite
        time.sleep(0.1)
    
    conn.commit()
    
    # Final stats
    cur.execute("SELECT COUNT(*) FROM papers WHERE title IS NOT NULL AND title != ''")
    total_with_title = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM papers WHERE authors IS NOT NULL AND authors != ''")
    total_with_authors = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM papers")
    total_papers = cur.fetchone()[0]
    
    print(f"\n=== DONE ===")
    print(f"Success: {success}")
    print(f"Not found: {not_found}")
    print(f"Failed: {failed}")
    print(f"Total papers with title: {total_with_title}/{total_papers}")
    print(f"Total papers with authors: {total_with_authors}/{total_papers}")
    
    conn.close()


if __name__ == "__main__":
    main()
