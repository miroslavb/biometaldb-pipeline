#!/usr/bin/env python3
"""
First-pass enrichment: sample 100 complexes and attempt to extract
coordination metadata (donor atoms, coord number, geometry) from papers.

Methods:
1. CrossRef API — full metadata (abstract, references, journal)
2. Semantic Scholar API — structured paper data + citation context
3. Text mining from abstracts — regex/heuristic extraction of coord info
"""
import sqlite3
import json
import re
import time
import sys
import os
from urllib.request import Request, urlopen
from urllib.parse import quote
from urllib.error import HTTPError

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "biometaldb.sqlite")
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")

# ==================== CrossRef ====================

def fetch_crossref_full(doi):
    """Fetch full metadata from CrossRef including abstract."""
    url = f"https://api.crossref.org/works/{quote(doi, safe='')}"
    headers = {
        "User-Agent": "MetalCytoToxDB/1.0 (research; mailto:sereb@example.com)",
        "Accept": "application/json"
    }
    try:
        req = Request(url, headers=headers)
        with urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
            msg = data.get("message", {})
            
            title_list = msg.get("title", [])
            title = title_list[0] if title_list else None
            
            abstract = msg.get("abstract", None)
            # Clean JATS XML from abstract
            if abstract:
                abstract = re.sub(r'<[^>]+>', '', abstract)
                abstract = re.sub(r'\s+', ' ', abstract).strip()
            
            container = msg.get("container-title", [])
            journal = container[0] if container else None
            
            authors = []
            for a in msg.get("author", []):
                name = f"{a.get('family', '')}, {a.get('given', '')}"
                if name.strip(", "):
                    authors.append(name.strip(", "))
            
            year = None
            for df in ("published-print", "published-online", "created"):
                dp = msg.get(df, {}).get("date-parts", [[]])
                if dp and dp[0]:
                    year = dp[0][0]
                    break
            
            # Subject/keywords
            subjects = msg.get("subject", [])
            keywords = msg.get("keyword", []) if isinstance(msg.get("keyword"), list) else []
            
            return {
                "source": "crossref",
                "doi": doi,
                "title": title,
                "authors": "; ".join(authors) if authors else None,
                "journal": journal,
                "year": year,
                "abstract": abstract,
                "subjects": subjects,
                "keywords": keywords if isinstance(keywords, list) else []
            }
    except HTTPError as e:
        if e.code == 404:
            return {"source": "crossref", "doi": doi, "error": "not_found"}
        return {"source": "crossref", "doi": doi, "error": f"http_{e.code}"}
    except Exception as e:
        return {"source": "crossref", "doi": doi, "error": str(e)[:100]}


# ==================== Semantic Scholar ====================

def fetch_semantic_scholar(doi):
    """Fetch paper data from Semantic Scholar (free API, no key needed)."""
    url = f"https://api.semanticscholar.org/graph/v1/paper/DOI:{quote(doi, safe='')}?fields=title,abstract,authors,year,citationCount,referenceCount,tldr,fieldsOfStudy"
    headers = {"Accept": "application/json"}
    try:
        req = Request(url, headers=headers)
        with urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode())
            
            tldr = data.get("tldr", {})
            tldr_text = tldr.get("text") if tldr else None
            
            authors = []
            for a in data.get("authors", []):
                if a.get("name"):
                    authors.append(a["name"])
            
            fields = data.get("fieldsOfStudy", []) or []
            
            return {
                "source": "semantic_scholar",
                "doi": doi,
                "title": data.get("title"),
                "authors": "; ".join(authors) if authors else None,
                "year": data.get("year"),
                "abstract": data.get("abstract"),
                "tldr": tldr_text,
                "citation_count": data.get("citationCount"),
                "reference_count": data.get("referenceCount"),
                "fields_of_study": fields
            }
    except HTTPError as e:
        if e.code == 404:
            return {"source": "semantic_scholar", "doi": doi, "error": "not_found"}
        return {"source": "semantic_scholar", "doi": doi, "error": f"http_{e.code}"}
    except Exception as e:
        return {"source": "semantic_scholar", "doi": doi, "error": str(e)[:100]}


# ==================== Text Mining ====================

def extract_coord_info(text, metal, smiles):
    """Extract coordination chemistry info from abstract/text using heuristics."""
    if not text:
        return {}
    
    text_lower = text.lower()
    result = {}
    
    # Coordination number patterns
    cn_patterns = [
        r'(\d+)[\s-]*coordinat(?:ion|ed|e)',
        r'coordinat(?:ion|ed|e)[\s-]*(?:number|sphere)?[\s-]*(?:of[\s-]*)?(\d+)',
        r'(?:hexa|penta|tetra|tri|di|mono|octa)(?:dentate|coordinate)',
        r'(\d+)[\s-]*(?:membered|coordinate)',
    ]
    coord_words = {
        "hexadentate": 6, "pentadentate": 5, "tetradentate": 4,
        "tridentate": 3, "bidentate": 2, "monodentate": 1,
        "octahedral": 6, "tetrahedral": 4, "square planar": 4,
        "trigonal bipyramidal": 5, "square pyramidal": 5,
    }
    
    for word, cn in coord_words.items():
        if word in text_lower:
            result["coordination_number"] = cn
            result["coordination_geometry"] = word
            break
    
    for pat in cn_patterns:
        m = re.search(pat, text_lower)
        if m:
            try:
                cn = int(m.group(1))
                if 1 <= cn <= 12:
                    result["coordination_number"] = cn
            except (ValueError, IndexError):
                pass
    
    # Donor atom patterns
    donor_patterns = [
        r'(?:N,O|O,N)[\s-]*(?:donor|chelat|bind|coordinat)',
        r'(?:N,N|O,O|S,N|P,N)[\s-]*(?:donor|chelat|bind|coordinat)',
        r'(?:through|via|by)[\s-]*(?:the[\s-]*)?(N|O|S|P|Cl)[\s-]*atoms?',
        r'(?:bound|coordinat\w*)[\s-]*(?:through|via|to)[\s-]*(?:the[\s-]*)?(N|O|S|P|Cl)',
        r'(?:deprotonated|anionic)[\s-]*(N|O|S)',
    ]
    
    donors_found = set()
    for pat in donor_patterns:
        matches = re.findall(pat, text_lower)
        for m in matches:
            if isinstance(m, str) and len(m) <= 2:
                elem = m.upper()
                if elem in ("N", "O", "S", "P", "CL", "BR", "I"):
                    donors_found.add(elem if elem != "CL" else "Cl")
    
    if "n,o" in text_lower or "o,n" in text_lower:
        donors_found.update(["N", "O"])
    if "n,n" in text_lower:
        donors_found.add("N")
    if "o,o" in text_lower:
        donors_found.add("O")
    if "s,n" in text_lower or "n,s" in text_lower:
        donors_found.update(["S", "N"])
    
    if donors_found:
        result["donor_elements"] = sorted(donors_found)
    
    # Ligand type patterns
    ligand_keywords = []
    ligand_types = {
        "cyclometalated": "cyclometalated",
        "cyclometallic": "cyclometalated",
        "arene": "arene",
        "p-cymene": "p-cymene",
        "cyclopentadienyl": "Cp",
        "bipyridine": "bpy",
        "phenanthroline": "phen",
        "phenylpyridine": "ppy",
        "acetylacetonate": "acac",
        "picolinate": "pic",
        "tetrazole": "tetrazole",
        "tetrazolato": "tetrazolate",
        "imidazole": "imidazole",
        "pyridine": "pyridine",
        "thiolate": "thiolate",
        "carbene": "NHC",
        "terpyridine": "terpy",
    }
    
    for keyword, lig_type in ligand_types.items():
        if keyword in text_lower:
            ligand_keywords.append(lig_type)
    
    if ligand_keywords:
        result["ligand_types"] = ligand_keywords
    
    # Geometry patterns
    geom_patterns = [
        r'(octahedral|tetrahedral|square[\s-]*planar|trigonal[\s-]*bipyramidal|square[\s-]*pyramidal|linear|trigonal[\s-]*prismatic)',
    ]
    for pat in geom_patterns:
        m = re.search(pat, text_lower)
        if m:
            result["coordination_geometry"] = m.group(1).replace("-", " ")
            break
    
    return result


# ==================== Main ====================

def main():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    
    # Sample 100 complexes: mix of metals, prefer those with donor_atoms already
    cur.execute("""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, 
               c.charge_complex, c.donor_atoms, m.doi
        FROM complexes c
        JOIN measurements m ON m.complex_id = c.id
        WHERE m.doi IS NOT NULL AND m.doi != ''
        GROUP BY c.id
        ORDER BY RANDOM()
        LIMIT 100
    """)
    sample = cur.fetchall()
    print(f"Sample: {len(sample)} complexes")
    
    # Collect unique DOIs
    dois = list(set(r[6] for r in sample if r[6]))
    print(f"Unique DOIs to fetch: {len(dois)}")
    
    results = []
    
    for i, (cid, smiles, metal, ox_state, charge, donor_atoms, doi) in enumerate(sample):
        print(f"\n[{i+1}/100] Complex #{cid}: {metal}({ox_state}) DOI: {doi[:50]}...")
        sys.stdout.flush()
        
        # Method 1: CrossRef
        cr = fetch_crossref_full(doi)
        abstract = cr.get("abstract", "")
        
        # Method 2: Semantic Scholar
        ss = fetch_semantic_scholar(doi)
        tldr = ss.get("tldr", "") or ""
        
        # Method 3: Text mining
        combined_text = f"{abstract} {tldr}"
        coord_info = extract_coord_info(combined_text, metal, smiles)
        
        result = {
            "complex_id": cid,
            "metal": metal,
            "oxidation_state": ox_state,
            "smiles": smiles[:100],
            "existing_donor_atoms": donor_atoms,
            "doi": doi,
            "crossref": {
                "title": cr.get("title", "")[:100] if cr.get("title") else None,
                "has_abstract": bool(abstract),
                "journal": cr.get("journal"),
                "year": cr.get("year"),
            },
            "semantic_scholar": {
                "has_tldr": bool(tldr),
                "citation_count": ss.get("citation_count"),
                "fields": ss.get("fields_of_study", []),
            },
            "extracted": coord_info,
        }
        
        results.append(result)
        
        # Rate limiting
        time.sleep(0.15)  # ~6-7 req/s total (2 APIs per complex)
    
    conn.close()
    
    # Save results
    output_path = os.path.join(OUTPUT_DIR, "enrichment_sample_100.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2)
    print(f"\n\nResults saved to {output_path}")
    
    # Summary statistics
    print("\n" + "=" * 60)
    print("  ENRICHMENT SAMPLE STATISTICS (100 complexes)")
    print("=" * 60)
    
    n_abstracts = sum(1 for r in results if r["crossref"]["has_abstract"])
    n_tldr = sum(1 for r in results if r["semantic_scholar"]["has_tldr"])
    n_coord_num = sum(1 for r in results if r["extracted"].get("coordination_number"))
    n_donor_extracted = sum(1 for r in results if r["extracted"].get("donor_elements"))
    n_geom = sum(1 for r in results if r["extracted"].get("coordination_geometry"))
    n_ligand_types = sum(1 for r in results if r["extracted"].get("ligand_types"))
    n_crossref_ok = sum(1 for r in results if not r["crossref"].get("error"))
    n_ss_ok = sum(1 for r in results if not r["semantic_scholar"].get("error"))
    
    print(f"\nCrossRef success: {n_crossref_ok}/100")
    print(f"Semantic Scholar success: {n_ss_ok}/100")
    print(f"Abstracts available: {n_abstracts}/100")
    print(f"TLDR available: {n_tldr}/100")
    print(f"Coordination number extracted: {n_coord_num}/100")
    print(f"Donor elements extracted: {n_donor_extracted}/100")
    print(f"Geometry extracted: {n_geom}/100")
    print(f"Ligand types extracted: {n_ligand_types}/100")
    
    # Show some examples
    print("\n=== Examples with extracted coord info ===")
    for r in results:
        ext = r["extracted"]
        if ext.get("coordination_number") or ext.get("donor_elements"):
            print(f"\n  #{r['complex_id']} {r['metal']}({r['oxidation_state']})")
            print(f"  DOI: {r['doi']}")
            if ext.get("coordination_number"):
                print(f"  Coord number: {ext['coordination_number']}")
            if ext.get("coordination_geometry"):
                print(f"  Geometry: {ext['coordination_geometry']}")
            if ext.get("donor_elements"):
                print(f"  Donors: {ext['donor_elements']}")
            if ext.get("ligand_types"):
                print(f"  Ligands: {ext['ligand_types']}")


if __name__ == "__main__":
    main()
