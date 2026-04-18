#!/usr/bin/env python3
"""
Scrape predictions from coordinate.rdmetallics.net for all BiometalDB complexes.
Stores results as scoring_version='rdmetallics.net' in dmpnn_summary + dmpnn_results.
Handles multi-ligand complexes sequentially (single ligand per API call).
Uses async HTTP (aiohttp) for concurrency instead of subprocess curl.
"""

import sqlite3
import re
import json
import os
import sys
import time
import asyncio
import aiohttp
import ssl as ssl_mod
from threading import Lock

DB_PATH = os.path.expanduser('~/.hermes-agent2/biometaldb/data/biometaldb.sqlite')
API_URL = 'https://coordinate.rdmetallics.net/coordinate'
SCORING_VERSION = 'rdmetallics.net'

# SSL context that skips verification (expired cert on rdmetallics.net)
SSL_CTX = ssl_mod.create_default_context()
SSL_CTX.check_hostname = False
SSL_CTX.verify_mode = ssl_mod.CERT_NONE

# Thread-safe DB writes (sqlite3 connection is shared across threads via asyncio)
db_lock = Lock()


def parse_ligands(smiles_ligands, metal_symbol):
    """Split multi-ligand SMILES into individual ligands, skip counterions/solvents."""
    if not smiles_ligands:
        return []
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles_ligands, sanitize=False)
    if mol is None:
        return []
    frags = Chem.GetMolFrags(mol, asMols=True)
    ligands = []
    halides = {'Cl', 'Br', 'I', 'F'}
    counterions = {'Na', 'K', 'Cs', 'Li', 'Mg', 'Ca', 'Ba'}
    skip = halides | counterions
    for frag in frags:
        atoms = [a.GetSymbol() for a in frag.GetAtoms()]
        # Skip single-atom counterions/halides
        if len(atoms) == 1 and atoms[0] in skip:
            continue
        # Skip single atoms that aren't the metal
        if len(atoms) <= 1:
            continue
        # Skip common solvents (acetonitrile, water, methanol, DMSO, etc.)
        smi = Chem.MolToSmiles(frag)
        if smi in ('CC#N', 'O', 'CO', 'CS(C)=O', 'CC(=O)C', 'CCO', 'CCOC',
                    'C1CCOC1', 'OCC', 'CC=O', 'N', 'CC(=O)O'):
            continue
        ligands.append(smi)
    return ligands


async def call_rdmetallics_api(session, metal_smi, ligand_smi, timeout=120, max_retries=3):
    """Call coordinate.rdmetallics.net via aiohttp and parse scores from HTML."""
    data = {
        'metal_smi': metal_smi,
        'ligand_smi': ligand_smi,
        'complex_smi': '',
    }

    for attempt in range(max_retries):
        try:
            async with session.post(
                API_URL,
                data=data,
                timeout=aiohttp.ClientTimeout(total=timeout, connect=10),
                ssl=SSL_CTX,
            ) as resp:
                if resp.status != 200:
                    if attempt < max_retries - 1:
                        await asyncio.sleep(2 ** attempt)
                        continue
                    return None
                html = await resp.text()

            # Parse candidates from caption divs: "0.84 ± 0.08<br/>SMILES"
            pattern = r'class="caption">\s*([\d.]+)\s*±\s*([\d.]+)<br/>(.*?)</div>'
            matches = re.findall(pattern, html, re.DOTALL)

            candidates = []
            for mean_str, std_str, smi in matches:
                smi = smi.strip()
                if not smi:
                    continue
                candidates.append({
                    'candidate_smi': smi,
                    'score': float(mean_str),
                    'score_std': float(std_str),
                })

            # Sort by score descending
            candidates.sort(key=lambda x: -x['score'])
            return candidates

        except (aiohttp.ClientError, asyncio.TimeoutError, OSError) as e:
            err_msg = str(e) or type(e).__name__
            if attempt < max_retries - 1:
                wait = 2 ** attempt
                print(f"  API error (attempt {attempt+1}/{max_retries}): {err_msg}, retrying in {wait}s",
                      file=sys.stderr)
                await asyncio.sleep(wait)
            else:
                print(f"  API error after {max_retries} attempts: {err_msg}", file=sys.stderr)
                return None

    return None


async def score_complex_sequential(session, complex_id, metal_symbol, ox_state, smiles_ligands):
    """
    Score a complex through the author API.
    Returns list of (step, ligand_smi, metal_smi, candidates) tuples.
    """
    ligands = parse_ligands(smiles_ligands, metal_symbol)
    if not ligands:
        return None

    metal_smi = f'[{metal_symbol}+{ox_state}]' if ox_state else f'[{metal_symbol}+2]'
    all_results = []

    for step, lig_smi in enumerate(ligands):
        candidates = await call_rdmetallics_api(session, metal_smi, lig_smi)
        if candidates is None:
            # API failed, skip this ligand entirely (don't save empty result)
            continue

        all_results.append((step, lig_smi, metal_smi, candidates))

        # Feed top-1 as new metal for next ligand
        if candidates:
            top_smi = candidates[0]['candidate_smi']
            metal_smi = top_smi

    # Return None if no API call succeeded at all
    return all_results if all_results else None


def save_to_db(conn, complex_id, metal, ox_state, donor_atoms_gt, results):
    """Save author API results to dmpnn_summary and dmpnn_results."""
    if not results:
        return

    # Find overall top-1 (from last step's top candidate)
    overall_top1 = None
    overall_score = 0.0
    total_candidates = 0
    top1_donor_pred = None

    for step, lig_smi, metal_smi, candidates in results:
        total_candidates += len(candidates)
        if candidates and candidates[0]['score'] >= overall_score:
            overall_top1 = candidates[0]['candidate_smi']
            overall_score = candidates[0]['score']
            # Extract donor atoms from top-1
            from rdkit import Chem
            top_mol = Chem.MolFromSmiles(overall_top1, sanitize=False)
            if top_mol:
                # Find metal atom
                for atom in top_mol.GetAtoms():
                    if atom.GetSymbol() == metal:
                        donors = []
                        for neighbor in atom.GetNeighbors():
                            sym = neighbor.GetSymbol()
                            if sym != 'C':
                                donors.append(sym)
                        donor_counts = {}
                        for d in donors:
                            donor_counts[d] = donor_counts.get(d, 0) + 1
                        top1_donor_pred = json.dumps(donor_counts)
                        break

    # Insert summary
    with db_lock:
        conn.execute("""
            INSERT OR REPLACE INTO dmpnn_summary
            (complex_id, metal, oxidation_state, donor_atoms_gt, top1_donor_pred, match_score, n_candidates, scoring_version)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (complex_id, metal, ox_state, donor_atoms_gt, top1_donor_pred, overall_score, total_candidates, SCORING_VERSION))

        # Insert per-candidate results
        for step, lig_smi, metal_smi, candidates in results:
            for rank, cand in enumerate(candidates):
                conn.execute("""
                    INSERT OR REPLACE INTO dmpnn_results
                    (complex_id, idx, ligand_smi, metal_smi, candidate_smi, metal_ox, score, score_std, rank, is_top1, scoring_version)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """, (complex_id, rank, lig_smi, metal_smi, cand['candidate_smi'],
                      ox_state, cand['score'], cand.get('score_std', 0.0),
                      rank + 1, 1 if rank == 0 else 0, SCORING_VERSION))

        conn.commit()


def get_complexes(conn, limit=0, offset=0, skip_done=True):
    """Get complexes from DB, optionally skipping already-processed ones."""
    if skip_done:
        query = """
            SELECT c.id, c.metal, c.oxidation_state, c.smiles_ligands, c.donor_atoms
            FROM complexes c
            WHERE c.smiles_ligands IS NOT NULL AND c.smiles_ligands != ''
              AND c.id NOT IN (
                  SELECT complex_id FROM dmpnn_summary WHERE scoring_version = ?
              )
            ORDER BY c.id
        """
        params = [SCORING_VERSION]
    else:
        query = """
            SELECT c.id, c.metal, c.oxidation_state, c.smiles_ligands, c.donor_atoms
            FROM complexes c
            WHERE c.smiles_ligands IS NOT NULL AND c.smiles_ligands != ''
            ORDER BY c.id
        """
        params = []
    if limit > 0:
        query += f" LIMIT {limit} OFFSET {offset}"
    return conn.execute(query, params).fetchall()


async def process_one(session, semaphore, complex_id, metal, ox_state, smiles_ligands, donor_atoms):
    """Process a single complex with semaphore-controlled concurrency."""
    async with semaphore:
        results = await score_complex_sequential(session, complex_id, metal, ox_state, smiles_ligands)
        return complex_id, metal, ox_state, donor_atoms, results


async def process_complexes_async(conn, complexes, max_workers, log_path):
    """Process all complexes asynchronously with bounded concurrency."""
    semaphore = asyncio.Semaphore(max_workers)
    total = len(complexes)

    success = 0
    failed = 0
    processed = 0
    start_time = time.time()

    connector = aiohttp.TCPConnector(limit=max_workers, limit_per_host=max_workers)
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = []
        for c in complexes:
            cid, metal, ox, ligands, donor_atoms = c
            task = asyncio.create_task(
                process_one(session, semaphore, cid, metal, ox or 2, ligands, donor_atoms)
            )
            tasks.append(task)

        for coro in asyncio.as_completed(tasks):
            try:
                cid, metal, ox, donor_atoms, results = await coro
                processed += 1
                if results:
                    save_to_db(conn, cid, metal, ox, donor_atoms, results)
                    success += 1
                else:
                    failed += 1
                    print(f"  [{processed}/{total}] ID={cid}: no ligands parsed")
            except Exception as e:
                processed += 1
                failed += 1
                print(f"  [{processed}/{total}] Error: {e}")

            if processed % 10 == 0:
                elapsed = time.time() - start_time
                rate = processed / elapsed if elapsed > 0 else 0
                eta = (total - processed) / rate if rate > 0 else 0
                print(f"  [{processed}/{total}] ok={success} fail={failed} "
                      f"rate={rate:.1f}/s ETA={eta/60:.0f}min")

    return success, failed


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--limit', type=int, default=0, help='Max complexes (0=all)')
    parser.add_argument('--offset', type=int, default=0)
    parser.add_argument('--workers', type=int, default=5, help='Max concurrent requests')
    parser.add_argument('--test', action='store_true', help='Test mode (5 random complexes)')
    parser.add_argument('--no-resume', action='store_true', help='Reprocess all (not just skipped)')
    parser.add_argument('--log', type=str, default='/root/.hermes-agent2/biometaldb/rdmetallics_scraper.log', help='Log file path')
    args = parser.parse_args()

    conn = sqlite3.connect(DB_PATH)

    # Create dmpnn_results if not exists with proper schema
    conn.execute("""
        CREATE TABLE IF NOT EXISTS dmpnn_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            complex_id INTEGER,
            idx INTEGER,
            ligand_smi TEXT,
            metal_smi TEXT,
            candidate_smi TEXT,
            metal_ox INTEGER,
            score REAL,
            rank INTEGER,
            is_top1 INTEGER,
            score_std REAL DEFAULT 0.0,
            FOREIGN KEY (complex_id) REFERENCES complexes(id)
        )
    """)

    if args.test:
        complexes = conn.execute("""
            SELECT c.id, c.metal, c.oxidation_state, c.smiles_ligands, c.donor_atoms
            FROM complexes c
            WHERE c.smiles_ligands IS NOT NULL AND c.smiles_ligands != ''
            ORDER BY RANDOM() LIMIT 5
        """).fetchall()
    else:
        complexes = get_complexes(conn, args.limit, args.offset, skip_done=not args.no_resume)

    total = len(complexes)
    print(f"Processing {total} complexes with max {args.workers} concurrent requests...")

    success, failed = asyncio.run(process_complexes_async(conn, complexes, args.workers, args.log))

    print(f"\nDone: {success} success, {failed} failed")
    conn.close()


if __name__ == '__main__':
    main()
