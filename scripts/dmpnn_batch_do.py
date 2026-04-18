#!/usr/bin/env python3
"""
Parallel D-MPNN batch scorer for BiometalDB.
Uses multiprocessing with per-worker temp dirs to avoid file conflicts.
Designed for DO droplets (4+ vCPUs).

Usage:
    python dmpnn_batch_do.py [--db PATH] [--workers N] [--limit N] [--model-dir DIR]
"""
import sys, os, json, sqlite3, subprocess, time, argparse, tempfile, shutil, logging
from multiprocessing import Pool, cpu_count
from collections import Counter
import pandas as pd
import numpy as np
from rdkit import Chem

sys.path.insert(0, '/root/RDMetallics_coordinate')
import rdmetallics
from rdmetallics.functions_wrapped import generate_candidate_answers, generate_model_input

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(message)s',
    datefmt='%H:%M:%S',
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger('dmpnn')

COMPLEX_TIMEOUT = 120
FOLD_TIMEOUT = 90


def extract_donors_haptic(candidate_smi):
    mol = Chem.MolFromSmiles(candidate_smi)
    if mol is None:
        return Counter(), False
    metal_idx = None
    for a in mol.GetAtoms():
        if a.GetAtomicNum() >= 21:
            metal_idx = a.GetIdx()
            break
    if metal_idx is None:
        return Counter(), False
    donors = Counter()
    has_haptic = False
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if a1.GetIdx() == metal_idx:
            other = a2
        elif a2.GetIdx() == metal_idx:
            other = a1
        else:
            continue
        sym = other.GetSymbol()
        if sym == 'C':
            has_haptic = True
        else:
            donors[sym] += 1
    return donors, has_haptic


def match_score(gt, pred, pred_haptic=False):
    if not gt and not pred and not pred_haptic:
        return 1.0
    if pred_haptic:
        return 0.0
    if not gt or not pred:
        return 0.0
    elems = set(gt.keys()) | set(pred.keys())
    inter = sum(min(gt.get(e, 0), pred.get(e, 0)) for e in elems)
    union = sum(max(gt.get(e, 0), pred.get(e, 0)) for e in elems)
    return inter / union if union else 0.0


def parse_donor_atoms(donor_json):
    if not donor_json or donor_json in ('None', 'null', ''):
        return Counter()
    try:
        return Counter(json.loads(donor_json))
    except:
        return Counter()


def run_chemprop_folds(input_csv, desc_npz, model_dir, tmpdir, fold_timeout=90):
    """Run chemprop predict for both folds in parallel, using tmpdir for outputs."""
    procs = []
    for fold in (0, 1):
        preds_path = os.path.join(tmpdir, f'prediction_fold_{fold}.csv')
        model_path = os.path.join(model_dir, f'best_{fold}.pt')
        cmd = [
            'chemprop', 'predict',
            '--test-path', input_csv,
            '--descriptors-path', desc_npz,
            '--model-path', model_path,
            '--smiles-columns', 'ligands_smi', 'metals_smi', 'candidates_smi',
            '--preds-path', preds_path,
            '--accelerator', 'gpu', '--devices', '1'
        ]
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        procs.append((fold, p, preds_path))

    results = {}
    for fold, p, preds_path in procs:
        try:
            p.wait(timeout=fold_timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            p.wait()
            return None
        if p.returncode != 0:
            return None
        try:
            df = pd.read_csv(preds_path)
            results[fold] = df.iloc[:, -1].values
        except:
            return None

    return np.column_stack([results[0], results[1]]) if len(results) == 2 else None


def score_one(args):
    """Score a single complex. args = (cid, donor_json, ligand_smi, metal_smi_db, ox, metal_sym, model_dir)"""
    cid, donor_json, ligand_smi, metal_smi_db, ox, metal_sym, model_dir = args
    metal_ion_smi = f'[{metal_sym}+{ox}]' if ox else f'[{metal_sym}+2]'
    metal = Chem.MolFromSmiles(metal_ion_smi)
    ligand = Chem.MolFromSmiles(ligand_smi)
    if metal is None or ligand is None:
        return (cid, metal_sym, ox or 0, donor_json, '', 0.0, 0, 'parse')

    try:
        candidates = generate_candidate_answers(metal, ligand, ligand_smi)
    except Exception as e:
        return (cid, metal_sym, ox or 0, donor_json, '', 0.0, 0, f'cand:{str(e)[:40]}')

    if not candidates:
        return (cid, metal_sym, ox or 0, donor_json, '', 0.0, 0, 'no_cand')

    candidates = candidates[:30]

    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = os.path.join(tmpdir, 'input.csv')
        desc_npz = os.path.join(tmpdir, 'desc.npz')
        generate_model_input(candidates, path_input=input_csv, path_descriptor=desc_npz)

        all_scores = run_chemprop_folds(input_csv, desc_npz, model_dir, tmpdir, FOLD_TIMEOUT)
        if all_scores is None:
            return (cid, metal_sym, ox or 0, donor_json, '', 0.0, len(candidates), 'chemprop')

        test_df = pd.read_csv(input_csv)
        mean_scores = all_scores.mean(axis=1)
        sorted_indices = np.argsort(-mean_scores)

        top_smi = test_df.iloc[sorted_indices[0]]['candidates_smi']
        pd_h, ph = extract_donors_haptic(top_smi)
        pred_str = "C:6(eta)" if (ph and not pd_h) else json.dumps(dict(pd_h))
        gt = parse_donor_atoms(donor_json)
        ms = match_score(gt, pd_h, ph)

        return (cid, metal_sym, ox or 0, donor_json, pred_str, ms, len(candidates), '')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db', default=os.path.expanduser('~/.hermes-agent2/biometaldb/data/biometaldb.sqlite'))
    parser.add_argument('--workers', type=int, default=min(cpu_count(), 4))
    parser.add_argument('--offset', type=int, default=0, help='Skip first N unscored complexes')
    parser.add_argument('--limit', type=int, default=0, help='Max complexes (0=all)')
    parser.add_argument('--model-dir', default='/root/RDMetallics_coordinate/rdmetallics/model_random_split')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    log.info(f"DB: {args.db}")
    log.info(f"Workers: {args.workers}")
    log.info(f"Model dir: {args.model_dir}")

    conn = sqlite3.connect(args.db)
    cur = conn.cursor()

    cur.execute("""
        SELECT c.id, c.donor_atoms, c.smiles_ligands, c.metal_smiles, c.oxidation_state, c.metal
        FROM complexes c
        WHERE c.has_mol3=1
        AND c.id NOT IN (SELECT complex_id FROM dmpnn_summary)
        ORDER BY c.id
    """)
    all_rows = cur.fetchall()
    if args.offset > 0:
        all_rows = all_rows[args.offset:]
    if args.limit > 0:
        rows = all_rows[:args.limit]
    else:
        rows = all_rows
    conn.close()

    log.info(f"Unscored complexes: {len(rows)}")
    if args.dry_run:
        return

    work_items = [(cid, da, smi, m_smi, ox, m, args.model_dir) for cid, da, smi, m_smi, ox, m in rows]

    t0 = time.time()
    ok = 0
    err = 0
    batch_conn = sqlite3.connect(args.db)

    with Pool(args.workers) as pool:
        for i, result in enumerate(pool.imap_unordered(score_one, work_items, chunksize=5)):
            cid, metal, ox, donor_json, pred, match, n_cand, err_str = result

            if err_str:
                err += 1
            else:
                ok += 1

            batch_conn.execute("""INSERT OR REPLACE INTO dmpnn_summary
                (complex_id, metal, oxidation_state, donor_atoms_gt, top1_donor_pred, match_score, n_candidates)
                VALUES (?, ?, ?, ?, ?, ?, ?)""",
                (cid, metal, ox, donor_json, pred, match, n_cand))
            batch_conn.commit()

            if (i + 1) % 50 == 0 or (i + 1) == len(work_items):
                elapsed = time.time() - t0
                rate = (i + 1) / elapsed
                eta = (len(work_items) - i - 1) / rate if rate > 0 else 0
                log.info(f"[{i+1}/{len(work_items)}] ok={ok} err={err} | {rate:.1f}/s | ETA {eta/3600:.1f}h")

    batch_conn.close()
    elapsed = time.time() - t0
    log.info(f"DONE: {ok} ok, {err} errors, {elapsed/3600:.1f}h total")


if __name__ == '__main__':
    main()
