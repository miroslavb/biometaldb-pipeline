#!/usr/bin/env python3
"""Flask server for BiometalDB: MOL download + 2D structure viewer.
Runs on port 8502 alongside Datasette (8501).
"""
import os
import base64
from flask import Flask, abort, send_file, Response, request, render_template_string
import sqlite3
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MOL_DIR = os.path.join(BASE_DIR, "data", "mol3")
DB_PATH = os.path.join(BASE_DIR, "data", "biometaldb.sqlite")
VIEWER_DIR = os.path.join(BASE_DIR, "viewer")

app = Flask(__name__)

# Register enhanced D-MPNN routes
from dmpnn_routes import register_dmpnn_routes
register_dmpnn_routes(app)

# Register enhanced complexes routes (pagination, filters, 3D structures)
from complexes_routes import register_complexes_routes
register_complexes_routes(app)


# Complexes routes moved to complexes_routes.py (registered above)


@app.route("/mol3/<int:complex_id>")
def download_mol(complex_id):
    """Download MOL V3000 file by complex ID."""
    path = os.path.join(MOL_DIR, f"complex_{complex_id}.mol")
    if not os.path.exists(path):
        abort(404)
    return send_file(path, mimetype="chemical/x-mdl-molfile",
                     as_attachment=True,
                     download_name=f"complex_{complex_id}.mol")


@app.route("/mol3/<int:complex_id>/view")
def view_mol(complex_id):
    """View MOL V3000 file inline in browser."""
    path = os.path.join(MOL_DIR, f"complex_{complex_id}.mol")
    if not os.path.exists(path):
        abort(404)
    return send_file(path, mimetype="chemical/x-mdl-molfile",
                     as_attachment=False)


@app.route("/mol3/<int:complex_id>/smiles")
def complex_smiles(complex_id):
    """Return SMILES for a complex (from DB)."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("SELECT smiles_ligands, metal, oxidation_state FROM complexes WHERE id = ?", (complex_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        abort(404)
    return {
        "complex_id": complex_id,
        "smiles_ligands": row[0],
        "metal": row[1],
        "oxidation_state": row[2],
        "mol3_url": f"/mol3/{complex_id}",
    }


@app.route("/mol3/<int:complex_id>/candidate")
def view_candidate_mol(complex_id):
    """Generate MOL from top-1 D-MPNN candidate SMILES (with coordination bonds)."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    # Try sequential first, then single-ligand
    for ver in ['sequential', 'single-ligand']:
        cur.execute("""SELECT r.candidate_smi FROM dmpnn_results r
            JOIN dmpnn_summary s ON s.complex_id = r.complex_id
            WHERE r.complex_id=? AND r.is_top1=1 AND s.scoring_version=?
            LIMIT 1""", (complex_id, ver))
        row = cur.fetchone()
        if row:
            break
    conn.close()

    if not row:
        abort(404)

    smi = row[0]
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if not mol:
        abort(500)

    # Compute 2D coords first (required for bonds to be preserved)
    AllChem.Compute2DCoords(mol)
    mol_block = Chem.MolToMolBlock(mol, forceV3000=True)
    return Response(mol_block, mimetype="chemical/x-mdl-molfile")


@app.route("/api/mol3/batch", methods=["POST"])
def batch_mol_info():
    """Return MOL availability for multiple complex IDs."""
    ids = request.json.get("ids", [])
    if not ids or len(ids) > 1000:
        return {"error": "Provide up to 1000 ids"}, 400

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    placeholders = ",".join("?" * len(ids))
    cur.execute(f"""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, c.donor_atoms, c.has_mol3
        FROM complexes c WHERE c.id IN ({placeholders})
    """, ids)
    results = []
    for row in cur.fetchall():
        results.append({
            "complex_id": row[0],
            "smiles_ligands": row[1][:100] if row[1] else None,
            "metal": row[2],
            "oxidation_state": row[3],
            "donor_atoms": row[4],
            "has_mol3": bool(row[5]),
            "mol3_url": f"/mol3/{row[0]}",
        })
    conn.close()
    return {"results": results}


@app.route("/api/complexes")
def api_complexes():
    """Paginated, filterable complex list — pure SQLite, no filesystem calls."""
    offset = int(request.args.get("offset", 0))
    limit = min(int(request.args.get("limit", 50)), 200)
    metal = request.args.get("metal", "")
    has_mol = request.args.get("has_mol3", "")
    q = request.args.get("q", "").strip()

    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()

    where = []
    params = []
    if metal:
        where.append("c.metal = ?")
        params.append(metal)
    if has_mol == "1":
        where.append("c.has_mol3 = 1")
    if q:
        if q.startswith("#"):
            # Search by ID
            try:
                where.append("c.id = ?")
                params.append(int(q[1:]))
            except ValueError:
                pass
        else:
            where.append("(c.smiles_ligands LIKE ? OR c.metal LIKE ? OR CAST(c.id AS TEXT) LIKE ?)")
            params.extend([f"%{q}%", f"%{q}%", f"%{q}%"])

    where_sql = ("WHERE " + " AND ".join(where)) if where else ""

    # Count
    cur.execute(f"SELECT COUNT(*) FROM complexes c {where_sql}", params)
    total = cur.fetchone()[0]

    # Fetch page
    cur.execute(f"""
        SELECT c.id, c.smiles_ligands, c.metal, c.oxidation_state, c.donor_atoms, c.has_mol3
        FROM complexes c {where_sql}
        ORDER BY c.id
        LIMIT ? OFFSET ?
    """, params + [limit, offset])

    results = []
    for row in cur.fetchall():
        results.append({
            "complex_id": row[0],
            "smiles_ligands": (row[1][:80] + "...") if row[1] and len(row[1]) > 80 else row[1],
            "metal": row[2],
            "oxidation_state": row[3],
            "donor_atoms": row[4],
            "has_mol3": bool(row[5]),
        })

    # Get metal list for filter
    cur.execute("SELECT DISTINCT metal FROM complexes ORDER BY metal")
    metals = [r[0] for r in cur.fetchall()]

    conn.close()
    return {"results": results, "total": total, "offset": offset, "limit": limit, "metals": metals}


@app.route("/viewer/")
@app.route("/viewer/<path:filename>")
def serve_viewer(filename="index.html"):
    """Serve the 3D MOL viewer."""
    path = os.path.join(VIEWER_DIR, filename)
    if not os.path.exists(path):
        abort(404)
    return send_file(path)


# D-MPNN routes moved to dmpnn_routes.py (registered above via register_dmpnn_routes)


@app.route("/backup")
def backup_status():
    """NAS backup dashboard."""
    import subprocess
    import json as _json

    log_path = "/root/.hermes/logs/nas-backup.log"
    backup_script = "/root/.hermes/scripts/nas-backup.sh"

    # Read last log entries
    try:
        with open(log_path) as f:
            lines = f.readlines()
        log_lines = lines[-30:] if lines else ["No log yet"]
        # Find last "===" header
        last_run = ""
        for line in reversed(lines):
            if line.startswith("==="):
                last_run = line.strip()
                break
    except FileNotFoundError:
        log_lines = ["Log not found"]
        last_run = "Never"

    # Check if rsync is running
    try:
        r = subprocess.run(["pgrep", "-f", "rsync.*hermes"],
                           capture_output=True, text=True)
        running = r.returncode == 0
    except:
        running = False

    # NAS disk usage
    try:
        r = subprocess.run(
            ["sshpass", "-p", "Nhmft1378", "ssh", "-o",
             "StrictHostKeyChecking=no", "bereft@100.97.192.84",
             "du -sh /volume1/Backups/hermes/ && df -h /volume1 | tail -1"],
            capture_output=True, text=True, timeout=10)
        nas_info = r.stdout.strip() if r.returncode == 0 else "N/A"
    except:
        nas_info = "Connection timeout"

    # Backup dirs
    sources = [
        "/root/.hermes",
        "/root/.hermes-agent2",
        "/etc/nginx",
        "/etc/ssh",
        "/root/.ssh",
    ]
    dirs_html = ""
    for d in sources:
        try:
            sz = subprocess.run(["du", "-sh", d], capture_output=True,
                                text=True, timeout=5).stdout.split()[0]
        except:
            sz = "—"
        exists = os.path.isdir(d)
        color = "#16a34a" if exists else "#dc2626"
        dirs_html += f'<tr><td style="color:{color}">{"●" if exists else "○"}</td><td><code>{d}</code></td><td>{sz}</td></tr>\n'

    # Cron job status
    try:
        r = subprocess.run(
            ["bash", "-c", "crontab -l 2>/dev/null | grep -i backup || echo 'No crontab'"],
            capture_output=True, text=True)
        cron_info = r.stdout.strip()
    except:
        cron_info = "Unknown"

    status_color = "#f59e0b" if running else "#16a34a"
    status_text = "⏳ Syncing..." if running else "✅ Idle"

    log_html = "".join(
        f"<div style='color:{'#f87171' if 'FAIL' in l else '#94a3b8'}'>"
        f"{l.rstrip()}</div>" for l in log_lines)

    return f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<link rel="icon" href="/static/favicon.ico" type="image/x-icon">
<title>NAS Backup — Hermes</title>
<style>
body{{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#0f172a;color:#e2e8f0}}
.hdr{{background:linear-gradient(135deg,#1e293b,#0f172a);padding:1.5rem 2rem;border-bottom:3px solid #3b82f6}}
.hdr h1{{margin:0;font-size:1.4rem;color:#f8fafc}}
.hdr .sub{{color:#94a3b8;font-size:.85rem;margin-top:.2rem}}
.wrap{{max-width:900px;margin:1.5rem auto;padding:0 1.5rem}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:1rem;margin-bottom:1.5rem}}
.stat{{background:#1e293b;border-radius:8px;padding:1rem 1.2rem;text-align:center}}
.stat .v{{font-size:1.6rem;font-weight:700}}
.stat .l{{font-size:.8rem;color:#94a3b8;margin-top:.2rem}}
.card{{background:#1e293b;border-radius:8px;padding:1.2rem;margin-bottom:1rem}}
.card h2{{font-size:1rem;color:#60a5fa;margin:0 0 .8rem;border-bottom:1px solid #334155;padding-bottom:.5rem}}
table{{width:100%;border-collapse:collapse;font-size:.85rem}}
td{{padding:.3rem .5rem;border-bottom:1px solid #1e293b}}
.log{{background:#020617;border-radius:6px;padding:.8rem;font-family:monospace;font-size:.75rem;max-height:400px;overflow-y:auto;line-height:1.6}}
a{{color:#60a5fa;text-decoration:none}}a:hover{{text-decoration:underline}}
.btn{{display:inline-block;padding:.5rem 1rem;background:#3b82f6;color:#fff;border-radius:6px;font-size:.85rem;cursor:pointer;border:none}}
.btn:hover{{background:#2563eb}}
</style></head><body>
<div class="hdr">
<h1>💾 NAS Backup Dashboard</h1>
<div class="sub">Synology DS220j — {status_text}</div>
</div>
<div class="wrap">
<div class="grid">
<div class="stat"><div class="v" style="color:{status_color}">{status_text}</div><div class="l">Status</div></div>
<div class="stat"><div class="v">{last_run.replace("=== ","").replace(" ===","")}</div><div class="l">Last run</div></div>
<div class="stat"><div class="v">{nas_info.split(chr(10))[0] if nas_info != "N/A" else "N/A"}</div><div class="l">NAS usage</div></div>
</div>
<div class="card"><h2>📁 Directories</h2>
<table>{dirs_html}</table>
</div>
<div class="card"><h2>📋 Backup Log</h2>
<div class="log">{log_html}</div>
</div>
<div class="card"><h2>ℹ️ Info</h2>
<table>
<tr><td style="color:#94a3b8">NAS</td><td>100.97.192.84 (Tailscale)</td></tr>
<tr><td style="color:#94a3b8">Protocol</td><td>rsync over SSH</td></tr>
<tr><td style="color:#94a3b8">Schedule</td><td>Daily 06:00 Kyiv (03:00 UTC)</td></tr>
<tr><td style="color:#94a3b8">Script</td><td><code>{backup_script}</code></td></tr>
<tr><td style="color:#94a3b8">Log</td><td><code>{log_path}</code></td></tr>
</table>
</div>
<p style="text-align:center;margin-top:1rem"><a href="/">← Back to BiometalDB</a></p>
</div></body></html>"""


@app.route("/")
def index():
    n_mols = len([f for f in os.listdir(MOL_DIR) if f.endswith(".mol")])
    return f"""<html><head><title>BiometalDB</title>
    <link rel="icon" href="/static/favicon.ico" type="image/x-icon">
    <style>body{{font-family:sans-serif;max-width:600px;margin:40px auto;padding:0 1rem}}
    a{{color:#2563eb;text-decoration:none}}a:hover{{text-decoration:underline}}</style></head>
    <body><h2>BiometalDB MOL Server</h2><p>{n_mols} MOLfiles available.</p>
    <ul>
    <li><a href="/complexes">Browse complexes</a> (with 2D structures)</li>
    <li><a href="/dmpnn">📊 D-MPNN Scoring</a> (Ru donor prediction results)</li>
    <li><a href="/viewer/">🔬 3D MOL Viewer</a> (interactive 3D)</li>
    <li><a href="/backup">💾 NAS Backup</a> — status & logs</li>
    <li><code>GET /mol3/{{id}}</code> — download MOL</li>
    <li><code>GET /complexes/{{id}}</code> — detail view with structure</li>
    </ul>
    <p>Datasette: <a href="http://localhost:8501/biometaldb/v_compounds">localhost:8501</a></p>
    </body></html>"""


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8502, debug=False)
