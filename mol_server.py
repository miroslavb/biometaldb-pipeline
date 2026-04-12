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


def render_2d(cid):
    """Render MOL to 2D PNG, return base64 string."""
    mol_path = os.path.join(MOL_DIR, f"complex_{cid}.mol")
    if not os.path.exists(mol_path):
        return None
    try:
        mol = Chem.MolFromMolFile(mol_path, sanitize=False, removeHs=False)
        if mol is None:
            return None
        try:
            AllChem.Compute2DCoords(mol)
        except:
            pass
        drawer = Draw.MolDraw2DCairo(520, 400)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().bondLineWidth = 1.5
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode()
    except:
        return None


HTML_PAGE = """
<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Complex #{{ cid }} — BiometalDB</title>
<style>
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#f0f2f5}
.hdr{background:linear-gradient(135deg,#0f172a,#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid #3b82f6}
.hdr h1{margin:0;font-size:1.4rem}.hdr a{color:#60a5fa;text-decoration:none;margin-right:1rem;font-size:0.85rem}
.wrap{max-width:1100px;margin:1.5rem auto;display:grid;grid-template-columns:520px 1fr;gap:1.5rem;padding:0 1.5rem}
.card{background:#fff;border-radius:8px;padding:1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08)}
.card h2{font-size:1rem;color:#1e40af;margin:0 0 .8rem;border-bottom:2px solid #e2e8f0;padding-bottom:.4rem}
.str{text-align:center}.str img{max-width:100%;border:1px solid #e2e8f0;border-radius:6px}
.dl{margin-top:1rem;display:flex;gap:.75rem;flex-wrap:wrap}
.btn{display:inline-block;padding:.5rem 1rem;border-radius:6px;text-decoration:none;font-size:.85rem;font-weight:500}
.bp{background:#3b82f6;color:#fff}.bp:hover{background:#2563eb}
.bs{background:#e2e8f0;color:#334155}.bs:hover{background:#cbd5e1}
.mt table{width:100%;border-collapse:collapse;font-size:.85rem}
.mt td{padding:.35rem .5rem;border-bottom:1px solid #f1f5f9}
.mt td:first-child{font-weight:600;color:#475569;width:140px}
.tuc{font-family:'JetBrains Mono',monospace;font-size:.72rem;color:#64748b;word-break:break-all;background:#f8fafc;padding:.5rem;border-radius:4px;margin-top:.5rem}
</style></head><body>
<div class="hdr">
<div><a href="/">home</a><a href="/complexes">complexes</a><a href="http://38.19.202.28:8501/biometaldb/complexes">Datasette</a></div>
<h1>Complex #{{ cid }} — {{ metal }}({{ ox }})</h1></div>
<div class="wrap">
<div class="card str"><h2>2D Structure</h2>
{% if img %}<img src="data:image/png;base64,{{ img }}" alt="Structure">
{% else %}<p style="color:#94a3b8">Could not render.</p>{% endif %}
<div class="dl"><a class="btn bp" href="/mol3/{{ cid }}">⬇ MOL V3000</a>
<a class="btn bs" href="/mol3/{{ cid }}/view">Raw MOL</a>
<a class="btn bs" href="/mol3/{{ cid }}/smiles">SMILES JSON</a></div></div>
<div>
<div class="card mt"><h2>Info</h2><table>
<tr><td>ID</td><td>{{ cid }}</td></tr>
<tr><td>Metal</td><td>{{ metal }}</td></tr>
<tr><td>Oxidation state</td><td>{{ ox }}</td></tr>
<tr><td>Charge</td><td>{{ charge }}</td></tr>
<tr><td>Donor atoms</td><td>{{ donors or '—' }}</td></tr>
<tr><td>Abbreviation</td><td>{{ abbr or '—' }}</td></tr>
<tr><td>SMILES</td><td style="font-family:monospace;font-size:.72rem">{{ smi[:200] }}{% if smi|length>200 %}...{% endif %}</td></tr>
</table></div>
<div class="card mt" style="margin-top:1rem"><h2>TUCAN</h2>
<div class="tuc">{{ tucan or 'Not generated' }}</div></div>
{% if meas %}
<div class="card mt" style="margin-top:1rem"><h2>Measurements ({{ meas|length }})</h2>
<table><tr><td>Cell line</td><td>IC50 dark (µM)</td><td>IC50 light (µM)</td><td>DOI</td><td>Year</td></tr>
{% for m in meas %}<tr>
<td>{{ m[0] or '—' }}</td><td>{{ m[1] or '—' }}</td><td>{{ m[2] or '—' }}</td>
<td>{% if m[3] %}<a href="https://doi.org/{{ m[3] }}">{{ m[3][:30] }}{% if m[3]|length>30 %}...{% endif %}</a>{% else %}—{% endif %}</td>
<td>{{ m[4] or '—' }}</td></tr>{% endfor %}</table></div>{% endif %}
</div></div></body></html>"""


@app.route("/complexes")
def list_complexes():
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""SELECT c.id, c.metal, c.oxidation_state, c.donor_atoms,
        (SELECT COUNT(*) FROM measurements WHERE complex_id=c.id) as n
        FROM complexes c ORDER BY c.id DESC LIMIT 100""")
    rows = cur.fetchall()
    conn.close()
    h = '<!DOCTYPE html><html><head><title>Complexes</title><style>body{font-family:sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem}table{width:100%;border-collapse:collapse}th,td{padding:.4rem;text-align:left;border-bottom:1px solid #eee}th{background:#f1f5f9;font-size:.8rem}a{color:#2563eb;text-decoration:none}</style></head><body><h1>Last 100 complexes</h1><table><tr><th>ID</th><th>Metal</th><th>Ox</th><th>Donors</th><th>Meas.</th></tr>'
    for cid, m, o, d, n in rows:
        h += f'<tr><td><a href="/complexes/{cid}">{cid}</a></td><td>{m}</td><td>{o}</td><td>{d or "—"}</td><td>{n}</td></tr>'
    return h + '</table></body></html>'


@app.route("/complexes/<int:cid>")
def complex_detail(cid):
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""SELECT smiles_ligands, metal, oxidation_state, charge_complex,
        donor_atoms, tucan FROM complexes WHERE id=?""", (cid,))
    row = cur.fetchone()
    if not row:
        conn.close(); abort(404)
    smi, metal, ox, charge, donors, tucan = row
    cur.execute("SELECT abbreviation FROM measurements WHERE complex_id=? AND abbreviation IS NOT NULL LIMIT 1", (cid,))
    abbr = cur.fetchone()
    cur.execute("""SELECT cell_line, ic50_dark, ic50_light, doi, year
        FROM measurements WHERE complex_id=? ORDER BY ic50_dark ASC NULLS LAST""", (cid,))
    meas = cur.fetchall()
    conn.close()
    img = render_2d(cid)
    return render_template_string(HTML_PAGE, cid=cid, metal=metal, ox=ox,
        charge=charge, donors=donors, smi=smi, tucan=tucan,
        abbr=abbr[0] if abbr else None, img=img, meas=meas)


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


DMPNN_PAGE = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>D-MPNN Scoring — BiometalDB</title>
<style>
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#f0f2f5}
.hdr{background:linear-gradient(135deg,#0f172a,#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid #8b5cf6}
.hdr h1{margin:0;font-size:1.4rem}.hdr a{color:#a78bfa;text-decoration:none;margin-right:1rem;font-size:.85rem}
.wrap{max-width:1000px;margin:1.5rem auto;padding:0 1.5rem}
.stats{display:flex;gap:1rem;margin-bottom:1.5rem}
.stat{background:#fff;border-radius:8px;padding:1rem 1.5rem;box-shadow:0 1px 3px rgba(0,0,0,.08);flex:1;text-align:center}
.stat .n{font-size:1.8rem;font-weight:700}.stat .l{font-size:.8rem;color:#64748b;margin-top:.2rem}
.card{background:#fff;border-radius:8px;padding:1rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:.8rem}
.card a{text-decoration:none}
code{background:#f1f5f9;padding:.1rem .3rem;border-radius:3px;font-size:.82rem}
</style></head><body>
<div class="hdr"><div><a href="/">home</a><a href="/complexes">complexes</a><a href="/dmpnn">D-MPNN</a></div>
<h1>D-MPNN Donor Prediction Results</h1></div>
<div class="wrap">
<div class="stats">
<div class="stat"><div class="n" style="color:#16a34a">{{ n_match }}</div><div class="l">Exact match</div></div>
<div class="stat"><div class="n" style="color:#d97706">{{ n_partial }}</div><div class="l">Partial match</div></div>
<div class="stat"><div class="n" style="color:#dc2626">{{ n_fail }}</div><div class="l">No match</div></div>
<div class="stat"><div class="n">{{ "%.2f"|format(avg) }}</div><div class="l">Mean score</div></div>
</div>
<h3 style="color:#334155">{{ total }} Ru complexes scored</h3>
{{ cards|safe }}
</div></body></html>"""

DMPNN_DETAIL_PAGE = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Complex #{{ cid }} D-MPNN — BiometalDB</title>
<style>
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#f0f2f5}
.hdr{background:linear-gradient(135deg,#0f172a,#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid #8b5cf6}
.hdr h1{margin:0;font-size:1.3rem}.hdr a{color:#a78bfa;text-decoration:none;margin-right:1rem;font-size:.85rem}
.hdr .info{margin-top:.3rem;font-size:.85rem;color:#c4b5fd}
.wrap{max-width:1100px;margin:1.5rem auto;padding:0 1.5rem}
.meta{background:#fff;border-radius:8px;padding:1rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:1.5rem;display:flex;gap:2rem;align-items:center;flex-wrap:wrap}
.meta .item{text-align:center}.meta .val{font-size:1.4rem;font-weight:700}.meta .lab{font-size:.75rem;color:#64748b}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:1rem}
.cand{background:#fff;border-radius:8px;padding:.8rem;box-shadow:0 1px 3px rgba(0,0,0,.08)}
.cand img{width:100%;border-radius:4px}
.cand .smi{font-family:monospace;font-size:.68rem;color:#64748b;word-break:break-all;margin-top:.4rem;background:#f8fafc;padding:.3rem;border-radius:3px}
</style></head><body>
<div class="hdr"><div><a href="/">home</a><a href="/dmpnn">← all results</a><a href="/complexes/{{ cid }}">complex #{{ cid }}</a></div>
<h1>Complex #{{ cid }} — {{ metal }}({{ ox }})</h1>
<div class="info">GT: {{ gt }} → Top-1: {{ pred }}</div></div>
<div class="wrap">
<div class="meta">
<div class="item"><div class="val">{{ nc }}</div><div class="lab">Candidates</div></div>
<div class="item"><div class="val">{{ gt }}</div><div class="lab">Ground Truth</div></div>
<div class="item"><div class="val">{{ pred }}</div><div class="lab">Top-1 Prediction</div></div>
<div class="item"><div class="val" style="color:{% if match==1 %}#16a34a{% elif match>0 %}#d97706{% else %}#dc2626{% endif %}">{{ "%.2f"|format(match) }}</div><div class="lab">Match Score</div></div>
</div>
<div class="grid">{{ cand_html|safe }}</div>
</div></body></html>"""


@app.route("/dmpnn")
def dmpnn_summary():
    """D-MPNN scoring results overview."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""SELECT complex_id, metal, oxidation_state, donor_atoms_gt, 
        top1_donor_pred, match_score, n_candidates
        FROM dmpnn_summary ORDER BY match_score DESC, complex_id""")
    rows = cur.fetchall()
    conn.close()
    
    def match_color(m):
        if m == 1.0: return '#16a34a'
        elif m > 0: return '#d97706'
        return '#dc2626'
    
    def match_badge(m):
        c = match_color(m)
        return f'<span style="background:{c};color:#fff;padding:.2rem .6rem;border-radius:4px;font-weight:600">{m:.2f}</span>'
    
    cards = ''
    for cid, metal, ox, gt, pred, match, nc in rows:
        cards += f'''<div class="card" style="border-left:4px solid {match_color(match)}">
        <div style="display:flex;justify-content:space-between;align-items:center">
        <div>
        <a href="/complexes/{cid}"><b>#{cid}</b></a> {metal}({ox}) — {nc} candidates
        <div style="margin-top:.3rem;font-size:.85rem">
        <span style="color:#64748b">Ground truth:</span> <code>{gt}</code>
        <span style="color:#64748b;margin-left:.5rem">Top-1:</span> <code>{pred}</code>
        </div></div>
        <div>{match_badge(match)}</div></div>
        <a href="/dmpnn/{cid}" style="display:inline-block;margin-top:.5rem;font-size:.8rem;color:#3b82f6">View candidates →</a></div>'''
    
    n_match = sum(1 for r in rows if r[5] == 1.0)
    n_partial = sum(1 for r in rows if 0 < r[5] < 1.0)
    n_fail = sum(1 for r in rows if r[5] == 0.0)
    avg = sum(r[5] for r in rows) / len(rows) if rows else 0
    
    return render_template_string(DMPNN_PAGE, cards=cards, total=len(rows),
        n_match=n_match, n_partial=n_partial, n_fail=n_fail, avg=avg)


@app.route("/dmpnn/<int:cid>")
def dmpnn_detail(cid):
    """D-MPNN scoring detail: candidates with rendered structures."""
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("""SELECT metal, oxidation_state, donor_atoms_gt, top1_donor_pred, match_score, n_candidates
        FROM dmpnn_summary WHERE complex_id=?""", (cid,))
    summ = cur.fetchone()
    if not summ:
        conn.close(); abort(404)
    
    cur.execute("""SELECT rank, candidate_smi, metal_ox, is_top1 
        FROM dmpnn_results WHERE complex_id=? ORDER BY rank""", (cid,))
    candidates = cur.fetchall()
    conn.close()
    
    metal, ox, gt, pred, match, nc = summ
    
    # Render each candidate
    cand_html = ''
    for rank, smi, m_ox, is_top1 in candidates:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        img_b64 = None
        if mol:
            try:
                AllChem.Compute2DCoords(mol)
                drawer = Draw.MolDraw2DCairo(380, 280)
                drawer.drawOptions().addStereoAnnotation = True
                drawer.drawOptions().bondLineWidth = 2
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                img_b64 = base64.b64encode(drawer.GetDrawingText()).decode()
            except:
                pass
        
        border = '#16a34a' if is_top1 else '#e2e8f0'
        badge = '<span style="background:#16a34a;color:#fff;padding:.1rem .4rem;border-radius:3px;font-size:.75rem">TOP-1</span>' if is_top1 else ''
        
        cand_html += f'''<div class="cand" style="border:2px solid {border}">
        <div style="font-size:.8rem;color:#64748b;margin-bottom:.3rem">#{rank} {badge}</div>
        {'<img src="data:image/png;base64,' + img_b64 + '">' if img_b64 else '<div style="color:#94a3b8;text-align:center;padding:2rem">Render failed</div>'}
        <div class="smi">{smi[:80]}{'...' if len(smi)>80 else ''}</div></div>'''
    
    return render_template_string(DMPNN_DETAIL_PAGE, cid=cid, metal=metal, ox=ox,
        gt=gt, pred=pred, match=match, nc=nc, cand_html=cand_html)


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
