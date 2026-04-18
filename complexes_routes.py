#!/usr/bin/env python3
"""Enhanced complexes routes for BiometalDB Flask server.
Pagination, filters, sorting, 3D structure downloads.
"""
import os
import json
import sqlite3
import base64
from flask import request, render_template_string, abort, send_file, Response
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "biometaldb.sqlite")
MOL3_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mol3")
STRUCT3D_DIR = "/root/ir_docking_structures"

# ─── CSS (shared with D-MPNN) ────────────────────────────────────────────────
COMPLEXES_CSS = """
<link rel="icon" href="/static/favicon.ico" type="image/x-icon">
<style>
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#f0f2f5;color:#1e293b}
.hdr{background:linear-gradient(135deg,#0f172a,#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid #3b82f6}
.hdr h1{margin:0;font-size:1.4rem}.hdr a{color:#60a5fa;text-decoration:none;margin-right:1rem;font-size:.85rem}
.hdr .sub{margin-top:.3rem;font-size:.85rem;color:#93c5fd}
.wrap{max-width:1200px;margin:1.5rem auto;padding:0 1.5rem}
.stats{display:flex;gap:1rem;margin-bottom:1.5rem;flex-wrap:wrap}
.stat{background:#fff;border-radius:8px;padding:1rem 1.5rem;box-shadow:0 1px 3px rgba(0,0,0,.08);flex:1;text-align:center;min-width:120px}
.stat .n{font-size:1.8rem;font-weight:700}.stat .l{font-size:.8rem;color:#64748b;margin-top:.2rem}
.filters{background:#fff;border-radius:8px;padding:1rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:1.5rem}
.filters form{display:flex;gap:.75rem;flex-wrap:wrap;align-items:flex-end}
.fg{display:flex;flex-direction:column;gap:.2rem}
.fg label{font-size:.75rem;color:#64748b;font-weight:600;text-transform:uppercase}
.fg select,.fg input{padding:.4rem .6rem;border:1px solid #e2e8f0;border-radius:6px;font-size:.85rem;background:#fff}
.fg select:focus,.fg input:focus{border-color:#3b82f6;outline:none}
.btn{padding:.4rem 1rem;border-radius:6px;font-size:.85rem;font-weight:500;cursor:pointer;border:none}
.bp{background:#3b82f6;color:#fff}.bp:hover{background:#2563eb}
.bs{background:#e2e8f0;color:#334155}.bs:hover{background:#cbd5e1}
.card{background:#fff;border-radius:8px;padding:1rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:.8rem}
.card a{text-decoration:none}
code{background:#f1f5f9;padding:.1rem .3rem;border-radius:3px;font-size:.82rem}
.pager{display:flex;gap:.5rem;justify-content:center;margin:1.5rem 0;align-items:center;flex-wrap:wrap}
.pager a,.pager span{padding:.4rem .8rem;border-radius:6px;font-size:.85rem;text-decoration:none}
.pager a{background:#fff;color:#3b82f6;box-shadow:0 1px 2px rgba(0,0,0,.08)}.pager a:hover{background:#eff6ff}
.pager .cur{background:#3b82f6;color:#fff;font-weight:600}
.pager .dis{color:#94a3b8}
.dl3d{display:flex;gap:.4rem;flex-wrap:wrap;margin-top:.4rem}
.dl3d a{font-size:.75rem;padding:.2rem .5rem;border-radius:4px;text-decoration:none;border:1px solid #e2e8f0;color:#3b82f6}
.dl3d a:hover{background:#eff6ff}
.badge{display:inline-block;padding:.15rem .5rem;border-radius:4px;font-size:.7rem;font-weight:600}
.badge-ir{background:#dbeafe;color:#1e40af}
.badge-ru{background:#fce7f3;color:#9d174d}
.badge-rh{background:#fef3c7;color:#92400e}
.badge-os{background:#e0e7ff;color:#3730a3}
.badge-re{background:#d1fae5;color:#065f46}
.badge-other{background:#f1f5f9;color:#475569}
</style>
"""

# ─── List Page ────────────────────────────────────────────────────────────────
COMPLEXES_LIST = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Complexes — BiometalDB</title>
""" + COMPLEXES_CSS + """</head><body>
<div class="hdr"><div><a href="/">home</a><a href="/complexes" style="color:#93c5fd;font-weight:600">complexes</a><a href="/dmpnn">D-MPNN</a></div>
<h1>Coordination Complexes</h1>
<div class="sub">{{ total }} complexes in database</div></div>
<div class="wrap">
<!-- Stats -->
<div class="stats">
<div class="stat"><div class="n">{{ total }}</div><div class="l">Total</div></div>
{% for metal, cnt in metal_stats %}
<div class="stat"><div class="n">{{ cnt }}</div><div class="l">{{ metal }}</div></div>
{% endfor %}
</div>
<!-- Filters -->
<div class="filters"><form method="get">
<div class="fg"><label>Metal</label>
<select name="metal"><option value="">All</option>
{% for m in metals %}<option value="{{ m }}" {{ 'selected' if m==f_metal else '' }}>{{ m }}</option>{% endfor %}
</select></div>
<div class="fg"><label>Ox. State</label>
<select name="ox"><option value="">All</option>
{% for o in ox_states %}<option value="{{ o }}" {{ 'selected' if o|string==f_ox else '' }}>{{ o }}</option>{% endfor %}
</select></div>
<div class="fg"><label>Donor atoms</label>
<input name="donors" value="{{ f_donors }}" placeholder="e.g. N, O, P"></div>
<div class="fg"><label>Has 3D struct</label>
<select name="has_3d"><option value="">Any</option>
<option value="1" {{ 'selected' if f_has3d=='1' else '' }}>Yes</option>
</select></div>
<div class="fg"><label>Search (ID / SMILES)</label>
<input name="q" value="{{ f_q }}" placeholder="#123 or c1ccn..." style="width:200px"></div>
<div class="fg"><label>Sort</label>
<select name="sort">
<option value="id_asc" {{ 'selected' if sort=='id_asc' else '' }}>ID ↑</option>
<option value="id_desc" {{ 'selected' if sort=='id_desc' else '' }}>ID ↓</option>
<option value="meas_desc" {{ 'selected' if sort=='meas_desc' else '' }}>Measurements ↓</option>
</select></div>
<div class="fg"><label>Per page</label>
<select name="per_page">
{% for pp in [50,100,500] %}<option value="{{ pp }}" {{ 'selected' if pp==per_page else '' }}>{{ pp }}</option>{% endfor %}
</select></div>
<div class="fg"><label>&nbsp;</label><button class="btn bp" type="submit">Apply</button></div>
<div class="fg"><label>&nbsp;</label><a class="btn bs" href="/complexes">Reset</a></div>
</form></div>
<!-- Results -->
{% for cid,metal,ox,smi,donors,n_meas,has_mol3,has_3d in rows %}
<div class="card" style="border-left:4px solid {{ metal_color(metal) }}">
<div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:.5rem">
<div style="flex:1;min-width:300px">
<a href="/complexes/{{ cid }}"><b>#{{ cid }}</b></a>
<span class="badge badge-{{ metal.lower() if metal.lower() in ['ir','ru','rh','os','re'] else 'other' }}">{{ metal }}({{ ox }})</span>
{% if donors and donors != 'None' %}<span style="color:#64748b;font-size:.85rem"> donors: {{ donors }}</span>{% endif %}
<div style="margin-top:.3rem;font-size:.82rem;color:#475569;font-family:monospace">{{ smi[:120] }}{% if smi|length>120 %}...{% endif %}</div>
</div>
<div style="text-align:right">
<div style="font-size:.8rem;color:#64748b">{{ n_meas }} meas.</div>
<div class="dl3d">
{% if has_mol3 %}<a href="/mol3/{{ cid }}">MOL</a>{% endif %}
{% if has_3d %}<a href="/structures/{{ cid }}/sdf">SDF</a><a href="/structures/{{ cid }}/pdb">PDB</a><a href="/structures/{{ cid }}/xyz">XYZ</a>{% endif %}
<a href="/viewer/?id={{ cid }}">3D</a>
</div></div></div></div>
{% endfor %}
<!-- Pagination -->
<div class="pager">
{% if page > 0 %}<a href="{{ pager_url(0) }}">«</a><a href="{{ pager_url(page-1) }}">‹</a>{% endif %}
{% for p in page_range %}<a href="{{ pager_url(p) }}" class="{{ 'cur' if p==page else '' }}">{{ p+1 }}</a>{% endfor %}
{% if page < total_pages-1 %}<a href="{{ pager_url(page+1) }}">›</a><a href="{{ pager_url(total_pages-1) }}">»</a>{% endif %}
<span class="dis">{{ page+1 }}/{{ total_pages }} ({{ total }} total)</span>
</div>
</div></body></html>"""


# ─── Detail Page ──────────────────────────────────────────────────────────────
COMPLEXES_DETAIL = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Complex #{{ cid }} — BiometalDB</title>
""" + COMPLEXES_CSS + """
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
.side{display:grid;grid-template-columns:520px 1fr;gap:1.5rem}
.str{text-align:center}.str img{max-width:100%;border:1px solid #e2e8f0;border-radius:6px}
.mt table{width:100%;border-collapse:collapse;font-size:.85rem}
.mt td{padding:.35rem .5rem;border-bottom:1px solid #f1f5f9}
.mt td:first-child{font-weight:600;color:#475569;width:140px}
.tuc{font-family:monospace;font-size:.72rem;color:#64748b;word-break:break-all;background:#f8fafc;padding:.5rem;border-radius:4px;margin-top:.5rem}
.viewer3d{background:#0d1117;border-radius:8px;overflow:hidden;margin-top:1.5rem;box-shadow:0 1px 3px rgba(0,0,0,.12)}
.viewer3d .vhdr{background:#161b22;padding:.6rem 1rem;color:#c9d1d9;font-size:.85rem;display:flex;justify-content:space-between;align-items:center;border-bottom:1px solid #30363d}
.viewer3d .vhdr b{color:#58a6ff}
.viewer3d .vctrl{display:flex;gap:.4rem}
.viewer3d .vctrl button{padding:.25rem .6rem;background:#21262d;border:1px solid #30363d;border-radius:4px;color:#c9d1d9;cursor:pointer;font-size:.75rem}
.viewer3d .vctrl button:hover{background:#30363d}
.viewer3d .vctrl button.active{background:#1c3a5e;border-color:#58a6ff}
#mol3d{width:100%;height:420px;position:relative}
</style></head><body>
<div class="hdr"><div><a href="/">home</a><a href="/complexes">← all complexes</a>
<a href="/dmpnn/{{ cid }}">D-MPNN</a></div>
<h1>Complex #{{ cid }} — {{ metal }}({{ ox }}){% if abbr %} <span style="color:#93c5fd">[{{ abbr }}]</span>{% endif %}</h1></div>
<div class="wrap">
<div class="side">
<!-- Structure -->
<div class="card str"><h2 style="font-size:1rem;color:#1e40af;margin:0 0 .8rem;border-bottom:2px solid #e2e8f0;padding-bottom:.4rem">2D Structure</h2>
{% if img %}<img src="data:image/png;base64,{{ img }}" alt="Structure">
{% else %}<p style="color:#94a3b8">Could not render.</p>{% endif %}
<div style="margin-top:1rem;display:flex;gap:.5rem;flex-wrap:wrap">
{% if has_mol3 %}<a class="btn bp" href="/mol3/{{ cid }}">⬇ MOL V3000</a>{% endif %}
{% if has_3d %}
<a class="btn bp" href="/structures/{{ cid }}/sdf">⬇ SDF (3D)</a>
<a class="btn bs" href="/structures/{{ cid }}/pdb">⬇ PDB</a>
<a class="btn bs" href="/structures/{{ cid }}/xyz">⬇ XYZ (xTB)</a>
{% endif %}
</div></div>
<!-- Info -->
<div>
<div class="card mt"><h2 style="font-size:1rem;color:#1e40af;margin:0 0 .8rem;border-bottom:2px solid #e2e8f0;padding-bottom:.4rem">Info</h2>
<table>
<tr><td>ID</td><td>{{ cid }}</td></tr>
<tr><td>Metal</td><td>{{ metal }}</td></tr>
<tr><td>Oxidation state</td><td>{{ ox }}</td></tr>
<tr><td>Charge</td><td>{{ charge }}</td></tr>
<tr><td>Donor atoms</td><td>{{ donors or '—' }}</td></tr>
<tr><td>Abbreviation</td><td>{{ abbr or '—' }}</td></tr>
<tr><td>SMILES</td><td style="font-family:monospace;font-size:.72rem">{{ smi[:200] }}{% if smi|length>200 %}...{% endif %}</td></tr>
{% if has_3d %}<tr><td>3D Structure</td><td style="color:#16a34a;font-weight:600">✓ Available (RDKit + xTB GFN2)</td></tr>{% endif %}
</table></div>
{% if tucan %}
<div class="card mt" style="margin-top:1rem"><h2 style="font-size:1rem;color:#1e40af;margin:0 0 .8rem;border-bottom:2px solid #e2e8f0;padding-bottom:.4rem">TUCAN</h2>
<div class="tuc">{{ tucan }}</div></div>{% endif %}
{% if meas %}
<div class="card mt" style="margin-top:1rem"><h2 style="font-size:1rem;color:#1e40af;margin:0 0 .8rem;border-bottom:2px solid #e2e8f0;padding-bottom:.4rem">Measurements ({{ meas|length }})</h2>
<table><tr><td>Cell line</td><td>IC₅₀ dark (µM)</td><td>IC₅₀ light (µM)</td><td>DOI</td><td>Year</td></tr>
{% for m in meas %}<tr>
<td>{{ m[0] or '—' }}</td><td>{{ '%.2g'|format(m[1]) if m[1] else '—' }}</td><td>{{ '%.2g'|format(m[2]) if m[2] else '—' }}</td>
<td>{% if m[3] %}<a href="https://doi.org/{{ m[3] }}">{{ m[3][:30] }}{% if m[3]|length>30 %}...{% endif %}</a>{% else %}—{% endif %}</td>
<td>{{ m[4] or '—' }}</td></tr>{% endfor %}</table></div>{% endif %}
</div></div>
{% if has_3d %}
<div class="viewer3d">
<div class="vhdr"><b>🔬 3D Structure (xTB GFN2)</b>
<div class="vctrl">
<button onclick="setStyle('stick')" id="btn-stick">Stick</button>
<button onclick="setStyle('sphere')" id="btn-sphere">Sphere</button>
<button onclick="setStyle('ball')" id="btn-ball" class="active">Ball+Stick</button>
<button onclick="viewer3d.spin(!viewer3d.spinState())">Spin</button>
</div></div>
<div id="mol3d"></div>
</div>
<script>
let viewer3d = $3Dmol.createViewer("mol3d", {backgroundColor: "#0d1117"});
let curStyle = "ball";
fetch("/api/structures/{{ cid }}/xyz")
  .then(r => r.json())
  .then(d => {
    viewer3d.addModel(d.xyz, "xyz");
    applyStyle();
    viewer3d.zoomTo();
    viewer3d.render();
  })
  .catch(e => document.getElementById("mol3d").innerHTML = '<p style="color:#f87171;padding:2rem">Failed to load 3D structure</p>');
function applyStyle() {
  viewer3d.setStyle({}, {});
  if (curStyle === "stick") viewer3d.setStyle({}, {stick: {radius: 0.12, colorscheme: "Jmol"}});
  else if (curStyle === "sphere") viewer3d.setStyle({}, {sphere: {scale: 0.4, colorscheme: "Jmol"}});
  else viewer3d.setStyle({}, {stick: {radius: 0.1, colorscheme: "Jmol"}, sphere: {scale: 0.25, colorscheme: "Jmol"}});
  viewer3d.render();
}
function setStyle(s) {
  curStyle = s;
  document.querySelectorAll(".vctrl button").forEach(b => b.classList.remove("active"));
  document.getElementById("btn-" + s).classList.add("active");
  applyStyle();
}
</script>
{% endif %}
</div></body></html>"""


def metal_color(m):
    colors = {'Ir': '#3b82f6', 'Ru': '#ec4899', 'Rh': '#f59e0b', 'Os': '#6366f1', 'Re': '#10b981'}
    return colors.get(m, '#94a3b8')


def render_2d(cid):
    """Render MOL to 2D PNG, return base64 string."""
    mol_path = os.path.join(MOL3_DIR, f"complex_{cid}.mol")
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


def has_3d_structure(cid):
    """Check if 3D structure files exist for this complex."""
    return any(os.path.exists(os.path.join(STRUCT3D_DIR, f"ir_{cid}.{ext}"))
               for ext in ['sdf', 'pdb'])


def register_complexes_routes(app):
    """Register enhanced complexes routes on the Flask app."""

    @app.route("/complexes")
    def complexes_list():
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()

        # Params
        per_page = int(request.args.get('per_page', 50))
        page = int(request.args.get('page', 0))
        f_metal = request.args.get('metal', '')
        f_ox = request.args.get('ox', '')
        f_donors = request.args.get('donors', '')
        f_has3d = request.args.get('has_3d', '')
        f_q = request.args.get('q', '').strip()
        sort = request.args.get('sort', 'id_desc')

        # Build query
        where = []
        params = []

        if f_metal:
            where.append("c.metal = ?")
            params.append(f_metal)
        if f_ox:
            where.append("c.oxidation_state = ?")
            params.append(int(f_ox))
        if f_donors:
            for d in f_donors.split(','):
                d = d.strip()
                if d:
                    where.append("c.donor_atoms LIKE ?")
                    params.append(f'%"{d}"%')
        if f_has3d == '1':
            ids_3d = set()
            if os.path.isdir(STRUCT3D_DIR):
                for fname in os.listdir(STRUCT3D_DIR):
                    if fname.startswith('ir_') and (fname.endswith('.sdf') or fname.endswith('.pdb')):
                        try:
                            ids_3d.add(int(fname.split('_')[1].split('.')[0]))
                        except (ValueError, IndexError):
                            pass
            if ids_3d:
                placeholders = ','.join('?' * len(ids_3d))
                where.append(f"c.id IN ({placeholders})")
                params.extend(sorted(ids_3d))
            else:
                where.append("1=0")
        if f_q:
            if f_q.startswith('#'):
                try:
                    where.append("c.id = ?")
                    params.append(int(f_q[1:]))
                except ValueError:
                    pass
            else:
                where.append("(c.smiles_ligands LIKE ? OR c.metal LIKE ? OR CAST(c.id AS TEXT) LIKE ?)")
                params.extend([f'%{f_q}%', f'%{f_q}%', f'%{f_q}%'])

        where_sql = ("WHERE " + " AND ".join(where)) if where else ""

        # Sort
        order_map = {
            'id_asc': 'c.id ASC',
            'id_desc': 'c.id DESC',
            'meas_desc': 'n_meas DESC, c.id DESC',
        }
        order_sql = order_map.get(sort, 'c.id DESC')

        # Count
        cur.execute(f"""
            SELECT COUNT(*) FROM complexes c
            LEFT JOIN (SELECT complex_id, COUNT(*) as n FROM measurements GROUP BY complex_id) m
            ON m.complex_id = c.id
            {where_sql}
        """, params)
        total = cur.fetchone()[0]

        # Pagination
        total_pages = max(1, (total + per_page - 1) // per_page)
        page = min(page, total_pages - 1)
        offset = page * per_page

        # Fetch
        cur.execute(f"""
            SELECT c.id, c.metal, c.oxidation_state, c.smiles_ligands, c.donor_atoms,
                   COALESCE(m.n, 0) as n_meas, c.has_mol3
            FROM complexes c
            LEFT JOIN (SELECT complex_id, COUNT(*) as n FROM measurements GROUP BY complex_id) m
            ON m.complex_id = c.id
            {where_sql}
            ORDER BY {order_sql}
            LIMIT ? OFFSET ?
        """, params + [per_page, offset])
        raw_rows = cur.fetchall()

        # Add 3D availability flag
        rows = []
        for r in raw_rows:
            cid = r[0]
            has3d = has_3d_structure(cid)
            rows.append((*r, has3d))

        # Filter by has_3d is handled in WHERE clause above

        # Metals for filter
        cur.execute("SELECT DISTINCT metal FROM complexes ORDER BY metal")
        metals = [r[0] for r in cur.fetchall()]

        # Ox states
        cur.execute("SELECT DISTINCT oxidation_state FROM complexes WHERE oxidation_state IS NOT NULL ORDER BY oxidation_state")
        ox_states = [r[0] for r in cur.fetchall()]

        # Metal stats
        cur.execute("SELECT metal, COUNT(*) FROM complexes GROUP BY metal ORDER BY COUNT(*) DESC")
        metal_stats = cur.fetchall()

        conn.close()

        # Pagination URL builder
        def pager_url(p):
            parts = [f"page={p}", f"per_page={per_page}"]
            if f_metal: parts.append(f"metal={f_metal}")
            if f_ox: parts.append(f"ox={f_ox}")
            if f_donors: parts.append(f"donors={f_donors}")
            if f_has3d: parts.append(f"has_3d={f_has3d}")
            if f_q: parts.append(f"q={f_q}")
            if sort != 'id_desc': parts.append(f"sort={sort}")
            return "?" + "&".join(parts)

        page_range = list(range(max(0, page - 3), min(total_pages, page + 4)))

        return render_template_string(COMPLEXES_LIST,
            rows=rows, total=total, page=page, per_page=per_page,
            total_pages=total_pages, page_range=page_range, pager_url=pager_url,
            metals=metals, ox_states=ox_states, metal_stats=metal_stats,
            f_metal=f_metal, f_ox=f_ox, f_donors=f_donors, f_has3d=f_has3d,
            f_q=f_q, sort=sort, metal_color=metal_color)

    @app.route("/complexes/<int:cid>")
    def complexes_detail(cid):
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
        has_mol3 = os.path.exists(os.path.join(MOL3_DIR, f"complex_{cid}.mol"))
        has_3d = has_3d_structure(cid)

        return render_template_string(COMPLEXES_DETAIL,
            cid=cid, metal=metal, ox=ox, charge=charge, donors=donors,
            smi=smi, tucan=tucan, abbr=abbr[0] if abbr else None,
            img=img, meas=meas, has_mol3=has_mol3, has_3d=has_3d)

    # ─── 3D Structure download routes ─────────────────────────────────────
    @app.route("/structures/<int:cid>/sdf")
    def download_sdf(cid):
        path = os.path.join(STRUCT3D_DIR, f"ir_{cid}.sdf")
        if not os.path.exists(path):
            abort(404)
        return send_file(path, mimetype="chemical/x-mdl-sdfile",
                         as_attachment=True, download_name=f"ir_{cid}.sdf")

    @app.route("/structures/<int:cid>/pdb")
    def download_pdb(cid):
        path = os.path.join(STRUCT3D_DIR, f"ir_{cid}.pdb")
        if not os.path.exists(path):
            abort(404)
        return send_file(path, mimetype="chemical/x-pdb",
                         as_attachment=True, download_name=f"ir_{cid}.pdb")

    @app.route("/structures/<int:cid>/xyz")
    def download_xyz(cid):
        path = os.path.join(STRUCT3D_DIR, f"ir_{cid}_opt.xyz")
        if not os.path.exists(path):
            abort(404)
        return send_file(path, mimetype="chemical/x-xyz",
                         as_attachment=True, download_name=f"ir_{cid}_opt.xyz")

    @app.route("/structures/<int:cid>/initial")
    def download_xyz_initial(cid):
        path = os.path.join(STRUCT3D_DIR, f"ir_{cid}_initial.xyz")
        if not os.path.exists(path):
            abort(404)
        return send_file(path, mimetype="chemical/x-xyz",
                         as_attachment=True, download_name=f"ir_{cid}_initial.xyz")

    @app.route("/api/structures/<int:cid>/xyz")
    def api_xyz_data(cid):
        """Return XYZ content as JSON for inline 3D viewer."""
        import json as _json
        path = os.path.join(STRUCT3D_DIR, f"ir_{cid}_opt.xyz")
        if not os.path.exists(path):
            abort(404)
        with open(path) as f:
            xyz = f.read()
        return Response(_json.dumps({"cid": cid, "xyz": xyz}),
                        mimetype="application/json")
