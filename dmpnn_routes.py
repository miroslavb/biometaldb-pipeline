#!/usr/bin/env python3
"""Enhanced D-MPNN routes for BiometalDB Flask server.
Replaces old /dmpnn and /dmpnn/<id> routes with full-featured UI.
"""
import os
import json
import sqlite3
import base64
from flask import request, render_template_string, abort
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "biometaldb.sqlite")

# ─── CSS (shared) ───────────────────────────────────────────────────────────
DMPNN_CSS = """
<link rel="icon" href="/static/favicon.ico" type="image/x-icon">
<style>
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#f0f2f5;color:#1e293b}
.hdr{background:linear-gradient(135deg,#0f172a,#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid #8b5cf6}
.hdr h1{margin:0;font-size:1.4rem}.hdr a{color:#a78bfa;text-decoration:none;margin-right:1rem;font-size:.85rem}
.hdr .sub{margin-top:.3rem;font-size:.85rem;color:#c4b5fd}
.wrap{max-width:1200px;margin:1.5rem auto;padding:0 1.5rem}
.stats{display:flex;gap:1rem;margin-bottom:1.5rem;flex-wrap:wrap}
.stat{background:#fff;border-radius:8px;padding:1rem 1.5rem;box-shadow:0 1px 3px rgba(0,0,0,.08);flex:1;text-align:center;min-width:120px}
.stat .n{font-size:1.8rem;font-weight:700}.stat .l{font-size:.8rem;color:#64748b;margin-top:.2rem}
.filters{background:#fff;border-radius:8px;padding:1rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:1.5rem}
.filters form{display:flex;gap:.75rem;flex-wrap:wrap;align-items:flex-end}
.fg{display:flex;flex-direction:column;gap:.2rem}
.fg label{font-size:.75rem;color:#64748b;font-weight:600;text-transform:uppercase}
.fg select,.fg input{padding:.4rem .6rem;border:1px solid #e2e8f0;border-radius:6px;font-size:.85rem;background:#fff}
.fg select:focus,.fg input:focus{border-color:#8b5cf6;outline:none}
.btn{padding:.4rem 1rem;border-radius:6px;font-size:.85rem;font-weight:500;cursor:pointer;border:none}
.bp{background:#8b5cf6;color:#fff}.bp:hover{background:#7c3aed}
.bs{background:#e2e8f0;color:#334155}.bs:hover{background:#cbd5e1}
.card{background:#fff;border-radius:8px;padding:1rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:.8rem}
.card a{text-decoration:none}
code{background:#f1f5f9;padding:.1rem .3rem;border-radius:3px;font-size:.82rem}
.pager{display:flex;gap:.5rem;justify-content:center;margin:1.5rem 0;align-items:center;flex-wrap:wrap}
.pager a,.pager span{padding:.4rem .8rem;border-radius:6px;font-size:.85rem;text-decoration:none}
.pager a{background:#fff;color:#3b82f6;box-shadow:0 1px 2px rgba(0,0,0,.08)}.pager a:hover{background:#eff6ff}
.pager .cur{background:#3b82f6;color:#fff;font-weight:600}
.pager .dis{color:#94a3b8}
.sort-bar{display:flex;gap:.75rem;align-items:center;margin-bottom:1rem;flex-wrap:wrap}
.sort-bar label{font-size:.8rem;color:#64748b}
.vtabs{display:flex;gap:.25rem;margin-bottom:1rem}
.vtabs a{padding:.5rem 1rem;border-radius:6px 6px 0 0;font-size:.85rem;text-decoration:none;color:#64748b;background:#e2e8f0}
.vtabs a.active{background:#8b5cf6;color:#fff;font-weight:600}
.cand{background:#fff;border-radius:8px;padding:.8rem;box-shadow:0 1px 3px rgba(0,0,0,.08)}
.cand img{width:100%;border-radius:4px}
.cand .smi{font-family:monospace;font-size:.68rem;color:#64748b;word-break:break-all;margin-top:.4rem;background:#f8fafc;padding:.3rem;border-radius:3px}
.score-badge{display:inline-flex;align-items:center;gap:.3rem;padding:.2rem .6rem;border-radius:4px;font-weight:600;font-size:.85rem}
.links-bar{margin-top:.5rem;display:flex;gap:.5rem;flex-wrap:wrap}
.links-bar a{font-size:.78rem;color:#3b82f6;text-decoration:none;padding:.2rem .5rem;border:1px solid #e2e8f0;border-radius:4px}
.links-bar a:hover{background:#eff6ff}
</style>
"""

# ─── Summary Page ────────────────────────────────────────────────────────────
DMPNN_LIST = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>D-MPNN Scoring — BiometalDB</title>
""" + DMPNN_CSS + """</head><body>
<div class="hdr"><div><a href="/">home</a><a href="/complexes">complexes</a><a href="/dmpnn">D-MPNN</a></div>
<h1>D-MPNN Donor Prediction</h1>
<div class="sub">{{ ver_label }} — {{ total }} complexes</div></div>
<div class="wrap">
<!-- Version tabs -->
<div class="vtabs">
{% for v,label in versions %}<a href="?version={{ v }}{% if per_page != 50 %}&per_page={{ per_page }}{% endif %}" class="{{ 'active' if v==version else '' }}">{{ label }}</a>{% endfor %}
</div>
<!-- Stats -->
<div class="stats">
<div class="stat"><div class="n" style="color:#16a34a">{{ n_match }}</div><div class="l">Exact (1.0)</div></div>
<div class="stat"><div class="n" style="color:#d97706">{{ n_partial }}</div><div class="l">Partial (>0)</div></div>
<div class="stat"><div class="n" style="color:#dc2626">{{ n_fail }}</div><div class="l">No match</div></div>
<div class="stat"><div class="n">{{ "%.3f"|format(avg) }}</div><div class="l">Mean score</div></div>
</div>
<!-- Filters -->
<div class="filters"><form method="get">
<input type="hidden" name="version" value="{{ version }}">
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
<div class="fg"><label>Score min</label>
<input name="score_min" type="number" step="0.05" min="0" max="1" value="{{ f_score_min }}" style="width:80px"></div>
<div class="fg"><label>Score max</label>
<input name="score_max" type="number" step="0.05" min="0" max="1" value="{{ f_score_max }}" style="width:80px"></div>
<div class="fg"><label>Search (ID / SMILES)</label>
<input name="q" value="{{ f_q }}" placeholder="#123 or c1ccn..." style="width:200px"></div>
<div class="fg"><label>Sort</label>
<select name="sort">
<option value="score_desc" {{ 'selected' if sort=='score_desc' else '' }}>Score ↓</option>
<option value="score_asc" {{ 'selected' if sort=='score_asc' else '' }}>Score ↑</option>
<option value="id_asc" {{ 'selected' if sort=='id_asc' else '' }}>ID ↑</option>
<option value="id_desc" {{ 'selected' if sort=='id_desc' else '' }}>ID ↓</option>
</select></div>
<div class="fg"><label>Per page</label>
<select name="per_page">
{% for pp in [50,100,500] %}<option value="{{ pp }}" {{ 'selected' if pp==per_page else '' }}>{{ pp }}</option>{% endfor %}
</select></div>
<div class="fg"><label>&nbsp;</label><button class="btn bp" type="submit">Apply</button></div>
<div class="fg"><label>&nbsp;</label><a class="btn bs" href="/dmpnn?version={{ version }}">Reset</a></div>
</form></div>
<!-- Results -->
{% for cid,metal,ox,gt,pred,match,nc,abbrevs in rows %}
<div class="card" style="border-left:4px solid {{ match_color(match) }}">
<div style="display:flex;justify-content:space-between;align-items:center">
<div>
<a href="/dmpnn/{{ cid }}?version={{ version }}"><b>#{{ cid }}</b></a>
{% if abbrevs %}<span style="color:#8b5cf6;font-weight:600;margin-left:.3rem">{{ abbrevs }}</span>{% endif %}
{{ metal }}({{ ox }}) — {{ nc }} candidates
<div style="margin-top:.3rem;font-size:.85rem">
<span style="color:#64748b">GT:</span> <code>{{ gt }}</code>
<span style="color:#64748b;margin-left:.5rem">Pred:</span> <code>{{ pred }}</code>
</div></div>
<div><span class="score-badge" style="background:{{ match_color(match) }};color:#fff">{{ "%.2f"|format(match) }}</span></div></div>
<div class="links-bar">
<a href="/dmpnn/{{ cid }}?version={{ version }}">Candidates →</a>
<a href="/complexes/{{ cid }}">Complex detail</a>
<a href="/mol3/{{ cid }}/candidate">MOL</a>
<a href="/viewer/?id={{ cid }}">3D Viewer</a>
</div></div>{% endfor %}
<!-- Pagination -->
<div class="pager">
{% if page > 1 %}<a href="{{ pager_url(0) }}">«</a><a href="{{ pager_url(page-1) }}">‹</a>{% endif %}
{% for p in page_range %}<a href="{{ pager_url(p) }}" class="{{ 'cur' if p==page else '' }}">{{ p+1 }}</a>{% endfor %}
{% if page < total_pages-1 %}<a href="{{ pager_url(page+1) }}">›</a><a href="{{ pager_url(total_pages-1) }}">»</a>{% endif %}
<span class="dis">{{ page+1 }}/{{ total_pages }} ({{ total }} total)</span>
</div>
</div></body></html>"""


# ─── Detail Page ─────────────────────────────────────────────────────────────
DMPNN_DETAIL = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Complex #{{ cid }} D-MPNN — BiometalDB</title>
""" + DMPNN_CSS + """</head><body>
<div class="hdr"><div><a href="/">home</a><a href="/dmpnn?version={{ version }}">← all results</a>
<a href="/complexes/{{ cid }}">complex #{{ cid }}</a></div>
<h1>Complex #{{ cid }} — {{ metal }}({{ ox }}){% if abbrevs %} <span style="color:#c4b5fd">[{{ abbrevs }}]</span>{% endif %}</h1>
<div class="info">GT: {{ gt }} → Top-1: {{ pred }}</div></div>
<div class="wrap">
<!-- Version tabs -->
<div class="vtabs">
{% for v,label in versions %}{% if v != 'compare' %}<a href="?version={{ v }}" class="{{ 'active' if v==version else '' }}">{{ label }}</a>{% endif %}{% endfor %}
</div>
<!-- Meta -->
<div style="display:flex;gap:1rem;flex-wrap:wrap;margin-bottom:1.5rem">
<div class="stat"><div class="n">{{ nc }}</div><div class="l">Candidates</div></div>
<div class="stat"><div class="n">{{ gt }}</div><div class="l">Ground Truth</div></div>
<div class="stat"><div class="n">{{ pred }}</div><div class="l">Top-1 Prediction</div></div>
<div class="stat"><div class="n" style="color:{{ '#16a34a' if match==1 else ('#d97706' if match>0 else '#dc2626') }}">{{ "%.3f"|format(match) }}</div><div class="l">Match Score</div></div>
</div>
<!-- Links -->
<div class="links-bar" style="margin-bottom:1.5rem">
<a href="/mol3/{{ cid }}/candidate">⬇ MOL (coordinated)</a>
<a href="/mol3/{{ cid }}">Raw MOL</a>
<a href="/viewer/?id={{ cid }}">🔬 3D Viewer</a>
<a href="/complexes/{{ cid }}">Complex detail</a>
</div>
<!-- Candidates grid -->
{% if not has_cands %}<div style="background:#fef3c7;color:#92400e;padding:.8rem;border-radius:6px;margin-bottom:1rem;font-size:.85rem">
⚠ Per-candidate scores unavailable for <b>{{ version }}</b>. Showing {{ candidates|length }} candidates from single-ligand scoring.
For sequential results only top-1 prediction is stored.</div>{% endif %}
<div style="display:grid;grid-template-columns:repeat(auto-fill,minmax(360px,1fr));gap:1rem">
{% for rank,smi,is_top1,score,score_std in candidates %}
<div class="cand" style="border:2px solid {{ '#16a34a' if is_top1 else '#e2e8f0' }}">
<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:.4rem">
<span style="font-size:.8rem;color:#64748b">#{{ rank }}</span>
{% if is_top1 %}<span style="background:#16a34a;color:#fff;padding:.1rem .4rem;border-radius:3px;font-size:.75rem">TOP-1</span>{% endif %}
<span class="score-badge" style="background:{{ '#16a34a' if score>0.5 else ('#d97706' if score>0.1 else '#dc2626') }};color:#fff;font-size:.8rem">{{ "%.2f"|format(score) }} ± {{ "%.2f"|format(score_std) }}</span>
</div>
{% if cand_imgs[rank] %}<img src="data:image/png;base64,{{ cand_imgs[rank] }}">
{% else %}<div style="color:#94a3b8;text-align:center;padding:2rem">Render failed</div>{% endif %}
<div class="smi">{{ smi[:100] }}{% if smi|length>100 %}...{% endif %}</div>
</div>{% endfor %}
</div>
</div></body></html>"""


def match_color(m):
    if m >= 0.99: return '#16a34a'
    elif m > 0: return '#d97706'
    return '#dc2626'


COMPARE_PAGE = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Compare v1 vs v2 — BiometalDB</title>
""" + DMPNN_CSS + """
<style>
.diff-pos{color:#16a34a;font-weight:700}.diff-neg{color:#dc2626;font-weight:700}.diff-zero{color:#94a3b8}
.compare-table{width:100%;border-collapse:collapse;font-size:.85rem}
.compare-table th{background:#f1f5f9;padding:.5rem .6rem;text-align:left;font-size:.75rem;color:#64748b;text-transform:uppercase;border-bottom:2px solid #e2e8f0}
.compare-table td{padding:.5rem .6rem;border-bottom:1px solid #f1f5f9}
.compare-table tr:hover{background:#f8fafc}
</style>
</head><body>
<div class="hdr"><div><a href="/">home</a><a href="/dmpnn?v=v1">v1</a><a href="/dmpnn?v=v2_full">v2</a>
<a href="/dmpnn?v=compare" style="color:#fbbf24;font-weight:600">⇄ Compare</a></div>
<h1>Сравнение: v1 (один SMILES) vs v2 (sequential multi-ligand)</h1>
<div class="sub">{{ total }} комплексов с обоими методами. Сортировка: |Δ| по убыванию.</div></div>
<div class="wrap">
<!-- Stats -->
<div class="stats">
<div class="stat"><div class="n">{{ total }}</div><div class="l">Общих комплексов</div></div>
<div class="stat"><div class="n" style="color:#16a34a">{{ v2_better }}</div><div class="l">v2 лучше</div></div>
<div class="stat"><div class="n" style="color:#dc2626">{{ v1_better }}</div><div class="l">v1 лучше</div></div>
<div class="stat"><div class="n" style="color:#94a3b8">{{ equal }}</div><div class="l">Равны</div></div>
<div class="stat"><div class="n">{{ "%.3f"|format(v1_mean) }}</div><div class="l">v1 средний</div></div>
<div class="stat"><div class="n">{{ "%.3f"|format(v2_mean) }}</div><div class="l">v2 средний</div></div>
</div>
<!-- Table -->
<table class="compare-table">
<thead><tr>
<th>ID</th><th>Abbr</th><th>Metal</th><th>GT</th>
<th>v1 Pred</th><th>v1</th>
<th>v2 Pred</th><th>v2</th>
<th>Δ</th><th>Лучше</th><th></th>
</tr></thead>
<tbody>
{% for cid,metal,ox,v1_pred,v1_score,v2_pred,v2_score,gt,abbr in rows %}
<tr>
<td><a href="/complexes/{{ cid }}">#{{ cid }}</a></td>
<td style="color:#8b5cf6;font-weight:600">{{ abbr or '' }}</td>
<td>{{ metal }}({{ ox }})</td>
<td><code>{{ gt }}</code></td>
<td><code>{{ v1_pred }}</code></td>
<td><a href="/dmpnn/compare/{{ cid }}" style="font-weight:600">{{ "%.2f"|format(v1_score) }}</a></td>
<td><code>{{ v2_pred }}</code></td>
<td><a href="/dmpnn/compare/{{ cid }}" style="font-weight:600">{{ "%.2f"|format(v2_score) }}</a></td>
<td class="{{ 'diff-pos' if v2_score>v1_score else ('diff-neg' if v2_score<v1_score else 'diff-zero') }}">
{{ "%+.2f"|format(v2_score - v1_score) }}</td>
<td style="font-weight:600">{% if v2_score > v1_score %}<span class="diff-pos">v2</span>
{% elif v2_score < v1_score %}<span class="diff-neg">v1</span>
{% else %}<span class="diff-zero">=</span>{% endif %}</td>
<td><a href="/dmpnn/compare/{{ cid }}" style="font-size:.75rem">сравнить</a> <a href="/viewer/?id={{ cid }}" style="font-size:.75rem">3D</a></td>
</tr>{% endfor %}
</tbody></table>
<!-- Pagination -->
<div class="pager">
{% if page > 0 %}<a href=\"?version=compare&per_page={{ per_page }}&page=0\">«</a>
<a href=\"?version=compare&per_page={{ per_page }}&page={{ page-1 }}\">‹</a>{% endif %}
<span class=\"dis\">{{ page+1 }}/{{ total_pages }} ({{ total }})</span>
{% if page < total_pages-1 %}<a href=\"?version=compare&per_page={{ per_page }}&page={{ page+1 }}\">›</a>
<a href=\"?version=compare&per_page={{ per_page }}&page={{ total_pages-1 }}\">»</a>{% endif %}
</div>
</div></body></html>"""


def render_compare(all_rows, per_page, page):
    """Render comparison table."""
    total = len(all_rows)
    total_pages = max(1, (total + per_page - 1) // per_page)
    page = min(page, total_pages - 1)
    rows = all_rows[page * per_page : (page + 1) * per_page]

    v1_scores = [r[4] for r in all_rows]
    v2_scores = [r[6] for r in all_rows]
    v1_mean = sum(v1_scores) / len(v1_scores) if v1_scores else 0
    v2_mean = sum(v2_scores) / len(v2_scores) if v2_scores else 0
    v2_better = sum(1 for a, b in zip(v2_scores, v1_scores) if a > b)
    v1_better = sum(1 for a, b in zip(v2_scores, v1_scores) if a < b)
    equal = sum(1 for a, b in zip(v2_scores, v1_scores) if a == b)

    versions = [('v1', 'v1'), ('v2_full', 'v2'), ('compare', '⇄ Compare')]

    return render_template_string(COMPARE_PAGE,
        rows=rows, total=total, page=page, per_page=per_page,
        total_pages=total_pages, versions=versions,
        v1_mean=v1_mean, v2_mean=v2_mean,
        v2_better=v2_better, v1_better=v1_better, equal=equal)


COMPARE_DETAIL = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Compare #{{ cid }} — BiometalDB</title>
""" + DMPNN_CSS + """
<style>
.side-by-side{display:grid;grid-template-columns:1fr 1fr;gap:1.5rem}
.ver-header{padding:.6rem 1rem;border-radius:6px 6px 0 0;font-weight:600;font-size:.9rem}
.ver-sl{background:#fef3c7;color:#92400e}.ver-seq{background:#dbeafe;color:#1e40af}
.ver-body{background:#fff;border-radius:0 0 8px 8px;padding:1rem;box-shadow:0 1px 3px rgba(0,0,0,.08)}
.cand-row{display:flex;align-items:center;gap:.5rem;padding:.3rem 0;border-bottom:1px solid #f1f5f9}
.cand-row:last-child{border:none}
.cand-row .rank{font-size:.75rem;color:#64748b;width:20px;text-align:center}
.cand-row .smi{font-family:monospace;font-size:.7rem;color:#475569;flex:1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.cand-row .score{font-weight:600;font-size:.8rem;min-width:50px;text-align:right}
.cand-row.top1{background:#f0fdf4;border-radius:4px}
</style></head><body>
<div class="hdr"><div><a href="/dmpnn?version=compare">← сравнение</a>
<a href="/complexes/{{ cid }}">complex #{{ cid }}</a>
<a href="/viewer/?id={{ cid }}">3D Viewer</a></div>
<h1>Complex #{{ cid }} — {{ metal }}({{ ox }})</h1>
<div class="sub">GT: {{ gt }} · Δ = {{ "%+.2f"|format(seq_score - sl_score) }}</div></div>
<div class="wrap">
<div style="display:flex;gap:2rem;margin-bottom:1rem;flex-wrap:wrap">
<div class="stat"><div class="n">{{ gt }}</div><div class="l">Ground Truth</div></div>
<div class="stat"><div class="n" style="color:#92400e">{{ "%.3f"|format(sl_score) }}</div><div class="l">Single-lig</div></div>
<div class="stat"><div class="n" style="color:#1e40af">{{ "%.3f"|format(seq_score) }}</div><div class="l">Sequential</div></div>
<div class="stat"><div class="n {{ 'diff-pos' if seq_score>sl_score else ('diff-neg' if seq_score<sl_score else 'diff-zero') }}">{{ "%+.3f"|format(seq_score - sl_score) }}</div><div class="l">Δ Sequential</div></div>
</div>
<div class="side-by-side">
<!-- Single-ligand -->
<div>
<div class="ver-header ver-sl">Single-ligand (v1) — {{ sl_cands|length }} candidates</div>
<div class="ver-body">
<div style="margin-bottom:.5rem;font-size:.85rem">Pred: <code>{{ sl_pred }}</code></div>
{% for rank,smi,score in sl_cands %}
<div class="cand-row {{ 'top1' if rank==0 else '' }}">
<span class="rank">#{{ rank+1 }}</span>
<span class="smi">{{ smi[:80] }}{% if smi|length>80 %}...{% endif %}</span>
<span class="score" style="color:{{ '#16a34a' if score>0.5 else ('#d97706' if score>0.1 else '#dc2626') }}">{{ "%.3f"|format(score) }}</span>
</div>{% endfor %}
</div></div>
<!-- Sequential -->
<div>
<div class="ver-header ver-seq">Sequential multi-lig (v2) — {{ seq_cands|length }} candidates</div>
<div class="ver-body">
<div style="margin-bottom:.5rem;font-size:.85rem">Pred: <code>{{ seq_pred }}</code></div>
{% for rank,smi,score in seq_cands %}
<div class="cand-row {{ 'top1' if rank==0 else '' }}">
<span class="rank">#{{ rank+1 }}</span>
<span class="smi">{{ smi[:80] }}{% if smi|length>80 %}...{% endif %}</span>
<span class="score" style="color:{{ '#16a34a' if score>0.5 else ('#d97706' if score>0.1 else '#dc2626') }}">{{ "%.3f"|format(score) }}</span>
</div>{% endfor %}
</div></div>
</div>
</div></body></html>"""


def render_candidate_smi(smi):
    """Render candidate SMILES to base64 PNG."""
    try:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if not mol:
            return None
        AllChem.Compute2DCoords(mol)
        drawer = Draw.MolDraw2DCairo(380, 280)
        drawer.drawOptions().bondLineWidth = 2
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode()
    except:
        return None


def get_versions(conn):
    """Get available scoring versions from DB."""
    cur = conn.cursor()
    versions = []
    for v, cnt in cur.execute("SELECT scoring_version, COUNT(*) FROM dmpnn_summary GROUP BY scoring_version ORDER BY scoring_version"):
        if v == 'single-ligand':
            versions.append(('single-ligand', f'Single-ligand ({cnt})'))
        elif v == 'sequential':
            versions.append(('sequential', f'Sequential multi-lig ({cnt})'))
        else:
            versions.append((v, f'{v} ({cnt})'))
    versions.append(('compare', '⇄ Compare'))
    return versions


def get_compare_rows(conn):
    """Get complexes scored by both v1 and v2_full for comparison."""
    cur = conn.cursor()
    cur.execute("""
        SELECT v1.complex_id, v1.metal, v1.oxidation_state,
            v1.top1_donor_pred, v1.match_score,
            v2.top1_donor_pred, v2.match_score,
            v1.donor_atoms_gt,
            GROUP_CONCAT(DISTINCT m.abbreviation)
        FROM dmpnn_summary v1
        JOIN dmpnn_summary v2 ON v1.complex_id = v2.complex_id
        LEFT JOIN measurements m ON m.complex_id = v1.complex_id
        WHERE v1.scoring_version='single-ligand' AND v2.scoring_version='sequential'
        GROUP BY v1.complex_id
        ORDER BY ABS(v1.match_score - v2.match_score) DESC
    """)
    return cur.fetchall()


def register_dmpnn_routes(app):
    """Register enhanced D-MPNN routes on the Flask app."""

    @app.route("/dmpnn")
    def dmpnn_summary():
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()

        # Params
        version = request.args.get('version', 'sequential')
        per_page = int(request.args.get('per_page', 50))
        page = int(request.args.get('page', 0))
        f_metal = request.args.get('metal', '')
        f_ox = request.args.get('ox', '')
        f_donors = request.args.get('donors', '')
        f_score_min = request.args.get('score_min', '')
        f_score_max = request.args.get('score_max', '')
        f_q = request.args.get('q', '').strip()
        sort = request.args.get('sort', 'score_desc')

        # Available versions
        versions = get_versions(conn)

        # Compare view
        if version == 'compare':
            compare_rows = get_compare_rows(conn)
            conn.close()
            return render_compare(compare_rows, per_page, page)

        # Build query
        where = ["s.scoring_version = ?"]
        params = [version]

        if f_metal:
            where.append("s.metal = ?")
            params.append(f_metal)
        if f_ox:
            where.append("s.oxidation_state = ?")
            params.append(int(f_ox))
        if f_donors:
            for d in f_donors.split(','):
                d = d.strip()
                if d:
                    where.append("s.donor_atoms_gt LIKE ?")
                    params.append(f'%"{d}"%')
        if f_score_min:
            where.append("s.match_score >= ?")
            params.append(float(f_score_min))
        if f_score_max:
            where.append("s.match_score <= ?")
            params.append(float(f_score_max))
        if f_q:
            if f_q.startswith('#'):
                try:
                    where.append("s.complex_id = ?")
                    params.append(int(f_q[1:]))
                except ValueError:
                    pass
            else:
                where.append("""(CAST(s.complex_id AS TEXT) LIKE ? OR s.top1_donor_pred LIKE ?
                    OR s.complex_id IN (SELECT m.complex_id FROM measurements m WHERE m.abbreviation LIKE ?))""")
                params.extend([f'%{f_q}%', f'%{f_q}%', f'%{f_q}%'])

        where_sql = "WHERE " + " AND ".join(where)

        # Sort
        order_map = {
            'score_desc': 's.match_score DESC, s.complex_id',
            'score_asc': 's.match_score ASC, s.complex_id',
            'id_asc': 's.complex_id ASC',
            'id_desc': 's.complex_id DESC',
        }
        order_sql = order_map.get(sort, 's.match_score DESC, s.complex_id')

        # Build full query with JOIN for abbreviations
        join_sql = "LEFT JOIN measurements m ON m.complex_id = s.complex_id"

        # Count
        cur.execute(f"SELECT COUNT(DISTINCT s.complex_id) FROM dmpnn_summary s {join_sql} {where_sql}", params)
        total = cur.fetchone()[0]

        # Page
        total_pages = max(1, (total + per_page - 1) // per_page)
        page = min(page, total_pages - 1)
        offset = page * per_page

        # Fetch with abbreviations
        cur.execute(f"""SELECT s.complex_id, s.metal, s.oxidation_state, s.donor_atoms_gt,
            s.top1_donor_pred, s.match_score, s.n_candidates,
            GROUP_CONCAT(DISTINCT m.abbreviation) as abbrevs
            FROM dmpnn_summary s {join_sql}
            {where_sql}
            GROUP BY s.complex_id
            ORDER BY {order_sql}
            LIMIT ? OFFSET ?""",
            params + [per_page, offset])
        rows = cur.fetchall()

        # Stats — use subquery to avoid JOIN duplication
        base_where = [c.replace('s.', '') for c in where]
        base_where[0] = "scoring_version = ?"  # first condition
        base_where_sql = "WHERE " + " AND ".join(base_where)
        cur.execute(f"""SELECT COUNT(*), AVG(match_score),
            SUM(CASE WHEN match_score >= 0.99 THEN 1 ELSE 0 END),
            SUM(CASE WHEN match_score > 0 AND match_score < 0.99 THEN 1 ELSE 0 END),
            SUM(CASE WHEN match_score = 0 THEN 1 ELSE 0 END)
            FROM dmpnn_summary {base_where_sql}""", params)
        s = cur.fetchone()
        total_all, avg, n_match, n_partial, n_fail = s[0] or 0, s[1] or 0, s[2] or 0, s[3] or 0, s[4] or 0

        # Filter options
        cur.execute("SELECT DISTINCT metal FROM dmpnn_summary WHERE scoring_version=? ORDER BY metal", [version])
        metals = [r[0] for r in cur.fetchall()]
        cur.execute("SELECT DISTINCT oxidation_state FROM dmpnn_summary WHERE scoring_version=? ORDER BY oxidation_state", [version])
        ox_states = [r[0] for r in cur.fetchall()]

        conn.close()

        # Ver label
        ver_label = dict(versions).get(version, version)

        def pager_url(p):
            parts = [f'version={version}', f'per_page={per_page}', f'page={p}', f'sort={sort}']
            if f_metal: parts.append(f'metal={f_metal}')
            if f_ox: parts.append(f'ox={f_ox}')
            if f_donors: parts.append(f'donors={f_donors}')
            if f_score_min: parts.append(f'score_min={f_score_min}')
            if f_score_max: parts.append(f'score_max={f_score_max}')
            if f_q: parts.append(f'q={f_q}')
            return '?' + '&'.join(parts)

        # Page range (show max 10 page links)
        start_p = max(0, page - 5)
        end_p = min(total_pages, start_p + 10)
        page_range = range(start_p, end_p)

        return render_template_string(DMPNN_LIST,
            version=version, versions=versions, ver_label=ver_label,
            total=total_all, avg=avg, n_match=n_match, n_partial=n_partial, n_fail=n_fail,
            rows=rows, metals=metals, ox_states=ox_states,
            f_metal=f_metal, f_ox=f_ox, f_donors=f_donors,
            f_score_min=f_score_min, f_score_max=f_score_max, f_q=f_q,
            sort=sort, per_page=per_page, page=page,
            total_pages=total_pages, page_range=page_range,
            pager_url=pager_url, match_color=match_color)

    @app.route("/dmpnn/<int:cid>")
    def dmpnn_detail(cid):
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()

        version = request.args.get('version', 'sequential')
        versions = get_versions(conn)

        # Summary
        cur.execute("""SELECT s.metal, s.oxidation_state, s.donor_atoms_gt, s.top1_donor_pred,
            s.match_score, s.n_candidates,
            GROUP_CONCAT(DISTINCT m.abbreviation)
            FROM dmpnn_summary s
            LEFT JOIN measurements m ON m.complex_id = s.complex_id
            WHERE s.complex_id=? AND s.scoring_version=?
            GROUP BY s.complex_id""", (cid, version))
        summ = cur.fetchone()

        # Fallback: if not found in requested version, try single-ligand
        if not summ and version != 'single-ligand':
            cur.execute("""SELECT metal, oxidation_state, donor_atoms_gt, top1_donor_pred,
                match_score, n_candidates, GROUP_CONCAT(DISTINCT m.abbreviation)
                FROM dmpnn_summary s LEFT JOIN measurements m ON m.complex_id=s.complex_id
                WHERE s.complex_id=? AND s.scoring_version='single-ligand'
                GROUP BY s.complex_id""", (cid,))
            summ = cur.fetchone()
            if summ:
                version = 'single-ligand'

        if not summ:
            conn.close()
            abort(404)

        metal, ox, gt, pred, match, nc, abbrevs = summ

        # Candidates — filter by scoring_version
        cur.execute("""SELECT rank, candidate_smi, is_top1, score, COALESCE(score_std, 0)
            FROM dmpnn_results WHERE complex_id=? AND scoring_version=?
            ORDER BY rank""", (cid, version))
        candidates = cur.fetchall()

        # Fallback: if no per-candidate data for this version, show what exists
        has_cands = len(candidates) > 0
        if not has_cands:
            cur.execute("""SELECT rank, candidate_smi, is_top1, score, COALESCE(score_std, 0)
                FROM dmpnn_results WHERE complex_id=? ORDER BY rank LIMIT 30""", (cid,))
            candidates = cur.fetchall()

        # Render structures
        cand_imgs = {}
        for rank, smi, is_top1, score, score_std in candidates:
            cand_imgs[rank] = render_candidate_smi(smi)

        conn.close()

        return render_template_string(DMPNN_DETAIL,
            cid=cid, metal=metal, ox=ox, gt=gt, pred=pred,
            match=match, nc=nc, abbrevs=abbrevs, version=version, versions=versions,
            candidates=candidates, cand_imgs=cand_imgs, has_cands=has_cands)

    @app.route("/dmpnn/compare/<int:cid>")
    def dmpnn_compare(cid):
        """Side-by-side comparison of single-ligand vs sequential for one complex."""
        conn = sqlite3.connect(DB_PATH)
        cur = conn.cursor()

        # Get both versions
        cur.execute("""SELECT scoring_version, metal, oxidation_state, donor_atoms_gt,
            top1_donor_pred, match_score
            FROM dmpnn_summary WHERE complex_id=? AND scoring_version IN ('single-ligand','sequential')""", (cid,))
        versions_data = {r[0]: r for r in cur.fetchall()}

        if not versions_data:
            conn.close()
            abort(404)

        # Use whichever version has data for meta
        meta = versions_data.get('sequential') or versions_data.get('single-ligand')
        metal, ox, gt = meta[1], meta[2], meta[3]

        sl = versions_data.get('single-ligand')
        seq = versions_data.get('sequential')
        sl_pred = sl[4] if sl else 'N/A'
        sl_score = sl[5] if sl else 0.0
        seq_pred = seq[4] if seq else 'N/A'
        seq_score = seq[5] if seq else 0.0

        # Get candidates — need to check which version's results are in dmpnn_results
        # For now, load all candidates (dmpnn_results doesn't have scoring_version yet)
        cur.execute("""SELECT rank, candidate_smi, score
            FROM dmpnn_results WHERE complex_id=? ORDER BY rank""", (cid,))
        all_cands = cur.fetchall()

        # The dmpnn_results currently has one version's candidates.
        # Show them as the sequential version since that's the main one.
        seq_cands = [(r[0], r[1], r[2]) for r in all_cands]
        # For single-ligand, we don't have per-candidate data stored separately.
        # Show summary prediction instead.
        sl_cands = seq_cands  # Same candidates, different ranking would need separate storage

        conn.close()

        return render_template_string(COMPARE_DETAIL,
            cid=cid, metal=metal, ox=ox, gt=gt,
            sl_pred=sl_pred, sl_score=sl_score, sl_cands=sl_cands,
            seq_pred=seq_pred, seq_score=seq_score, seq_cands=seq_cands)
