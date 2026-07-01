#!/usr/bin/env python3
"""Conformer Browser routes for BiometalDB Flask server.
Shows GFN2-xTB-optimized (correct per-complex charge/spin), RMSD+energy-deduped distinct conformers for 427 Ir(III) complexes.
"""
import os, json
from flask import request, render_template_string, send_file, abort, Response, redirect

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ZENITH_JSON = os.path.join(BASE_DIR, "data", "zenith_conformers.json")
CONFORMER_DIR = "/root/conformer-generation/results/conformers"
CROSSVAL_REPORT = os.path.join(BASE_DIR, "conformers", "crest_crossval_report.md")

# Nav HTML (matches complexes_routes.py)
NAV_HTML = (
    '<a href="/">home</a>'
    '<a href="/complexes">complexes</a>'
    '<a href="/ir_cn_families">Ir C^N families</a>'
    '<a href="/conformers">conformers</a>'
    '<a href="/conformers/validation">CREST validation</a>'
    '<a href="/complexes/cell-death-heatmap">cell death map</a>'
    '<a href="/dmpnn">D-MPNN</a>'
)

ZENITH_PAGE = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Ir(III) Conformer Browser — BiometalDB</title>
<link rel="icon" href="/static/favicon.ico" type="image/x-icon">
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
:root {
  --bg: #f0f2f5;
  --card: #fff;
  --text: #1e293b;
  --muted: #64748b;
  --accent: #3b82f6;
  --accent2: #1e40af;
  --dark: #0f172a;
  --good: #dcfce7;
  --good-text: #166534;
  --warn: #fef3c7;
  --warn-text: #92400e;
}
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:var(--bg);color:var(--text)}
.hdr{background:linear-gradient(135deg,var(--dark),#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid var(--accent)}
.hdr h1{margin:0;font-size:1.4rem}.hdr a{color:#60a5fa;text-decoration:none;margin-right:1rem;font-size:.85rem}
.hdr .sub{margin-top:.3rem;font-size:.85rem;color:#93c5fd}
.wrap{max-width:1300px;margin:1.5rem auto;padding:0 1.5rem}
.stats{display:flex;gap:1rem;margin-bottom:1.5rem;flex-wrap:wrap}
.stat{background:var(--card);border-radius:8px;padding:1rem 1.5rem;box-shadow:0 1px 3px rgba(0,0,0,.08);flex:1;text-align:center;min-width:120px}
.stat .n{font-size:1.8rem;font-weight:700}.stat .l{font-size:.8rem;color:var(--muted);margin-top:.2rem}
.toolbar{background:var(--card);border-radius:8px;padding:.8rem 1.2rem;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:1.2rem;display:flex;gap:1rem;flex-wrap:wrap;align-items:center}
.toolbar label{font-size:.78rem;color:var(--muted);font-weight:600;display:flex;align-items:center;gap:.3rem}
.toolbar select,.toolbar input{padding:.3rem .5rem;border:1px solid #e2e8f0;border-radius:5px;font-size:.82rem}
.toolbar select:focus,.toolbar input:focus{border-color:var(--accent);outline:none}
.btn-sm{padding:.3rem .7rem;border-radius:5px;font-size:.78rem;cursor:pointer;border:1px solid #e2e8f0;background:var(--card);color:var(--accent)}
.btn-sm:hover{background:#eff6ff}
.btn-sm.active{background:var(--accent);color:#fff;border-color:var(--accent)}
.count{margin-left:auto;font-size:.82rem;color:var(--muted)}

/* Complex card */
.card{background:var(--card);border-radius:8px;padding:0;box-shadow:0 1px 3px rgba(0,0,0,.08);margin-bottom:.6rem;overflow:hidden}
.card-hdr{padding:.7rem 1rem;cursor:pointer;display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:.5rem;border-left:4px solid var(--accent);transition:background .15s}
.card-hdr:hover{background:#f8fafc}
.card-hdr .cid{font-weight:700;color:var(--accent2);font-size:.95rem}
.card-hdr .cid a{color:var(--accent2);text-decoration:none}
.card-hdr .meta{font-size:.75rem;color:var(--muted)}
.card-hdr .nconf{font-size:.82rem;font-weight:600;color:var(--accent);white-space:nowrap}
.card-hdr .erange{font-size:.72rem;color:var(--muted)}
.card-body{display:none;padding:.5rem 1rem 1rem;border-top:1px solid #e2e8f0}
.card-body.open{display:block}
.card-body .summary{font-size:.78rem;color:var(--muted);margin-bottom:.8rem}

/* Energy table */
.etbl{width:100%;border-collapse:collapse;font-size:.78rem;margin-bottom:.8rem}
.etbl th{background:#f1f5f9;padding:.35rem .5rem;text-align:left;font-weight:600;color:var(--muted);font-size:.7rem;text-transform:uppercase;position:sticky;top:0}
.etbl td{padding:.3rem .5rem;border-bottom:1px solid #f1f5f9;font-family:monospace;font-size:.75rem}
.etbl tr:hover td{background:#f8fafc}
.etbl tr.sel td{background:#eff6ff;font-weight:600}
.etbl .de{color:#94a3b8}
.etbl .best{color:var(--good-text)}

/* 3D viewer */
.viewer-wrap{display:flex;gap:1rem;flex-wrap:wrap}
.viewer-panel{flex:1;min-width:350px}
.cx3d{width:100%;height:360px;position:relative;background:var(--dark);border-radius:6px}
.table-panel{flex:1;min-width:300px;max-height:400px;overflow-y:auto}

/* SDF source label */
.sdf-src{font-size:.65rem;color:var(--muted);margin-top:.3rem}

/* Jump bar */
.jump{display:flex;gap:.4rem;flex-wrap:wrap;margin-bottom:1rem;max-height:8rem;overflow-y:auto}
.jump a{background:var(--card);border:1px solid #e2e8f0;border-radius:4px;padding:.2rem .5rem;font-size:.7rem;text-decoration:none;color:var(--accent2);font-family:monospace}
.jump a:hover{background:#eff6ff}

/* Pagination */
.pager{display:flex;gap:.3rem;justify-content:center;margin:1.5rem 0;flex-wrap:wrap}
.pager a,.pager span{padding:.3rem .7rem;border-radius:5px;font-size:.8rem;text-decoration:none}
.pager a{background:var(--card);color:var(--accent);box-shadow:0 1px 2px rgba(0,0,0,.08)}
.pager a:hover{background:#eff6ff}
.pager .cur{background:var(--accent);color:#fff;font-weight:600}

/* Conformer selector buttons */
.conf-sel{display:flex;gap:.2rem;flex-wrap:wrap;margin-bottom:.5rem}
.conf-sel button{padding:.2rem .45rem;border-radius:3px;font-size:.68rem;cursor:pointer;border:1px solid #e2e8f0;background:var(--card);color:var(--text);font-family:monospace}
.conf-sel button:hover{background:#eff6ff}
.conf-sel button.active{background:var(--accent);color:#fff;border-color:var(--accent)}
</style>
</head><body>
<div class="hdr"><div>""" + NAV_HTML + """</div>
<h1>Ir(III) Conformer Browser</h1>
<div class="sub">{{ total_compounds }} Ir(III) complexes · {{ total_conformers }} distinct conformers · Uniconf torsional seeds (MK≤3000) → GFN2-xTB geometry optimization at the correct per-complex charge/spin → RMSD (0.5&#8202;Å) + energy (0.1&#8202;kcal/mol) dedup · <a href="/conformers/validation" style="color:#93c5fd">CREST cross-validated ✓</a></div></div>
<div class="wrap">
<div class="stats">
<div class="stat"><div class="n">{{ total_compounds }}</div><div class="l">Complexes</div></div>
<div class="stat"><div class="n">{{ total_conformers }}</div><div class="l">Total Conformers</div></div>
<div class="stat"><div class="n">{{ '{:.1f}'.format(avg_conf) }}</div><div class="l">Avg / Complex</div></div>
<div class="stat"><div class="n">{{ with_multiple }}</div><div class="l">With ≥3 Conformers</div></div>
<div class="stat"><div class="n">{{ rigid_single }}</div><div class="l">Rigid (1 conformer)</div></div>
</div>

<div class="toolbar">
<label>Search ID: <input type="text" id="search_id" placeholder="#4836" style="width:80px" onkeyup="applyFilters()"></label>
<label>Min conformers: <input type="number" id="min_conf" value="0" min="0" max="60" style="width:60px" onchange="applyFilters()"></label>
<label>Sort:
<select id="sort_sel" onchange="applyFilters()">
<option value="id_asc">ID ↑</option>
<option value="id_desc">ID ↓</option>
<option value="nconf_desc">Conformers ↓</option>
<option value="nconf_asc">Conformers ↑</option>
<option value="mw_asc">Molar mass ↑</option>
<option value="mw_desc">Molar mass ↓</option>
<option value="erange_desc">Energy range ↓</option>
</select>
</label>
<span class="count" id="result_count"></span>
</div>

<div class="jump" id="jump_bar"></div>
<div id="cards_container"></div>
<div class="pager" id="pager"></div>
<p style="font-size:.75rem;color:var(--muted);margin-top:2rem">
Conformer generation for Ir(III) C,N-chelate complexes.
Method: Uniconf v1.0.1 (UFF/NLOPT) seeds. 
Then GFN2-xTB geometry optimization at the correct per-complex charge/spin, with RMSD+energy deduplication to distinct minima.
SDF sorted by GFN2-xTB energy (lowest first). 
Regenerate with <code>scripts/build_zenith_conformers.py</code>.</p>
</div>

<script>
const RAW_DATA = {{ data_json|safe }};
const RECORDS = RAW_DATA.records;
const PER_PAGE = 50;

// Viewer state
const viewers = {};  // compound_id -> { viewer, models }
let currentPage = 0;
let activeCompound = null;

function applyFilters() {
  const searchId = document.getElementById('search_id').value.trim();
  const minConf = parseInt(document.getElementById('min_conf').value) || 0;
  const sortVal = document.getElementById('sort_sel').value;

  let filtered = RECORDS.filter(r => {
    if (searchId && !String(r.id).includes(searchId)) return false;
    if (r.conformer_count < minConf) return false;
    return true;
  });

  const sorters = {
    'id_asc': (a,b) => a.id - b.id,
    'id_desc': (a,b) => b.id - a.id,
    'nconf_desc': (a,b) => b.conformer_count - a.conformer_count,
    'nconf_asc': (a,b) => a.conformer_count - b.conformer_count,
    'mw_asc': (a,b) => (a.mw||0) - (b.mw||0),
    'mw_desc': (a,b) => (b.mw||0) - (a.mw||0),
    'erange_desc': (a,b) => (b.energy_range_kcal||0) - (a.energy_range_kcal||0),
  };
  filtered.sort(sorters[sortVal] || sorters['id_asc']);

  document.getElementById('result_count').textContent = filtered.length + ' complexes';
  renderCards(filtered);
  renderJump(filtered);
}

function renderJump(filtered) {
  // Quick-jump bar: IDs grouped by 100s
  const groups = {};
  for (const r of filtered) {
    const g = Math.floor(r.id / 100) * 100;
    if (!groups[g]) groups[g] = r.id;
  }
  let html = '';
  for (const g of Object.keys(groups).sort((a,b) => a-b)) {
    html += `<a href="#card-${groups[g]}">${g}-${parseInt(g)+99}</a>`;
  }
  document.getElementById('jump_bar').innerHTML = html;
}

function renderCards(filtered) {
  const totalPages = Math.ceil(filtered.length / PER_PAGE);
  if (currentPage >= totalPages) currentPage = totalPages - 1;
  if (currentPage < 0) currentPage = 0;

  const start = currentPage * PER_PAGE;
  const page = filtered.slice(start, start + PER_PAGE);

  let html = '';
  for (const r of page) {
    const nconf = r.conformer_count;
    const erange = r.energy_range_kcal ? r.energy_range_kcal.toFixed(2) : '—';
    const emin = r.energy_min_kcal ? r.energy_min_kcal.toFixed(1) : '';
    const mwStr = r.mw ? `MW: ${r.mw.toFixed(1)}` : '';
    const ic50 = r.min_ic50_dark ? `IC50: ${r.min_ic50_dark.toFixed(1)} µM` : '';
    const metal = r.metal || 'Ir';
    const charge = r.charge != null ? (r.charge > 0 ? '+' + r.charge : r.charge) : '';

    html += `<div class="card" id="card-${r.id}">`;
    html += `<div class="card-hdr" onclick="toggleCard(${r.id})">`;
    html += `<span><span class="cid">#${r.id}</span> <span class="meta">${metal}(${r.oxidation_state||'?'})${charge}</span></span>`;
    html += `<span style="display:flex;gap:1rem;align-items:center">`;
    html += `<span class="erange">ΔE: ${erange} kcal/mol</span>`;
    html += `<span class="nconf">${nconf} conformer${nconf!==1?'s':''}</span>`;
    html += `</span></div>`;
    html += `<div class="card-body" id="body-${r.id}">`;
    html += `<div class="summary">`;
    if (mwStr) html += `<span>${mwStr}</span> · `;
    if (ic50) html += `<span>${ic50}</span> · `;
    html += `<span>${r.n_rotatable_bonds} rotatable bonds</span>`;
    if (r.rmsd_max > 0) html += ` · <span>RMSD max: ${r.rmsd_max.toFixed(2)} Å</span>`;
    html += `</div>`;
    html += `<div class="viewer-wrap">`;
    html += `<div class="viewer-panel"><div class="cx3d" id="viewer-${r.id}"><div style="color:#64748b;position:absolute;top:50%;left:50%;transform:translate(-50%,-50%)">Click to load 3D</div></div>`;
    html += `<div class="sdf-src">SDF: conformers_opt_sorted.sdf (GFN2-xTB optimized)</div></div>`;
    html += `<div class="table-panel">`;
    html += `<table class="etbl"><thead><tr><th>#</th><th>E (Hartree)</th><th>E (kcal/mol)</th><th>ΔE (kcal/mol)</th></tr></thead><tbody>`;
    if (r.energies_kcal && r.energies_kcal.length > 0) {
      const eminVal = r.energies_kcal[0]; // sorted ascending
      for (let i = 0; i < r.energies_kcal.length; i++) {
        const de = r.energies_kcal[i] - eminVal;
        const cls = i === 0 ? 'best' : (de < 1.0 ? '' : 'de');
        html += `<tr class="${cls}" onclick="showConformer(${r.id}, ${i})" style="cursor:pointer" title="Click to view conformer ${i+1}">`;
        html += `<td>${i+1}</td>`;
        html += `<td>${r.energies_hartree ? r.energies_hartree[i].toFixed(6) : '—'}</td>`;
        html += `<td>${r.energies_kcal[i].toFixed(2)}</td>`;
        html += `<td>${i===0 ? '0.00 (best)' : de.toFixed(2)}</td>`;
        html += `</tr>`;
      }
    } else {
      html += `<tr><td colspan="4">No energy data</td></tr>`;
    }
    html += `</tbody></table></div></div>`;
    html += `</div></div>`;
  }
  document.getElementById('cards_container').innerHTML = html;

  // Pagination
  let pagerHtml = '';
  if (totalPages > 1) {
    pagerHtml += `<span>Page ${currentPage+1} of ${totalPages}</span>`;
    for (let i = 0; i < totalPages; i++) {
      if (i === currentPage) {
        pagerHtml += `<span class="cur">${i+1}</span>`;
      } else {
        pagerHtml += `<a href="#" onclick="goPage(${i});return false">${i+1}</a>`;
      }
    }
  }
  document.getElementById('pager').innerHTML = pagerHtml;
}

function goPage(n) {
  currentPage = n;
  applyFilters();
  window.scrollTo(0, 0);
}

let loadedCompounds = {};

function toggleCard(cid) {
  const body = document.getElementById('body-' + cid);
  const isOpen = body.classList.contains('open');
  if (isOpen) {
    body.classList.remove('open');
  } else {
    body.classList.add('open');
    if (!loadedCompounds[cid]) {
      loadViewer(cid);
    }
  }
}

function loadViewer(cid) {
  const containerId = 'viewer-' + cid;
  const container = document.getElementById(containerId);
  if (!container) return;
  container.innerHTML = '<div style="color:#64748b;position:absolute;top:50%;left:50%;transform:translate(-50%,-50%)">Loading...</div>';

  const sdfUrl = '/conformers/sdf/' + cid;
  fetch(sdfUrl)
    .then(r => {
      if (!r.ok) throw new Error('SDF not found');
      return r.text();
    })
    .then(sdfText => {
      // Count conformers
      const parts = sdfText.split('$$$$\n').filter(p => p.trim());
      const nModels = parts.length;

      let element = document.getElementById(containerId);
      // Container might have been replaced
      if (!element) { element = document.getElementById(containerId); }
      if (!element) return;
      element.innerHTML = '';

      let viewer = $3Dmol.createViewer(element, {
        backgroundColor: '0x0f172a',
        antialias: true,
      });

      // Add all models
      for (let i = 0; i < nModels; i++) {
        let chunk = parts[i];
        if (!chunk.includes('$$$$')) chunk += '\n$$$$';
        viewer.addModel(chunk, 'sdf');
      }

      viewer.setStyle({}, {stick: {radius: 0.15, colorscheme: 'Jmol'}});
      viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.2, colorscheme: 'Jmol'});
      viewer.zoomTo();
      viewer.render();

      // Show first model, hide rest
      for (let i = 1; i < nModels; i++) {
        viewer.setStyle({model: i}, {stick: {hidden: true}, sphere: {hidden: true}});
      }
      viewer.render();

      loadedCompounds[cid] = { viewer, nModels, currentModel: 0 };
    })
    .catch(err => {
      const element = document.getElementById(containerId);
      if (element) {
        element.innerHTML = '<div style="color:#ef4444;position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);font-size:.85rem">3D structure not available</div>';
      }
    });
}

function showConformer(cid, modelIdx) {
  if (!loadedCompounds[cid]) {
    // Open the card first
    const body = document.getElementById('body-' + cid);
    if (body && !body.classList.contains('open')) {
      body.classList.add('open');
    }
    loadViewer(cid);
    // Retry after load
    setTimeout(() => showConformer(cid, modelIdx), 500);
    return;
  }

  const state = loadedCompounds[cid];
  const v = state.viewer;

  // Show selected model, hide others
  for (let i = 0; i < state.nModels; i++) {
    if (i === modelIdx) {
      v.setStyle({model: i}, {stick: {radius: 0.15, colorscheme: 'Jmol'}});
      v.addSurface($3Dmol.SurfaceType.VDW, {opacity: 0.2, colorscheme: 'Jmol'}, {model: i});
    } else {
      v.setStyle({model: i}, {stick: {hidden: true}, sphere: {hidden: true}});
    }
  }
  v.zoomTo();
  v.render();
  state.currentModel = modelIdx;

  // Highlight table row
  const table = document.getElementById('body-' + cid).querySelector('table');
  if (table) {
    table.querySelectorAll('tr').forEach((tr, idx) => {
      tr.classList.toggle('sel', idx === modelIdx + 1); // +1 for header
    });
  }

  // Scroll viewer into view
  const viewerEl = document.getElementById('viewer-' + cid);
  if (viewerEl) viewerEl.scrollIntoView({behavior: 'smooth', block: 'center'});
}

// Init
applyFilters();
</script>
</body></html>"""


def load_zenith_data():
    """Load the zenith conformers JSON."""
    if not os.path.exists(ZENITH_JSON):
        return None
    with open(ZENITH_JSON) as f:
        return json.load(f)


def register_conformer_routes(app):
    """Register conformer browser routes on the Flask app."""

    @app.route("/conformers")
    def conformer_browser():
        data = load_zenith_data()
        if not data:
            return Response("Conformer data unavailable. Run scripts/build_zenith_conformers.py first.", status=503)

        total = data['total_compounds']
        total_conf = data['total_conformers']
        avg_conf = total_conf / total if total > 0 else 0

        return render_template_string(
            ZENITH_PAGE,
            total_compounds=total,
            total_conformers=total_conf,
            with_multiple=data['with_multiple'],
            rigid_single=data['rigid_single'],
            avg_conf=avg_conf,
            data_json=json.dumps(data, ensure_ascii=False),
        )

    @app.route("/conformers/sdf/<int:cid>")
    def conformer_sdf(cid):
        """Serve the GFN2-xTB-optimized multi-conformer SDF (falls back to the original)."""
        base = os.path.join(CONFORMER_DIR, f"complex_{cid}")
        sdf_path = os.path.join(base, "opt", "conformers_opt_sorted.sdf")
        if not os.path.exists(sdf_path):
            sdf_path = os.path.join(base, "conformers_sorted.sdf")
        if not os.path.exists(sdf_path):
            abort(404)
        return send_file(sdf_path, mimetype="chemical/x-mdl-sdfile")

    @app.route("/conformers/validation")
    def conformer_validation():
        """Render the CREST cross-validation + remediation report (markdown -> HTML)."""
        if not os.path.exists(CROSSVAL_REPORT):
            return Response("Validation report unavailable.", status=404)
        try:
            import markdown
            body = markdown.markdown(open(CROSSVAL_REPORT).read(),
                                     extensions=["tables", "fenced_code", "sane_lists"])
        except Exception:
            # graceful fallback: serve raw markdown in a <pre>
            import html as _html
            body = "<pre>" + _html.escape(open(CROSSVAL_REPORT).read()) + "</pre>"
        page = """<!DOCTYPE html><html><head><meta charset="utf-8">
<title>CREST cross-validation — Ir(III) Conformer Browser</title>
<link rel="icon" href="/static/favicon.ico" type="image/x-icon">
<style>
body{font-family:'Inter',-apple-system,sans-serif;margin:0;background:#f0f2f5;color:#1e293b}
.hdr{background:linear-gradient(135deg,#0f172a,#1e3a5f);color:#fff;padding:1rem 2rem;border-bottom:3px solid #3b82f6}
.hdr a{color:#60a5fa;text-decoration:none;margin-right:1rem;font-size:.85rem}
.wrap{max-width:900px;margin:1.5rem auto;padding:2rem 2.5rem;background:#fff;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,.08)}
.wrap h1{font-size:1.5rem;margin-top:0}.wrap h2{font-size:1.15rem;margin-top:1.6rem;color:#1e40af;border-bottom:1px solid #e2e8f0;padding-bottom:.3rem}
.wrap h3{font-size:1rem;color:#334155}
.wrap table{border-collapse:collapse;width:100%;margin:1rem 0;font-size:.85rem}
.wrap th,.wrap td{border:1px solid #e2e8f0;padding:.35rem .6rem;text-align:left}
.wrap th{background:#f1f5f9;font-weight:600}
.wrap code{background:#f1f5f9;padding:.1rem .3rem;border-radius:3px;font-size:.85em}
.wrap pre{background:#0f172a;color:#e2e8f0;padding:1rem;border-radius:6px;overflow-x:auto;font-size:.8rem}
.wrap a{color:#2563eb}.wrap hr{border:none;border-top:1px solid #e2e8f0;margin:2rem 0}
</style></head><body>
<div class="hdr"><div>""" + NAV_HTML + """</div></div>
<div class="wrap">""" + body + """</div></body></html>"""
        return Response(page, mimetype="text/html")

    # Backwards-compatible 301 redirects from the old /zenith route.
    # "Zenith" was the mission orchestrator, not the scientific result.
    @app.route("/zenith")
    def zenith_browser_redirect():
        return redirect("/conformers", code=301)

    @app.route("/zenith/sdf/<int:cid>")
    def zenith_sdf_redirect(cid):
        return redirect(f"/conformers/sdf/{cid}", code=301)
