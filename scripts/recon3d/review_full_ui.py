"""Paginated, data-driven review UI for the FULL BiometalDB 3D reconstruction.

Unlike review_ui.py (one self-contained page per ~100 complexes — inlines every
card and every 3Dmol viewer div, which does not scale past a few hundred), this
generator emits:

  * index.json  — one compact record per complex (assignment summary, ligand
                  table, per-isomer file refs + validation), ~a few MB for 9.4k.
  * review.html — a small shell that fetches index.json, offers metal/status/id
                  filters, paginates (PAGE_SIZE per view), and renders only the
                  current page's cards + lazy 3Dmol viewers. Viewers are created
                  on scroll-in and disposed on scroll-out / page change so the
                  browser never exceeds its WebGL-context limit (~16).

Structures (.mol2) and depictions (.png) are fetched on demand from struct/ and
img/ — never embedded — so the shell stays tiny regardless of DB size.
"""
from __future__ import annotations
import argparse, json, os


def _light_record(r):
    """Trim a full manifest record down to what the review card needs."""
    ligs = []
    for L in r.get("ligands", []):
        o = L.get("oracle") or {}
        ligs.append({
            "smi": (L.get("smiles") or "")[:60],
            "role": L.get("role"),
            "via": L.get("via"),
            "sites": L.get("sites"),
            "comp": L.get("comp", {}),
            "prob": o.get("prob"),
            "hemi": bool(L.get("hemilabile")),
            "hemi_prob": o.get("hemi_prob"),
        })
    isos = []
    for iso in r.get("isomers", []):
        g = iso.get("gates", {})
        f = iso.get("files", {})
        isos.append({
            "label": iso.get("label"),
            "valid": iso.get("valid"),
            "energy": iso.get("energy"),
            "q": iso.get("total_charge"),
            "mult": iso.get("mult"),
            "n_atoms": iso.get("n_atoms"),
            "method": iso.get("method"),
            "ent": bool(iso.get("enantiomer")),
            "gates": {"no_clash": g.get("no_clash"), "bonds": g.get("bond_lengths_ok")},
            "mol2": f.get("mol2"),
            "xyz": f.get("xyz"),
            "h5": f.get("trexio"),
            "txt": f.get("trexio_txt"),
        })
    rm = r.get("retry_meta") or {}
    return {
        "id": r.get("id"),
        "metal": r.get("metal"),
        "ox": r.get("ox"),
        "charge": r.get("charge_complex"),
        "geom": r.get("geometry"),
        "cn": r.get("cn"),
        "conf": r.get("assign_conf"),
        "hemi": bool(r.get("hemilabile")),
        "status": r.get("status"),
        "png": r.get("png"),
        "archive": r.get("archive"),
        "n_valid": sum(1 for i in r.get("isomers", []) if i.get("valid")),
        "ligands": ligs,
        "isomers": isos,
        # failure-reason surfacing (set by retry_tail for unbuildable complexes)
        "fail_reason": r.get("fail_reason"),
        "fail_label": r.get("fail_label"),
        "fail_hint": r.get("fail_hint"),
        "fail_flags": r.get("fail_flags") or [],
        "recovered_by": rm.get("recovered_by"),
    }


def build(out_dir, manifest_path=None):
    man = json.load(open(manifest_path or os.path.join(out_dir, "manifest.json")))
    recs = [_light_record(r) for r in man["records"]]
    recs.sort(key=lambda r: (r["metal"] or "", r["id"] or 0))

    # summary by metal/status for the filter chips + headline
    import collections
    by_metal = collections.Counter(r["metal"] for r in recs)
    by_status = collections.Counter(r["status"] for r in recs)
    by_fail = collections.Counter(r["fail_reason"] for r in recs
                                  if r["status"] != "ok" and r.get("fail_reason"))
    n_ok = by_status.get("ok", 0)
    n_struct = sum(len(r["isomers"]) for r in recs)
    n_valid = sum(r["n_valid"] for r in recs)
    n_recovered = sum(1 for r in recs if r.get("recovered_by"))
    summary = {
        "n": len(recs), "n_ok": n_ok, "n_struct": n_struct, "n_valid": n_valid,
        "n_recovered": n_recovered,
        "by_metal": dict(by_metal.most_common()), "by_status": dict(by_status),
        "by_fail_reason": dict(by_fail.most_common()),
    }

    index_path = os.path.join(out_dir, "index.json")
    json.dump({"summary": summary, "records": recs}, open(index_path, "w"))

    metals = "".join(
        f'<button class="chip" data-metal="{m}">{m} <i>{c}</i></button>'
        for m, c in by_metal.most_common())
    statuses = "".join(
        f'<button class="chip schip" data-status="{s}">{s} <i>{c}</i></button>'
        for s, c in by_status.most_common())
    fails = "".join(
        f'<button class="chip fchip" data-fail="{s}">{s} <i>{c}</i></button>'
        for s, c in by_fail.most_common())

    html = """<!doctype html><html><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BiometalDB — 3D reconstruction review (full DB)</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body{font-family:system-ui,sans-serif;margin:0;background:#0f172a;color:#e2e8f0}
.top{position:sticky;top:0;background:#1e293b;padding:10px 16px;border-bottom:2px solid #6366f1;z-index:10}
.top b{font-size:1.05rem} .hl{color:#a5b4fc}
.bar{margin-top:8px;display:flex;flex-wrap:wrap;gap:6px;align-items:center}
.chip{background:#0f172a;color:#cbd5e1;border:1px solid #334155;border-radius:14px;padding:3px 10px;font-size:.78rem;cursor:pointer}
.chip i{color:#64748b;font-style:normal;font-size:.72rem}
.chip.on{background:#4338ca;color:#fff;border-color:#6366f1}
.chip.schip.on{background:#0e7490;border-color:#22d3ee}
.chip.fchip{border-color:#7c2d12;color:#fca5a5}.chip.fchip.on{background:#b91c1c;color:#fff;border-color:#f87171}
.top input{background:#0f172a;color:#e2e8f0;border:1px solid #334155;border-radius:5px;padding:3px 8px}
.top button.act{background:#16a34a;color:#fff;border:0;padding:4px 10px;border-radius:5px;cursor:pointer;font-size:.78rem}
.top button.gho{background:transparent;color:#94a3b8;border:1px solid #334155;padding:4px 10px;border-radius:5px;cursor:pointer;font-size:.78rem}
.nav{display:flex;gap:8px;align-items:center;margin-left:auto;font-size:.8rem}
.nav button{background:#334155;color:#e2e8f0;border:0;border-radius:5px;padding:3px 9px;cursor:pointer}
.nav button:disabled{opacity:.4;cursor:default}
.card{background:#1e293b;margin:12px;border-radius:8px;padding:10px 14px;border-left:4px solid #475569}
.card.ok{border-left-color:#16a34a}.card.no_structure{border-left-color:#b45309}
.card.crash_quarantine{border-left-color:#dc2626}.card.timeout_or_crash{border-left-color:#dc2626}.card.exception{border-left-color:#dc2626}
.h{font-size:1rem;margin-bottom:6px} .body{display:flex;gap:14px;align-items:flex-start;flex-wrap:wrap}
.left{flex:0 0 280px} .left img{width:280px;background:#fff;border-radius:4px;min-height:60px}
.right{flex:1;display:flex;gap:12px;flex-wrap:wrap}
.iso{background:#0f172a;border-radius:6px;padding:6px;width:540px;max-width:100%}
.ih{font-size:.82rem;color:#cbd5e1;margin-bottom:4px}
.v{width:528px;max-width:100%;height:440px;position:relative;background:#020617;border-radius:4px;cursor:pointer;display:flex;align-items:center;justify-content:center}
.v .hint{color:#64748b;font-size:1.3rem;letter-spacing:1px}
.lt{font-size:.7rem;border-collapse:collapse;margin-top:6px;width:100%}
.lt th,.lt td{border:1px solid #334155;padding:2px 4px;text-align:left}
.smi{font-family:monospace;color:#93c5fd}
.b{font-size:.68rem;padding:1px 6px;border-radius:3px;margin-left:4px}
.b.ok{background:#166534}.b.no{background:#7f1d1d}.b.hemi{background:#7c2d12}
.st{font-size:.7rem;padding:1px 6px;border-radius:3px;background:#334155}.st.ok{background:#166534}
.ctrl{margin-top:8px;font-size:.8rem}.ctrl button{margin:0 4px;cursor:pointer}
.acc{background:#16a34a;color:#fff;border:0;padding:3px 8px;border-radius:4px}
.rej{background:#dc2626;color:#fff;border:0;padding:3px 8px;border-radius:4px}
.flg{background:#d97706;color:#fff;border:0;padding:3px 8px;border-radius:4px}
.note{width:300px;background:#0f172a;color:#e2e8f0;border:1px solid #334155;border-radius:4px;padding:2px 6px}
.rv{font-weight:700;margin:0 8px}
.b.rec{background:#0e7490;color:#fff}
.reason{flex:1;min-width:300px;background:#1c1410;border:1px solid #b45309;border-left:4px solid #f59e0b;border-radius:6px;padding:10px 14px}
.reason.hard{background:#1f1112;border-color:#dc2626;border-left-color:#ef4444}
.reason .rt{font-weight:700;color:#fbbf24;font-size:.92rem}
.reason.hard .rt{color:#fca5a5}
.reason .rc{font-family:monospace;font-size:.68rem;color:#94a3b8;margin-left:6px}
.reason .rh{color:#cbd5e1;font-size:.82rem;margin-top:5px;line-height:1.45}
.reason .rf{margin-top:6px;font-size:.68rem;color:#64748b;font-family:monospace;word-break:break-all}
a{color:#93c5fd}
#empty{padding:40px;text-align:center;color:#64748b}
</style></head><body>
<div class="top">
<b>BiometalDB — genuine metal-coordinated 3D · <span class="hl" id="hcount">loading…</span></b>
<div class="bar">
<span style="font-size:.78rem;color:#94a3b8">metal:</span> __METALS__
<button class="chip" id="metalAll">all</button>
</div>
<div class="bar">
<span style="font-size:.78rem;color:#94a3b8">why-failed:</span> __FAILS__
<button class="chip fchip" id="failAll">all</button>
</div>
<div class="bar">
<span style="font-size:.78rem;color:#94a3b8">status:</span> __STATUSES__
<button class="chip schip" id="statusAll">all</button>
<input id="q" placeholder="#id or substring..." style="margin-left:8px">
<button class="act" id="acceptValid">✓ validate all passing (page)</button>
<button class="gho" id="exp">⬇ export reviews</button>
<button class="gho" id="clr">↺ clear</button>
<span class="nav"><button id="prev">‹ prev</button><span id="pageinfo"></span><button id="next">next ›</button></span>
</div>
</div>
<div id="list"></div>
<div id="empty" style="display:none">no complexes match the current filter</div>
<script>
const PAGE_SIZE = __PAGE_SIZE__;
const LS = 'biometaldb_full_review';
let DATA = [], VIEW = {}, page = 0, fMetal = null, fStatus = null, fFail = null, fQ = '';
let R = JSON.parse(localStorage.getItem(LS) || '{}');

const esc = s => (s==null?'':(''+s)).replace(/[&<>"]/g, c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
function save(){ localStorage.setItem(LS, JSON.stringify(R)); }

function filtered(){
  return DATA.filter(r=>{
    if(fMetal && r.metal!==fMetal) return false;
    if(fStatus && r.status!==fStatus) return false;
    if(fFail && r.fail_reason!==fFail) return false;
    if(fQ){ const q=fQ.toLowerCase();
      if(!((''+r.id)===fQ.replace('#','') || (''+r.id).includes(fQ.replace('#','')) ||
           JSON.stringify(r.ligands).toLowerCase().includes(q))) return false; }
    return true;
  });
}

function badge(t, ok){ return `<span class="b ${ok?'ok':'no'}">${t}</span>`; }

const HARD_FAILS = {native_runtime_failure:1, corrupt_source_donors:1, unparseable_smiles:1};
function reasonHTML(r){
  const code = r.fail_reason || 'no_structure';
  const label = r.fail_label || 'No 3D structure built';
  const hint = r.fail_hint || 'Architector produced no conformer for this complex.';
  const hard = HARD_FAILS[code] ? ' hard' : '';
  const flags = (r.fail_flags && r.fail_flags.length)
    ? `<div class="rf">assignment flags: ${esc(r.fail_flags.join(' · '))}</div>` : '';
  return `<div class="reason${hard}"><div class="rt">⚠ ${esc(label)}`+
    `<span class="rc">[${esc(code)}]</span></div>`+
    `<div class="rh">${esc(hint)}</div>${flags}</div>`;
}
function cardHTML(r){
  const hemi = r.hemi ? '<span class="b hemi">hemilabile</span>' : '';
  const rec = r.recovered_by ? `<span class="b rec">↻ recovered: ${esc(r.recovered_by)}</span>` : '';
  const arc = r.archive ? ` · <a href="archive/${esc(r.archive)}" download>⬇ .zip</a>` : '';
  const head = `<div class="h"><b>#${r.id}</b> ${esc(r.metal)}(${esc(r.ox)}) · `+
    `${esc(r.geom||'?')} CN${esc(r.cn||'?')} · conf=${esc(r.conf||'?')} ${hemi} `+
    `<span class="st ${esc(r.status)}">${esc(r.status)}</span>${rec}${arc}</div>`;
  let ligrows='';
  for(const L of (r.ligands||[])){
    ligrows += `<tr><td class="smi">${esc(L.smi)}</td><td>${esc(L.role)}</td>`+
      `<td>${esc(L.via)}</td><td>${esc(L.sites)}</td>`+
      `<td>${esc(JSON.stringify(L.comp||{}))}</td>`+
      `<td>${L.hemi?('hemi p='+esc(L.hemi_prob)+' '):''}${L.prob!=null?('p='+esc(L.prob)):''}</td></tr>`;
  }
  const ligtbl = `<table class="lt"><tr><th>ligand</th><th>role</th><th>via</th><th>sites</th><th>donors</th><th>pydentate</th></tr>${ligrows}</table>`;
  let isos='';
  for(const iso of (r.isomers||[])){
    const key = r.id+'__'+iso.label;
    const vb = badge('valid',iso.valid)+badge('no-clash',iso.gates.no_clash)+badge('bonds',iso.gates.bonds);
    const dl = (iso.txt?` · <a href="struct/${esc(iso.txt)}" target="_blank">TREXIO·txt</a>`:'')+
      (iso.h5?` · <a href="struct/${esc(iso.h5)}" download>·h5</a>`:'')+
      (iso.xyz?` · <a href="struct/${esc(iso.xyz)}" download>xyz</a>`:'');
    const has3d = iso.mol2 ? `data-mol2="${esc(iso.mol2)}" onclick="show3d('${esc(key)}')"` : '';
    isos += `<div class="iso"><div class="ih">${esc(iso.label)} · E=${esc(iso.energy)} · q=${esc(iso.q)} · mult=${esc(iso.mult)} ${vb}${dl}</div>`+
      `<div class="v" id="v_${esc(key)}" data-key="${esc(key)}" ${has3d}>`+
      `<span class="hint">${iso.mol2?'⬢ 3D':'no structure'}</span></div></div>`;
  }
  const ctrl = `<div class="ctrl" data-id="${r.id}">review: `+
    `<button class="acc" onclick="setv(${r.id},'accept')">✓ accept</button>`+
    `<button class="rej" onclick="setv(${r.id},'reject')">✗ reject</button>`+
    `<button class="flg" onclick="setv(${r.id},'flag')">⚑ flag</button>`+
    `<span class="rv" id="rv_${r.id}"></span>`+
    `<input class="note" id="nt_${r.id}" placeholder="notes..." oninput="note(${r.id})"></div>`;
  const img = r.png ? `<img loading="lazy" src="img/${esc(r.png)}"/>` : '<div style="color:#64748b">no depiction</div>';
  const right = isos || reasonHTML(r);
  return `<div class="card ${esc(r.status)}" id="card_${r.id}">${head}`+
    `<div class="body"><div class="left">${img}${ligtbl}</div><div class="right">${right}</div></div>${ctrl}</div>`;
}

let io = new IntersectionObserver(ents=>{ents.forEach(e=>{
  const k=e.target.dataset.key; if(!k) return;
  if(e.isIntersecting) show3d(k); else dispose(k);
});},{rootMargin:'300px 0px'});

function show3d(key){
  const el=document.getElementById('v_'+key); if(!el || el.dataset.inited) return;
  const mol2=el.dataset.mol2; if(!mol2) return;
  el.dataset.inited='1'; el.innerHTML='<span class="hint">loading 3D…</span>';
  fetch('struct/'+mol2).then(r=>r.text()).then(txt=>{
    el.innerHTML='';
    const v=$3Dmol.createViewer(el,{backgroundColor:'#020617'});
    v.addModel(txt,'mol2');
    v.setStyle({},{stick:{radius:0.14},sphere:{scale:0.27}});
    v.setStyle({elem:['Ir','Ru','Os','Rh','Re','Au']},{sphere:{scale:0.5,color:'gold'}});
    v.zoomTo(); v.render(); VIEW[key]=v;
  }).catch(()=>{el.dataset.inited='';el.innerHTML='<span class="hint">3D load failed</span>';});
}
function dispose(key){
  if(VIEW[key]){try{VIEW[key].clear();}catch(e){} delete VIEW[key];}
  const el=document.getElementById('v_'+key);
  if(el){el.dataset.inited='';el.innerHTML='<span class="hint">⬢ 3D</span>';}
}
function disposeAll(){ Object.keys(VIEW).forEach(dispose); }

function render(){
  disposeAll();
  const list=document.getElementById('list');
  const f=filtered();
  const pages=Math.max(1,Math.ceil(f.length/PAGE_SIZE));
  if(page>=pages) page=pages-1; if(page<0) page=0;
  const slice=f.slice(page*PAGE_SIZE,(page+1)*PAGE_SIZE);
  document.getElementById('empty').style.display=f.length?'none':'block';
  list.innerHTML=slice.map(cardHTML).join('');
  document.getElementById('pageinfo').textContent=` ${f.length?page+1:0}/${pages} (${f.length}) `;
  document.getElementById('prev').disabled=page<=0;
  document.getElementById('next').disabled=page>=pages-1;
  slice.forEach(r=>{paint(r.id); const n=R[r.id]&&R[r.id].note; if(n){const el=document.getElementById('nt_'+r.id); if(el)el.value=n;}});
  document.querySelectorAll('.v').forEach(el=>io.observe(el));
  window.scrollTo(0,0);
}

function paint(id){ const s=R[id]&&R[id].v||''; const e=document.getElementById('rv_'+id);
  if(!e) return; e.textContent=s?('→ '+s):''; e.style.color=s=='accept'?'#22c55e':s=='reject'?'#ef4444':'#f59e0b'; }
function setv(id,v){ R[id]=R[id]||{}; R[id].v=v; save(); paint(id); count(); }
function note(id){ R[id]=R[id]||{}; R[id].note=document.getElementById('nt_'+id).value; save(); }
function count(){ let a=0,r=0,f=0; for(const k in R){if(R[k].v=='accept')a++;if(R[k].v=='reject')r++;if(R[k].v=='flag')f++;}
  document.getElementById('hcount').textContent=`${DATA.length} complexes · ✓${a} ✗${r} ⚑${f}`; }
function exp(){ const b=new Blob([JSON.stringify(R,null,2)],{type:'application/json'});
  const a=document.createElement('a'); a.href=URL.createObjectURL(b); a.download='biometaldb_reviews.json'; a.click(); }

document.getElementById('prev').onclick=()=>{page--;render();};
document.getElementById('next').onclick=()=>{page++;render();};
document.getElementById('q').oninput=e=>{fQ=e.target.value.trim();page=0;render();};
document.getElementById('acceptValid').onclick=()=>{ filtered().slice(page*PAGE_SIZE,(page+1)*PAGE_SIZE)
  .filter(r=>r.status=='ok'&&r.isomers.length&&r.isomers.every(i=>i.valid))
  .forEach(r=>{R[r.id]=R[r.id]||{};R[r.id].v='accept';}); save(); render(); count(); };
document.getElementById('clr').onclick=()=>{ if(!confirm('Clear all review marks & notes?'))return; R={}; save(); render(); count(); };
function wireChips(sel,attr,setter){ document.querySelectorAll(sel).forEach(c=>{ c.onclick=()=>{
  const v=c.dataset[attr]; setter(v===undefined?null:v);
  document.querySelectorAll(sel).forEach(x=>x.classList.remove('on')); if(v!==undefined)c.classList.add('on');
  page=0; render(); };});}

fetch('index.json').then(r=>r.json()).then(j=>{
  DATA=j.records; const s=j.summary;
  document.getElementById('hcount').textContent=`${s.n} complexes · ${s.n_ok} built · ${s.n_struct} structures · ${s.n_valid} valid`+(s.n_recovered?` · ↻${s.n_recovered} recovered`:'');
  wireChips('.chip[data-metal]','metal',v=>fMetal=v);
  wireChips('.chip[data-status]','status',v=>fStatus=v);
  wireChips('.chip[data-fail]','fail',v=>fFail=v);
  document.getElementById('metalAll').onclick=()=>{fMetal=null;document.querySelectorAll('.chip[data-metal]').forEach(x=>x.classList.remove('on'));page=0;render();};
  document.getElementById('statusAll').onclick=()=>{fStatus=null;document.querySelectorAll('.chip[data-status]').forEach(x=>x.classList.remove('on'));page=0;render();};
  document.getElementById('failAll').onclick=()=>{fFail=null;document.querySelectorAll('.chip[data-fail]').forEach(x=>x.classList.remove('on'));page=0;render();};
  count(); render();
});
</script></body></html>"""
    html = (html.replace("__METALS__", metals).replace("__STATUSES__", statuses)
                .replace("__FAILS__", fails or '<i style="color:#475569;font-size:.72rem">none</i>')
                .replace("__PAGE_SIZE__", "48"))
    out = os.path.join(out_dir, "review.html")
    open(out, "w").write(html)
    print(f"wrote {out} ({len(html)//1024} KB shell) + {index_path} "
          f"({os.path.getsize(index_path)//1024} KB, {len(recs)} records)")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/full")
    ap.add_argument("--manifest", default="")
    a = ap.parse_args()
    build(a.out, a.manifest or None)
