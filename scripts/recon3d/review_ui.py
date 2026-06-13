"""Generate a self-contained manual-review HTML interface from an ir100 manifest.

Per complex: 2D depiction, each stereoisomer in a lazy-loaded 3Dmol.js viewer
(structure embedded inline -> works from file://), pydentate coordination
prediction + hemilability per ligand, validation badges, TREXIO download, and
Accept/Reject/Flag + notes controls (persisted to localStorage, exportable).
"""
from __future__ import annotations
import argparse, base64, json, os


def _b64png(path):
    try:
        return "data:image/png;base64," + base64.b64encode(open(path, "rb").read()).decode()
    except Exception:
        return ""


def build(out_dir, manifest_path=None):
    man = json.load(open(manifest_path or os.path.join(out_dir, "manifest.json")))
    recs = man["records"]
    struct, imgs = {}, {}
    for r in recs:
        if r.get("png"):
            imgs[r["id"]] = _b64png(os.path.join(out_dir, "img", r["png"]))
        for iso in r.get("isomers", []):
            mol2 = iso.get("files", {}).get("mol2")
            if mol2:
                try:
                    struct[f"{r['id']}_{iso['label']}"] = open(os.path.join(out_dir, "struct", mol2)).read()
                except Exception:
                    pass

    cards = []
    for r in recs:
        rid = r["id"]
        badge = lambda t, ok: f'<span class="b {"ok" if ok else "no"}">{t}</span>'
        hemi = '<span class="b hemi">hemilabile</span>' if r.get("hemilabile") else ""
        head = (f'<div class="h"><b>#{rid}</b> {r.get("metal")}({r.get("ox")}) · '
                f'{r.get("geometry","?")} CN{r.get("cn","?")} · conf={r.get("assign_conf","?")} '
                f'{hemi} <span class="st {r.get("status")}">{r.get("status")}</span></div>')
        # ligand table
        lig_rows = ""
        for L in r.get("ligands", []):
            o = L.get("oracle") or {}
            alt = ""
            lig_rows += (f'<tr><td class="smi">{L["smiles"][:48]}</td><td>{L["role"]}</td>'
                         f'<td>{L.get("via")}</td><td>{L.get("sites")}</td>'
                         f'<td>{json.dumps(L.get("comp",{}))}</td>'
                         f'<td>{("hemi p="+str(o.get("hemi_prob")) ) if L.get("hemilabile") else ""} '
                         f'{("p="+str(o.get("prob"))) if o else ""}</td></tr>')
        ligtbl = (f'<table class="lt"><tr><th>ligand</th><th>role</th><th>via</th>'
                  f'<th>sites</th><th>donors</th><th>pydentate</th></tr>{lig_rows}</table>')
        # isomers
        isos = ""
        for iso in r.get("isomers", []):
            key = f'{rid}_{iso["label"]}'
            g = iso.get("gates", {})
            vb = (badge("valid", iso.get("valid")) + badge("no-clash", g.get("no_clash")) +
                  badge("bonds", g.get("bond_lengths_ok")))
            tx = iso.get("files", {}).get("trexio", "")
            xyz = iso.get("files", {}).get("xyz", "")
            dl = (f' · <a href="struct/{tx}" download>TREXIO</a>' if tx else "") + \
                 (f' · <a href="struct/{xyz}" download>xyz</a>' if xyz else "")
            isos += (f'<div class="iso"><div class="ih">{iso["label"]} '
                     f'· E={iso.get("energy")} · q={iso.get("total_charge")} '
                     f'· mult={iso.get("mult")} {vb}{dl}</div>'
                     f'<div class="v" id="v_{key}"></div>'
                     f'<button onclick="show3d(\'{key}\')">⬢ 3D</button></div>')
        ctrl = (f'<div class="ctrl" data-id="{rid}">review: '
                f'<button class="acc" onclick="setv({rid},\'accept\')">✓ accept</button>'
                f'<button class="rej" onclick="setv({rid},\'reject\')">✗ reject</button>'
                f'<button class="flg" onclick="setv({rid},\'flag\')">⚑ flag</button>'
                f'<span class="rv" id="rv_{rid}"></span>'
                f'<input class="note" id="nt_{rid}" placeholder="notes..." '
                f'oninput="note({rid})"></div>')
        cards.append(f'<div class="card" id="card_{rid}">{head}'
                     f'<div class="body"><div class="left"><img src="{imgs.get(rid,"")}"/>'
                     f'{ligtbl}</div><div class="right">{isos}</div></div>{ctrl}</div>')

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>BiometalDB Ir — 3D reconstruction review</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body{{font-family:system-ui,sans-serif;margin:0;background:#0f172a;color:#e2e8f0}}
.top{{position:sticky;top:0;background:#1e293b;padding:10px 16px;border-bottom:2px solid #6366f1;z-index:10}}
.top b{{font-size:1.1rem}} .top button{{margin-left:8px}}
.card{{background:#1e293b;margin:12px;border-radius:8px;padding:10px 14px;border-left:4px solid #475569}}
.h{{font-size:1rem;margin-bottom:6px}} .body{{display:flex;gap:14px}}
.left{{flex:0 0 540px}} .left img{{width:540px;background:#fff;border-radius:4px}}
.right{{flex:1;display:flex;gap:10px;flex-wrap:wrap}}
.iso{{background:#0f172a;border-radius:6px;padding:6px;width:330px}}
.ih{{font-size:.78rem;color:#cbd5e1;margin-bottom:4px}}
.v{{width:316px;height:300px;position:relative;background:#020617;border-radius:4px}}
.lt{{font-size:.7rem;border-collapse:collapse;margin-top:6px;width:100%}}
.lt th,.lt td{{border:1px solid #334155;padding:2px 4px;text-align:left}}
.smi{{font-family:monospace;color:#93c5fd}}
.b{{font-size:.68rem;padding:1px 6px;border-radius:3px;margin-left:4px}}
.b.ok{{background:#166534}} .b.no{{background:#7f1d1d}} .b.hemi{{background:#7c2d12}}
.st{{font-size:.7rem;padding:1px 6px;border-radius:3px}} .st.ok{{background:#166534}}
.ctrl{{margin-top:8px;font-size:.8rem}} .ctrl button{{margin:0 4px;cursor:pointer}}
.acc{{background:#16a34a;color:#fff;border:0;padding:3px 8px;border-radius:4px}}
.rej{{background:#dc2626;color:#fff;border:0;padding:3px 8px;border-radius:4px}}
.flg{{background:#d97706;color:#fff;border:0;padding:3px 8px;border-radius:4px}}
.note{{width:300px;background:#0f172a;color:#e2e8f0;border:1px solid #334155;border-radius:4px;padding:2px 6px}}
.rv{{font-weight:700;margin:0 8px}}
button{{cursor:pointer}}
</style></head><body>
<div class="top"><b>BiometalDB Ir(III) — 3D reconstruction · manual review</b>
&nbsp; {man['n']} complexes · {man['n_isomers_total']} structures · {man['n_valid_struct']} valid
<button onclick="exp()">⬇ export reviews</button>
<span id="cnt"></span></div>
<div id="list">{''.join(cards)}</div>
<script>
const STRUCT={json.dumps(struct)};
const inited={{}};
function show3d(key){{
  if(inited[key])return; inited[key]=1;
  let el=document.getElementById('v_'+key);
  let v=$3Dmol.createViewer(el,{{backgroundColor:'#020617'}});
  v.addModel(STRUCT[key],'mol2');
  v.setStyle({{}},{{stick:{{radius:0.13}},sphere:{{scale:0.25}}}});
  v.setStyle({{elem:['Ir','Ru','Os','Rh','Re','Au']}},{{sphere:{{scale:0.45,color:'gold'}}}});
  v.zoomTo(); v.render();
}}
let R=JSON.parse(localStorage.getItem('ir100_review')||'{{}}');
function paint(id){{let s=R[id]&&R[id].v||''; let e=document.getElementById('rv_'+id);
  e.textContent=s?('→ '+s):''; e.style.color=s=='accept'?'#22c55e':s=='reject'?'#ef4444':'#f59e0b';
  let n=document.getElementById('nt_'+id); if(R[id]&&R[id].note)n.value=R[id].note;}}
function setv(id,v){{R[id]=R[id]||{{}}; R[id].v=v; save(); paint(id); count();}}
function note(id){{R[id]=R[id]||{{}}; R[id].note=document.getElementById('nt_'+id).value; save();}}
function save(){{localStorage.setItem('ir100_review',JSON.stringify(R));}}
function count(){{let a=0,r=0,f=0; for(k in R){{if(R[k].v=='accept')a++;if(R[k].v=='reject')r++;if(R[k].v=='flag')f++;}}
  document.getElementById('cnt').textContent=` · ✓${{a}} ✗${{r}} ⚑${{f}}`;}}
function exp(){{let b=new Blob([JSON.stringify(R,null,2)],{{type:'application/json'}});
  let a=document.createElement('a'); a.href=URL.createObjectURL(b); a.download='ir100_reviews.json'; a.click();}}
for(let id in {json.dumps({r["id"]: 1 for r in recs})}) paint(id);
count();
</script></body></html>"""
    out = os.path.join(out_dir, "review.html")
    open(out, "w").write(html)
    print(f"wrote {out} ({len(html)//1024} KB, {len(struct)} structures embedded)")
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/ir100")
    ap.add_argument("--manifest", default="")
    a = ap.parse_args()
    build(a.out, a.manifest or None)
