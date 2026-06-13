"""Build the ligand-library web interface (data-driven) from ligands.json."""
import argparse, os

HTML = """<!doctype html><html><head><meta charset="utf-8">
<title>BiometalDB — Ligand Library</title>
<style>
body{font-family:system-ui,sans-serif;margin:0;background:#0f172a;color:#e2e8f0}
.top{position:sticky;top:0;background:#1e293b;padding:10px 16px;border-bottom:2px solid #6366f1;z-index:10;display:flex;gap:12px;align-items:center;flex-wrap:wrap}
.top b{font-size:1.05rem} input,select{background:#0f172a;color:#e2e8f0;border:1px solid #334155;border-radius:5px;padding:5px 8px}
input#q{width:280px}
details{margin:10px 12px;background:#1e293b;border-radius:8px;border-left:4px solid #475569}
summary{cursor:pointer;padding:10px 14px;font-size:1rem;font-weight:600;list-style:none}
summary::-webkit-details-marker{display:none} summary:before{content:"▸ ";color:#6366f1}
details[open] summary:before{content:"▾ "}
summary .n{color:#94a3b8;font-weight:400;font-size:.85rem;margin-left:6px}
.grid{display:flex;flex-wrap:wrap;gap:10px;padding:0 14px 14px}
.card{background:#0f172a;border:1px solid #334155;border-radius:6px;width:210px;padding:6px}
.card img{width:198px;height:150px;object-fit:contain;background:#fff;border-radius:4px}
.smi{font-family:monospace;font-size:.62rem;color:#93c5fd;word-break:break-all;max-height:2.4em;overflow:hidden;margin:3px 0}
.meta{font-size:.66rem;color:#cbd5e1} .b{font-size:.6rem;padding:1px 5px;border-radius:3px;margin-right:3px}
.b.cnt{background:#1e3a8a} .b.dent{background:#374151} .b.hemi{background:#7c2d12} .b.met{background:#334155}
a{color:#60a5fa;text-decoration:none} .muted{color:#64748b}
</style></head><body>
<div class="top"><b>🧬 BiometalDB Ligand Library</b>
<span id="stat" class="muted"></span>
<input id="q" placeholder="filter by SMILES…" oninput="render()">
<select id="metal" onchange="render()"><option value="">all metals</option></select>
<select id="grp" onchange="jump()"><option value="">jump to group…</option></select>
<button onclick="toggleAll(true)">expand all</button><button onclick="toggleAll(false)">collapse all</button>
</div>
<div id="list" class="muted" style="padding:20px">loading…</div>
<script>
let DATA=null;
fetch('ligands.json').then(r=>r.json()).then(d=>{DATA=d;init();});
function init(){
  let mets=new Set(); DATA.ligands.forEach(l=>l.metals.forEach(m=>mets.add(m)));
  let ms=document.getElementById('metal'); [...mets].sort().forEach(m=>ms.add(new Option(m,m)));
  let gs=document.getElementById('grp'); DATA.group_order.forEach(g=>{if(DATA.group_counts[g])gs.add(new Option(g+' ('+DATA.group_counts[g]+')',g));});
  render();
}
function card(l){
  let img=l.png?`<img loading="lazy" src="img/${l.png}">`:'<div class="muted" style="height:150px">no 2D</div>';
  let dent=l.denticity!=null?`<span class="b dent">κ${l.denticity} ${(l.donors||[]).join('')}</span>`:'';
  let hemi=l.hemilabile?'<span class="b hemi">hemilabile</span>':'';
  return `<div class="card">${img}<div class="smi" title="${l.smiles}">${l.smiles}</div>
   <div class="meta"><span class="b cnt">×${l.count}</span>${dent}${hemi}
   <span class="b met">${l.metals.join(',')}</span>
   <a href="/complexes/${l.example_id}" target="_blank">ex #${l.example_id}</a></div></div>`;
}
function render(){
  let q=document.getElementById('q').value.toLowerCase(), mf=document.getElementById('metal').value;
  let byg={}; let shown=0;
  for(let l of DATA.ligands){
    if(q && !l.smiles.toLowerCase().includes(q)) continue;
    if(mf && !l.metals.includes(mf)) continue;
    (byg[l.group]=byg[l.group]||[]).push(l); shown++;
  }
  let html='';
  for(let g of DATA.group_order){
    let arr=byg[g]; if(!arr||!arr.length) continue;
    let open=(q||mf)?'open':'';
    html+=`<details ${open}><summary>${g}<span class="n">${arr.length} ligand${arr.length>1?'s':''}</span></summary>
      <div class="grid">${arr.map(card).join('')}</div></details>`;
  }
  document.getElementById('list').innerHTML=html||'<div class="muted" style="padding:20px">no matches</div>';
  document.getElementById('stat').textContent=`${DATA.n_ligands} unique ligands · ${shown} shown`;
}
function jump(){let g=document.getElementById('grp').value;if(!g)return;
  let d=[...document.querySelectorAll('details')].find(x=>x.querySelector('summary').textContent.startsWith(g));
  if(d){d.open=true;d.scrollIntoView();}}
function toggleAll(o){document.querySelectorAll('details').forEach(d=>d.open=o);}
</script></body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="out/ligands")
    a = ap.parse_args()
    open(os.path.join(a.out, "index.html"), "w").write(HTML)
    print("wrote", os.path.join(a.out, "index.html"))


if __name__ == "__main__":
    main()
