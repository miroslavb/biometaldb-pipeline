"""Build a per-compound .zip bundling all its artifacts, link it in the manifest.

Each complex_<id>.zip contains every isomer/enantiomer's geometry (xyz, mol2,
sdf), the TREXIO records (.h5 binary + .trexio.txt text), the 2D depiction, and
a metadata.json (assignment + per-isomer validation). Served from viewer/.
"""
import json, os, sys, zipfile

OUT = sys.argv[1] if len(sys.argv) > 1 else "out/ir100"
STRUCT = os.path.join(OUT, "struct")
IMG = os.path.join(OUT, "img")
ARC = os.path.join(OUT, "archive")
os.makedirs(ARC, exist_ok=True)

man = json.load(open(os.path.join(OUT, "manifest.json")))
n = 0
for r in man["records"]:
    cid = r["id"]
    root = f"complex_{cid}"
    files = []
    for iso in r.get("isomers", []):
        for fn in iso.get("files", {}).values():
            if fn and ("struct", fn) not in files:
                files.append(("struct", fn))
    if r.get("png"):
        files.append(("img", r["png"]))
    zpath = os.path.join(ARC, f"{root}.zip")
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        for kind, fn in files:
            src = os.path.join(STRUCT if kind == "struct" else IMG, fn)
            if os.path.exists(src):
                z.write(src, arcname=f"{root}/{fn}")
        z.writestr(f"{root}/metadata.json", json.dumps(r, indent=2))
        readme = (f"BiometalDB complex #{cid} — {r.get('metal')}({r.get('ox')}) "
                  f"{r.get('geometry')}, CN {r.get('cn')}\n"
                  f"isomers: {', '.join(i['label'] for i in r.get('isomers', [])) or 'none'}\n"
                  f"files per isomer: .xyz (coords) · .mol2/.sdf (with bonds) · "
                  f".h5 (TREXIO binary) · .trexio.txt (TREXIO text)\n"
                  f"status: {r.get('status')}; reconstructed with Architector+GFN2-xTB, "
                  f"coordinating atoms via pydentate.\n")
        z.writestr(f"{root}/README.txt", readme)
    r["archive"] = f"{root}.zip"
    n += 1
json.dump(man, open(os.path.join(OUT, "manifest.json"), "w"), indent=2)
print(f"wrote {n} per-compound archives -> {ARC}")
