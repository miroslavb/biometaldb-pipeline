"""Write a consolidated text TREXIO record (.trexio.txt) for every structure."""
import json, os, sys
sys.path.insert(0, "recon3d")
import trexio_writer as TX

OUT = sys.argv[1] if len(sys.argv) > 1 else "out/ir100"
STRUCT = os.path.join(OUT, "struct")


def read_xyz(p):
    L = open(p).read().splitlines(); n = int(L[0]); sym, C = [], []
    for ln in L[2:2 + n]:
        q = ln.split(); sym.append(q[0]); C.append([float(x) for x in q[1:4]])
    return sym, C


man = json.load(open(os.path.join(OUT, "manifest.json")))
cnt = 0
for r in man["records"]:
    for iso in r.get("isomers", []):
        if iso.get("files", {}).get("trexio_txt"):
            continue  # idempotent: already written
        xyzf = iso.get("files", {}).get("xyz")
        if not xyzf or not os.path.exists(os.path.join(STRUCT, xyzf)):
            continue
        sym, C = read_xyz(os.path.join(STRUCT, xyzf))
        stem = xyzf[:-4]
        txt = stem + ".trexio.txt"
        TX.write_trexio_text(os.path.join(STRUCT, txt), sym, C,
                             iso.get("total_charge"), (iso.get("mult", 1) - 1),
                             title=f"BiometalDB #{r['id']} {iso['label']}")
        iso["files"]["trexio_txt"] = txt
        cnt += 1
json.dump(man, open(os.path.join(OUT, "manifest.json"), "w"), indent=2)
print(f"wrote {cnt} TREXIO text records")
