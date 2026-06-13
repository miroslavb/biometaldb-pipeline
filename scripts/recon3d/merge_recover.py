"""Merge recovered complex records into the main ir100 manifest + assets."""
import json, os, shutil, sys, collections

MAIN = "out/ir100"
REC = sys.argv[1] if len(sys.argv) > 1 else "out/recover"

main = json.load(open(f"{MAIN}/manifest.json"))
rec = json.load(open(f"{REC}/manifest.json"))
by_id = {r["id"]: r for r in main["records"]}

merged = 0
for r in rec["records"]:
    if r.get("status") == "ok" and r.get("isomers"):
        by_id[r["id"]] = r
        for iso in r["isomers"]:
            for kind, fn in iso.get("files", {}).items():
                src = os.path.join(REC, "struct", fn)
                if os.path.exists(src):
                    shutil.copy(src, os.path.join(MAIN, "struct", fn))
        if r.get("png"):
            p = os.path.join(REC, "img", r["png"])
            if os.path.exists(p):
                shutil.copy(p, os.path.join(MAIN, "img", r["png"]))
        merged += 1

recs = sorted(by_id.values(), key=lambda r: r["id"])
main["records"] = recs
main["n_ok"] = sum(1 for r in recs if r.get("status") == "ok")
main["n_isomers_total"] = sum(len(r.get("isomers", [])) for r in recs)
main["n_valid_struct"] = sum(1 for r in recs for i in r.get("isomers", []) if i.get("valid"))
json.dump(main, open(f"{MAIN}/manifest.json", "w"), indent=2)
d = collections.Counter(len(r.get("isomers", [])) for r in recs)
print(f"merged {merged} | ok={main['n_ok']}/{main['n']} structures={main['n_isomers_total']} "
      f"valid={main['n_valid_struct']} isomers/complex={dict(d)} "
      f"not-ok={[r['id'] for r in recs if r.get('status')!='ok']}")
