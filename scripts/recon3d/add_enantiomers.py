"""Add the mirror-image enantiomer for every chiral-at-metal complex.

The build/stereo step keeps diastereomers but collapses enantiomers (Λ/Δ share a
trans-pair signature). For a chiral complex the mirror image is a genuinely
distinct 3D structure (same energy, same connectivity), so we generate it by
reflecting coordinates — no extra xtb needed. Achiral complexes (mirror
superimposable by a proper rotation) are left with a single structure.
"""
import json, os, sys, glob
import numpy as np
sys.path.insert(0, "recon3d")
import trexio_writer as TX

OUT = sys.argv[1] if len(sys.argv) > 1 else "out/ir100"
STRUCT = os.path.join(OUT, "struct")
CHIRAL_RMSD = 0.4  # Å; mirror not superimposable by proper rotation => chiral


def read_xyz(path):
    L = open(path).read().splitlines()
    n = int(L[0]); sym, C = [], []
    for ln in L[2:2 + n]:
        p = ln.split(); sym.append(p[0]); C.append([float(x) for x in p[1:4]])
    return sym, np.array(C), L[1]


def rmsd_proper(A, B):
    A0 = A - A.mean(0); B0 = B - B.mean(0)
    H = B0.T @ A0; U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(((A0 - (B0 @ R.T)) ** 2).sum() / len(A0)))


def mirror_mol2(src, dst):
    out = []
    in_atoms = False
    for ln in open(src):
        s = ln.strip()
        if s.startswith("@<TRIPOS>"):
            in_atoms = s.startswith("@<TRIPOS>ATOM")
            out.append(ln); continue
        if in_atoms and s and not s.startswith("@"):
            f = ln.split()
            if len(f) >= 5:
                try:
                    f[2] = f"{-float(f[2]):.4f}"
                    out.append("      " + " ".join(f) + "\n"); continue
                except ValueError:
                    pass
        out.append(ln)
    open(dst, "w").write("".join(out))


def main():
    man = json.load(open(os.path.join(OUT, "manifest.json")))
    added = 0
    for r in man["records"]:
        if r.get("status") != "ok":
            continue
        # idempotent: skip records already enantiomer-enriched (Λ/Δ or flagged)
        if any(i.get("enantiomer") or str(i.get("label", "")).endswith(("Λ", "Δ"))
               for i in r.get("isomers", [])):
            continue
        extra = []
        for iso in list(r["isomers"]):
            if iso["label"].endswith("Δ") or iso.get("enantiomer"):
                continue
            xyzf = iso.get("files", {}).get("xyz")
            if not xyzf:
                continue
            sym, C, comment = read_xyz(os.path.join(STRUCT, xyzf))
            M = C.copy(); M[:, 0] *= -1
            if rmsd_proper(C, M) < CHIRAL_RMSD:
                continue  # achiral
            stem = xyzf[:-4]                # complex_<id>_<label>
            ent_stem = stem + "_ent"
            # xyz
            with open(os.path.join(STRUCT, ent_stem + ".xyz"), "w") as f:
                f.write(f"{len(sym)}\n{comment} | enantiomer (mirror)\n")
                for s, p in zip(sym, M):
                    f.write(f"{s} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
            files = {"xyz": ent_stem + ".xyz"}
            # mol2 (mirror x of original)
            if iso["files"].get("mol2") and os.path.exists(os.path.join(STRUCT, iso["files"]["mol2"])):
                mirror_mol2(os.path.join(STRUCT, iso["files"]["mol2"]),
                            os.path.join(STRUCT, ent_stem + ".mol2"))
                files["mol2"] = ent_stem + ".mol2"
            # trexio
            try:
                TX.write_trexio(os.path.join(STRUCT, ent_stem + ".h5"), sym, M,
                                iso.get("total_charge"), (iso.get("mult", 1) - 1),
                                title=f"#{r['id']} enantiomer")
                files["trexio"] = ent_stem + ".h5"
            except Exception:
                pass
            # relabel original as Λ, add Δ enantiomer
            base = iso["label"]
            iso["label"] = (base + "·Λ") if base not in ("only",) else "Λ"
            ent = dict(iso)
            ent["label"] = (base + "·Δ") if base not in ("only",) else "Δ"
            ent["files"] = files
            ent["enantiomer"] = True
            ent["note"] = "enantiomer (mirror image; identical energy)"
            extra.append(ent)
            added += 1
        r["isomers"].extend(extra)
    man["n_isomers_total"] = sum(len(r.get("isomers", [])) for r in man["records"])
    man["n_valid_struct"] = sum(1 for r in man["records"] for i in r.get("isomers", []) if i.get("valid"))
    json.dump(man, open(os.path.join(OUT, "manifest.json"), "w"), indent=2)
    import collections
    d = collections.Counter(len(r.get("isomers", [])) for r in man["records"])
    print(f"added {added} enantiomers | structures now {man['n_isomers_total']} | isomers/complex {dict(d)}")


if __name__ == "__main__":
    main()
