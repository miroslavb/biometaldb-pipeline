"""Write a TREXIO ("T-Rex") record for a reconstructed complex.

TREXIO (TREX-CoE quantum-chemistry data format) record carrying the geometry
(nucleus coords in Bohr, atomic numbers, labels) and electron counts (up/down
from net charge + multiplicity) — a standard, QM-code-readable starting record.
HDF5 single-file back-end, with a text-directory fallback.
"""
from __future__ import annotations
import os
import numpy as np
import trexio
from ase.data import atomic_numbers

ANG2BOHR = 1.8897259886


def write_trexio(path, symbols, positions_ang, total_charge, n_unpaired,
                 title="", extra: dict | None = None):
    """Write a TREXIO file. Returns the path written."""
    Z = [int(atomic_numbers[s]) for s in symbols]
    nelec = int(sum(Z)) - int(total_charge or 0)
    nunp = int(n_unpaired or 0)
    if (nelec - nunp) % 2 != 0:          # keep parity consistent
        nunp = nelec % 2
    nup = (nelec + nunp) // 2
    ndn = nelec - nup
    coord_bohr = (np.asarray(positions_ang, dtype=float) * ANG2BOHR).tolist()

    if os.path.exists(path):
        os.remove(path)
    backend = trexio.TREXIO_HDF5
    try:
        tf = trexio.File(path, "w", back_end=backend)
    except Exception:
        backend = trexio.TREXIO_TEXT
        if os.path.exists(path):
            os.remove(path)
        tf = trexio.File(path, "w", back_end=backend)
    try:
        trexio.write_nucleus_num(tf, len(Z))
        trexio.write_nucleus_charge(tf, [float(z) for z in Z])
        trexio.write_nucleus_label(tf, list(symbols))
        trexio.write_nucleus_coord(tf, coord_bohr)
        trexio.write_electron_num(tf, nelec)
        trexio.write_electron_up_num(tf, nup)
        trexio.write_electron_dn_num(tf, ndn)
        try:
            trexio.write_metadata_code_num(tf, 1)
            trexio.write_metadata_code(tf, ["recon3d/Architector+GFN2-xTB"])
            if title:
                trexio.write_metadata_description(tf, title[:2040])
        except Exception:
            pass
    finally:
        tf.close()
    return path


def write_trexio_text(txt_path, symbols, positions_ang, total_charge, n_unpaired, title=""):
    """Write the TREXIO record in TEXT form as a single readable .trexio.txt.

    Uses TREXIO's native TEXT back-end (a directory of nucleus/electron/metadata
    text files), then consolidates them into one annotated text file (with the
    Cartesian coordinates also echoed in Angstrom for human reading).
    """
    import shutil
    Z = [int(atomic_numbers[s]) for s in symbols]
    nelec = int(sum(Z)) - int(total_charge or 0)
    nunp = int(n_unpaired or 0)
    if (nelec - nunp) % 2 != 0:
        nunp = nelec % 2
    nup = (nelec + nunp) // 2
    ndn = nelec - nup
    coord_bohr = (np.asarray(positions_ang, dtype=float) * ANG2BOHR)

    d = txt_path + ".trexdir"
    if os.path.exists(d):
        shutil.rmtree(d)
    tf = trexio.File(d, "w", back_end=trexio.TREXIO_TEXT)
    try:
        trexio.write_nucleus_num(tf, len(Z))
        trexio.write_nucleus_charge(tf, [float(z) for z in Z])
        trexio.write_nucleus_label(tf, list(symbols))
        trexio.write_nucleus_coord(tf, coord_bohr.tolist())
        trexio.write_electron_num(tf, nelec)
        trexio.write_electron_up_num(tf, nup)
        trexio.write_electron_dn_num(tf, ndn)
        try:
            trexio.write_metadata_code_num(tf, 1)
            trexio.write_metadata_code(tf, ["recon3d/Architector+GFN2-xTB"])
            if title:
                trexio.write_metadata_description(tf, title[:2040])
        except Exception:
            pass
    finally:
        tf.close()

    lines = [f"# TREXIO record (text) — {title}",
             f"# charge={int(total_charge or 0)} multiplicity={nunp + 1} "
             f"electrons={nelec} (up={nup} dn={ndn})  nuclei={len(Z)}",
             "# coordinates stored in Bohr (TREXIO standard); Angstrom echo below.", ""]
    for name in ("metadata.txt", "nucleus.txt", "electron.txt"):
        p = os.path.join(d, name)
        if os.path.exists(p):
            lines.append(f"### {name} ###")
            lines.append(open(p).read().rstrip("\n"))
            lines.append("")
    lines.append("### geometry (Angstrom, human-readable) ###")
    lines.append(f"{len(symbols)}")
    for s, p in zip(symbols, positions_ang):
        lines.append(f"{s:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}")
    open(txt_path, "w").write("\n".join(lines) + "\n")
    shutil.rmtree(d, ignore_errors=True)
    return txt_path


def trexio_summary(path):
    """Read back a few fields for display/verification."""
    backend = trexio.TREXIO_HDF5 if path.endswith((".h5", ".trexio")) else trexio.TREXIO_TEXT
    with trexio.File(path, "r", back_end=backend) as tf:
        return {
            "nucleus_num": trexio.read_nucleus_num(tf),
            "electron_num": trexio.read_electron_num(tf),
            "electron_up_num": trexio.read_electron_up_num(tf),
            "electron_dn_num": trexio.read_electron_dn_num(tf),
        }
