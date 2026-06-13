"""Donor-site assignment: DB row -> per-ligand coordinating atoms + geometry.

Given metal, oxidation state, the dot-separated ligand SMILES, and the
(imperfect) aggregate donor_atoms composition, decide which atom of which
ligand coordinates the metal, the denticity of each ligand, the coordination
number / polyhedron, and a confidence flag. This is the inference that the
original pipeline never did (it dropped a disconnected metal into a 2D layout).
"""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass, field

from rdkit import Chem

import priors


@dataclass
class Unit:
    smiles: str
    mol: "Chem.Mol"
    donors: list                      # [(idx, element, strength)]
    pocket: list = field(default_factory=list)   # selected coordinating atom idxs
    pocket_comp: Counter = field(default_factory=Counter)
    role: str = "ligand"              # ligand | halide | counterion | haptic
    charge: int = 0
    ligType: str = None               # mono | bi_cis | ... (None -> Architector auto)
    sites: int = 0                    # coordination sites occupied (haptic -> 3)
    hkind: str = None                 # 'cp' | 'arene' for haptic ligands
    oracle: dict = None               # pydentate prediction (if used)


@dataclass
class Assignment:
    metal: str
    ox: int
    target_comp: Counter
    target_cn: int
    pref_cn: int
    geometry: str
    spin: int
    units: list                       # selected coordinating Units (with pocket set)
    final_comp: Counter
    final_cn: int
    confidence: str
    flags: list
    spectators: list                  # counterion/unbound fragment smiles
    needs_xtb: bool = False           # haptic present -> UFF insufficient, force xtb
    hemilabile: bool = False          # any coordinating ligand predicted hemilabile


def _mol(smi):
    m = Chem.MolFromSmiles(smi, sanitize=True)
    if m is None:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is not None:
            try:
                Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
                                 ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                                 ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
            except Exception:
                pass
    return m


def _bond_path_len(mol, a, b):
    try:
        path = Chem.GetShortestPath(mol, a, b)
        return len(path) - 1 if path else 999
    except Exception:
        return 999


def _natural_pocket(mol, donors):
    """Chelate-aware pocket selection.

    Chelate edge = two donors whose bond path is 2..3 (i.e. a 5- or 6-membered
    metallacycle, by far the most common). Donors are grouped into components
    over those edges; the component with the highest total donor strength wins
    (so a strong N^N beats a distal weak O...O on the same fused ring). Finally,
    if the pocket contains a strong donor (P/S/strong-C, >=85) we drop any much
    weaker tag-along (ether/carbonyl/nitrile, <=45) — e.g. a furan O next to a
    phosphine P — without disturbing genuine weak-donor chelates like acac
    (whose O's are themselves moderate-strength and 1,3-related)."""
    if not donors:
        return []
    idxs = [d[0] for d in donors]
    strength = {d[0]: d[2] for d in donors}
    parent = {i: i for i in idxs}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(x, y):
        parent[find(x)] = find(y)
    for i in range(len(idxs)):
        for j in range(i + 1, len(idxs)):
            if 2 <= _bond_path_len(mol, idxs[i], idxs[j]) <= 4:
                union(idxs[i], idxs[j])
    comps = {}
    for i in idxs:
        comps.setdefault(find(i), []).append(i)
    best = max(comps.values(), key=lambda c: sum(strength[i] for i in c))
    if best and max(strength[i] for i in best) >= 85:
        strong = [i for i in best if strength[i] > 45]
        if strong:
            best = strong
    return sorted(best)


def _frag_charge(mol):
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def build_units(metal, ox, smiles_ligands, oracle=None):
    units = []
    oracle = oracle or {}
    for frag in smiles_ligands.split("."):
        frag = frag.strip()
        if not frag:
            continue
        mol0 = _mol(frag)
        if mol0 is None:
            units.append(Unit(frag, None, [], role="unparsed"))
            continue
        canon = Chem.MolToSmiles(mol0)
        # Re-parse canonical SMILES so atom indices == SMILES string order,
        # which is the order Architector/OpenBabel assigns when given `canon`.
        mol = _mol(canon) or mol0
        if _mol(canon) is None:
            canon = frag
        donors = priors.find_candidate_donors(mol)
        u = Unit(canon, mol, donors, charge=_frag_charge(mol))
        hap = priors.detect_haptic(mol)
        if canon in priors.COUNTERIONS:
            u.role = "counterion"
        elif canon in priors.HALIDES_CANON:
            u.role = "halide"
        elif hap:
            kind, ring, nsites = hap
            u.role = "haptic"; u.hkind = kind
            u.pocket = ring
            u.pocket_comp = Counter({"C": len(ring)})
            u.sites = nsites
        elif canon in oracle and isinstance(oracle[canon], dict) \
                and oracle[canon].get("coord_atoms") and "error" not in oracle[canon]:
            # pydentate-predicted coordinating atoms (authoritative for hetero donors)
            pred = oracle[canon]
            valid = [i for i in pred["coord_atoms"] if i < mol.GetNumAtoms()]
            if valid:
                u.oracle = pred
                u.pocket = sorted(valid)
                u.pocket_comp = Counter(mol.GetAtomWithIdx(i).GetSymbol() for i in u.pocket)
            elif not donors:
                u.role = "counterion"
        elif not donors:
            u.role = "counterion"   # no donor atoms -> outer sphere
        units.append(u)
    return units


def assign(metal, ox, smiles_ligands, donor_atoms, oracle=None):
    """Return an Assignment.

    donor_atoms: dict like {'N':5,'P':1} (may be None) — composition cross-check.
    oracle: optional {canonical_smiles: pydentate_prediction} — authoritative
    coordinating-atom assignment for heteroatom-donor ligands.
    """
    flags = []
    metal = metal.strip()
    ox = int(ox) if ox is not None else None
    pref_cn = priors.preferred_cn(metal, ox)
    spin = priors.low_spin_unpaired(metal, ox)

    target_comp = Counter({k: int(v) for k, v in (donor_atoms or {}).items() if v})
    cn_from_donor = sum(target_comp.values())
    # CN driven by metal+oxidation state (reliable); donor_atoms is only a
    # composition hint (it came from a fragile regex and often miscounts).
    target_cn = pref_cn if pref_cn else (cn_from_donor or 6)
    donor_prior_sane = bool(target_comp) and (cn_from_donor == target_cn)
    if not target_comp:
        flags.append("no_donor_atoms_prior")
    elif not donor_prior_sane:
        flags.append(f"donor_prior_cn{cn_from_donor}_ne_pref{target_cn}")

    units = build_units(metal, ox, smiles_ligands, oracle=oracle)

    # natural pocket (SMARTS heuristic) only where the oracle didn't assign one
    for u in units:
        if u.role in ("ligand", "halide") and u.donors and not u.oracle and not u.pocket:
            u.pocket = _natural_pocket(u.mol, u.donors)
            el = {d[0]: d[1] for d in u.donors}
            u.pocket_comp = Counter(el[i] for i in u.pocket)

    haptics = [u for u in units if u.role == "haptic"]
    ligands = [u for u in units if u.role == "ligand" and u.pocket]
    halides = [u for u in units if u.role == "halide" and u.pocket]
    spectators = [u.smiles for u in units if u.role in ("counterion", "unparsed")]

    def smap(u):  # atom idx -> (element, strength)
        return {d[0]: (d[1], d[2]) for d in u.donors}

    # order ligands: chelators first, then strongest donor, then larger pocket
    def lig_key(u):
        s = smap(u)
        maxs = max((s.get(i, ("", 50))[1] for i in u.pocket), default=0)
        return (-(len(u.pocket) >= 2), -maxs, -len(u.pocket))
    ligands.sort(key=lig_key)

    # Fill EXACTLY target_cn: take whole pockets that fit; trim the overshooting
    # one to its strongest donors; top up with halides (monodentate).
    # geometric denticity cap: linear/CN2 metals (Au(I)) cannot chelate -> each
    # ligand binds monodentate via its strongest donor.
    dcap = 1 if target_cn <= 2 else 999
    selected, slots = [], target_cn
    # haptic ligands first (definite; each η5-Cp / η6-arene occupies 3 facial sites)
    for u in haptics:
        if slots < u.sites:
            flags.append("haptic_no_room")
            continue
        selected.append(u)
        slots -= u.sites
    for u in ligands:
        if slots <= 0:
            break
        s = smap(u)
        pk = sorted(u.pocket, key=lambda i: -s.get(i, ("", 50))[1])
        take = pk[:min(len(pk), slots, dcap)]
        u.pocket = sorted(take)
        u.pocket_comp = Counter(u.mol.GetAtomWithIdx(i).GetSymbol() for i in take)
        u.sites = len(take)
        u.ligType = {1: "mono", 2: "bi_cis"}.get(len(take))
        selected.append(u)
        slots -= len(take)
    for u in sorted(halides, key=lambda x: -max((d[2] for d in x.donors), default=0)):
        if slots <= 0:
            break
        u.pocket = [u.donors[0][0]]
        u.pocket_comp = Counter([u.donors[0][1]])
        u.sites = 1
        u.ligType = "mono"
        selected.append(u)
        slots -= 1
    if slots > 0:
        flags.append(f"cn_unfilled_by_{slots}")
    unused = len(ligands) - len([u for u in selected if u.role == "ligand"])
    if unused > 0:
        flags.append(f"unused_ligands_{unused}")

    final_units = [u for u in selected if u.pocket]
    final_comp = Counter()
    for u in final_units:
        final_comp += u.pocket_comp
    final_cn = sum(u.sites for u in final_units)   # site count (haptic ring = 3)
    needs_xtb = bool(haptics)
    hemilabile = any(u.oracle and u.oracle.get("hemilabile") for u in final_units)

    comp_match = bool(target_comp) and (final_comp == target_comp)
    hard_flag = any(f.startswith(("cn_unfilled", "unused", "haptic_no_room")) for f in flags)
    if final_cn == target_cn and comp_match and not hard_flag:
        conf = "high"
    elif final_cn == target_cn and not hard_flag:
        conf = "medium"
        if target_comp and not comp_match:
            flags.append(f"comp_got{dict(final_comp)}_want{dict(target_comp)}")
    elif final_cn == target_cn:
        conf = "medium"
    else:
        conf = "low"

    geometry = priors.geometry_for(metal, ox, final_cn if final_cn > 0 else pref_cn)

    return Assignment(
        metal=metal, ox=ox, target_comp=target_comp, target_cn=target_cn,
        pref_cn=pref_cn, geometry=geometry, spin=spin,
        units=final_units, final_comp=final_comp, final_cn=final_cn,
        confidence=conf, flags=flags, spectators=spectators, needs_xtb=needs_xtb,
        hemilabile=hemilabile,
    )
