"""Chemistry priors for coordination-complex 3D reconstruction.

Maps (metal, oxidation_state) -> preferred coordination number, polyhedron,
and low-spin unpaired-electron count, plus counterion recognition and donor-atom
SMARTS. Tuned for the BiometalDB antiproliferative set (Ru, Au, Ir, Rh, Os, Re;
4d/5d, overwhelmingly low-spin).
"""
from __future__ import annotations
from rdkit import Chem

# Group number (number of valence d+s electrons of the neutral atom block).
METAL_GROUP = {
    "Sc": 3, "Ti": 4, "V": 5, "Cr": 6, "Mn": 7, "Fe": 8, "Co": 9, "Ni": 10,
    "Cu": 11, "Zn": 12,
    "Y": 3, "Zr": 4, "Nb": 5, "Mo": 6, "Tc": 7, "Ru": 8, "Rh": 9, "Pd": 10,
    "Ag": 11, "Cd": 12,
    "Hf": 4, "Ta": 5, "W": 6, "Re": 7, "Os": 8, "Ir": 9, "Pt": 10, "Au": 11,
    "Hg": 12,
}


def d_count(metal: str, ox: int) -> int:
    """d-electron count = group - oxidation state (clamped 0..10)."""
    g = METAL_GROUP.get(metal)
    if g is None or ox is None:
        return -1
    return max(0, min(10, g - int(ox)))


def low_spin_unpaired(metal: str, ox: int) -> int:
    """Unpaired electrons assuming low spin (valid for 4d/5d strong-field).

    Low-spin parity rule: even d-count -> 0 unpaired, odd -> 1.
    (Square-planar d8 and octahedral d6 -> 0; d5/d7 low-spin -> 1.)
    """
    d = d_count(metal, ox)
    if d < 0:
        return 0
    return d % 2


# Preferred coordination number by (metal, ox). Falls back to CN_BY_DCOUNT.
PREFERRED_CN = {
    ("Au", 1): 2, ("Au", 3): 4, ("Au", 2): 4,
    ("Ru", 2): 6, ("Ru", 3): 6, ("Ru", 4): 6,
    ("Os", 2): 6, ("Os", 3): 6, ("Os", 4): 6, ("Os", 6): 6,
    ("Ir", 3): 6, ("Ir", 1): 4, ("Ir", 4): 6,
    ("Rh", 3): 6, ("Rh", 1): 4, ("Rh", 2): 6,
    ("Re", 1): 6, ("Re", 2): 6, ("Re", 3): 6, ("Re", 4): 6, ("Re", 5): 6, ("Re", 7): 4,
    ("Pt", 2): 4, ("Pt", 4): 6, ("Pd", 2): 4,
    ("Fe", 2): 6, ("Fe", 3): 6, ("Co", 3): 6, ("Co", 2): 6, ("Cu", 2): 4,
}

# Fallback CN purely by d-count when (metal,ox) not tabulated.
CN_BY_DCOUNT = {0: 6, 1: 6, 2: 6, 3: 6, 4: 6, 5: 6, 6: 6, 7: 5, 8: 4, 9: 4, 10: 2}


def preferred_cn(metal: str, ox: int) -> int:
    if (metal, ox) in PREFERRED_CN:
        return PREFERRED_CN[(metal, ox)]
    return CN_BY_DCOUNT.get(d_count(metal, ox), 6)


def geometry_for(metal: str, ox: int, cn: int) -> str:
    """Architector coreType for a given CN, disambiguating CN4 by d-count."""
    if cn == 2:
        return "linear"
    if cn == 3:
        return "trigonal_planar"
    if cn == 4:
        # d8 (Au III, Pt II, Pd II, Ir I, Rh I, Ni II) -> square planar; else tetrahedral
        return "square_planar" if d_count(metal, ox) == 8 else "tetrahedral"
    if cn == 5:
        return "square_pyramidal"
    if cn == 6:
        return "octahedral"
    if cn == 7:
        return "pentagonal_bipyramidal"
    if cn == 8:
        return "square_antiprismatic"
    return "octahedral"


# --- Counterion / spectator recognition (matched on canonical SMILES) ---
_COUNTERION_SMILES_RAW = [
    "F[B-](F)(F)F",                       # BF4-
    "F[P-](F)(F)(F)(F)F",                 # PF6-
    "F[As-](F)(F)(F)(F)F",                # AsF6-
    "F[Sb-](F)(F)(F)(F)F",                # SbF6-
    "[O-]Cl(=O)(=O)=O",                   # ClO4-
    "[O-][Cl+3]([O-])([O-])[O-]",         # ClO4- alt
    "O=S(=O)([O-])C(F)(F)F",              # triflate OTf-
    "O=S(=O)([O-])[O-]",                  # SO4 2-
    "[O-]S(=O)(=O)O",                     # HSO4-
    "O=[N+]([O-])[O-]",                   # NO3-
    "[O-][N+](=O)[O-]",                   # NO3- alt
    "[O-]C(=O)C(F)(F)F",                  # trifluoroacetate (often counterion)
    "FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F",  # NTf2-
    "[B-](c1ccccc1)(c1ccccc1)(c1ccccc1)c1ccccc1",  # BPh4-
    "[Na+]", "[K+]", "[Li+]", "[Cs+]", "[Rb+]",
    "[NH4+]", "[OH-]", "O",               # water (outer sphere unless donor budget needs it)
    "[PF6-]", "[BF4-]",
]
_HALIDES = {"[Cl-]", "[Br-]", "[I-]", "[F-]"}


def _canon(smi: str) -> str | None:
    try:
        m = Chem.MolFromSmiles(smi, sanitize=True)
        if m is None:
            m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None


COUNTERIONS = set(filter(None, (_canon(s) for s in _COUNTERION_SMILES_RAW)))
HALIDES_CANON = set(filter(None, (_canon(s) for s in _HALIDES)))


# --- Donor-atom SMARTS. (pattern, element, strength). Higher = binds preferentially. ---
# strength only breaks ties when more candidate donors than CN slots exist.
DONOR_SMARTS = [
    ("[c-]",                              "C", 100),  # cyclometalated aromatic carbanion (C^N)
    ("[C-]",                              "C", 100),  # carbanion / carbene-like
    ("[#6;X2;$([#6](=[#7])[#7]),$([#6]([#7])=[#7])]", "C", 98),  # NHC carbene C between 2 N
    ("[P;X3;!$([P]=O);!$([P+])]",         "P", 92),   # phosphine
    ("[As;X3]",                           "As", 88),
    ("[S-]",                              "S", 86),   # thiolate
    ("[#7;X3;-]",                         "N", 82),   # amide / deprotonated N-
    ("[O-]",                              "O", 78),   # carboxylate/phenolate/alkoxide O-
    ("[#7;a;X2;!$([n+]);!$([nH])]",       "N", 72),   # pyridine-type aromatic N
    ("[#7;X2]=[#6]",                      "N", 70),   # imine / azomethine N
    ("[Se]",                              "Se", 68),
    ("[#16;X2;!$([#16]=O)]",              "S", 64),   # thioether
    ("[#7;X3;!$([#7][#6]=[O,S]);!$([n]);!+;!$([#7][#7])]", "N", 60),  # amine (not amide)
    ("[Cl-]",                             "Cl", 52),
    ("[Br-]",                             "Br", 50),
    ("[I-]",                              "I", 48),
    ("[#8;X1]=[#6]",                      "O", 42),   # carbonyl O
    ("[OX2H]",                            "O", 40),   # hydroxyl / phenol / water
    ("[#8;X2;!$([#8][#6]=O);!$([O-])]",   "O", 38),   # ether / alkoxy
    ("[F-]",                              "F", 36),
    ("[#7;X1]#[#6]",                      "N", 34),   # nitrile N
]
_COMPILED = [(Chem.MolFromSmarts(p), el, s) for p, el, s in DONOR_SMARTS]


def detect_haptic(mol: "Chem.Mol"):
    """Detect a haptic (η5-Cp/Cp* or η6-arene) ligand.

    Returns (kind, sorted_ring_atom_idxs, n_sites) or None. A haptic ring
    occupies 3 facial coordination sites (piano-stool). All ring carbons are
    passed to Architector as the coordList; xtb relaxation then pulls the ring
    into η5/η6 (UFF alone leaves it floating).

    Guards: η6-arene only when the fragment carries NO heteroatom donor
    (so PPh3 / backbone phenyls are never mistaken for an arene ligand);
    pyridine/fused-N rings are not all-carbon so never match.
    """
    has_hetero_donor = any(a.GetSymbol() in ("N", "O", "S", "P", "Se", "As")
                           for a in mol.GetAtoms())
    for ring in mol.GetRingInfo().AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(a.GetSymbol() == "C" for a in atoms):
            continue
        aromatic = all(a.GetIsAromatic() for a in atoms)
        anionic = any(a.GetFormalCharge() < 0 for a in atoms)
        if len(ring) == 5 and (aromatic or anionic):
            return ("cp", sorted(ring), 3)
        if len(ring) == 6 and aromatic and not has_hetero_donor:
            return ("arene", sorted(ring), 3)
    return None


def find_candidate_donors(mol: "Chem.Mol") -> list[tuple[int, str, int]]:
    """Return [(atom_idx, element, strength)] of plausible donor atoms.

    Each atom keeps its highest-strength match only.
    """
    best: dict[int, tuple[str, int]] = {}
    for patt, el, strength in _COMPILED:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            idx = match[0]
            if idx not in best or strength > best[idx][1]:
                best[idx] = (el, strength)
    return [(i, el, s) for i, (el, s) in best.items()]
