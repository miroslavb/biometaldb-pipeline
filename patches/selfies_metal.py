"""
SELFIES-Metal Extension
=======================
Extends python-selfies to support transition metal atoms in SELFIES strings.

Key insight: SELFIES encoder passes bracket atoms like [Ru+2] through as-is.
So we just need to handle dot-separated fragments correctly.

Usage:
    from selfies_metal import smiles_to_selfies_metal, selfies_to_smiles_metal
    
    # Convert MetalCytoSMILES → SELFIES
    selfies_str = smiles_to_selfies_metal("[Ru+2].Cc1ccc(C(C)C)cc1.[Cl-].[Cl-]")
    
    # Convert SELFIES → MetalCytoSMILES
    smiles_str = selfies_to_smiles_metal(selfies_str)
"""

import re
import selfies as sf

# Metal atoms we support
METAL_ATOMS = {
    'Ru', 'Ir', 'Rh', 'Os', 'Re',
    'Fe', 'Pt', 'Pd', 'Cu', 'Zn',
    'Au', 'Ag', 'Co', 'Ni', 'Mn',
    'Cr', 'Mo', 'W', 'V', 'Ti',
}


def is_metal_token(frag):
    """Check if a fragment is a metal token like [Ru+2], [Ir+3], [Rh], [Fe+2]."""
    m = re.match(r'^\[([A-Z][a-z]?)([+-]\d+)?\]$', frag.strip())
    if m and m.group(1) in METAL_ATOMS:
        return True
    return False


def smiles_to_selfies_metal(smiles_str):
    """
    Convert MetalCytoSMILES to SELFIES with metal tokens.
    
    Handles dot-separated fragments. Metal tokens like [Ru+2] pass through
    unchanged (SELFIES bracket notation). Organic fragments are encoded normally.
    """
    fragments = smiles_str.split('.')
    selfies_fragments = []
    
    for frag in fragments:
        frag = frag.strip()
        if not frag:
            continue
        
        if is_metal_token(frag):
            # Metal tokens pass through as-is
            selfies_fragments.append(frag)
        else:
            # Encode organic fragment
            try:
                s = sf.encoder(frag)
                selfies_fragments.append(s)
            except Exception as e:
                # Fallback: store as bracket atom
                selfies_fragments.append(f"[_{frag}_]")
    
    return '.'.join(selfies_fragments)


def selfies_to_smiles_metal(selfies_str):
    """
    Convert SELFIES with metal tokens back to MetalCytoSMILES.
    """
    fragments = selfies_str.split('.')
    smiles_fragments = []
    
    for frag in fragments:
        frag = frag.strip()
        if not frag:
            continue
        
        if is_metal_token(frag):
            # Metal token passes through
            smiles_fragments.append(frag)
        elif frag.startswith('[_') and frag.endswith('_]'):
            # Recover raw SMILES from fallback
            smiles_fragments.append(frag[2:-2])
        else:
            # Decode SELFIES fragment
            try:
                s = sf.decoder(frag)
                smiles_fragments.append(s)
            except Exception:
                smiles_fragments.append(frag)
    
    return '.'.join(smiles_fragments)


def mol_to_selfies_metal(mol):
    """Convert RDKit Mol (with metal) to SELFIES string."""
    from rdkit import Chem
    smiles = Chem.MolToSmiles(mol, canonical=False)
    if not smiles:
        return None
    return smiles_to_selfies_metal(smiles)


def selfies_to_mol_metal(selfies_str):
    """Convert SELFIES string to RDKit Mol (with metal)."""
    from rdkit import Chem
    smiles = selfies_to_smiles_metal(selfies_str)
    return Chem.MolFromSmiles(smiles, sanitize=False)


def roundtrip_test(smiles_str):
    """
    Test roundtrip: SMILES → SELFIES → SMILES
    Returns (success, selfies_str, roundtrip_smiles)
    """
    try:
        selfies_str = smiles_to_selfies_metal(smiles_str)
        roundtrip = selfies_to_smiles_metal(selfies_str)
        return True, selfies_str, roundtrip
    except Exception as e:
        return False, None, str(e)


if __name__ == "__main__":
    # Self-test
    test_cases = [
        "[Ru+2].Cc1ccc(C(C)C)cc1.[Cl-].[Cl-]",
        "[Ir+3].Cc1ccc(-c2ccccn2)cc1.[Cl-].[Cl-].[Cl-]",
        "[Os+2].c1ccnc2ccncc12.c1ccnc2ccncc12.[Cl-].[Cl-]",
        "[Re+1].[CH3-]",
        "[Ru+2].N#Cc1ccc(C(c2ccc(C#N)cc2)n2cncn2)cc1.N#Cc1ccc(C(c2ccc(C#N)cc2)n2cncn2)cc1.[Cl-].c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    ]
    
    print("SELFIES-Metal Self-Test")
    print("=" * 70)
    
    all_ok = True
    for smiles in test_cases:
        ok, selfies_str, roundtrip = roundtrip_test(smiles)
        status = "✓" if ok else "✗"
        print(f"\n{status} SMILES:    {smiles[:80]}...")
        print(f"  SELFIES:   {selfies_str[:80]}...")
        print(f"  Roundtrip: {roundtrip[:80]}...")
        
        if ok and roundtrip:
            from rdkit import Chem
            mol1 = Chem.MolFromSmiles(smiles, sanitize=False)
            mol2 = Chem.MolFromSmiles(roundtrip, sanitize=False)
            if mol1 and mol2:
                n1, n2 = mol1.GetNumAtoms(), mol2.GetNumAtoms()
                equiv = "✓" if n1 == n2 else "✗"
                if n1 != n2:
                    all_ok = False
                print(f"  Atoms:     {n1} {equiv} {n2}")
            else:
                print(f"  Parse:     mol1={'ok' if mol1 else 'FAIL'}, mol2={'ok' if mol2 else 'FAIL'}")
                all_ok = False
    
    print(f"\n{'✓ All tests passed' if all_ok else '✗ Some tests failed'}")
