"""
RDKit SELFIES-Metal Bridge
==========================
Unified API for converting between RDKit Mol, SMILES, and SELFIES
with support for transition metal atoms.

Usage:
    from rdkit_selfies_metal import mol_to_selfies, selfies_to_mol, smiles_to_selfies
    
    # RDKit Mol → SELFIES
    selfies_str = mol_to_selfies(mol)
    
    # SELFIES → RDKit Mol
    mol = selfies_to_mol(selfies_str)
    
    # MetalCytoSMILES → SELFIES
    selfies_str = smiles_to_selfies("[Ru+2].Cc1ccc(C(C)C)cc1.[Cl-].[Cl-]")
    
    # SELFIES → MetalCytoSMILES
    smiles_str = selfies_to_smiles(selfies_str)
"""

import sys
import os

# Add patches directory to path
sys.path.insert(0, os.path.dirname(__file__))

from selfies_metal import (
    smiles_to_selfies_metal,
    selfies_to_smiles_metal,
    mol_to_selfies_metal,
    selfies_to_mol_metal,
    is_metal_token,
    METAL_ATOMS,
)


def smiles_to_selfies(smiles_str):
    """Convert MetalCytoSMILES to SELFIES."""
    return smiles_to_selfies_metal(smiles_str)


def selfies_to_smiles(selfies_str):
    """Convert SELFIES to MetalCytoSMILES."""
    return selfies_to_smiles_metal(selfies_str)


def mol_to_selfies(mol):
    """Convert RDKit Mol to SELFIES."""
    return mol_to_selfies_metal(mol)


def selfies_to_mol(selfies_str):
    """Convert SELFIES to RDKit Mol."""
    return selfies_to_mol_metal(selfies_str)


# Re-export for convenience
__all__ = [
    'smiles_to_selfies',
    'selfies_to_smiles',
    'mol_to_selfies',
    'selfies_to_mol',
    'is_metal_token',
    'METAL_ATOMS',
]
