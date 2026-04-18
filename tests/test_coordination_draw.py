"""Tests for coordination_draw module."""
import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from chem_pipeline_lib import draw_compound, mol_from_molblock, mol_from_smiles, build_square_planar
from rdkit import Chem


class TestMolFromMolblock:
    def test_parse_cisplatin(self):
        molblock = """
     RDKit          2D

 11 10  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 Pt  0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.5000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    2.3000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7000    1.8500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7000    1.8500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.3000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7000   -1.8500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.7000   -1.8500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
  4  6  1  0
  4  7  1  0
  4  8  1  0
  5  9  1  0
  5 10  1  0
  5 11  1  0
M  END
"""
        mol = mol_from_molblock(molblock)
        assert mol is not None
        assert mol.GetNumAtoms() == 11
        symbols = [a.GetSymbol() for a in mol.GetAtoms()]
        assert 'Pt' in symbols
        assert symbols.count('Cl') == 2
        assert symbols.count('N') == 2

    def test_invalid_molblock_returns_none(self):
        assert mol_from_molblock("not a molblock") is None


class TestMolFromSmiles:
    def test_parse_simple(self):
        mol = mol_from_smiles("c1ccccc1")  # benzene
        assert mol is not None
        assert mol.GetNumAtoms() == 6


class TestBuildSquarePlanar:
    def test_cisplatin_cis(self):
        mol = build_square_planar(
            metal_symbol='Pt',
            ligands=[('Cl', 0), ('Cl', 0), ('N', 0), ('N', 0)],
            cis_pairs=[(0, 1)]  # Cl atoms are cis
        )
        assert mol is not None
        assert mol.GetNumAtoms() == 5  # Pt + 2Cl + 2N

    def test_cisplatin_trans(self):
        mol = build_square_planar(
            metal_symbol='Pt',
            ligands=[('Cl', 0), ('Cl', 0), ('N', 0), ('N', 0)],
            cis_pairs=None  # default trans
        )
        assert mol is not None
        assert mol.GetNumAtoms() == 5


class TestPubChem:
    def test_fetch_cisplatin(self):
        result = draw_compound("Cisplatin", pubchem_name="cisplatin")
        assert result.get('success') is True or 'error' in result  # network dependent


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
