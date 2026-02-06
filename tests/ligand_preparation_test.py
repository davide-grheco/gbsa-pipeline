from __future__ import annotations
from src.gbsa_pipeline.ligand_preparation import load_ligand_sdf, ligand_standardizer


def test_hydrogens_added() -> None:
    """Test standardization adds correct number of hydrogens"""
    mol = load_ligand_sdf("data/complex3/complex3.sdf")
    mol = ligand_standardizer(mol)

    n_atoms = mol.GetNumAtoms()
    assert n_atoms == 72
