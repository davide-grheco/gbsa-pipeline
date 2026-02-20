from __future__ import annotations

import BioSimSpace as BSS

from gbsa_pipeline.parametrization import load_protein_pdb, parameterise_protein_amber
from gbsa_pipeline.solvation_box import run_solvation


def test_solvation() -> None:
    # Load raw protein
    mol = load_protein_pdb("tests/testdata/test1.pdb")

    # Parameterise -> returns Molecule
    mol_param = parameterise_protein_amber(mol, ff="ff14SB")

    # Create System explicitly
    system = BSS._SireWrappers.System(mol_param)
    n_atoms_before = system.nAtoms()

    # Solvate (assumed in-place)
    solvated = run_solvation(system=system, water_model="tip3p", box_size=8)

    # Check that solvent was added
    assert solvated.nAtoms() > n_atoms_before
