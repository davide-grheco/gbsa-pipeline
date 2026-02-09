import BioSimSpace as BSS
from src.gbsa_pipeline.parametrization import load_protein_pdb
from src.gbsa_pipeline.parametrization import parameterise_protein_amber
from src.gbsa_pipeline.solvation_box import run_solvation


def test_solvation() -> None:
    mol_1 = load_protein_pdb("tests/testdata/test1.pdb")
    mol_2 = parameterise_protein_amber(mol_1)
    run_solvation(system=mol_2, water_model="tip3p", box_size=8)

    assert not mol_1.nAtoms() == mol_2.nAtoms()

    from src.gbsa_pipeline.ligand_preparation import load_ligand_sdf
