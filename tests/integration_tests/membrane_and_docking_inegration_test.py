"""Integration test: full [membrane] docking -> GBSA-prep chain."""

from __future__ import annotations

import shutil
from pathlib import Path

import BioSimSpace as BSS
import pytest

from gbsa_pipeline.docking import (
    DockingBox,
    DockingRequest,
    VinaEngine,
    convert_receptor_pdb_to_pdbqt,
    export_pdbqt_to_sdf,
    load_first_sdf_molecule,
    prepare_ligand_with_meeko,
)
from gbsa_pipeline.membrane import extract_receptor_pdb
from gbsa_pipeline.parametrization import parameterise_ligand_gaff2
from gbsa_pipeline.solvation_box import SolvationParams, WaterModel, solvate_membrane

TESTDATA = Path(__file__).parent.parent / "testdata" / "membrane" / "2rh1"

# centered on the co-crystallised ligand.sdf; centroid, max span + 10 A padding
BOX = DockingBox(center=(60.7, 52.5, 74.3), size=(17.8, 17.8, 17.8))


@pytest.mark.integration
def test_membrane_docking_to_gbsa_prep_chain(tmp_path: Path) -> None:
    """Dock into a protein-only 2rh1 receptor, then prep for GBSA chain in the full membrane system."""
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    dock_dir = tmp_path / "docking"
    dock_dir.mkdir()

    receptor_pdb = extract_receptor_pdb(TESTDATA / "system_unsolvated.gro", dock_dir / "receptor.pdb")
    receptor_pdbqt = convert_receptor_pdb_to_pdbqt(receptor_pdb)

    ligand_sdf = TESTDATA / "ligand.sdf"
    ligand_mol = load_first_sdf_molecule(ligand_sdf, remove_hs=False)
    ligand_pdbqt = dock_dir / "ligand.pdbqt"
    prepare_ligand_with_meeko(ligand_mol, ligand_pdbqt, name="LIG")

    engine = VinaEngine()
    request = DockingRequest(
        receptor=receptor_pdbqt,
        ligands=[ligand_pdbqt],
        box=BOX,
        workdir=dock_dir,
        parameters={"exhaustiveness": 4, "num_modes": 3},
    )
    result = engine.dock(request)
    assert result.poses

    docked_sdf = dock_dir / "docked_ligand.sdf"
    export_pdbqt_to_sdf(
        result.poses[0].pose_path,
        docked_sdf,
        template_mol=ligand_mol,
        add_hydrogens_after_template=True,
    )
    assert docked_sdf.exists()

    system = BSS.IO.readMolecules(
        [
            str(TESTDATA / "system_unsolvated.gro"),
            str(TESTDATA / "system_unsolvated.top"),
        ],
        make_whole=True,
    )
    n_atoms_before = system.nAtoms()

    ligand_mol = BSS.IO.readMolecules(str(docked_sdf)).getMolecules()[0]
    ligand = parameterise_ligand_gaff2(ligand_mol, net_charge=0, work_dir=tmp_path)
    system.addMolecules(ligand)
    assert system.nAtoms() > n_atoms_before

    solvated = solvate_membrane(
        system=system,
        params=SolvationParams(water_model=WaterModel.TIP3P, ion_concentration=0.15, neutralize=True),
        z_padding_nm=1.5,
        work_dir=tmp_path,
    )
    assert solvated.getWaterMolecules().nMolecules() > 0
