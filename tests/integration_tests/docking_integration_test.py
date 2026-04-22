"""Integration tests for the lightweight docking adapter layer."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolAlign

from gbsa_pipeline.docking import (
    DockingBox,
    DockingRequest,
    VinaEngine,
    export_pdbqt_to_sdf,
    prepare_ligand_with_meeko,
)

TESTDATA = Path(__file__).parents[1] / "testdata"
DOCKING_TESTDATA = TESTDATA / "docking"
DOCKLIGAND_SDF = DOCKING_TESTDATA / "dockligand.sdf"
DOCKPROTEIN_PDB = DOCKING_TESTDATA / "dockprotein.pdb"

ROUNDTRIP_PDBQT = DOCKING_TESTDATA / "dockligand_roundtrip.pdbqt"
ROUNDTRIP_SDF = DOCKING_TESTDATA / "dockligand_roundtrip_exported.sdf"

DOCKPROTEIN_BOX = DockingBox(
    center=(-3.245, 29.915, 53.639),
    size=(10.0, 10.0, 10.0),
)


def _read_first_sdf_molecule(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecule = supplier[0]

    if molecule is None:
        raise ValueError(f"Could not read first molecule from SDF: {path}")

    return molecule


@pytest.mark.integration
def test_pdbqt_to_sdf_roundtrip_preserves_pose() -> None:
    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    input_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)

    prepare_ligand_with_meeko(input_molecule, ROUNDTRIP_PDBQT, name="DOCKLIG")
    export_pdbqt_to_sdf(ROUNDTRIP_PDBQT, ROUNDTRIP_SDF)

    assert ROUNDTRIP_PDBQT.exists()
    assert ROUNDTRIP_SDF.exists()

    exported_molecule = _read_first_sdf_molecule(ROUNDTRIP_SDF)

    input_heavy = Chem.RemoveHs(Chem.Mol(input_molecule))
    exported_heavy = Chem.RemoveHs(Chem.Mol(exported_molecule))

    assert input_heavy.GetNumAtoms() > 0
    assert exported_heavy.GetNumAtoms() > 0
    assert input_heavy.GetNumAtoms() == exported_heavy.GetNumAtoms()

    heavy_atom_rmsd = rdMolAlign.CalcRMS(exported_heavy, input_heavy)

    assert heavy_atom_rmsd < 1.0


@pytest.mark.integration
def test_vina_binary_smoke(tmp_path: Path) -> None:
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    if not DOCKPROTEIN_PDB.exists():
        pytest.skip(f"missing receptor test file: {DOCKPROTEIN_PDB}")

    engine = VinaEngine(binary="vina")

    ligand = tmp_path / "dockligand.pdbqt"
    prepare_ligand_with_meeko(
        _read_first_sdf_molecule(DOCKLIGAND_SDF),
        ligand,
        name="DOCKLIG",
    )

    request = DockingRequest(
        receptor=DOCKPROTEIN_PDB,
        ligands=[ligand],
        box=DOCKPROTEIN_BOX,
        workdir=tmp_path,
    )

    result = engine.dock(request=request)

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0
    assert result.poses[0].score is not None
