"""Integration tests for the lightweight docking adapter layer."""

from __future__ import annotations

import math
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

DOCKPROTEIN_BOX = DockingBox(
    center=(10.115, 39.148, 53.112),
    size=(10.0, 10.0, 10.0),
)


def _read_first_sdf_molecule(path: Path) -> Chem.Mol:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecule = supplier[0]

    if molecule is None:
        raise ValueError(f"Could not read first molecule from SDF: {path}")

    return molecule


def _centroid(molecule: Chem.Mol) -> tuple[float, float, float]:
    if molecule.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformer.")

    conformer = molecule.GetConformer()
    atom_count = molecule.GetNumAtoms()

    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0

    for atom_index in range(atom_count):
        position = conformer.GetAtomPosition(atom_index)
        x_sum += float(position.x)
        y_sum += float(position.y)
        z_sum += float(position.z)

    return (
        x_sum / atom_count,
        y_sum / atom_count,
        z_sum / atom_count,
    )


def _distance(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
) -> float:
    return math.sqrt((left[0] - right[0]) ** 2 + (left[1] - right[1]) ** 2 + (left[2] - right[2]) ** 2)


def _heavy_molecule(molecule: Chem.Mol) -> Chem.Mol:
    return Chem.RemoveHs(Chem.Mol(molecule))


def _bond_order_sum(molecule: Chem.Mol) -> float:
    heavy = _heavy_molecule(molecule)
    return sum(bond.GetBondTypeAsDouble() for bond in heavy.GetBonds())


def _aromatic_bond_count(molecule: Chem.Mol) -> int:
    heavy = _heavy_molecule(molecule)
    return sum(1 for bond in heavy.GetBonds() if bond.GetIsAromatic())


def test_meeko_smiles_to_pdbqt(tmp_path: Path) -> None:
    output = tmp_path / "ligand.pdbqt"

    path = prepare_ligand_with_meeko("CCO", output, name="ETH")

    assert path == output
    assert output.exists()

    content = output.read_text(encoding="utf-8")
    assert "ROOT" in content
    assert "ATOM" in content


def test_vina_build_command(tmp_path: Path) -> None:
    engine = VinaEngine(binary="vina")

    receptor = tmp_path / "receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    output = tmp_path / "out.pdbqt"

    receptor.write_text("", encoding="utf-8")
    ligand.write_text("", encoding="utf-8")

    cmd = engine._build_command(
        receptor=receptor,
        ligand=ligand,
        output=output,
        box=DOCKPROTEIN_BOX,
        seed=42,
        num_modes=5,
        exhaustiveness=3,
        energy_range=4.5,
        extra_flags={"--cpu": 2},
    )

    assert cmd[:2] == ["vina", "--receptor"]
    assert "--seed" in cmd
    assert "42" in cmd
    assert "--cpu" in cmd
    assert "2" in cmd
    assert "--num_modes" in cmd
    assert "5" in cmd
    assert "--exhaustiveness" in cmd
    assert "3" in cmd
    assert "--energy_range" in cmd
    assert "4.5" in cmd


@pytest.mark.integration
def test_pdbqt_to_sdf_roundtrip_preserves_pose(tmp_path: Path) -> None:
    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    roundtrip_pdbqt = tmp_path / "dockligand_roundtrip.pdbqt"
    roundtrip_raw_sdf = tmp_path / "dockligand_roundtrip_raw.sdf"

    input_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)

    prepare_ligand_with_meeko(input_molecule, roundtrip_pdbqt, name="DOCKLIG")
    export_pdbqt_to_sdf(roundtrip_pdbqt, roundtrip_raw_sdf)

    assert roundtrip_pdbqt.exists()
    assert roundtrip_raw_sdf.exists()

    exported_molecule = _read_first_sdf_molecule(roundtrip_raw_sdf)

    input_heavy = _heavy_molecule(input_molecule)
    exported_heavy = _heavy_molecule(exported_molecule)

    assert input_heavy.GetNumAtoms() > 0
    assert exported_heavy.GetNumAtoms() > 0
    assert input_heavy.GetNumAtoms() == exported_heavy.GetNumAtoms()

    heavy_atom_rmsd = rdMolAlign.CalcRMS(exported_heavy, input_heavy)

    assert heavy_atom_rmsd < 1.0


@pytest.mark.integration
def test_pdbqt_to_sdf_template_reconstruction_restores_chemistry_without_pose_drift(
    tmp_path: Path,
) -> None:
    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    roundtrip_pdbqt = tmp_path / "dockligand_roundtrip.pdbqt"
    roundtrip_raw_sdf = tmp_path / "dockligand_roundtrip_raw.sdf"
    roundtrip_rebuilt_sdf = tmp_path / "dockligand_roundtrip_rebuilt.sdf"

    template_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)

    prepare_ligand_with_meeko(template_molecule, roundtrip_pdbqt, name="DOCKLIG")

    raw_export_path = export_pdbqt_to_sdf(
        roundtrip_pdbqt,
        roundtrip_raw_sdf,
        template_bond_orders=False,
    )
    rebuilt_export_path = export_pdbqt_to_sdf(
        roundtrip_pdbqt,
        roundtrip_rebuilt_sdf,
        template_mol=template_molecule,
        template_bond_orders=True,
        add_hydrogens_after_template=True,
    )

    assert raw_export_path.exists()
    assert rebuilt_export_path.exists()

    raw_molecule = _read_first_sdf_molecule(raw_export_path)
    rebuilt_molecule = _read_first_sdf_molecule(rebuilt_export_path)

    raw_heavy = _heavy_molecule(raw_molecule)
    rebuilt_heavy = _heavy_molecule(rebuilt_molecule)
    template_heavy = _heavy_molecule(template_molecule)

    assert raw_heavy.GetNumAtoms() > 0
    assert rebuilt_heavy.GetNumAtoms() > 0
    assert template_heavy.GetNumAtoms() > 0

    assert raw_heavy.GetNumAtoms() == rebuilt_heavy.GetNumAtoms()
    assert rebuilt_heavy.GetNumAtoms() == template_heavy.GetNumAtoms()

    raw_centroid = _centroid(raw_heavy)
    rebuilt_centroid = _centroid(rebuilt_heavy)
    centroid_shift = _distance(raw_centroid, rebuilt_centroid)

    raw_to_rebuilt_rmsd = rdMolAlign.CalcRMS(rebuilt_heavy, raw_heavy)

    assert centroid_shift < 0.5
    assert raw_to_rebuilt_rmsd < 1.0

    raw_bond_order_delta = abs(_bond_order_sum(raw_molecule) - _bond_order_sum(template_molecule))
    rebuilt_bond_order_delta = abs(_bond_order_sum(rebuilt_molecule) - _bond_order_sum(template_molecule))

    raw_aromatic_delta = abs(_aromatic_bond_count(raw_molecule) - _aromatic_bond_count(template_molecule))
    rebuilt_aromatic_delta = abs(_aromatic_bond_count(rebuilt_molecule) - _aromatic_bond_count(template_molecule))

    assert rebuilt_bond_order_delta <= raw_bond_order_delta
    assert rebuilt_aromatic_delta <= raw_aromatic_delta
    assert rebuilt_bond_order_delta == 0.0
    assert rebuilt_aromatic_delta == 0


@pytest.mark.integration
def test_vina_binary_smoke(tmp_path: Path) -> None:
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    if not DOCKPROTEIN_PDB.exists():
        pytest.skip(f"missing receptor test file: {DOCKPROTEIN_PDB}")

    engine = VinaEngine(binary="vina")
    docking_input_pdbqt = tmp_path / "dockligand_for_docking.pdbqt"

    prepare_ligand_with_meeko(
        _read_first_sdf_molecule(DOCKLIGAND_SDF),
        docking_input_pdbqt,
        name="DOCKLIG",
    )

    request = DockingRequest(
        receptor=DOCKPROTEIN_PDB,
        ligands=[docking_input_pdbqt],
        box=DOCKPROTEIN_BOX,
        workdir=tmp_path,
    )

    result = engine.dock(request=request)

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0
    assert result.poses[0].score is not None


@pytest.mark.integration
def test_docked_pose_reconstruction_restores_chemistry_and_keeps_docked_position(
    tmp_path: Path,
) -> None:
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKPROTEIN_PDB.exists():
        pytest.skip(f"missing receptor test file: {DOCKPROTEIN_PDB}")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    template_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)
    docking_input_pdbqt = tmp_path / "dockligand_for_docking.pdbqt"
    docking_output_pdbqt = tmp_path / "dockligand_for_docking_vina_out.pdbqt"
    docking_output_raw_sdf = tmp_path / "dockligand_for_docking_vina_out_raw.sdf"
    docking_output_rebuilt_sdf = tmp_path / "dockligand_for_docking_vina_out_rebuilt.sdf"

    prepare_ligand_with_meeko(
        template_molecule,
        docking_input_pdbqt,
        name="DOCKLIG",
    )

    engine = VinaEngine(binary="vina")
    request = DockingRequest(
        receptor=DOCKPROTEIN_PDB,
        ligands=[docking_input_pdbqt],
        box=DOCKPROTEIN_BOX,
        workdir=tmp_path,
    )

    result = engine.dock(request=request)

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path == docking_output_pdbqt
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0

    raw_docked_sdf = export_pdbqt_to_sdf(
        docking_output_pdbqt,
        docking_output_raw_sdf,
        template_bond_orders=False,
    )
    rebuilt_docked_sdf = export_pdbqt_to_sdf(
        docking_output_pdbqt,
        docking_output_rebuilt_sdf,
        template_mol=template_molecule,
        template_bond_orders=True,
        add_hydrogens_after_template=True,
    )

    assert raw_docked_sdf.exists()
    assert rebuilt_docked_sdf.exists()

    raw_docked_molecule = _read_first_sdf_molecule(raw_docked_sdf)
    rebuilt_docked_molecule = _read_first_sdf_molecule(rebuilt_docked_sdf)

    raw_docked_heavy = _heavy_molecule(raw_docked_molecule)
    rebuilt_docked_heavy = _heavy_molecule(rebuilt_docked_molecule)
    template_heavy = _heavy_molecule(template_molecule)

    assert raw_docked_heavy.GetNumAtoms() > 0
    assert rebuilt_docked_heavy.GetNumAtoms() > 0
    assert template_heavy.GetNumAtoms() > 0

    assert raw_docked_heavy.GetNumAtoms() == rebuilt_docked_heavy.GetNumAtoms()
    assert rebuilt_docked_heavy.GetNumAtoms() == template_heavy.GetNumAtoms()

    raw_centroid = _centroid(raw_docked_heavy)
    rebuilt_centroid = _centroid(rebuilt_docked_heavy)
    centroid_shift = _distance(raw_centroid, rebuilt_centroid)

    raw_to_rebuilt_rmsd = rdMolAlign.CalcRMS(
        rebuilt_docked_heavy,
        raw_docked_heavy,
    )

    assert centroid_shift < 0.5
    assert raw_to_rebuilt_rmsd < 1.0

    rebuilt_bond_order_delta = abs(_bond_order_sum(rebuilt_docked_molecule) - _bond_order_sum(template_molecule))
    rebuilt_aromatic_delta = abs(
        _aromatic_bond_count(rebuilt_docked_molecule) - _aromatic_bond_count(template_molecule)
    )

    assert rebuilt_bond_order_delta == 0.0
    assert rebuilt_aromatic_delta == 0
