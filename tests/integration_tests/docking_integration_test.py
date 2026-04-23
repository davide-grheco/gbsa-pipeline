"""Integration tests for the lightweight docking adapter layer.

This module keeps only integration-level checks.
Unit-style checks such as command construction and simple Meeko output
generation should live in tests/unit_tests/, not here.
The goal here is to verify external-tool workflows end to end.
"""

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

ROUNDTRIP_PDBQT = DOCKING_TESTDATA / "dockligand_roundtrip.pdbqt"
ROUNDTRIP_RAW_SDF = DOCKING_TESTDATA / "dockligand_roundtrip_raw.sdf"
ROUNDTRIP_REBUILT_SDF = DOCKING_TESTDATA / "dockligand_roundtrip_rebuilt.sdf"

DOCKING_INPUT_PDBQT = DOCKING_TESTDATA / "dockligand_for_docking.pdbqt"
DOCKING_OUTPUT_PDBQT = DOCKING_TESTDATA / "dockligand_for_docking_vina_out.pdbqt"
DOCKING_OUTPUT_RAW_SDF = DOCKING_TESTDATA / "dockligand_for_docking_vina_out_raw.sdf"
DOCKING_OUTPUT_REBUILT_SDF = DOCKING_TESTDATA / "dockligand_for_docking_vina_out_rebuilt.sdf"

DOCKPROTEIN_BOX = DockingBox(
    center=(10.115, 39.148, 53.112),
    size=(10.0, 10.0, 10.0),
)


def _read_first_sdf_molecule(path: Path) -> Chem.Mol:
    """Read the first molecule from an SDF file.

    This helper exists so every integration test uses the same SDF-loading
    behavior and fails in the same way when the fixture is broken.
    The `path` parameter is required because these tests use several derived
    SDF files written during roundtrip and reconstruction steps.
    We are currently checking that the requested SDF exists structurally as a
    usable RDKit molecule, not just as a text file on disk.
    """
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    molecule = supplier[0]

    if molecule is None:
        raise ValueError(f"Could not read first molecule from SDF: {path}")

    return molecule


def _centroid(molecule: Chem.Mol) -> tuple[float, float, float]:
    """Compute the geometric centroid of all atoms in one conformer.

    This helper is used as a coarse pose-preservation check when we compare a
    raw exported ligand to its reconstructed version.
    The `molecule` parameter is required because we want to evaluate the
    coordinates of whichever intermediate object the test is currently studying.
    We are currently checking that reconstruction does not shift the whole
    ligand to a different part of space even if bonding is rewritten.
    """
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
    """Return Euclidean distance between two 3D points.

    This helper exists so centroid-shift checks stay readable inside the tests.
    The `left` and `right` parameters are needed because the tests compare
    centroids from raw and rebuilt molecules, or other pose-like summaries.
    We are currently checking whether two objects occupy almost the same region
    of space after export and chemistry reconstruction.
    """
    return math.sqrt((left[0] - right[0]) ** 2 + (left[1] - right[1]) ** 2 + (left[2] - right[2]) ** 2)


def _heavy_molecule(molecule: Chem.Mol) -> Chem.Mol:
    """Return a copy of a molecule with hydrogens removed.

    We use heavy-atom-only comparisons because RDKit's own documentation warns
    that hydrogen-rich symmetry matching can lead to combinatorial explosion in
    RMS calculations, which is exactly the wrong failure mode for pose tests.
    The `molecule` parameter is required because every test stage may generate
    a different molecule object that still needs a fair, comparable metric.
    We are currently checking pose preservation and chemistry on the stable
    heavy-atom scaffold rather than on hydrogen placement details.

    Reference:
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
    """
    return Chem.RemoveHs(Chem.Mol(molecule))


def _bond_order_sum(molecule: Chem.Mol) -> float:
    """Return the sum of bond orders over the heavy-atom graph.

    This is a deliberately simple chemistry summary used to compare whether a
    reconstructed ligand moved closer to the template chemistry than the raw
    export did.
    The `molecule` parameter is needed because we compare raw export, rebuilt
    export, and template against the same scalar metric.
    We are currently checking coarse bond-order restoration, not a full graph
    isomorphism proof of chemical correctness.
    """
    heavy = _heavy_molecule(molecule)
    return sum(bond.GetBondTypeAsDouble() for bond in heavy.GetBonds())


def _aromatic_bond_count(molecule: Chem.Mol) -> int:
    """Count aromatic bonds on the heavy-atom graph.

    This helper complements `_bond_order_sum()` because aromatic perception is
    one of the places where raw export and reconstructed chemistry may differ.
    The `molecule` parameter is required because we compare template, raw, and
    rebuilt forms under the same aromaticity summary.
    We are currently checking whether reconstruction restores at least this
    basic aromatic feature set without disturbing the coordinates.
    """
    heavy = _heavy_molecule(molecule)
    return sum(1 for bond in heavy.GetBonds() if bond.GetIsAromatic())


@pytest.mark.integration
def test_pdbqt_to_sdf_roundtrip_preserves_pose() -> None:
    """Check that raw PDBQT-to-SDF export preserves the ligand pose reasonably.

    This test isolates export behavior without docking so that failures can be
    attributed to Meeko preparation or mk_export-based roundtrip handling
    rather than to Vina sampling.
    The outputs are intentionally written into tests/testdata/docking during
    review so the produced intermediate files can be inspected directly.
    We are currently checking that a ligand prepared to PDBQT and exported back
    to SDF stays close in heavy-atom geometry to the original template ligand.

    Reference:
    RDKit documents `CalcRMS()` as an in-place RMS measure that is useful for
    docking-pose comparisons and does not pre-align the probe molecule:
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
    """
    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    input_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)

    prepare_ligand_with_meeko(input_molecule, ROUNDTRIP_PDBQT, name="DOCKLIG")
    export_pdbqt_to_sdf(ROUNDTRIP_PDBQT, ROUNDTRIP_RAW_SDF)

    assert ROUNDTRIP_PDBQT.exists()
    assert ROUNDTRIP_RAW_SDF.exists()

    exported_molecule = _read_first_sdf_molecule(ROUNDTRIP_RAW_SDF)

    input_heavy = _heavy_molecule(input_molecule)
    exported_heavy = _heavy_molecule(exported_molecule)

    assert input_heavy.GetNumAtoms() > 0
    assert exported_heavy.GetNumAtoms() > 0
    assert input_heavy.GetNumAtoms() == exported_heavy.GetNumAtoms()

    heavy_atom_rmsd = rdMolAlign.CalcRMS(exported_heavy, input_heavy)

    assert heavy_atom_rmsd < 1.0


@pytest.mark.integration
def test_pdbqt_to_sdf_template_reconstruction_restores_chemistry_without_pose_drift() -> None:
    """Check isolated chemistry reconstruction on a non-docked roundtrip ligand.

    This test is intentionally narrower than the end-to-end docking test and
    exists so chemistry-rebuild failures can be debugged without involving Vina.
    The outputs are intentionally written into tests/testdata/docking during
    review so the raw and rebuilt SDF files can be inspected directly.
    We are currently checking two things at once: the rebuilt molecule should
    remain where the raw exported molecule is, and its chemistry should match
    the original template at least as well as the raw export does.

    Reference:
    RDKit `AssignBondOrdersFromTemplate()` is the core mechanism used for
    restoring bond orders from a trusted template molecule:
    https://www.rdkit.org/docs/source/rdkit.Chem.AllChem.html
    """
    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    template_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)

    prepare_ligand_with_meeko(template_molecule, ROUNDTRIP_PDBQT, name="DOCKLIG")

    raw_export_path = export_pdbqt_to_sdf(
        ROUNDTRIP_PDBQT,
        ROUNDTRIP_RAW_SDF,
        template_bond_orders=False,
    )
    rebuilt_export_path = export_pdbqt_to_sdf(
        ROUNDTRIP_PDBQT,
        ROUNDTRIP_REBUILT_SDF,
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
def test_vina_binary_smoke() -> None:
    """Check that the external Vina workflow runs and returns one pose.

    This is intentionally just a smoke test and should not be interpreted as a
    scientific validation of the docking setup or box choice.
    The prepared ligand PDBQT and Vina pose are intentionally written into
    tests/testdata/docking during review so the produced files can be inspected.
    We are currently checking only that the adapter can prepare input, call
    Vina, and produce one readable pose file with a parsed score.

    Reference:
    The Vina basic workflow requires a receptor, ligand, box center, and box
    size for a standard docking run:
    https://autodock-vina.readthedocs.io/en/latest/docking_basic.html
    """
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    if not DOCKPROTEIN_PDB.exists():
        pytest.skip(f"missing receptor test file: {DOCKPROTEIN_PDB}")

    engine = VinaEngine(binary="vina")

    prepare_ligand_with_meeko(
        _read_first_sdf_molecule(DOCKLIGAND_SDF),
        DOCKING_INPUT_PDBQT,
        name="DOCKLIG",
    )

    request = DockingRequest(
        receptor=DOCKPROTEIN_PDB,
        ligands=[DOCKING_INPUT_PDBQT],
        box=DOCKPROTEIN_BOX,
        workdir=DOCKING_TESTDATA,
    )

    result = engine.dock(request=request)

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path == DOCKING_OUTPUT_PDBQT
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0
    assert result.poses[0].score is not None


# ┌─────────────────────────────────────────────────────────────────────────────┐
# │ TEMPORARY PULL-REQUEST REVIEW NOTE                                         │
# ├─────────────────────────────────────────────────────────────────────────────┤
# │ This end-to-end test intentionally writes generated docking workflow files  │
# │ into tests/testdata/docking/ so the full path can be inspected manually    │
# │ during review. This is temporary and should be cleaned up after the PR is  │
# │ accepted as correct.                                                       │
# │                                                                             │
# │ Static inputs from tests/testdata/docking/                                  │
# │   - dockprotein.pdb                                                         │
# │   - dockligand.sdf                                                          │
# │                                                                             │
# │ Generated review artifacts in tests/testdata/docking/                       │
# │   - dockligand_for_docking.pdbqt                                            │
# │   - dockligand_for_docking_vina_out.pdbqt                                   │
# │   - dockligand_for_docking_vina_out_raw.sdf                                 │
# │   - dockligand_for_docking_vina_out_rebuilt.sdf                             │
# │                                                                             │
# │ Scientific intent                                                           │
# │   - template ligand SDF provides chemistry                                  │
# │   - docked Vina pose PDBQT provides coordinates                             │
# │   - rebuilt final SDF must keep the docked position while recovering        │
# │     chemistry from the original ligand template                             │
# └─────────────────────────────────────────────────────────────────────────────┘
@pytest.mark.integration
def test_docked_pose_reconstruction_restores_chemistry_and_keeps_docked_position() -> None:
    """Run the full ligand-prep, docking, export, and reconstruction workflow.

    This is the real end-to-end test for the current docking segment and it is
    the one that matters most scientifically for the reconstruction logic.
    The generated docking artifacts are intentionally written into
    tests/testdata/docking during review so the full workflow can be inspected.
    We are currently checking the core rule of the pipeline: the final rebuilt
    ligand must inherit chemistry from the original template while inheriting
    coordinates from the docked pose, not from the pre-docking template.

    Reference:
    RDKit `CalcRMS()` is used here specifically because its documentation notes
    that it computes RMS in place and is useful for docking-pose comparisons:
    https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html
    """
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    if shutil.which("mk_export.py") is None:
        pytest.skip("mk_export.py not available in PATH")

    if not DOCKPROTEIN_PDB.exists():
        pytest.skip(f"missing receptor test file: {DOCKPROTEIN_PDB}")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    template_molecule = _read_first_sdf_molecule(DOCKLIGAND_SDF)

    prepare_ligand_with_meeko(
        template_molecule,
        DOCKING_INPUT_PDBQT,
        name="DOCKLIG",
    )

    engine = VinaEngine(binary="vina")
    request = DockingRequest(
        receptor=DOCKPROTEIN_PDB,
        ligands=[DOCKING_INPUT_PDBQT],
        box=DOCKPROTEIN_BOX,
        workdir=DOCKING_TESTDATA,
    )

    result = engine.dock(request=request)

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path == DOCKING_OUTPUT_PDBQT
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0

    raw_docked_sdf = export_pdbqt_to_sdf(
        DOCKING_OUTPUT_PDBQT,
        DOCKING_OUTPUT_RAW_SDF,
        template_bond_orders=False,
    )
    rebuilt_docked_sdf = export_pdbqt_to_sdf(
        DOCKING_OUTPUT_PDBQT,
        DOCKING_OUTPUT_REBUILT_SDF,
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
