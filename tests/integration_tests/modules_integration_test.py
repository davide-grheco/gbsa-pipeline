# tests/integration_tests/modules_integration_test.py

"""Integration tests for visually inspectable module-level workflows.

This module intentionally keeps generated files on disk because the current
development step requires manual inspection of docking and later MD artefacts.
That is an exception to the normal test policy where runtime files should use
``tmp_path`` and be discarded after the test run. The output directory is kept
under ``tests/integration_tests/module_runs/`` so it is easy to inspect and easy
to exclude from commits. Once the workflows are stable, these tests should move
back to disposable runtime directories or become explicit manual smoke tests.
"""

from __future__ import annotations

import shutil
from pathlib import Path

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

TESTDATA = Path(__file__).parents[1] / "testdata"
DOCKING_TESTDATA = TESTDATA / "docking"
DOCKLIGAND_SDF = DOCKING_TESTDATA / "dockligand.sdf"
DOCKPROTEIN_PDB = DOCKING_TESTDATA / "dockprotein.pdb"

MODULE_RUNS = Path(__file__).parent / "module_runs"

MEEKO_RECEPTOR_BINARY = "mk_prepare_receptor.py"

DOCKPROTEIN_BOX = DockingBox(
    center=(10.115, 39.148, 53.112),
    size=(10.0, 10.0, 10.0),
)


@pytest.fixture
def visual_run_dir(request: pytest.FixtureRequest) -> Path:
    """Return a persistent output directory for one integration test.

    This fixture is intentionally not based on ``tmp_path`` because the current
    module-integration work needs files to remain available after the test run.
    The pytest node name is used as the directory name so each test has a stable
    and inspectable output location. Existing files are not removed here,
    because preserving artefacts is the point of this temporary test module.
    The directory should be ignored by git and should not be committed.
    """
    run_dir = MODULE_RUNS / request.node.name
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


@pytest.mark.integration
def test_prepare_inputs_and_run_docking_keeps_outputs(visual_run_dir: Path) -> None:
    """Prepare docking inputs and keep the Vina output files for inspection.

    This test is a workflow smoke test, not a detailed unit test of the docking
    module. The ligand starts from SDF and is converted to PDBQT through the
    public ligand-preparation helper, while the receptor starts from PDB and is
    converted to PDBQT through the public receptor-preparation helper. The Vina
    run then uses those prepared files and writes all artefacts into a persistent
    module-run directory so they can be inspected manually. We currently check
    only that the workflow produces the expected files, a successful return
    code, and a parsed docking score.
    """
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    if shutil.which(MEEKO_RECEPTOR_BINARY) is None:
        pytest.skip(f"{MEEKO_RECEPTOR_BINARY} not available in PATH")

    if not DOCKLIGAND_SDF.exists():
        pytest.skip(f"missing ligand test file: {DOCKLIGAND_SDF}")

    if not DOCKPROTEIN_PDB.exists():
        pytest.skip(f"missing receptor test file: {DOCKPROTEIN_PDB}")

    ligand_pdbqt = visual_run_dir / "dockligand.pdbqt"
    receptor_pdbqt = visual_run_dir / "dockprotein.pdbqt"

    ligand_molecule = load_first_sdf_molecule(DOCKLIGAND_SDF, remove_hs=False)

    prepared_ligand = prepare_ligand_with_meeko(
        ligand_molecule,
        ligand_pdbqt,
        name="DOCKLIG",
    )
    prepared_receptor = convert_receptor_pdb_to_pdbqt(
        DOCKPROTEIN_PDB,
        output_path=receptor_pdbqt,
        mk_prepare_receptor_binary=MEEKO_RECEPTOR_BINARY,
    )

    assert prepared_ligand == ligand_pdbqt
    assert prepared_receptor == receptor_pdbqt
    assert ligand_pdbqt.exists()
    assert receptor_pdbqt.exists()

    engine = VinaEngine(
        binary="vina",
        mk_prepare_receptor_binary=MEEKO_RECEPTOR_BINARY,
    )
    request = DockingRequest(
        receptor=receptor_pdbqt,
        ligands=[ligand_pdbqt],
        box=DOCKPROTEIN_BOX,
        workdir=visual_run_dir,
    )

    result = engine.dock(request=request)

    docking_output = visual_run_dir / "dockligand_vina_out.pdbqt"
    docking_log = visual_run_dir / "dockligand_vina.log"
    docking_output_sdf = visual_run_dir / "dockligand_vina_out.sdf"

    exported_sdf = export_pdbqt_to_sdf(
        docking_output,
        docking_output_sdf,
        template_mol=ligand_molecule,
        add_hydrogens_after_template=True,
    )

    assert exported_sdf == docking_output_sdf
    assert docking_output_sdf.exists()

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path == docking_output
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0
    assert result.poses[0].score is not None
    assert docking_output.exists()
    assert docking_log.exists()
