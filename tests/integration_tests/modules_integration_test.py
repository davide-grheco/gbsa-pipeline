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
from gbsa_pipeline.parametrization import (
    ParametrizationInput,
    parametrize,
)
from gbsa_pipeline.solvation_box import SolvationParams
from gbsa_pipeline.solvation_openmm import solvate_openmm

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
def test_prepare_inputs_run_docking_parametrize_and_solvate_keeps_outputs(
    visual_run_dir: Path,
) -> None:
    """Run docking, parametrize the docked complex, and solvate it.

    This test is a workflow smoke test, not a detailed unit test of the docking,
    parametrization, solvation, or BioSimSpace MD modules. The ligand starts
    from SDF, the receptor starts from PDB, and the docking output is exported
    back to SDF before it is passed into the parametrization entry point. The
    parametrized GROMACS coordinate and topology files are then passed into the
    OpenMM/ParmEd solvation bridge, where retained crystallographic waters are
    restored before bulk solvent is added. We currently check only that each
    bridge produces the expected files and that the solvated system can be
    handed to the later MD helpers.
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
    docked_sdf = visual_run_dir / "dockligand_vina_out.sdf"

    parametrization_dir = visual_run_dir / "parametrization"
    crystal_waters_pdb = parametrization_dir / "crystal_waters.pdb"

    solvation_dir = visual_run_dir / "solvation"
    restored_crystal_waters_pdb = solvation_dir / "restored_crystal_waters.pdb"
    solvation_dir.mkdir(parents=True, exist_ok=True)

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

    docking_result = engine.dock(request=request)

    docking_output = visual_run_dir / "dockligand_vina_out.pdbqt"
    docking_log = visual_run_dir / "dockligand_vina.log"

    assert docking_result.engine == "vina"
    assert len(docking_result.poses) == 1
    assert docking_result.poses[0].pose_path == docking_output
    assert docking_result.poses[0].pose_path.exists()
    assert docking_result.poses[0].metadata["returncode"] == 0
    assert docking_result.poses[0].score is not None
    assert docking_output.exists()
    assert docking_log.exists()

    exported_sdf = export_pdbqt_to_sdf(
        docking_output,
        docked_sdf,
        template_mol=ligand_molecule,
        add_hydrogens_after_template=True,
    )

    assert exported_sdf == docked_sdf
    assert docked_sdf.exists()

    parametrized = parametrize(
        ParametrizationInput(
            protein_pdb=DOCKPROTEIN_PDB,
            ligand_sdf=docked_sdf,
            work_dir=parametrization_dir,
        )
    )

    assert crystal_waters_pdb.exists()
    assert crystal_waters_pdb.read_text(encoding="utf-8").strip()
    assert parametrized.crystal_waters_pdb == crystal_waters_pdb
    assert parametrized.gro_file == parametrization_dir / "complex.gro"
    assert parametrized.top_file == parametrization_dir / "complex.top"
    assert parametrized.gro_file.exists()
    assert parametrized.top_file.exists()

    solvated = solvate_openmm(
        parametrized=parametrized,
        params=SolvationParams(
            water_model="tip3p",
            shape="truncated_octahedron",
            padding=1.0,
            ion_concentration=0.15,
            neutralize=True,
        ),
        output_gro=solvation_dir / "solvated.gro",
        output_top=solvation_dir / "solvated.top",
    )

    assert solvated.gro_file == solvation_dir / "solvated.gro"
    assert solvated.top_file == solvation_dir / "solvated.top"
    assert solvated.gro_file.exists()
    assert solvated.top_file.exists()
    assert restored_crystal_waters_pdb.exists()
    assert restored_crystal_waters_pdb.read_text(encoding="utf-8").strip()

    bss_system = solvated.load_bss()

    assert bss_system is not None
