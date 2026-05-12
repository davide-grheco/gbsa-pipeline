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

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Any

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
from gbsa_pipeline.md import (
    run_heating,
    run_minimization,
    run_npt_equilibration,
    run_production,
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

SD_PARAMS: dict[str, Any] = {
    "integrator": "steep",
    "nsteps": 50000,
    "emtol": 1000.0,
    "emstep": 0.001,
    "cutoff_scheme": "Verlet",
    "nstlist": 20,
    "pbc": "xyz",
    "rlist": 1.26,
    "coulombtype": "PME",
    "rcoulomb": 1.2,
    "fourierspacing": 0.0,
    "pme_order": 4,
    "ewald_rtol": 1e-5,
    "vdwtype": "Cut-off",
    "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0,
    "rvdw": 1.2,
    "constraints": "none",
    "tcoupl": "no",
    "pcoupl": "no",
    "gen_vel": "no",
    "nstlog": 500,
    "nstenergy": 500,
    "nstxout_compressed": 0,
}

CG_PARAMS: dict[str, Any] = {
    "integrator": "cg",
    "nsteps": 5000,
    "emtol": 100.0,
    "emstep": 0.001,
    "nstcgsteep": 50,
    "cutoff_scheme": "Verlet",
    "nstlist": 20,
    "pbc": "xyz",
    "rlist": 1.221,
    "coulombtype": "PME",
    "rcoulomb": 1.2,
    "fourierspacing": 0.16,
    "pme_order": 4,
    "ewald_rtol": 1e-5,
    "vdwtype": "Cut-off",
    "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0,
    "rvdw": 1.2,
    "constraints": "none",
    "tcoupl": "no",
    "pcoupl": "no",
    "gen_vel": "no",
    "nstlog": 500,
    "nstenergy": 500,
    "nstxout_compressed": 0,
}

HEATING_PARAMS: dict[str, Any] = {
    "integrator": "md",
    "dt": 0.001,
    "nsteps": 100000,
    "cutoff_scheme": "Verlet",
    "nstlist": 20,
    "pbc": "xyz",
    "rlist": 1.221,
    "coulombtype": "PME",
    "rcoulomb": 1.2,
    "fourierspacing": 0.16,
    "pme_order": 4,
    "ewald_rtol": 1e-5,
    "vdwtype": "Cut-off",
    "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0,
    "rvdw": 1.2,
    "constraints": "h-bonds",
    "constraint_algorithm": "LINCS",
    "lincs_order": 4,
    "lincs_warnangle": 30.0,
    "continuation": "no",
    "gen_vel": "yes",
    "gen_temp": 20.0,
    "tcoupl": "V-rescale",
    "tc_grps": "System",
    "tau_t": 0.1,
    "ref_t": 300.0,
    "pcoupl": "no",
    "annealing": "single",
    "annealing_npoints": 8,
    "annealing_time": "0 40.0 50.0 60.0 70.0 80.0 90.0 100.0",
    "annealing_temp": "20.0 50.0 50.0 100.0 150.0 200.0 250.0 300.0 ",
    "define": "-DPOSRES",
    "nstlog": 500,
    "nstenergy": 500,
    "nstcalcenergy": 100,
    "nstxout_compressed": 1000,
}

NPT_RESTRAINED_PARAMS: dict[str, Any] = {
    "integrator": "md",
    "dt": 0.002,
    "nsteps": 50000,
    "cutoff_scheme": "Verlet",
    "nstlist": 20,
    "pbc": "xyz",
    "rlist": 1.221,
    "coulombtype": "PME",
    "rcoulomb": 1.2,
    "fourierspacing": 0.16,
    "pme_order": 4,
    "ewald_rtol": 1e-5,
    "vdwtype": "Cut-off",
    "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0,
    "rvdw": 1.2,
    "constraints": "h-bonds",
    "constraint_algorithm": "LINCS",
    "lincs_order": 4,
    "lincs_warnangle": 30.0,
    "continuation": "yes",
    "gen_vel": "no",
    "tcoupl": "V-rescale",
    "tc_grps": "System",
    "tau_t": 0.1,
    "ref_t": 300.0,
    "pcoupl": "C-rescale",
    "pcoupltype": "isotropic",
    "tau_p": 5.0,
    "ref_p": 1.0,
    "compressibility": 4.5e-5,
    "refcoord_scaling": "com",
    "define": "-DPOSRES",
    "nstlog": 500,
    "nstenergy": 500,
    "nstcalcenergy": 100,
    "nstxout_compressed": 1000,
}

NPT_PARAMS: dict[str, Any] = {
    "integrator": "md",
    "dt": 0.002,
    "nsteps": 50000,
    "cutoff_scheme": "Verlet",
    "nstlist": 20,
    "pbc": "xyz",
    "rlist": 1.221,
    "coulombtype": "PME",
    "rcoulomb": 1.2,
    "fourierspacing": 0.16,
    "pme_order": 4,
    "ewald_rtol": 1e-5,
    "vdwtype": "Cut-off",
    "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0,
    "rvdw": 1.2,
    "constraints": "h-bonds",
    "constraint_algorithm": "LINCS",
    "lincs_order": 4,
    "lincs_warnangle": 30.0,
    "continuation": "yes",
    "gen_vel": "no",
    "tcoupl": "V-rescale",
    "tc_grps": "System",
    "tau_t": 0.1,
    "ref_t": 300.0,
    "pcoupl": "C-rescale",
    "pcoupltype": "isotropic",
    "tau_p": 5.0,
    "ref_p": 1.0,
    "compressibility": 4.5e-5,
    "define": "",
    "nstlog": 500,
    "nstenergy": 500,
    "nstcalcenergy": 100,
    "nstxout_compressed": 1000,
}

PRODUCTION_PARAMS: dict[str, Any] = {
    "integrator": "md",
    "dt": 0.002,
    "nsteps": 250000,
    "cutoff_scheme": "Verlet",
    "nstlist": 20,
    "pbc": "xyz",
    "rlist": 1.221,
    "coulombtype": "PME",
    "rcoulomb": 1.2,
    "fourierspacing": 0.16,
    "pme_order": 4,
    "ewald_rtol": 1e-5,
    "vdwtype": "Cut-off",
    "vdw_modifier": "Force-switch",
    "rvdw_switch": 1.0,
    "rvdw": 1.2,
    "constraints": "h-bonds",
    "constraint_algorithm": "LINCS",
    "lincs_order": 4,
    "lincs_warnangle": 30.0,
    "continuation": "yes",
    "gen_vel": "no",
    "tcoupl": "V-rescale",
    "tc_grps": "System",
    "tau_t": 0.1,
    "ref_t": 300.0,
    "pcoupl": "C-rescale",
    "pcoupltype": "isotropic",
    "tau_p": 5.0,
    "ref_p": 1.0,
    "compressibility": 4.5e-5,
    "define": "",
    "nstlog": 500,
    "nstenergy": 500,
    "nstcalcenergy": 100,
    "nstxout_compressed": 1000,
}


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
    """Run docking, solvation, minimization, equilibration, and production.

    This test is a workflow smoke test, not a detailed unit test of the docking,
    parametrization, solvation, or MD helper modules. The ligand starts from SDF,
    the receptor starts from PDB, and the docking output is exported back to SDF
    before it is passed into the parametrization entry point. The parametrized
    complex is solvated with retained crystallographic waters restored before
    bulk solvent placement, then loaded into BioSimSpace. The MD bridge runs the
    exact six stage parameter blocks copied from the combined direct-GROMACS MD
    runner: SD, CG, restrained NVT heating, restrained NPT, unrestrained NPT,
    and production.
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

    sd_minimization_dir = visual_run_dir / "sd"
    cg_minimization_dir = visual_run_dir / "cg"
    nvt_restrained_dir = visual_run_dir / "nvt_res"
    npt_restrained_dir = visual_run_dir / "npt_res"
    npt_unrestrained_dir = visual_run_dir / "npt"
    production_dir = visual_run_dir / "production"

    sd_minimization_dir.mkdir(parents=True, exist_ok=True)
    cg_minimization_dir.mkdir(parents=True, exist_ok=True)
    nvt_restrained_dir.mkdir(parents=True, exist_ok=True)
    npt_restrained_dir.mkdir(parents=True, exist_ok=True)
    npt_unrestrained_dir.mkdir(parents=True, exist_ok=True)
    production_dir.mkdir(parents=True, exist_ok=True)

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

    sd_minimized = run_minimization(
        bss_system,
        work_dir=sd_minimization_dir,
        params=SD_PARAMS,
    )

    assert sd_minimized is not None
    assert any(sd_minimization_dir.iterdir())

    cg_minimized = run_minimization(
        sd_minimized,
        work_dir=cg_minimization_dir,
        params=CG_PARAMS,
    )

    assert cg_minimized is not None
    assert any(cg_minimization_dir.iterdir())

    nvt_restrained = run_heating(
        50 * BSS.Units.Time.picosecond,
        cg_minimized,
        work_dir=nvt_restrained_dir,
        params=HEATING_PARAMS,
        temperature_start=50 * BSS.Units.Temperature.kelvin,
        temperature_end=300 * BSS.Units.Temperature.kelvin,
        restraint="backbone",
    )

    assert nvt_restrained is not None
    assert any(nvt_restrained_dir.iterdir())

    npt_restrained = run_npt_equilibration(
        100 * BSS.Units.Time.picosecond,
        nvt_restrained,
        work_dir=npt_restrained_dir,
        params=NPT_RESTRAINED_PARAMS,
        restraint="backbone",
    )

    assert npt_restrained is not None
    assert any(npt_restrained_dir.iterdir())

    npt_unrestrained = run_npt_equilibration(
        100 * BSS.Units.Time.picosecond,
        npt_restrained,
        work_dir=npt_unrestrained_dir,
        params=NPT_PARAMS,
        restraint=None,
    )

    assert npt_unrestrained is not None
    assert any(npt_unrestrained_dir.iterdir())

    production = run_production(
        500 * BSS.Units.Time.picosecond,
        npt_unrestrained,
        work_dir=production_dir,
        params=PRODUCTION_PARAMS,
    )

    assert production is not None
    assert any(production_dir.iterdir())
