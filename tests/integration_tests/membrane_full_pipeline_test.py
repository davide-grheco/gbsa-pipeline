"""Full run_pipeline() integration test for a real membrane-protein system.

Exercises every stage (1-9) end-to-end -- ligand parametrization + merge,
membrane solvation, SD/CG minimization, NVT, NPT (restrained/unrestrained),
production MD, and gmx_MMPBSA -- using short simulation times. This proves
the pipeline *plumbing* is correct (file names, molecule ordering, index
groups, stage wiring), not production-quality sampling. Docking itself is
already covered by membrane_and_docking_inegration_test.py; this test starts
from the already-posed ligand.sdf to keep scope focused on stages 1-9.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest

from gbsa_pipeline.config import (
    EquilibrationConfig,
    MembraneConfig,
    MinimizationConfig,
    NptConfig,
    RunConfig,
    SystemConfig,
)
from gbsa_pipeline.mdp import GromacsParams
from gbsa_pipeline.pipeline import run_pipeline

TESTDATA = Path(__file__).resolve().parents[1] / "testdata" / "membrane" / "2rh1"


@pytest.mark.integration
def test_run_pipeline_membrane_end_to_end(tmp_path: Path) -> None:
    """Parametrize -> merge -> solvate -> minimize -> NVT -> NPT -> production -> GBSA."""
    if shutil.which("gmx_MMPBSA") is None:
        pytest.skip("gmx_MMPBSA not available in PATH")

    config = RunConfig(
        system=SystemConfig(
            gro_file=TESTDATA / "system_unsolvated.gro",
            top_file=TESTDATA / "system_unsolvated.top",
            ligand=TESTDATA / "ligand.sdf",
            net_charge=0,
            membrane=True,
        ),
        membrane=MembraneConfig(z_padding_nm=1.5),
        minimization=MinimizationConfig(nsteps=10000),
        equilibration=EquilibrationConfig(simulation_time_ps=2.0),
        npt_equilibration=NptConfig(simulation_time_ps=2.0),
        md=GromacsParams(nsteps=100, dt=0.001, nstxout_compressed=10),
    )

    output_dir = tmp_path / "run"
    run_pipeline(config, output_dir)

    for stage_label in (
        "01_parametrize",
        "02_solvated",
        "03_sd",
        "04_cg",
        "05_nvt_res",
        "06_npt_res",
        "07_npt",
        "08_production",
        "09_mmbsa",
    ):
        stage_dir = output_dir / stage_label
        assert stage_dir.exists(), f"missing stage output directory: {stage_label}"
        assert any(stage_dir.iterdir()), f"stage output directory is empty: {stage_label}"

    # gmx_MMPBSA's own success artifact -- stronger than "any file exists",
    # since run_gmx_mmpbsa_from_gromacs uses check=False and _stage_mmbsa's
    # subprocess result isn't propagated up through run_pipeline().
    final_results = output_dir / "09_mmbsa" / "FINAL_RESULTS_MMPBSA.dat"
    assert final_results.exists()

    # Membrane PBSA used to copy the final-frame PDB as a one-frame "trajectory",
    # so FINAL_RESULTS always reported "Calculations performed using 1 complex
    # frames" regardless of production length. After the trjconv-based
    # extraction, the N-frame production trajectory must flow through to the
    # per-frame sander calls so that user sampling budget is honoured.
    text = final_results.read_text()
    match = re.search(r"Calculations performed using (\d+) complex frames", text)
    assert match is not None, "FINAL_RESULTS_MMPBSA.dat is missing the frame-count line"
    n_frames = int(match.group(1))
    assert n_frames > 1, (
        f"gmx_MMPBSA ran on only {n_frames} frame -- _stage_mmbsa is still using "
        "the single-frame shortcut instead of the full production trajectory"
    )
