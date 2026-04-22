from pathlib import Path

import BioSimSpace as BSS
import pytest

from gbsa_pipeline.minimization import run_minimization


@pytest.mark.integration
def test_run_minimization() -> None:
    testdata = Path(__file__).resolve().parents[1] / "testdata" / "minimization"
    system = BSS.IO.readMolecules([str(testdata / "solvated.gro"), str(testdata / "solvated.top")])
    minimized = run_minimization(
        nsteps=500,
        system=system,
        ignore_warnings=True,
    )
    assert minimized is not None
