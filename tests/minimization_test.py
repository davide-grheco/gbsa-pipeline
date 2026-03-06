from __future__ import annotations

from pathlib import Path

import BioSimSpace as BSS
import pytest

from gbsa_pipeline.minimization import run_minimization


@pytest.mark.integration
def test_run_minimization() -> None:
    testdata = Path(__file__).resolve().parent / "testdata" / "minimization"
    system = BSS.IO.readMolecules([str(testdata / "solvated.gro"), str(testdata / "solvated.top")])
    minimized = run_minimization(nsteps=500, system=system)
    assert minimized
