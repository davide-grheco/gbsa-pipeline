# tests/integration_tests/solvation_test.py

from __future__ import annotations

import pickle
from pathlib import Path

import BioSimSpace as BSS
import pytest

# importiere hier deine echten Funktionen / Klassen wie bisher
# from gbsa_pipeline.solvation_box import run_solvation
# from gbsa_pipeline.solvation_openmm import ...
# usw.


def _solvation_testdata_dir() -> Path:
    """Return the shared solvation testdata directory."""
    return Path(__file__).resolve().parents[1] / "testdata" / "solvation"


@pytest.mark.integration
def test_openmm_preparametrized(tmp_path: Path) -> None:
    """Load pre-parametrized complex from pickle testdata."""
    testdata_dir = _solvation_testdata_dir()
    pickle_file = testdata_dir / "complex.pickle"

    assert pickle_file.exists(), f"Missing testdata file: {pickle_file}"

    with pickle_file.open("rb") as f:
        complex_obj = pickle.load(f)

    assert complex_obj is not None


@pytest.mark.integration
def test_bss_solvation() -> None:
    """Read prebuilt GROMACS solvated system from testdata."""
    testdata_dir = _solvation_testdata_dir()
    gro_file = testdata_dir / "complex.gro"
    top_file = testdata_dir / "complex.top"

    assert gro_file.exists(), f"Missing testdata file: {gro_file}"
    assert top_file.exists(), f"Missing testdata file: {top_file}"

    system = BSS.IO.readMolecules([str(gro_file), str(top_file)])
    assert system is not None
