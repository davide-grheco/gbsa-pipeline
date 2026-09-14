"""Unit tests for gbsa_pipeline.membrane.

estimates membrane geometry is checked againts a real MemProtMD system.
"""

from __future__ import annotations

from pathlib import Path

from gbsa_pipeline.membrane import estimate_membrane_geometry

TESTDATA = Path(__file__).resolve().parents[1] / "testdata" / "membrane" / "1py6"

STRUCTURE = TESTDATA / "atomistic-system.pdb"


def test_estimate_membrane_geometry_matches_testdata() -> None:
    """Checl the 209/DPPC biilayer match hand checked values."""
    geometry = estimate_membrane_geometry(
        STRUCTURE,
        lipid_resnames=["DPP"],
    )

    assert geometry.n_phosphates == 209
    assert 35.0 < geometry.mthick < 45.0  # Based on crystal structure 39.4
    assert 0.0 < geometry.mctrdz < 97.334  # Based on crystal structure
