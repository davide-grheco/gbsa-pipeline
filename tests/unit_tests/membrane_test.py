"""Unit tests for gbsa_pipeline.membrane.

estimates membrane geometry is checked against a real MemProtMD system.
"""

from __future__ import annotations

from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.core.universe import Merge

from gbsa_pipeline.membrane import (
    MembraneGeometry,
    estimate_membrane_geometry,
)

TESTDATA = Path(__file__).resolve().parents[1] / "testdata" / "membrane" / "1py6"

STRUCTURE = TESTDATA / "atomistic-system.pdb"


def test_estimate_membrane_geometry_matches_testdata() -> None:
    """Check the 209/DPPC bilayer matches hand-checked values."""
    universe = mda.Universe(str(STRUCTURE))
    geometry = estimate_membrane_geometry(
        universe,
        lipid_resnames=["DPP"],
    )

    assert geometry.n_phosphates == 209
    assert 35.0 < geometry.mthick < 45.0  # Based on crystal structure 39.4
    assert 0.0 < geometry.mctrdz < 97.334  # Based on crystal structure


def test_estimate_membrane_geometry_ignores_contaminant_p_residue() -> None:
    """Contaminant lipid-like residues must not leak into the phosphate count.

    A residue starting with 'P' and containing a real phosphorus atom (e.g. PLM,
    palmitic acid -- a common crystallization additive) must not leak into the
    count unless its name is in lipid_resnames.
    """
    universe = mda.Universe(str(STRUCTURE))

    contaminant = mda.Universe.empty(n_atoms=1, n_residues=1, atom_resindex=[0], trajectory=True)
    contaminant.add_TopologyAttr("resnames", ["PLM"])
    contaminant.add_TopologyAttr("names", ["P1"])
    contaminant.atoms.positions = [[0.0, 0.0, 0.0]]

    merged = Merge(universe.atoms, contaminant.atoms)
    merged.trajectory.ts.positions = np.concatenate(
        [universe.atoms.positions, contaminant.atoms.positions],
    )

    geometry = estimate_membrane_geometry(merged, lipid_resnames=["DPP"])

    assert geometry.n_phosphates == 209


def test_estimate_membrane_geometry_raises_when_no_phosphates_found() -> None:
    universe = mda.Universe(str(STRUCTURE))
    with pytest.raises(ValueError, match="No phosphate atoms"):
        estimate_membrane_geometry(universe, lipid_resnames=["NOTALIPID"])


def test_membrane_geometry_pb_params() -> None:
    """pb_params() bridges a measured geometry into gmx_MMPBSA's membrane PBParams."""
    geometry = MembraneGeometry(mctrdz=50.0, mthick=39.4, n_phosphates=209)

    params = geometry.pb_params()

    assert params.memopt == 1
    assert params.mctrdz == 50.0
    assert params.mthick == 39.4
    assert params.eneopt == 1
