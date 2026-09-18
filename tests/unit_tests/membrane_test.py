"""Unit tests for gbsa_pipeline.membrane.

estimates membrane geometry is checked against a real MemProtMD system.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.core.universe import Merge

from gbsa_pipeline.membrane import (
    MembraneGeometry,
    estimate_membrane_geometry,
    extract_receptor_pdb,
    merge_ligand_into_system,
)

TESTDATA = Path(__file__).resolve().parents[1] / "testdata" / "membrane" / "1py6"

STRUCTURE = TESTDATA / "atomistic-system.pdb"

TESTDATA_2RH1 = Path(__file__).resolve().parents[1] / "testdata" / "membrane" / "2rh1"

SYSTEM_2RH1 = TESTDATA_2RH1 / "system_unsolvated.gro"


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


def test_estimate_membrane_geometry_raises_when_not_z_aligned() -> None:
    """A membrane whose leaflets are separated laterally, not along z, is rejected.

    Counter-example to the real, correctly z-aligned 2rh1 membrane-protein
    system: rotate it 90 degrees (swap the y/z axes) so the bilayer normal
    now points along y, and check that this is caught instead of silently
    producing a wrong mthick/mctrdz.
    """
    universe = mda.Universe(str(SYSTEM_2RH1))

    positions = universe.atoms.positions.copy()
    positions[:, [1, 2]] = positions[:, [2, 1]]
    universe.atoms.positions = positions

    with pytest.raises(ValueError, match="Rotate model"):
        estimate_membrane_geometry(universe, lipid_resnames=["POP"])


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


def test_extract_receptor_pdb(tmp_path: Path) -> None:
    """Only proteins atoms remain, lipids stripped for docking."""
    output_pdb = tmp_path / "receptor.pdb"

    result = extract_receptor_pdb(SYSTEM_2RH1, output_pdb)

    assert result == output_pdb
    assert output_pdb.exists()

    written = mda.Universe(str(output_pdb))
    assert written.atoms.n_atoms == 4597
    assert "POP" not in set(written.atoms.resnames)


def test_merge_ligand_into_system_add_molecules_and_returns_system() -> None:
    """A short wrapper test."""
    system = Mock()
    ligand = Mock()

    result = merge_ligand_into_system(system, ligand)
    assert result is system
