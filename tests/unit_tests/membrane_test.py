"""Unit tests for gbsa_pipeline.membrane.

estimates membrane geometry is checked against a real MemProtMD system.
"""

from __future__ import annotations

from pathlib import Path

import gemmi
import pytest

from gbsa_pipeline.membrane import (
    MembraneGeometry,
    _is_phosphate_atom,
    estimate_membrane_geometry,
)

TESTDATA = Path(__file__).resolve().parents[1] / "testdata" / "membrane" / "1py6"

STRUCTURE = TESTDATA / "atomistic-system.pdb"


def test_estimate_membrane_geometry_matches_testdata() -> None:
    """Check the 209/DPPC bilayer matches hand-checked values."""
    struct = gemmi.read_structure(str(STRUCTURE))
    geometry = estimate_membrane_geometry(
        struct,
        lipid_resnames=["DPP"],
    )

    assert geometry.n_phosphates == 209
    assert 35.0 < geometry.mthick < 45.0  # Based on crystal structure 39.4
    assert 0.0 < geometry.mctrdz < 97.334  # Based on crystal structure


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [("P", True), ("Pt", False), ("Pb", False), ("Pd", False), ("Po", False)],
)
def test_is_phosphate_atom_excludes_other_p_elements(symbol: str, expected: bool) -> None:
    """Element symbols that merely start with 'P' (Pt, Pb, Pd, Po, ...) aren't phosphorus."""
    atom = gemmi.Atom()
    atom.element = gemmi.Element(symbol)

    assert _is_phosphate_atom(atom) is expected


def test_estimate_membrane_geometry_ignores_contaminant_p_residue() -> None:
    """Contaminant lipid-like residues must not leak into the phosphate count.

    A residue starting with 'P' and containing a real phosphorus atom (e.g. PLM,
    palmitic acid -- a common crystallization additive) must not leak into the
    count unless its name is in lipid_resnames.
    """
    struct = gemmi.read_structure(str(STRUCTURE))
    chain = struct[0][0]

    contaminant = gemmi.Residue()
    contaminant.name = "PLM"
    contaminant.seqid = gemmi.SeqId(99999, " ")
    atom = gemmi.Atom()
    atom.name = "P1"
    atom.element = gemmi.Element("P")
    atom.pos = gemmi.Position(0.0, 0.0, 0.0)
    contaminant.add_atom(atom)
    chain.add_residue(contaminant)

    geometry = estimate_membrane_geometry(struct, lipid_resnames=["DPP"])

    assert geometry.n_phosphates == 209


def test_estimate_membrane_geometry_raises_when_no_phosphates_found() -> None:
    struct = gemmi.read_structure(str(STRUCTURE))
    with pytest.raises(ValueError, match="No phosphate atoms"):
        estimate_membrane_geometry(struct, lipid_resnames=["NOTALIPID"])


def test_membrane_geometry_pb_params() -> None:
    """pb_params() bridges a measured geometry into gmx_MMPBSA's membrane PBParams."""
    geometry = MembraneGeometry(mctrdz=50.0, mthick=39.4, n_phosphates=209)

    params = geometry.pb_params()

    assert params.memopt == 1
    assert params.mctrdz == 50.0
    assert params.mthick == 39.4
    assert params.eneopt == 1
