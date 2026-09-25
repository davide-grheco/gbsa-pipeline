"""Tests for the shared residue-name constants."""

from __future__ import annotations

from gbsa_pipeline._constants import ION_RESIDUE_NAMES, SOLVENT_RESIDUE_NAMES, WATER_RESIDUE_NAMES


def test_solvent_residue_names_is_union_of_water_and_ions() -> None:
    assert SOLVENT_RESIDUE_NAMES == WATER_RESIDUE_NAMES | ION_RESIDUE_NAMES


def test_water_and_ion_names_are_disjoint() -> None:
    assert not WATER_RESIDUE_NAMES & ION_RESIDUE_NAMES
