"""Tests for the shared gemmi structure helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import gemmi

from gbsa_pipeline._gemmi_utils import append_chains, filter_residues, iter_chains, residue_is_water, write_pdb

if TYPE_CHECKING:
    from pathlib import Path

_TWO_CHAIN_PDB = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
HETATM    3  O   HOH B   1      10.000  10.000  10.000  1.00  0.00           O
HETATM    4  ZN  ZN  B   2       5.000   5.000   5.000  1.00  0.00          ZN
END
"""


def _residue_names(st: gemmi.Structure) -> list[str]:
    return [res.name for model in st for chain in model for res in chain]


def test_iter_chains_yields_chains_across_models() -> None:
    st = gemmi.read_pdb_string(_TWO_CHAIN_PDB)

    assert [chain.name for chain in iter_chains(st)] == ["A", "B"]


def test_residue_is_water_matches_known_names_case_insensitively() -> None:
    st = gemmi.read_pdb_string(_TWO_CHAIN_PDB)
    flags = {res.name: residue_is_water(res) for model in st for chain in model for res in chain}

    assert flags == {"ALA": False, "HOH": True, "ZN": False}


def test_filter_residues_keeps_only_matching() -> None:
    st = gemmi.read_pdb_string(_TWO_CHAIN_PDB)

    result = filter_residues(st, residue_is_water)

    assert result is st
    assert _residue_names(st) == ["HOH"]


def test_filter_residues_keep_all_is_noop() -> None:
    st = gemmi.read_pdb_string(_TWO_CHAIN_PDB)

    filter_residues(st, lambda _res: True)

    assert _residue_names(st) == ["ALA", "HOH", "ZN"]


def test_append_chains_merges_and_renames_conflicts() -> None:
    base = gemmi.read_pdb_string(_TWO_CHAIN_PDB)
    donor = gemmi.read_pdb_string(_TWO_CHAIN_PDB)
    n_base_chains = len(base[0])

    append_chains(base, donor)

    assert len(base[0]) == n_base_chains + len(donor[0])
    chain_names = [chain.name for chain in base[0]]
    assert len(chain_names) == len(set(chain_names))
    assert _residue_names(base) == ["ALA", "HOH", "ZN", "ALA", "HOH", "ZN"]


def test_write_pdb_creates_parent_directories(tmp_path: Path) -> None:
    st = gemmi.read_pdb_string(_TWO_CHAIN_PDB)
    out = tmp_path / "nested" / "out.pdb"

    result = write_pdb(st, out)

    assert result == out
    assert "ALA" in out.read_text(encoding="utf-8")
