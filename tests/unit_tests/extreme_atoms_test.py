"""Tests for the between-stage extreme-coordinate scanner."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from gbsa_pipeline.md_diagnostics import find_extreme_atoms

if TYPE_CHECKING:
    from pathlib import Path


def _gro_line(res_num: int, res_name: str, atom_name: str, idx: int, x: float, y: float, z: float) -> str:
    return f"{res_num:5d}{res_name:<5s}{atom_name:>5s}{idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


_GRO = "\n".join(
    [
        "extreme test",
        "    2",
        _gro_line(1, "ALA", "CA", 1, 1.000, 1.000, 1.000),
        _gro_line(2, "SOL", "OW", 2, 1.000, -12.500, 1.000),
        "   5.00000   5.00000   5.00000",
        "",
    ]
)


def test_find_extreme_atoms_flags_out_of_range_coordinates(tmp_path: Path) -> None:
    gro = tmp_path / "system.gro"
    gro.write_text(_GRO, encoding="utf-8")

    extreme = find_extreme_atoms(gro, threshold_nm=10.0)

    assert len(extreme) == 1
    atom_id, res_name, atom_name, x, y, z = extreme[0]
    assert (atom_id, res_name, atom_name) == (2, "SOL", "OW")
    assert (x, y, z) == pytest.approx((1.0, -12.5, 1.0))


def test_find_extreme_atoms_missing_file_returns_empty(tmp_path: Path) -> None:
    assert find_extreme_atoms(tmp_path / "missing.gro") == []
