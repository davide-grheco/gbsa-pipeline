"""Tests for GRO clash-removal I/O."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from gbsa_pipeline._gro_io import _parse_gro, _write_cleaned_gro

if TYPE_CHECKING:
    from pathlib import Path


def _gro_line(res_num: int, res_name: str, atom_name: str, idx: int, x: float, y: float, z: float) -> str:
    return f"{res_num:5d}{res_name:<5s}{atom_name:>5s}{idx:5d}{x:8.3f}{y:8.3f}{z:8.3f}"


_GRO = "\n".join(
    [
        "test system",
        "    3",
        _gro_line(1, "LIG", "C1", 1, 1.000, 1.000, 1.000),
        _gro_line(2, "SOL", "OW", 2, 1.100, 1.000, 1.000),  # 0.1 nm from solute → clash
        _gro_line(3, "SOL", "OW", 3, 3.000, 3.000, 3.000),  # far away → kept
        "   5.00000   5.00000   5.00000",
        "",
    ]
)


def test_write_cleaned_gro_removes_clashing_water_and_renumbers(tmp_path: Path) -> None:
    input_gro = tmp_path / "in.gro"
    output_gro = tmp_path / "out.gro"
    input_gro.write_text(_GRO, encoding="utf-8")

    removed = _write_cleaned_gro(input_gro, output_gro, cutoff_nm=0.25, water_resnames={"SOL"})

    assert removed == {"SOL": 1}
    atoms = _parse_gro(output_gro)
    assert [(a.atom_idx, a.res_num, a.res_name) for a in atoms] == [(1, 1, "LIG"), (2, 3, "SOL")]
    assert atoms[1].x == pytest.approx(3.0)


def test_write_cleaned_gro_without_clashes_keeps_everything(tmp_path: Path) -> None:
    input_gro = tmp_path / "in.gro"
    output_gro = tmp_path / "out.gro"
    input_gro.write_text(_GRO, encoding="utf-8")

    removed = _write_cleaned_gro(input_gro, output_gro, cutoff_nm=0.05, water_resnames={"SOL"})

    assert removed == {}
    assert len(_parse_gro(output_gro)) == 3
