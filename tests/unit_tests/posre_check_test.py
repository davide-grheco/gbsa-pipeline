"""Tests for the position-restraint consistency check."""

from __future__ import annotations

from typing import TYPE_CHECKING

from gbsa_pipeline.md_diagnostics import _parse_posre_indices, check_posre_consistency

if TYPE_CHECKING:
    from pathlib import Path


def _gro_line(res_num: int, res_name: str, atom_name: str, idx: int) -> str:
    return f"{res_num:5d}{res_name:<5s}{atom_name:>5s}{idx:5d}{1.0:8.3f}{1.0:8.3f}{1.0:8.3f}"


_GRO = "\n".join(
    [
        "posre test",
        "    3",
        _gro_line(1, "ALA", "CA", 1),
        _gro_line(1, "ALA", "CB", 2),
        _gro_line(2, "SOL", "OW", 3),
        "   5.00000   5.00000   5.00000",
        "",
    ]
)

_POSRE = """\
[ position_restraints ]
; atom  functype  fx  fy  fz
   1     1  1000 1000 1000
   3     1  1000 1000 1000
"""


def test_parse_posre_indices_reads_restraint_lines(tmp_path: Path) -> None:
    posre = tmp_path / "posre.itp"
    posre.write_text(_POSRE, encoding="utf-8")

    assert _parse_posre_indices(posre) == [1, 3]


def test_check_posre_consistency_flags_solvent_and_nonbackbone(tmp_path: Path) -> None:
    gro = tmp_path / "system.gro"
    posre = tmp_path / "posre.itp"
    gro.write_text(_GRO, encoding="utf-8")
    posre.write_text(_POSRE, encoding="utf-8")

    result = check_posre_consistency(gro, posre)

    assert not result.ok
    assert result.n_restrained == 2
    assert result.unexpected == [(3, "SOL", "OW")]  # backbone CA passes, restrained water fails
    assert result.missing_indices == []


def test_check_posre_consistency_passes_for_backbone_only(tmp_path: Path) -> None:
    gro = tmp_path / "system.gro"
    posre = tmp_path / "posre.itp"
    gro.write_text(_GRO, encoding="utf-8")
    posre.write_text("   1     1  1000 1000 1000\n", encoding="utf-8")

    result = check_posre_consistency(gro, posre)

    assert result.ok
    assert result.n_restrained == 1


def test_check_posre_consistency_flags_index_without_atom(tmp_path: Path) -> None:
    gro = tmp_path / "system.gro"
    posre = tmp_path / "posre.itp"
    gro.write_text(_GRO, encoding="utf-8")
    posre.write_text("  99     1  1000 1000 1000\n", encoding="utf-8")

    result = check_posre_consistency(gro, posre)

    assert not result.ok
    assert result.missing_indices == [99]
    assert result.unexpected == []


def test_check_posre_consistency_missing_gro_is_not_ok(tmp_path: Path) -> None:
    posre = tmp_path / "posre.itp"
    posre.write_text(_POSRE, encoding="utf-8")

    result = check_posre_consistency(tmp_path / "missing.gro", posre)

    assert not result.ok
    assert result.n_restrained == 0
    assert result.error is not None
