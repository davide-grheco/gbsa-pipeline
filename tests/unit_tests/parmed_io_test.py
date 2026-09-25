"""Tests for the shared ParmEd I/O helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING
from unittest.mock import MagicMock

from gbsa_pipeline._parmed_io import export_parmed_gromacs

if TYPE_CHECKING:
    from pathlib import Path


def test_export_parmed_gromacs_saves_both_files_with_overwrite(tmp_path: Path) -> None:
    structure = MagicMock()

    gro_file, top_file = export_parmed_gromacs(structure, tmp_path)

    assert gro_file == tmp_path / "complex.gro"
    assert top_file == tmp_path / "complex.top"
    structure.save.assert_any_call(str(top_file), format="gromacs", overwrite=True)
    structure.save.assert_any_call(str(gro_file), overwrite=True)


def test_export_parmed_gromacs_honours_custom_stem(tmp_path: Path) -> None:
    gro_file, top_file = export_parmed_gromacs(MagicMock(), tmp_path, stem="system")

    assert gro_file.name == "system.gro"
    assert top_file.name == "system.top"
