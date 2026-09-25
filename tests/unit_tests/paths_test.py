"""Tests for the shared filesystem helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from gbsa_pipeline._paths import require_file, resolve_work_dir

if TYPE_CHECKING:
    from pathlib import Path


def test_require_file_returns_resolved_path(tmp_path: Path) -> None:
    target = tmp_path / "input.gro"
    target.write_text("data")

    assert require_file(target) == target.resolve()


def test_require_file_missing_raises_with_label(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError, match="GROMACS coordinate file not found"):
        require_file(tmp_path / "missing.gro", "GROMACS coordinate file")


def test_require_file_directory_raises_value_error(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="path is not a file"):
        require_file(tmp_path, "Input")


def test_resolve_work_dir_creates_explicit_directory(tmp_path: Path) -> None:
    explicit = tmp_path / "nested" / "work"

    assert resolve_work_dir(explicit, prefix="gbsa_test_") == explicit
    assert explicit.is_dir()


def test_resolve_work_dir_creates_temporary_directory() -> None:
    work_dir = resolve_work_dir(None, prefix="gbsa_test_")

    assert work_dir.is_dir()
    assert work_dir.name.startswith("gbsa_test_")
