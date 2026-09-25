"""Shared ParmEd I/O helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

    import parmed as pmd


def export_parmed_gromacs(structure: pmd.Structure, work_dir: Path, stem: str = "complex") -> tuple[Path, Path]:
    """Write ``structure`` as GROMACS ``<stem>.gro``/``<stem>.top`` in ``work_dir``.

    Existing files are overwritten. Returns ``(gro_file, top_file)``.
    """
    gro_file = work_dir / f"{stem}.gro"
    top_file = work_dir / f"{stem}.top"
    structure.save(str(top_file), format="gromacs", overwrite=True)
    structure.save(str(gro_file), overwrite=True)
    return gro_file, top_file
