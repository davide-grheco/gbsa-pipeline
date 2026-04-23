"""BioSimSpace file I/O utilities shared across MD stages."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import BioSimSpace as BSS

if TYPE_CHECKING:
    from BioSimSpace._SireWrappers import System


def export_gromacs_top_gro(system: System, prefix: str) -> list[Path]:
    """Export GROMACS .gro and .top files from a BSS System."""
    out_gro = Path(f"{prefix}.gro")
    out_top = Path(f"{prefix}.top")
    BSS.IO.saveMolecules(str(out_gro), system, fileformat="gro87")
    BSS.IO.saveMolecules(str(out_top), system, fileformat="grotop")
    return [out_gro, out_top]


def _run_bss_process(process: BSS.Process.Gromacs) -> System:
    """Start a BSS Gromacs process, wait for completion, and return the system."""
    process.start()
    process.wait()
    return process.getSystem(block=True)
