"""BioSimSpace/GROMACS IO helpers for workflow stage boundaries.

This module contains small file-based helpers for moving BioSimSpace/Sire
systems between independent workflow steps. The functions are intentionally
workflow-engine agnostic: they do not import grubicy, signac, row, or any
project-specific orchestration code. They only load and save GROMACS coordinate
and topology files through BioSimSpace. This keeps ``gbsa-pipeline`` usable as a
library while allowing external workflow layers to persist intermediate MD
states between separate processes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import BioSimSpace as BSS

if TYPE_CHECKING:
    from pathlib import Path

    import sire


def load_bss_system_from_gromacs(gro_file: Path, top_file: Path) -> sire.System:
    """Load a BioSimSpace/Sire system from GROMACS ``.gro`` and ``.top`` files.

    This helper is intended for workflow boundaries where one process has
    written a completed system and another process needs to continue from it.
    The caller must provide the coordinate file and the matching topology file,
    because BioSimSpace needs both to reconstruct the system for later MD
    stages. The function validates file existence before calling BioSimSpace so
    missing products fail with clear file paths instead of lower-level loader
    errors. It does not infer parent directories, stage names, or workflow
    metadata; those decisions belong to the external orchestration layer.
    """
    gro_file = gro_file.resolve()
    top_file = top_file.resolve()

    if not gro_file.exists():
        raise FileNotFoundError(f"GROMACS coordinate file not found: {gro_file}")

    if not gro_file.is_file():
        raise ValueError(f"GROMACS coordinate path is not a file: {gro_file}")

    if not top_file.exists():
        raise FileNotFoundError(f"GROMACS topology file not found: {top_file}")

    if not top_file.is_file():
        raise ValueError(f"GROMACS topology path is not a file: {top_file}")

    return BSS.IO.readMolecules([str(gro_file), str(top_file)])


def save_bss_system_to_gromacs(
    system: sire.System,
    output_prefix: Path,
) -> tuple[Path, Path]:
    """Save a BioSimSpace/Sire system as GROMACS ``.gro`` and ``.top`` files.

    Workflow engines usually run each stage in a separate Python process, so
    in-memory BioSimSpace objects cannot be passed directly between MD steps.
    This helper writes a coordinate file and a topology file using a shared
    output prefix and returns both resolved paths. The function creates the
    output directory before writing so callers can use fresh stage directories
    without preparing them separately. It does not collect trajectories,
    energies, checkpoints, or process logs; this helper only persists the system
    needed by the next stage.
    """
    output_prefix = output_prefix.resolve()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    gro_file = output_prefix.with_suffix(".gro")
    top_file = output_prefix.with_suffix(".top")

    BSS.IO.saveMolecules(str(gro_file), system, fileformat="gro87")
    BSS.IO.saveMolecules(str(top_file), system, fileformat="grotop")

    return gro_file, top_file
