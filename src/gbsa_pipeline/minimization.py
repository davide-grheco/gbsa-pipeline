"""Preparing minimization protocol and running."""

from __future__ import annotations

from typing import TYPE_CHECKING

import BioSimSpace as BSS

if TYPE_CHECKING:
    from pathlib import Path

    import sire


def run_minimization(nsteps: int, system: sire.System, work_dir: Path | None = None) -> BSS.System:
    """Run energy minimization via GROMACS.

    Args:
        nsteps: Maximum number of minimization steps.
        system: Solvated system to minimize.
        work_dir: Optional working directory for GROMACS files.

    Returns:
        Minimized BioSimSpace system.
    """
    protocol = BSS.Protocol.Minimisation(steps=nsteps)
    kwargs = {"work_dir": str(work_dir)} if work_dir else {}
    process = BSS.Process.Gromacs(system, protocol, name="min", **kwargs)
    process.start()
    return process.getSystem(block=True)
