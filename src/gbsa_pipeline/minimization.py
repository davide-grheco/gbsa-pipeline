"""Preparing minimization protocol and running."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import BioSimSpace as BSS

if TYPE_CHECKING:
    from pathlib import Path

    import sire


def run_minimization(
    nsteps: int,
    system: sire.System,
    work_dir: Path | None = None,
    **kwargs: Any,
) -> BSS.System:
    """Run energy minimization via GROMACS.

    Args:
        nsteps: Maximum number of minimization steps.
        system: Solvated system to minimize.
        work_dir: Optional working directory for GROMACS files.
        **kwargs: Extra keyword arguments forwarded to ``BSS.Process.Gromacs``.

    Returns:
        Minimized BioSimSpace system.
    """
    protocol = BSS.Protocol.Minimisation(steps=nsteps)

    gromacs_kwargs: dict[str, Any] = dict(kwargs)
    if work_dir is not None:
        gromacs_kwargs["work_dir"] = str(work_dir)

    process = BSS.Process.Gromacs(system, protocol, name="min", **gromacs_kwargs)
    process.start()
    return process.getSystem(block=True)
