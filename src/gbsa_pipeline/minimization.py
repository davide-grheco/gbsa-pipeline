"""Preparing minimization protocol and running."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import BioSimSpace as BSS

if TYPE_CHECKING:
    from pathlib import Path

    import sire

    from gbsa_pipeline.config import MinimizationConfig


def run_minimization(
    system: sire.System,
    config: MinimizationConfig,
    work_dir: Path | None = None,
    **kwargs: Any,
) -> BSS.System:
    """Run energy minimization via GROMACS.

    Args:
        system: Solvated system to minimize.
        config: Minimization settings (nsteps, emtol).
        work_dir: Optional working directory for GROMACS files.
        **kwargs: Extra keyword arguments forwarded to ``BSS.Process.Gromacs``.

    Returns:
        Minimized BioSimSpace system.
    """
    protocol = BSS.Protocol.Minimisation(steps=config.nsteps)

    gromacs_kwargs: dict[str, Any] = dict(kwargs)
    if work_dir is not None:
        gromacs_kwargs["work_dir"] = str(work_dir)

    process = BSS.Process.Gromacs(system, protocol, name="min", **gromacs_kwargs)
    process.start()
    return process.getSystem(block=True)
