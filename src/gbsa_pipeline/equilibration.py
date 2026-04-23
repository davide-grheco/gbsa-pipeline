"""Preparing and running a heating procedure."""

from __future__ import annotations

from typing import TYPE_CHECKING

import BioSimSpace as BSS

if TYPE_CHECKING:
    from pathlib import Path

    import sire

    from gbsa_pipeline.config import EquilibrationConfig


def run_heating(
    system: sire.System,
    config: EquilibrationConfig,
    work_dir: Path | None = None,
) -> sire.System:
    """Run NVT heating from 0 K to 300 K via GROMACS."""
    simulation_time = config.simulation_time_ps * BSS.Units.Time.picosecond
    heating_protocol = BSS.Protocol.Equilibration(
        runtime=simulation_time,
        temperature_start=0 * BSS.Units.Temperature.kelvin,
        temperature_end=300 * BSS.Units.Temperature.kelvin,
        restraint="backbone",
    )
    kwargs = {"work_dir": str(work_dir)} if work_dir else {}
    heating_process = BSS.Process.Gromacs(protocol=heating_protocol, system=system, **kwargs)

    heating_process.start()
    heating_process.wait()
    return heating_process.getSystem(block=True)
