# src/gbsa_pipeline/md.py

"""Small BioSimSpace MD protocol helpers for gbsa-pipeline.

This module currently contains the first isolated BioSimSpace building blocks
for the MD stage. The helpers assume that the caller already provides a valid
parametrized BioSimSpace/Sire system, because parametrization, solvation, file
conversion, and workflow orchestration belong to later commits. Keeping these
functions small makes the first MD protocol changes easy to review and easy to
replace if the final combined MD stage needs a different internal structure.
The functions use BioSimSpace protocol objects directly by default, while the
staged workflow can pass explicit GROMACS parameter blocks through the existing
helpers when tighter control over MDP settings is required.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import BioSimSpace as BSS

from gbsa_pipeline.change_defaults import GromacsParams, run_gro_custom

if TYPE_CHECKING:
    from collections.abc import Mapping
    from pathlib import Path

    import sire


def run_minimization(
    system: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
) -> sire.System:
    """Run a BioSimSpace/GROMACS energy minimization.

    This helper is the first small MD-stage building block and only handles
    minimization of an already prepared BioSimSpace/Sire system. The input
    ``system`` is expected to be parametrized already, because this function
    does not assign force-field parameters, solvate, neutralize, or prepare
    topology files. When ``params`` is not provided, the helper keeps the
    original conservative behaviour and uses BioSimSpace's standard
    ``Minimisation`` protocol. When ``params`` is provided, it is passed through
    the custom GROMACS protocol helper so staged workflows can use explicit MDP
    blocks such as steepest-descent and conjugate-gradient minimization.
    """
    if params is not None:
        minimized, _protocol = run_gro_custom(
            parameters=params,
            system=system,
            work_dir=work_dir,
        )
        return minimized

    minimization_protocol = BSS.Protocol.Minimisation()

    kwargs = {"work_dir": str(work_dir)} if work_dir else {}
    minimization_process = BSS.Process.Gromacs(
        protocol=minimization_protocol,
        system=system,
        **kwargs,
    )

    minimization_process.start()
    minimization_process.wait()

    minimized = minimization_process.getSystem(block=True)

    return minimized


def run_heating(
    simulation_time: BSS.Types.Time,
    minimized: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
) -> sire.System:
    """Run a BioSimSpace/GROMACS NVT heating procedure.

    This helper creates a heating stage for an already minimized system. The
    default path uses ``BSS.Protocol.Equilibration`` from 0 K to 300 K with
    backbone restraints, preserving the original BioSimSpace-driven behaviour.
    When ``params`` is provided, the helper routes the system through the custom
    GROMACS protocol runner instead, so callers can provide exact MDP settings
    such as annealing, velocity generation, LINCS options, and ``define =
    -DPOSRES``. The ``simulation_time`` argument remains part of the public
    signature for the default BioSimSpace path, but the custom parameter path
    uses the runtime encoded by ``dt`` and ``nsteps`` in the supplied MDP block.
    """
    if params is not None:
        heated, _protocol = run_gro_custom(
            parameters=params,
            system=minimized,
            work_dir=work_dir,
        )
        return heated

    heating_protocol = BSS.Protocol.Equilibration(
        runtime=simulation_time,
        temperature_start=0 * BSS.Units.Temperature.kelvin,
        temperature_end=300 * BSS.Units.Temperature.kelvin,
        restraint="backbone",
    )

    kwargs = {"work_dir": str(work_dir)} if work_dir else {}
    heating_process = BSS.Process.Gromacs(
        protocol=heating_protocol,
        system=minimized,
        **kwargs,
    )

    heating_process.start()
    heating_process.wait()

    equilibrated = heating_process.getSystem(block=True)

    return equilibrated


def run_npt_equilibration(
    simulation_time: BSS.Types.Time,
    heated: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
) -> sire.System:
    """Run a BioSimSpace/GROMACS NPT equilibration procedure.

    This helper creates an NPT equilibration stage for an already heated system.
    The default path uses ``BSS.Protocol.Equilibration`` at 300 K and 1 atm with
    backbone restraints, preserving the original BioSimSpace-driven behaviour.
    When ``params`` is provided, the helper routes the system through the custom
    GROMACS protocol runner instead, so callers can provide exact MDP settings
    such as C-rescale pressure coupling, reference-coordinate scaling, LINCS
    options, and ``define = -DPOSRES``. The ``simulation_time`` argument remains
    part of the public signature for the default BioSimSpace path, but the
    custom parameter path uses the runtime encoded by ``dt`` and ``nsteps`` in
    the supplied MDP block.
    """
    if params is not None:
        equilibrated, _protocol = run_gro_custom(
            parameters=params,
            system=heated,
            work_dir=work_dir,
        )
        return equilibrated

    equilibration_protocol = BSS.Protocol.Equilibration(
        runtime=simulation_time,
        temperature=300 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
        restraint="backbone",
    )

    kwargs = {"work_dir": str(work_dir)} if work_dir else {}
    equilibration_process = BSS.Process.Gromacs(
        protocol=equilibration_protocol,
        system=heated,
        **kwargs,
    )

    equilibration_process.start()
    equilibration_process.wait()

    equilibrated = equilibration_process.getSystem(block=True)

    return equilibrated


def run_production(
    simulation_time: BSS.Types.Time,
    equilibrated: sire.System,
    work_dir: Path | None = None,
) -> sire.System:
    """Run a BioSimSpace production MD procedure.

    This helper runs the first isolated production step for an already
    equilibrated BioSimSpace/Sire system. The ``simulation_time`` argument
    controls the length of the production stage and should be passed as a
    BioSimSpace time object by the caller. The ``equilibrated`` system is
    assumed to come from previous minimization, heating, and equilibration
    stages, because this function does not check temperature, pressure, density,
    or prior convergence. The helper intentionally uses BioSimSpace's standard
    production protocol without adding output collection or custom GROMACS
    ``.mdp`` handling yet, because those concerns belong to the later pipeline
    orchestration layer.
    """
    production_protocol = BSS.Protocol.Production(
        runtime=simulation_time,
        temperature=300 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
    )

    kwargs = {"work_dir": str(work_dir)} if work_dir else {}
    production_process = BSS.Process.Gromacs(
        protocol=production_protocol,
        system=equilibrated,
        **kwargs,
    )

    production_process.start()
    production_process.wait()

    production = production_process.getSystem(block=True)

    return production
