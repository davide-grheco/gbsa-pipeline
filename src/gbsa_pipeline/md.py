"""Small BioSimSpace MD protocol helpers for gbsa-pipeline.

This module currently contains the first isolated BioSimSpace building blocks
for the MD stage. The helpers assume that the caller already provides a valid
parametrized BioSimSpace/Sire system, because parametrization, solvation, file
conversion, and workflow orchestration belong to later commits. Keeping these
functions small makes the first MD protocol changes easy to review and easy to
replace if the final combined MD stage needs a different internal structure.
The functions use BioSimSpace protocol objects directly instead of manually
writing GROMACS input files at this stage.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import BioSimSpace as BSS

if TYPE_CHECKING:
    from pathlib import Path

    import sire


def run_minimization(
    system: sire.System,
    work_dir: Path | None = None,
) -> sire.System:
    """Run a BioSimSpace energy minimization.

    This helper is the first small MD-stage building block and only handles
    minimization of an already prepared BioSimSpace/Sire system. The input
    ``system`` is expected to be parametrized already, because this function
    does not assign force-field parameters, solvate, neutralize, or prepare
    topology files. The optional ``work_dir`` is passed to the BioSimSpace
    GROMACS process so the generated intermediate files stay in a predictable
    directory when the caller wants that. The function intentionally uses
    BioSimSpace's standard ``Minimisation`` protocol rather than manually
    writing GROMACS ``.mdp`` files, because this first step should stay close
    to the BioSimSpace public API.
    """
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
) -> sire.System:
    """Run a BioSimSpace NVT heating procedure.

    This helper creates a BioSimSpace equilibration protocol that heats an
    already minimized system from 0 K to 300 K. The ``simulation_time`` argument
    controls the runtime of the heating stage and is expected to be a
    BioSimSpace time object, for example a value built with
    ``BSS.Units.Time.picosecond``. The ``minimized`` system is assumed to be the
    output of a previous minimization step, because this function does not check
    energy convergence or repair unstable structures. Backbone restraints are
    kept as the first conservative default, matching the current base protocol
    before a later restraint-selection helper is added.
    """
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
) -> sire.System:
    """Run a BioSimSpace NPT equilibration procedure.

    This helper creates a BioSimSpace equilibration protocol for an already
    heated system at 300 K and 1 atm. The ``simulation_time`` argument controls
    only the length of this equilibration stage and should be provided as a
    BioSimSpace time object by the caller. The ``heated`` system is assumed to
    be structurally stable enough to start pressure coupling, because this
    function does not run minimization, diagnose LINCS problems, or inspect
    density convergence. Backbone restraints remain enabled as the conservative
    first default so this helper matches the early heating protocol before a
    later configurable restraint helper is introduced.
    """
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


def run_staged_bss_workflow(
    system: sire.System,
    heating_time: BSS.Types.Time,
    npt_equilibration_time: BSS.Types.Time,
    work_dir: Path | None = None,
) -> sire.System:
    """Run the current staged BioSimSpace MD workflow.

    This function composes the small helpers currently available in this module:
    minimization, NVT heating, and NPT equilibration. The input ``system`` is
    expected to be an already parametrized BioSimSpace/Sire system because this
    workflow does not yet prepare files, solvate, or assign force-field
    parameters. The optional ``work_dir`` is used as a parent directory for the
    individual process directories when provided, so intermediate files stay
    separated by stage without adding a larger output model yet. The returned
    system is the final system produced by the NPT equilibration helper.
    """
    minimized = run_minimization(
        system=system,
        work_dir=work_dir / "minimization" if work_dir else None,
    )
    heated = run_heating(
        simulation_time=heating_time,
        minimized=minimized,
        work_dir=work_dir / "heating" if work_dir else None,
    )
    equilibrated = run_npt_equilibration(
        simulation_time=npt_equilibration_time,
        heated=heated,
        work_dir=work_dir / "npt_equilibration" if work_dir else None,
    )

    return equilibrated
