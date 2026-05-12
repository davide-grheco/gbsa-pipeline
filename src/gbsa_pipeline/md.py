"""Small BioSimSpace MD protocol helpers for gbsa-pipeline.

This module currently contains the first isolated BioSimSpace building blocks
for the MD stage. The helpers assume that the caller already provides a valid
parametrized BioSimSpace/Sire system, because parametrization, solvation, file
conversion, and workflow orchestration belong to later commits. Keeping these
functions small makes the first MD protocol changes easy to review and easy to
replace if the final combined MD stage needs a different internal structure.
The functions create normal BioSimSpace protocol objects first and can then
apply explicit GROMACS parameter overrides to the generated MDP config when
tighter control over test-stage settings is required.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import BioSimSpace as BSS

from gbsa_pipeline.change_defaults import GromacsParams
from gbsa_pipeline.change_params import set_mdp_key

if TYPE_CHECKING:
    from pathlib import Path

    import Any
    import Mapping
    import sire


def _tail_text_file(path: Path, max_lines: int = 80) -> str:
    """Return the final lines of a process text file.

    GROMACS failures are usually diagnosed from the end of ``gromacs.out``,
    ``gromacs.log``, or the generated ``.mdp`` file. This helper keeps the
    failure-report formatting local to the MD helper module and avoids adding a
    broader logging abstraction. Missing files are represented explicitly so a
    failed BioSimSpace setup still gives useful context. Files are decoded with
    replacement enabled because external tools can occasionally write unusual
    characters to log files.
    """
    if not path.exists():
        return f"<missing: {path}>"

    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    if len(lines) <= max_lines:
        return text

    return "\n".join(lines[-max_lines:])


def _format_gromacs_failure_report(work_dir: Path | None, stage_name: str) -> str:
    """Build a compact diagnostic report for a failed GROMACS process.

    BioSimSpace can return ``None`` from ``getSystem`` after an ``mdrun``
    failure, which otherwise makes callers fail later with an unhelpful
    assertion error. This report includes the working directory and the tails of
    the most relevant files produced by the process. It intentionally does not
    try to parse or classify the chemistry problem because this helper only
    owns process execution. The goal is to expose the original GROMACS failure
    quickly enough for the integration test to be actionable.
    """
    if work_dir is None:
        return (
            f"{stage_name} finished without a readable output system. "
            "No work_dir was supplied, so no process logs can be summarized."
        )

    report_parts = [
        f"{stage_name} finished without a readable output system.",
        f"Work directory: {work_dir}",
    ]

    for filename in (
        "gromacs.out",
        "gromacs.log",
        "gromacs.err",
        "gromacs.mdp",
        "gromacs.out.mdp",
    ):
        path = work_dir / filename
        if path.exists():
            report_parts.append(f"\n--- tail: {path} ---\n{_tail_text_file(path)}")

    lincs_pdbs = sorted(work_dir.glob("step*.pdb"))
    if lincs_pdbs:
        report_parts.append("\nLINCS diagnostic PDB files:")
        report_parts.extend(str(path) for path in lincs_pdbs[-10:])

    return "\n".join(report_parts)


def _normalise_mdp_key(key: str) -> str:
    """Return a comparison-safe GROMACS MDP key.

    BioSimSpace and local parameter models do not always write MDP keys with
    identical spelling. For example, one source can write ``DispCorr`` while
    another writes ``dispcorr``; similarly, Python-facing code may use
    underscores while GROMACS normally uses hyphens. GROMACS treats those as the
    same logical parameter often enough that exact string matching is unsafe for
    config merging. This helper normalises only for comparison and does not
    decide how the final key should be written.
    """
    return key.strip().lower().replace("_", "-")


def _mdp_key_from_line(line: str) -> str | None:
    """Extract a normalised MDP key from one config line.

    GROMACS MDP files use ``key = value`` assignments and may include comments.
    This parser is intentionally small because it only needs to identify
    existing assignment lines before local overrides are applied. Blank lines,
    comment-only lines, and non-assignment lines are left untouched by returning
    ``None``. Inline comments are ignored so a line such as ``dt = 0.001 ;
    timestep`` is still recognised as the ``dt`` key.
    """
    body = line.split(";", 1)[0].strip()

    if not body or "=" not in body:
        return None

    key, _value = body.split("=", 1)
    return _normalise_mdp_key(key)


def _remove_existing_mdp_key(config: list[str], key: str) -> list[str]:
    """Remove existing assignments for one MDP key.

    This keeps BioSimSpace-generated configs valid when explicit overrides are
    applied afterwards. ``set_mdp_key`` can update exact key matches, but it may
    not catch spelling variants such as ``DispCorr`` versus ``dispcorr`` or
    underscore versus hyphen forms. Removing all equivalent assignments first
    guarantees that the later write produces a single final value for the key.
    Comments and unrelated config lines are preserved.
    """
    target_key = _normalise_mdp_key(key)

    return [line for line in config if _mdp_key_from_line(line) != target_key]


def _apply_gromacs_params_to_config(
    config: list[str],
    params: GromacsParams | Mapping[str, Any],
) -> list[str]:
    """Apply validated GROMACS parameters to a generated BioSimSpace config.

    BioSimSpace generates the initial MDP from the selected protocol object,
    which keeps stage-specific setup such as restraint topology generation tied
    to the normal BSS code path. This helper then updates or appends individual
    MDP keys using the serialized ``GromacsParams`` mapping. Before each key is
    written, any existing equivalent assignment is removed using a
    case-insensitive and hyphen/underscore-tolerant comparison. This prevents
    GROMACS errors such as ``Parameter "dispcorr" doubly defined`` when BSS and
    the local parameter model spell the same MDP key differently.
    """
    final_params: GromacsParams
    if isinstance(params, GromacsParams):
        final_params = params
    else:
        GromacsParams.from_mapping(params)

    updated_config = list(config)

    for key, value in final_params.to_mapping().items():
        updated_config = _remove_existing_mdp_key(updated_config, key)
        set_mdp_key(updated_config, key, value, inplace=True)

    return updated_config


def _run_bss_protocol(
    system: sire.System,
    protocol: Any,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
    *,
    ignore_warnings: bool = True,
    stage_name: str = "GROMACS stage",
) -> sire.System:
    """Run a BioSimSpace GROMACS protocol with optional MDP overrides.

    The process is created from a normal BioSimSpace protocol first, rather than
    from a standalone custom MDP. This matters for restrained stages because BSS
    can prepare the topology and reference-coordinate machinery associated with
    the protocol. If ``params`` is provided, the generated MDP config is then
    modified key-by-key and passed back into the process before execution. A
    failed process raises a diagnostic error with log tails instead of returning
    ``None`` to the caller.
    """
    kwargs: dict[str, object] = {"ignore_warnings": ignore_warnings}
    if work_dir:
        kwargs["work_dir"] = str(work_dir)

    process = BSS.Process.Gromacs(
        protocol=protocol,
        system=system,
        **kwargs,
    )

    if params is not None:
        config = _apply_gromacs_params_to_config(
            config=list(process.getConfig()),
            params=params,
        )
        process.setConfig(config)

    process.start()
    process.wait()

    result = process.getSystem(block=True)
    if result is None:
        raise RuntimeError(_format_gromacs_failure_report(work_dir, stage_name))

    return result


def run_minimization(
    system: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
    *,
    ignore_warnings: bool = True,
) -> sire.System:
    """Run a BioSimSpace energy minimization.

    This helper is the first small MD-stage building block and only handles
    minimization of an already prepared BioSimSpace/Sire system. The input
    ``system`` is expected to be parametrized already, because this function
    does not assign force-field parameters, solvate, neutralize, or prepare
    topology files. A normal ``BSS.Protocol.Minimisation`` object is always
    created first so BioSimSpace owns the base GROMACS setup. When ``params`` is
    provided, the generated MDP config is modified before execution so staged
    workflows can still request explicit steepest-descent or conjugate-gradient
    settings.
    """
    minimization_protocol = BSS.Protocol.Minimisation()

    return _run_bss_protocol(
        system=system,
        protocol=minimization_protocol,
        work_dir=work_dir,
        params=params,
        ignore_warnings=ignore_warnings,
        stage_name="GROMACS minimization",
    )


def run_heating(
    simulation_time: BSS.Types.Time,
    minimized: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
    temperature_start: BSS.Types.Temperature = 50 * BSS.Units.Temperature.kelvin,
    temperature_end: BSS.Types.Temperature = 300 * BSS.Units.Temperature.kelvin,
    restraint: str | None = "backbone",
    *,
    ignore_warnings: bool = True,
) -> sire.System:
    """Run a BioSimSpace/GROMACS restrained or unrestrained NVT stage.

    This helper creates an equilibration protocol for an already minimized or
    equilibrated system. If no explicit parameters are supplied, the protocol
    uses a 1 fs timestep and the requested temperature ramp. If explicit
    ``GromacsParams`` are supplied, BioSimSpace still creates the base protocol
    first and the generated MDP is then overwritten with the requested settings.
    The ``restraint`` argument is retained so tests can create a restrained BSS
    protocol before applying exact GROMACS parameters, which is safer than using
    a standalone custom MDP when position restraints are needed.
    """
    if params is None:
        heating_protocol = BSS.Protocol.Equilibration(
            timestep=1 * BSS.Units.Time.femtosecond,
            runtime=simulation_time,
            temperature_start=temperature_start,
            temperature_end=temperature_end,
            restraint=restraint,
        )
    else:
        heating_protocol = BSS.Protocol.Equilibration(
            timestep=1 * BSS.Units.Time.femtosecond,
            runtime=simulation_time,
            temperature=temperature_end,
            restraint=restraint,
        )

    return _run_bss_protocol(
        system=minimized,
        protocol=heating_protocol,
        work_dir=work_dir,
        params=params,
        ignore_warnings=ignore_warnings,
        stage_name="GROMACS NVT heating",
    )


def run_npt_equilibration(
    simulation_time: BSS.Types.Time,
    heated: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
    restraint: str | None = "backbone",
    *,
    ignore_warnings: bool = True,
) -> sire.System:
    """Run a BioSimSpace/GROMACS NPT equilibration procedure.

    This helper creates a pressure-equilibration stage for an already heated
    system. BioSimSpace creates the base NPT protocol first, including the
    restraint setup when requested. Explicit ``GromacsParams`` can then modify
    the generated MDP config before execution so the integration test can use
    exact short-stage settings without bypassing the normal BSS protocol path.
    The default timestep remains conservative at 1 fs because these stages are
    used directly after system preparation in the inspectable integration test.
    """
    equilibration_protocol = BSS.Protocol.Equilibration(
        timestep=1 * BSS.Units.Time.femtosecond,
        runtime=simulation_time,
        temperature=300 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
        restraint=restraint,
    )

    return _run_bss_protocol(
        system=heated,
        protocol=equilibration_protocol,
        work_dir=work_dir,
        params=params,
        ignore_warnings=ignore_warnings,
        stage_name="GROMACS NPT equilibration",
    )


def run_production(
    simulation_time: BSS.Types.Time,
    equilibrated: sire.System,
    work_dir: Path | None = None,
    params: GromacsParams | Mapping[str, Any] | None = None,
    *,
    ignore_warnings: bool = True,
) -> sire.System:
    """Run a BioSimSpace production MD procedure.

    This helper runs an isolated production step for an already equilibrated
    BioSimSpace/Sire system. A normal ``BSS.Protocol.Production`` object is
    created first so BSS controls the base setup. If explicit parameters are
    supplied, they are applied to the generated MDP config before the process is
    started. The function does not collect native output artefacts yet because
    stable output bundling belongs to a later pipeline orchestration change.
    """
    production_protocol = BSS.Protocol.Production(
        runtime=simulation_time,
        temperature=300 * BSS.Units.Temperature.kelvin,
        pressure=1 * BSS.Units.Pressure.atm,
    )

    return _run_bss_protocol(
        system=equilibrated,
        protocol=production_protocol,
        work_dir=work_dir,
        params=params,
        ignore_warnings=ignore_warnings,
        stage_name="GROMACS production",
    )
