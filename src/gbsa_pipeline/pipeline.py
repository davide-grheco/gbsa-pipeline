"""Functional pipeline runner — orchestrates all MD simulation stages."""

from __future__ import annotations

import logging
import time
from typing import TYPE_CHECKING, Any, Callable, TypeVar

import BioSimSpace as BSS

from gbsa_pipeline.md import (
    npt_barostat_overrides,
    remove_clashing_solvent_waters,
    run_heating,
    run_minimization,
    run_npt_equilibration,
    run_production,
    run_solvent_relaxation,
)
from gbsa_pipeline.md_io import save_bss_system_to_gromacs
from gbsa_pipeline.membrane import merge_ligand_into_system, parametrize_ligand_only
from gbsa_pipeline.parametrization import parametrize
from gbsa_pipeline.solvation_box import solvate_membrane
from gbsa_pipeline.solvation_bss import solvate_bss

if TYPE_CHECKING:
    from pathlib import Path

    from gbsa_pipeline.config import MembraneSystemConfig, RunConfig
    from gbsa_pipeline.parametrization import ParametrisedComplex

logger = logging.getLogger(__name__)

_T = TypeVar("_T")


# ---------------------------------------------------------------------------
# Stage runner
# ---------------------------------------------------------------------------


def _run_stage(name: str, fn: Callable[[], _T]) -> _T:
    """Run a named pipeline stage with logging and elapsed-time reporting."""
    logger.info("  [%s] starting …", name)
    t0 = time.perf_counter()
    try:
        result = fn()
    except Exception:
        elapsed = time.perf_counter() - t0
        logger.exception("  [%s] failed after %.1f s", name, elapsed)
        raise
    elapsed = time.perf_counter() - t0
    logger.info("  [%s] completed in %.1f s", name, elapsed)
    return result


def _run_md_stage(
    title: str,
    name: str,
    label: str,
    output_dir: Path,
    fn: Callable[[Path], Any],
) -> Any:
    """Run one MD stage: log banner, mkdir, run, save gro/top, return system."""
    logger.info("─── %s ───", title)
    stage_dir = output_dir / label
    stage_dir.mkdir(parents=True, exist_ok=True)
    system = _run_stage(name, lambda: fn(stage_dir))
    save_bss_system_to_gromacs(system, stage_dir / "system")
    logger.info("  Saved → %s/system.gro / .top", label)
    return system


# ---------------------------------------------------------------------------
# Individual stage helpers — pure functions over validated inputs
# ---------------------------------------------------------------------------


def _stage_parametrize(config: RunConfig, stage_dir: Path) -> ParametrisedComplex:
    """Assign force field parameters to the protein-ligand complex."""
    logger.info(
        "  protein_ff=%s  ligand_ff=%s  charge_method=%s",
        config.forcefield.protein_ff,
        config.forcefield.ligand_ff,
        config.forcefield.charge_method,
    )
    return parametrize(config.to_parametrization_input(stage_dir))


def _stage_solvate(
    config: RunConfig,
    parametrized: ParametrisedComplex,
    stage_dir: Path,
) -> Any:
    """Solvate with BSS.Solvent (gmx solvate + gmx genion) and return loaded BSS system."""
    sol = config.solvation
    box_desc = f"padding={sol.padding} nm" if sol.padding is not None else f"box_size={sol.box_size} nm"
    logger.info(
        "  water_model=%s  shape=%s  %s  ion_conc=%s mol/L",
        sol.water_model,
        sol.shape,
        box_desc,
        sol.ion_concentration,
    )
    solvated = solvate_bss(
        parametrized=parametrized,
        params=sol,
        output_gro=stage_dir / "solvated.gro",
        output_top=stage_dir / "solvated.top",
    )
    logger.info("  Saved → %s / %s", solvated.gro_file.name, solvated.top_file.name)

    logger.info("  Loading solvated system into BSS …")
    system = solvated.load_bss()
    logger.info("  Loaded %d molecules (%d atoms)", system.nMolecules(), system.nAtoms())
    return system


def _stage_parametrize_membrane(config: RunConfig, stage_dir: Path) -> Any:
    """Load a pre-built [membrane] system and merge with parametrised ligand."""
    membrane: MembraneSystemConfig | None = config.membrane
    if membrane is None:
        raise ValueError("stage parametrize_membrane requires membrane to be set.")

    logger.info(
        "gro_file=%s  top_file=%s ligand=%s net_charge=%s",
        membrane.gro_file.name,
        membrane.top_file.name,
        membrane.ligand.name,
        membrane.net_charge,
    )
    system = BSS.IO.readMolecules(
        [str(membrane.gro_file), str(membrane.top_file)],
        make_whole=True,
    )
    ligand = parametrize_ligand_only(
        membrane.ligand,
        net_charge=membrane.net_charge,
        work_dir=stage_dir,
    )
    return merge_ligand_into_system(system, ligand)


def _stage_solvate_membrane(config: RunConfig, system: Any, stage_dir: Path) -> Any:
    """Solvate a membrane system, or pass it through unchanged if already solvated."""
    membrane = config.membrane
    if membrane is None:
        raise ValueError("_stage_solvate_membrane requires [membrane] to be set.")

    if not membrane.solvate:
        logger.info("solvate=False - system is already solvated, skipping.")
        return system

    logger.info("z_padding_nm=%.2f water_models=%s.", membrane.z_padding_nm, config.solvation.water_model)
    return solvate_membrane(
        system=system,
        params=config.solvation,
        z_padding_nm=membrane.z_padding_nm,
        work_dir=stage_dir,
    )


def _stage_minimize_sd(config: RunConfig, system: Any, stage_dir: Path) -> Any:
    """Steepest-descent energy minimization."""
    logger.info(
        "  nsteps=%d  emtol=%.1f kJ/mol/nm",
        config.minimization.nsteps,
        config.minimization.emtol,
    )
    return run_minimization(
        system,
        work_dir=stage_dir,
        params={
            "nsteps": config.minimization.nsteps,
            "emtol": config.minimization.emtol,
        },
    )


def _stage_minimize_cg(system: Any, stage_dir: Path) -> Any:
    """Conjugate-gradient energy minimization."""
    return run_minimization(system, work_dir=stage_dir, params={"integrator": "cg"})


def _stage_nvt_restrained(config: RunConfig, system: Any, stage_dir: Path) -> Any:
    """Water clash removal → short solvent relax → NVT heating 50→300 K with backbone restraints."""
    logger.info("  NVT heating over %.1f ps", config.equilibration.simulation_time_ps)

    system = remove_clashing_solvent_waters(system, work_dir=stage_dir / "water_cleanup")
    system = run_solvent_relaxation(system, work_dir=stage_dir / "solvent_relax")

    equil_time = config.equilibration.simulation_time_ps * BSS.Units.Time.picosecond
    return run_heating(
        equil_time,
        system,
        work_dir=stage_dir,
        temperature_start=50 * BSS.Units.Temperature.kelvin,
        temperature_end=300 * BSS.Units.Temperature.kelvin,
        restraint="backbone",
    )


def _stage_npt(config: RunConfig, system: Any, stage_dir: Path, *, restraint: str | None = None) -> Any:
    """NPT equilibration, optionally with backbone restraints.

    Uses the same barostat as the [md] section so a memprot configured with
    pcouple = semiisotropic gets consistent, not isotropic, values during
    equilibration.
    """
    logger.info(
        "  %.1f ps  restraint=%s  pcoupltype=%s",
        config.npt_equilibration.simulation_time_ps,
        restraint or "none",
        config.md.pcoupltype,
    )

    npt_time = config.npt_equilibration.simulation_time_ps * BSS.Units.Time.picosecond
    return run_npt_equilibration(
        npt_time,
        system,
        work_dir=stage_dir,
        restraint=restraint,
        params=npt_barostat_overrides(config.md),
    )


def _stage_production(config: RunConfig, system: Any, stage_dir: Path) -> Any:
    """Production MD."""
    sim_time = config.md.nsteps * config.md.dt * BSS.Units.Time.picosecond
    logger.info(
        "  nsteps=%d  dt=%s ps  sim_time=%.1f ps  tcoupl=%s  pcoupl=%s",
        config.md.nsteps,
        config.md.dt,
        config.md.nsteps * config.md.dt,
        config.md.tcoupl,
        config.md.pcoupl,
    )
    return run_production(sim_time, system, work_dir=stage_dir, params=config.md)


# ---------------------------------------------------------------------------
# Pipeline entry point
# ---------------------------------------------------------------------------


def run_pipeline(config: RunConfig, output_dir: Path) -> None:
    """Run the full GBSA pipeline from a validated :class:`~gbsa_pipeline.config.RunConfig`.

    Stages (each writes output to a numbered subdirectory):

    1. **Parametrize** — assign force field parameters to protein + ligand.
    2. **Solvate** — add water box and counter-ions via BSS.Solvent.
    3. **SD Minimization** — steepest-descent energy minimization.
    4. **CG Minimization** — conjugate-gradient energy minimization.
    5. **NVT Restrained** — water cleanup, solvent relax, NVT heating 50→300 K.
    6. **NPT Restrained** — NPT equilibration with backbone restraints.
    7. **NPT** — NPT equilibration without restraints.
    8. **Production MD** — NpT simulation driven by ``[md]`` section params.

    Parameters
    ----------
    config:
        Validated run configuration (usually loaded via
        :meth:`~gbsa_pipeline.config.RunConfig.from_toml`).
    output_dir:
        Root directory for all output. Created if it does not exist.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    _log_config(config, output_dir)

    if config.membrane is not None:
        # Stage 1: Parametrize ligand + merge into pre-built membrane system
        logger.info("─── Stage 1/8: Ligand parametrization + membrane merge ───")
        param_dir = output_dir / "01_parametrize"
        param_dir.mkdir(parents=True, exist_ok=True)
        system = _run_stage(
            "parametrize_membrane",
            lambda: _stage_parametrize_membrane(config, param_dir),
        )

        # Stage 2: Solvate (membrane-aware, or skip if already solvated)
        logger.info("─── Stage 2/8: Membrane solvation ───")
        sol_dir = output_dir / "02_solvated"
        sol_dir.mkdir(parents=True, exist_ok=True)
        system = _run_stage(
            "solvate_membrane",
            lambda: _stage_solvate_membrane(config, system, sol_dir),
        )
    else:
        # Stage 1: Parametrize
        logger.info("─── Stage 1/8: Parametrization ───")
        param_dir = output_dir / "01_parametrize"
        parametrized = _run_stage("parametrize", lambda: _stage_parametrize(config, param_dir))
        logger.info("  Done → %s, %s", parametrized.gro_file.name, parametrized.top_file.name)

        # Stage 2: Solvate
        logger.info("─── Stage 2/8: Solvation ───")
        sol_dir = output_dir / "02_solvated"
        system = _run_stage("solvation", lambda: _stage_solvate(config, parametrized, sol_dir))

    system = _run_md_stage(
        "Stage 3/8: SD Minimization",
        "sd_minimization",
        "03_sd",
        output_dir,
        lambda d: _stage_minimize_sd(config, system, d),
    )
    system = _run_md_stage(
        "Stage 4/8: CG Minimization",
        "cg_minimization",
        "04_cg",
        output_dir,
        lambda d: _stage_minimize_cg(system, d),
    )
    system = _run_md_stage(
        "Stage 5/8: NVT Restrained Heating",
        "nvt_restrained",
        "05_nvt_res",
        output_dir,
        lambda d: _stage_nvt_restrained(config, system, d),
    )
    system = _run_md_stage(
        "Stage 6/8: NPT Restrained Equilibration",
        "npt_restrained",
        "06_npt_res",
        output_dir,
        lambda d: _stage_npt(config, system, d, restraint="backbone"),
    )
    system = _run_md_stage(
        "Stage 7/8: NPT Equilibration",
        "npt",
        "07_npt",
        output_dir,
        lambda d: _stage_npt(config, system, d),
    )
    system = _run_md_stage(
        "Stage 8/8: Production MD",
        "production_md",
        "08_production",
        output_dir,
        lambda d: _stage_production(config, system, d),
    )

    logger.info("Pipeline complete. Output written to %s", output_dir)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _log_config(config: RunConfig, output_dir: Path) -> None:
    """Write a JSON snapshot of the resolved config to ``output_dir/run_config.json``."""
    config_path = output_dir / "run_config.json"
    config_path.write_text(config.model_dump_json(indent=2))
    logger.info("Config written to %s", config_path)
