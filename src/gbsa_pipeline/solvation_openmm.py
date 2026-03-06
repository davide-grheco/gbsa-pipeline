"""OpenMM/ParmEd solvation that bypasses BioSimSpace IO."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import parmed as pmd
from openmm import Vec3
from openmm import unit as mm_unit
from openmm.app import ForceField, Modeller, NoCutoff

from gbsa_pipeline.solvation_box import BoxShape, SolvationParams, WaterModel

if TYPE_CHECKING:
    from pathlib import Path

    from gbsa_pipeline.parametrization import ParametrisedComplex

logger = logging.getLogger(__name__)

# OpenMM force-field XML files for each water model.
_WATER_XML: dict[WaterModel, str] = {
    WaterModel.TIP3P: "amber14/tip3p.xml",
    WaterModel.TIP4P: "amber14/tip4pew.xml",
    WaterModel.SPC: "amber14/spce.xml",
    WaterModel.SPCE: "amber14/spce.xml",
    WaterModel.TIP5P: "tip5p.xml",
}

# String name accepted by Modeller.addSolvent(model=...).
_WATER_NAME: dict[WaterModel, str] = {
    WaterModel.TIP3P: "tip3p",
    WaterModel.TIP4P: "tip4pew",
    WaterModel.SPC: "spce",
    WaterModel.SPCE: "spce",
    WaterModel.TIP5P: "tip5p",
}


def solvate_openmm(
    parametrized: ParametrisedComplex,
    params: SolvationParams,
    output_gro: Path,
    output_top: Path,
) -> tuple[Path, Path]:
    """Solvate a parametrised complex with OpenMM + ParmEd (no BSS IO).

    Reuses the ``ForceField`` and ParmEd ``Structure`` carried by
    *parametrized* — built in the parametrization stage — so that AM1-BCC
    charge assignment and disk I/O are not repeated.

    Parameters
    ----------
    parametrized:
        Output of :func:`~gbsa_pipeline.parametrization.parametrize`.
        Must carry ``forcefield`` and ``parmed_structure`` (populated by the
        OpenMM parametrization path).
    params:
        Solvation settings (water model, box shape/size, ion concentration).
    output_gro:
        Destination path for the solvated GROMACS coordinate file.
    output_top:
        Destination path for the solvated GROMACS topology file.

    Returns:
    -------
    tuple[Path, Path]
        ``(output_gro, output_top)`` — the paths that were written.

    Raises:
    ------
    ValueError
        If *parametrized* was not produced by the OpenMM parametrization path
        (i.e. ``forcefield`` or ``parmed_structure`` is ``None``).
    """
    if parametrized.forcefield is None or parametrized.parmed_structure is None:
        raise ValueError(
            "solvate_openmm requires parametrized.forcefield and "
            "parametrized.parmed_structure to be set. "
            "Use parametrize() (not load_amber_complex) before calling this function."
        )

    output_gro.parent.mkdir(parents=True, exist_ok=True)

    water_model = (
        WaterModel(str(params.water_model).lower())
        if not isinstance(params.water_model, WaterModel)
        else params.water_model
    )
    box_shape = BoxShape(str(params.shape).lower()) if not isinstance(params.shape, BoxShape) else params.shape

    # ------------------------------------------------------------------
    # 1. Reuse the in-memory ParmEd structure from parametrization
    #    (all protein+ligand force field parameters already assigned)
    # ------------------------------------------------------------------
    existing: Any = parametrized.parmed_structure
    n_orig = len(existing.atoms)
    logger.debug("Using in-memory ParmEd structure (%d atoms).", n_orig)

    # ------------------------------------------------------------------
    # 2. Extend the existing ForceField with water templates
    #    GAFF is already registered with pre-assigned ligand charges,
    #    so AM1-BCC will not re-run.
    # ------------------------------------------------------------------
    water_xml = _WATER_XML[water_model]
    ff: ForceField = parametrized.forcefield
    logger.debug("Loading water FF (%s) into existing ForceField …", water_xml)
    ff.loadFile(water_xml)
    # ff.loadFile("amber14/ions.xml")

    # ------------------------------------------------------------------
    # 3. Add solvent (water + ions) via OpenMM Modeller
    #    Full ForceField ensures correct clash detection for all residues.
    # ------------------------------------------------------------------
    logger.debug("Creating Modeller from existing topology …")
    modeller = Modeller(existing.topology, existing.positions)

    kwargs: dict[str, Any] = {
        "model": _WATER_NAME[water_model],
        "neutralize": params.is_neutral,
    }
    if params.ion_concentration is not None:
        kwargs["ionicStrength"] = params.ion_concentration * mm_unit.molar
    if params.padding is not None:
        kwargs["padding"] = params.padding * mm_unit.nanometer
    else:
        kwargs["boxSize"] = boxSize = Vec3(params.box_size, params.box_size, params.box_size) * mm_unit.nanometer  # type: ignore[operator]
    if box_shape is BoxShape.TRUNCATED_OCTAHEDRON:
        kwargs["boxShape"] = "octahedron"

    logger.debug(
        "Adding solvent (model=%s, %s, box_shape=%s) …",
        kwargs["model"],
        f"padding={params.padding} nm" if params.padding is not None else f"box_size={params.box_size} nm",
        box_shape,
    )
    modeller.addSolvent(parametrized.forcefield, **kwargs)
    logger.debug("Solvated topology: %d atoms.", modeller.topology.getNumAtoms())

    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)

    # ------------------------------------------------------------------
    # 7. Write GROMACS gro + top
    # ------------------------------------------------------------------
    logger.debug("Writing GROMACS topology → %s …", output_top)
    structure.save(str(output_top), format="gromacs", overwrite=True)
    logger.debug("Writing GROMACS coordinates → %s …", output_gro)
    structure.save(str(output_gro), overwrite=True)
    logger.debug("Solvated files written.")

    return output_gro, output_top
