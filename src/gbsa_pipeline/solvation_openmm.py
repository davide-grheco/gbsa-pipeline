"""OpenMM/ParmEd solvation that bypasses BioSimSpace IO."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
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


@dataclass(frozen=True)
class SolvatedComplex:
    """Solvated protein-ligand complex produced by :func:`solvate_openmm`.

    Carries both the on-disk GROMACS files (written for inspection and
    checkpoint purposes) and the in-memory ParmEd structure so that
    downstream steps can use either without repeating disk I/O.

    Attributes:
    ----------
    gro_file:
        GROMACS coordinate file (``.gro``) of the solvated system.
    top_file:
        GROMACS topology file (``.top``) of the solvated system.
    parmed_structure:
        In-memory ParmEd ``Structure`` with all force-field parameters.
        ``None`` when the object is constructed from existing files rather
        than from a live parametrization run.
    """

    gro_file: Path
    top_file: Path
    parmed_structure: Any = field(default=None, hash=False, compare=False, repr=False)

    def load_bss(self) -> Any:
        """Load this complex as a BioSimSpace System for MD stages.

        Reads from the GROMACS files already written to disk.  The returned
        system is ready for minimization, equilibration, and production MD.

        Returns:
        -------
        BioSimSpace._SireWrappers.System
            The solvated system loaded into BioSimSpace.
        """
        if not self.gro_file.exists() or not self.top_file.exists():
            raise FileNotFoundError(f"SolvatedComplex files not found: {self.gro_file}, {self.top_file}.")
        import BioSimSpace as BSS  # noqa: PLC0415

        return BSS.IO.readMolecules([str(self.gro_file), str(self.top_file)])


def solvate_openmm(
    parametrized: ParametrisedComplex,
    params: SolvationParams,
    output_gro: Path,
    output_top: Path,
) -> SolvatedComplex:
    """Solvate a parametrised complex with OpenMM + ParmEd (no BSS IO).

    Reuses the ``ForceField`` and ParmEd ``Structure`` carried by
    *parametrized* — built in the parametrization stage — so that AM1-BCC
    charge assignment is not repeated.  Writes GROMACS ``.gro`` and ``.top``
    files and returns a :class:`SolvatedComplex` that carries both the file
    paths and the in-memory ParmEd structure.

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
    SolvatedComplex
        Dataclass holding the written file paths and the in-memory ParmEd
        structure of the solvated system.

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
    # 2. Extend the existing ForceField with water templates.
    #    GAFF is already registered with pre-assigned ligand charges,
    #    so AM1-BCC will not re-run.
    # ------------------------------------------------------------------
    water_xml = _WATER_XML[water_model]
    ff: ForceField = parametrized.forcefield
    logger.debug("Loading water FF (%s) into existing ForceField …", water_xml)
    ff.loadFile(water_xml)

    # ------------------------------------------------------------------
    # 3. Add solvent (water + ions) via OpenMM Modeller.
    #    Full ForceField ensures correct clash detection for all residues.
    # ------------------------------------------------------------------
    logger.debug("Creating Modeller from existing topology …")
    modeller = Modeller(existing.topology, existing.positions)

    kwargs: dict[str, Any] = {
        "model": _WATER_NAME[water_model],
        "neutralize": params.neutralize,
    }
    if params.ion_concentration is not None:
        kwargs["ionicStrength"] = params.ion_concentration * mm_unit.molar
    if params.padding is not None:
        kwargs["padding"] = params.padding * mm_unit.nanometer
    else:
        kwargs["boxSize"] = Vec3(params.box_size, params.box_size, params.box_size) * mm_unit.nanometer
    if box_shape is BoxShape.TRUNCATED_OCTAHEDRON:
        kwargs["boxShape"] = "octahedron"

    logger.debug(
        "Adding solvent (model=%s, %s, box_shape=%s) …",
        kwargs["model"],
        f"padding={params.padding} nm" if params.padding is not None else f"box_size={params.box_size} nm",
        box_shape,
    )
    modeller.addSolvent(ff, **kwargs)
    logger.debug("Solvated topology: %d atoms.", modeller.topology.getNumAtoms())

    # ------------------------------------------------------------------
    # 4. Build full OpenMM system → ParmEd structure.
    #    rigidWater=False keeps O-H bonds in HarmonicBondForce so ParmEd
    #    can resolve SETTLE geometry when writing GROMACS topology.
    # ------------------------------------------------------------------
    logger.debug("Creating solvated OpenMM system …")
    system = ff.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None, rigidWater=False)
    logger.debug("Converting to ParmEd structure …")
    structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)

    # ------------------------------------------------------------------
    # 5. Write GROMACS gro + top
    # ------------------------------------------------------------------
    logger.debug("Writing GROMACS topology → %s …", output_top)
    structure.save(str(output_top), format="gromacs", overwrite=True)
    logger.debug("Writing GROMACS coordinates → %s …", output_gro)
    structure.save(str(output_gro), overwrite=True)
    logger.debug("Solvated files written.")

    return SolvatedComplex(
        gro_file=output_gro,
        top_file=output_top,
        parmed_structure=structure,
    )
