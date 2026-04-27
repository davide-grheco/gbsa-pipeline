# src/gbsa_pipeline/solvation_openmm.py

"""OpenMM/ParmEd solvation that bypasses BioSimSpace IO."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import parmed as pmd
from openmm import Vec3
from openmm import unit as mm_unit
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile

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

    Carries both the on-disk GROMACS files written for inspection/checkpointing
    and the in-memory ParmEd structure so downstream stages can use either path
    without repeating disk I/O. The object is intentionally small because MD
    orchestration belongs in a later pipeline layer, not in this solvation
    helper. The optional in-memory structure is useful when the caller continues
    in Python immediately, while the files remain the stable interface for
    BioSimSpace loading and visual inspection.
    """

    gro_file: Path
    top_file: Path
    parmed_structure: Any = field(default=None, hash=False, compare=False, repr=False)

    def load_bss(self) -> Any:
        """Load this complex as a BioSimSpace System for MD stages.

        Reads from the GROMACS files already written to disk. The returned
        system is ready for later minimization, equilibration, and production MD
        helpers. This method intentionally performs only loading and does not
        start any simulation stage. Missing files are reported explicitly because
        this object is often used after long integration-test runs.
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
    """Solvate a parametrised complex with OpenMM + ParmEd.

    Reuses the ``ForceField`` and ParmEd ``Structure`` carried by
    *parametrized*, so ligand charges and protein-ligand parameters are not
    regenerated. If parametrization extracted crystallographic waters, those
    waters are restored into the OpenMM modeller before bulk solvent is added.
    This makes retained waters part of the pre-solvation system, so newly
    generated solvent is placed around protein, ligand, and crystal waters
    together. The function writes GROMACS ``.gro`` and ``.top`` files and
    returns a small object carrying the paths and in-memory ParmEd structure.
    """
    if parametrized.forcefield is None or parametrized.parmed_structure is None:
        raise ValueError(
            "solvate_openmm requires parametrized.forcefield and "
            "parametrized.parmed_structure to be set. "
            "Use parametrize() (not load_amber_complex) before calling this function."
        )

    output_gro.parent.mkdir(parents=True, exist_ok=True)

    water_model = _coerce_water_model(params.water_model)
    box_shape = _coerce_box_shape(params.shape)

    # ------------------------------------------------------------------
    # 1. Reuse the in-memory ParmEd structure from parametrization.
    #    This is the dry protein-ligand complex, with all non-water
    #    protein and ligand parameters already assigned.
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
    # 3. Build a complete pre-solvation modeller.
    #    Crystal waters are restored before addSolvent so the bulk solvent
    #    placement sees them as existing atoms and avoids overlaps.
    # ------------------------------------------------------------------
    logger.debug("Creating Modeller from existing topology …")
    modeller = Modeller(existing.topology, existing.positions)

    restored_waters_pdb = _restore_crystal_waters_before_solvation(
        modeller=modeller,
        forcefield=ff,
        crystal_waters_pdb=parametrized.crystal_waters_pdb,
        output_pdb=output_gro.parent / "restored_crystal_waters.pdb",
    )
    if restored_waters_pdb is not None:
        logger.debug("Restored crystallographic waters from %s.", restored_waters_pdb)

    # ------------------------------------------------------------------
    # 4. Add bulk solvent and ions via OpenMM Modeller.
    # ------------------------------------------------------------------
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
    # 5. Build full OpenMM system → ParmEd structure.
    #    rigidWater=False keeps O-H bonds in HarmonicBondForce so ParmEd
    #    can resolve water geometry when writing the GROMACS topology.
    # ------------------------------------------------------------------
    logger.debug("Creating solvated OpenMM system …")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
        rigidWater=False,
    )
    logger.debug("Converting to ParmEd structure …")
    structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)

    # ------------------------------------------------------------------
    # 6. Write GROMACS gro + top.
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


def _restore_crystal_waters_before_solvation(
    *,
    modeller: Modeller,
    forcefield: ForceField,
    crystal_waters_pdb: Path | None,
    output_pdb: Path,
) -> Path | None:
    """Restore extracted crystallographic waters before bulk solvation.

    The parametrization stage stores crystal waters separately because the
    protein-ligand system is parameterized dry. This helper uses OpenMM's own
    ``Modeller`` machinery to add missing hydrogens to those waters, writes a
    small inspection PDB, and adds the completed waters to the pre-solvation
    modeller before ``addSolvent`` is called. The added waters therefore become
    part of the pre-solvation system and are included in clash avoidance when
    bulk water is placed. ``None`` is returned when no retained-water file is
    available.
    """
    if crystal_waters_pdb is None:
        _remove_stale_file(output_pdb)
        return None

    if not crystal_waters_pdb.exists():
        _remove_stale_file(output_pdb)
        logger.debug("Crystal water file not found: %s.", crystal_waters_pdb)
        return None

    water_pdb = PDBFile(str(crystal_waters_pdb))
    water_modeller = Modeller(water_pdb.topology, water_pdb.positions)

    atoms_before = water_modeller.topology.getNumAtoms()
    water_modeller.addHydrogens(forcefield)
    atoms_after = water_modeller.topology.getNumAtoms()

    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    with output_pdb.open("w", encoding="utf-8") as handle:
        PDBFile.writeFile(
            water_modeller.topology,
            water_modeller.positions,
            handle,
            keepIds=True,
        )

    modeller.add(water_modeller.topology, water_modeller.positions)

    logger.debug(
        "Restored crystal waters with OpenMM hydrogens (%d -> %d atoms).",
        atoms_before,
        atoms_after,
    )
    return output_pdb


def _coerce_box_shape(shape: BoxShape | str) -> BoxShape:
    """Return a validated box-shape enum.

    ``SolvationParams`` normally validates this already, but this helper keeps
    ``solvate_openmm`` robust for direct callers and tests that pass strings.
    It intentionally mirrors the small coercion pattern used in
    ``solvation_box.py`` instead of introducing a shared abstraction. Invalid
    values should fail with the normal ``BoxShape`` error message.
    """
    if isinstance(shape, BoxShape):
        return shape

    return BoxShape(str(shape).lower().strip())


def _coerce_water_model(model: WaterModel | str) -> WaterModel:
    """Return a validated water-model enum.

    ``SolvationParams`` normally validates this already, but keeping the
    conversion local makes the OpenMM helper safe for direct use. The resulting
    enum is used both to select the OpenMM XML file and to pass the model name
    accepted by ``Modeller.addSolvent``. Unsupported values fail through the
    ``WaterModel`` enum constructor. No model-specific chemistry is implemented
    in this helper.
    """
    if isinstance(model, WaterModel):
        return model

    return WaterModel(str(model).lower().strip())


def _remove_stale_file(path: Path) -> None:
    """Remove a generated file if it exists.

    Persistent integration-test directories can otherwise keep stale artefacts
    from an earlier run, which makes visual inspection misleading. The helper is
    intentionally small and silent because the absence of a generated water file
    is valid when no crystal waters were retained. Only ``FileNotFoundError`` is
    suppressed; other filesystem errors should still surface. This keeps stale
    cleanup local to the generated restore artefact.
    """
    try:
        path.unlink()
    except FileNotFoundError:
        return
