# /home/grheco/repositorios/gbsa-pipeline/src/gbsa_pipeline/parametrization.py

"""Parametrize protein-ligand complexes."""

from __future__ import annotations

import contextlib
import logging
import pickle
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

import BioSimSpace as BSS
import parmed as pmd
from openff.toolkit.topology import Molecule
from openmm.app import ForceField, Modeller, NoCutoff, PDBFile
from openmmforcefields.generators import GAFFTemplateGenerator
from pydantic import BaseModel, ConfigDict, Field, field_validator

if TYPE_CHECKING:
    from BioSimSpace._SireWrappers import System

from gbsa_pipeline.parametrization_enum import ChargeMethod, LigandFF, ProteinFF

PathLike = Union[str, Path]

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Force field configuration
# ---------------------------------------------------------------------------


class ParametrizationConfig(BaseModel):
    """Force field and charge method choices for a parametrization run.

    Defaults to AMBER ff14SB + GAFF2 + AM1-BCC.
    Use the class-method presets for the most common combinations, or
    construct directly to override individual axes.
    """

    model_config = ConfigDict(frozen=True, extra="forbid", validate_default=True)

    protein_ff: ProteinFF = ProteinFF.FF14SB
    ligand_ff: LigandFF = LigandFF.GAFF2
    charge_method: ChargeMethod = ChargeMethod.AM1BCC
    extra_ff_files: tuple[Path, ...] = ()

    @field_validator("extra_ff_files", mode="before")
    @classmethod
    def _check_extra_ff_files(cls, paths: Any) -> tuple[Path, ...]:
        result = tuple(Path(p) for p in paths)
        missing = [p for p in result if not p.exists()]
        if missing:
            raise ValueError("Extra force field files not found: " + ", ".join(str(p) for p in missing))
        return result

    @classmethod
    def amber14_gaff2(cls) -> ParametrizationConfig:
        """AMBER ff14SB + GAFF2 + AM1-BCC charges (default)."""
        return cls(protein_ff=ProteinFF.FF14SB, charge_method=ChargeMethod.AM1BCC)

    @classmethod
    def amber19_gaff2(cls) -> ParametrizationConfig:
        """AMBER ff19SB + GAFF2 + AM1-BCC charges."""
        return cls(protein_ff=ProteinFF.FF19SB, charge_method=ChargeMethod.AM1BCC)

    @classmethod
    def amber14_gaff2_nagl(cls) -> ParametrizationConfig:
        """AMBER ff14SB + GAFF2 + NAGL charges."""
        return cls(charge_method=ChargeMethod.NAGL)


# ---------------------------------------------------------------------------
# User-facing input model
# ---------------------------------------------------------------------------


class ParametrizationInput(BaseModel):
    """Validated inputs for a parametrization run."""

    model_config = ConfigDict(frozen=True, extra="forbid", validate_default=True)

    protein_pdb: Path
    ligand_sdf: Path
    config: ParametrizationConfig = Field(default_factory=ParametrizationConfig)
    net_charge: int | None = None
    work_dir: Path | None = None

    @field_validator("protein_pdb", "ligand_sdf")
    @classmethod
    def _check_exists(cls, path: Path) -> Path:
        if not path.exists():
            raise ValueError(f"File not found: {path}")
        return path


# ---------------------------------------------------------------------------
# Output type
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ParametrisedComplex:
    """Parametrised protein-ligand complex ready for solvation and MD."""

    gro_file: Path
    top_file: Path
    config: ParametrizationConfig
    forcefield: Any = field(default=None, hash=False, compare=False, repr=False)
    parmed_structure: Any = field(default=None, hash=False, compare=False, repr=False)


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def parametrize(inp: ParametrizationInput) -> ParametrisedComplex:
    """Parametrize a protein-ligand complex without running tleap."""
    return _parametrize_openmm(inp)


# ---------------------------------------------------------------------------
# OpenMM implementation
# ---------------------------------------------------------------------------

_PROTEIN_FF_XML: dict[ProteinFF, list[str]] = {
    ProteinFF.FF14SB: ["amber14-all.xml"],
    ProteinFF.FF19SB: ["amber/protein.ff19SB.xml"],
    ProteinFF.FF99SB: ["amber/protein.ff99SBildn.xml"],
}

_GAFF_FF_VERSION: dict[LigandFF, str] = {
    LigandFF.GAFF: "gaff-1.81",
    LigandFF.GAFF2: "gaff-2.11",
}


def _load_single_ligand_from_sdf(ligand_sdf: Path) -> Molecule:
    """Load one ligand from SDF and normalize multi-molecule returns.

    Exported docking SDFs may contain multiple poses/molecules.
    For the current pipeline stage, we take the first molecule.
    """
    logger.debug("Loading ligand SDF: %s …", ligand_sdf)

    loaded = Molecule.from_file(str(ligand_sdf))

    if isinstance(loaded, list):
        if len(loaded) == 0:
            raise ValueError(f"No ligand molecules found in SDF: {ligand_sdf}")

        if len(loaded) > 1:
            logger.warning(
                "Ligand SDF '%s' contained %d molecules/poses. Using the first one.",
                ligand_sdf,
                len(loaded),
            )

        ligand = loaded[0]
    else:
        ligand = loaded

    if not isinstance(ligand, Molecule):
        raise TypeError(f"Expected OpenFF Molecule from '{ligand_sdf}', got {type(ligand).__name__}")

    if not ligand.conformers:
        raise ValueError(
            f"Ligand SDF '{ligand_sdf}' contains no 3-D conformers. Provide an SDF file with embedded 3-D coordinates."
        )

    logger.debug(
        "Ligand loaded (%d atoms, %d conformers).",
        ligand.n_atoms,
        len(ligand.conformers),
    )

    return ligand


def _parametrize_openmm(inp: ParametrizationInput) -> ParametrisedComplex:
    work_dir = inp.work_dir or Path(tempfile.mkdtemp(prefix="gbsa_param_"))
    work_dir.mkdir(parents=True, exist_ok=True)

    # --- Protein -------------------------------------------------------
    logger.debug("Loading protein PDB: %s …", inp.protein_pdb)
    pdb = PDBFile(str(inp.protein_pdb))
    logger.debug("Protein PDB loaded (%d atoms).", pdb.topology.getNumAtoms())

    # --- Ligand --------------------------------------------------------
    ligand = _load_single_ligand_from_sdf(inp.ligand_sdf)

    logger.debug(
        "Assigning partial charges (method=%s) …",
        inp.config.charge_method.value,
    )

    if inp.net_charge is not None:
        logger.warning(
            "ParametrizationInput.net_charge=%s was provided, but the current "
            "OpenFF assign_partial_charges() path does not explicitly use it. "
            "Proceeding with toolkit-assigned charges for now.",
            inp.net_charge,
        )

    kwargs: dict[str, Any] = {
        "partial_charge_method": inp.config.charge_method.value,
        "normalize_partial_charges": True,
        "use_conformers": ligand.conformers,
    }

    ligand.assign_partial_charges(**kwargs)
    logger.debug("Partial charges assigned.")

    # --- Force field ---------------------------------------------------
    protein_xmls = _PROTEIN_FF_XML[inp.config.protein_ff]
    extra_xmls = [str(p) for p in inp.config.extra_ff_files]
    logger.debug(
        "Building force field (protein=%s, extra=%d files) …",
        inp.config.protein_ff.value,
        len(extra_xmls),
    )
    forcefield = ForceField(*protein_xmls, *extra_xmls)
    logger.debug("Force field built.")

    logger.debug(
        "Registering GAFF template generator (%s) …",
        _GAFF_FF_VERSION[inp.config.ligand_ff],
    )
    gaff = GAFFTemplateGenerator(
        molecules=[ligand],
        forcefield=_GAFF_FF_VERSION[inp.config.ligand_ff],
        cache=None,
    )
    forcefield.registerTemplateGenerator(gaff.generator)
    logger.debug("GAFF template generator registered.")

    logger.debug("Combining protein+ligand topology …")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0].to_openmm())
    logger.debug("Combined topology: %d atoms.", modeller.topology.getNumAtoms())

    logger.debug("Creating OpenMM system (may take a moment) …")
    system: System = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    logger.debug("OpenMM system created (%d particles).", system.getNumParticles())

    logger.debug("Converting OpenMM system to ParmEd structure …")
    structure = pmd.openmm.load_topology(modeller.topology, system, modeller.positions)
    logger.debug("ParmEd structure ready (%d atoms).", len(structure.atoms))

    gro_file = work_dir / "complex.gro"
    top_file = work_dir / "complex.top"
    logger.debug("Writing GROMACS topology → %s …", top_file)
    structure.save(str(top_file), format="gromacs")
    logger.debug("Writing GROMACS coordinates → %s …", gro_file)
    structure.save(str(gro_file))
    logger.debug("GROMACS files written.")

    complex_ = ParametrisedComplex(
        gro_file=gro_file,
        top_file=top_file,
        config=inp.config,
        forcefield=forcefield,
        parmed_structure=structure,
    )

    cache_file = work_dir / "complex.pickle"
    with contextlib.suppress(Exception):
        cache_file.write_bytes(pickle.dumps(complex_))

    return complex_


# ---------------------------------------------------------------------------
# Legacy BSS wrappers — kept for API compatibility
# ---------------------------------------------------------------------------


def load_protein_pdb(pdb_path: PathLike) -> BSS._SireWrappers.Molecule:
    """Load a protein from a PDB file and return the first molecule."""
    system = BSS.IO.readMolecules(str(pdb_path))
    mols = system.getMolecules()
    if not mols:
        raise ValueError(f"No molecules found in {pdb_path}")
    return mols[0]


def parameterise_protein_amber(
    protein: BSS._SireWrappers.Molecule,
    ff: str = "ff14SB",
    water_model: str | None = None,
    work_dir: PathLike | None = None,
) -> BSS._SireWrappers.Molecule:
    """Parameterize a protein via BioSimSpace/tleap."""
    ff = ff.lower().strip()
    kwargs: dict[str, Any] = {}
    if water_model is not None:
        kwargs["water_model"] = water_model
    if work_dir is not None:
        kwargs["work_dir"] = str(work_dir)

    if ff == "ff14sb":
        out = BSS.Parameters.ff14SB(protein, **kwargs)
    elif ff == "ff19sb":
        out = BSS.Parameters.ff19SB(protein, **kwargs)
    elif ff == "ff99sb":
        out = BSS.Parameters.ff99SB(protein, **kwargs)
    else:
        raise ValueError(f"Unsupported protein FF '{ff}'. Try ff14SB, ff19SB, ff99SB.")

    return _ensure_molecule(out)


def _ensure_molecule(x: Any) -> BSS._SireWrappers.Molecule:
    """Ensure that a parameterization output is returned as a BSS Molecule."""
    if hasattr(x, "getMolecule"):
        return x.getMolecule()
    return x


def parameterise_ligand_gaff2(
    ligand: BSS._SireWrappers.Molecule,
    net_charge: int | None = None,
    charge_method: str = "BCC",
    work_dir: PathLike | None = None,
) -> BSS._SireWrappers.Molecule:
    """Parameterise a ligand via BioSimSpace/antechamber (GAFF2)."""
    kwargs: dict[str, Any] = {
        "net_charge": net_charge,
        "charge_method": charge_method,
    }
    if work_dir is not None:
        kwargs["work_dir"] = str(work_dir)

    return _ensure_molecule(BSS.Parameters.gaff2(ligand, **kwargs))


def make_protein_ligand_system(
    protein: BSS._SireWrappers.Molecule,
    ligand: BSS._SireWrappers.Molecule,
) -> BSS._SireWrappers.System:
    """Combine a parametrised protein and ligand into a BSS System."""
    system = BSS._SireWrappers.System(protein)
    system.addMolecules(ligand)
    return system


def load_and_parameterise(
    protein_pdb: PathLike,
    ligand: PathLike | BSS._SireWrappers.Molecule,
    protein_ff: str = "ff14SB",
    ligand_net_charge: int | None = None,
    ligand_charge_method: str = "BCC",
    work_dir: PathLike | None = None,
) -> ParametrisedComplex:
    """Load and parametrize a protein-ligand complex.

    .. deprecated::
        Use :func:`parametrize` with :class:`ParametrizationInput` instead.
    """
    warnings.warn(
        "load_and_parameterise() is deprecated. Use parametrize() with ParametrizationInput instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    if not isinstance(ligand, (str, Path)):
        raise TypeError(
            "load_and_parameterise() now requires ligand as a file path. "
            "Save the molecule to SDF first, or use parametrize() directly."
        )
    return parametrize(
        ParametrizationInput(
            protein_pdb=Path(protein_pdb),
            ligand_sdf=Path(ligand),
            config=ParametrizationConfig(protein_ff=ProteinFF.from_str(protein_ff)),
            net_charge=ligand_net_charge,
            work_dir=Path(work_dir) if work_dir else None,
        )
    )


def export_gromacs_top_gro(
    system: System,
    prefix: str,
) -> list[Path]:
    """Export GROMACS .gro and .top files from a BSS System."""
    out_gro = Path(f"{prefix}.gro")
    out_top = Path(f"{prefix}.top")

    BSS.IO.saveMolecules(str(out_gro), system, fileformat="gro87")
    BSS.IO.saveMolecules(str(out_top), system, fileformat="grotop")

    return [out_gro, out_top]
