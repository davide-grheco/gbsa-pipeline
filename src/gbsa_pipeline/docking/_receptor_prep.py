"""Receptor preparation: PDB → PDBQT conversion and crystal-water merging."""

from __future__ import annotations

import logging
from pathlib import Path

import gemmi
from meeko import (
    MoleculePreparation,
    PDBQTWriterLegacy,
    Polymer,
    PolymerCreationError,
    ResidueChemTemplates,
)
from rdkit import Chem

from gbsa_pipeline._gemmi_utils import append_chains, filter_residues, write_pdb
from gbsa_pipeline._paths import require_file

LOGGER = logging.getLogger(__name__)


def _strip_hetatm(receptor_pdb: Path, dest: Path) -> Path:
    """Write a copy of receptor_pdb with HETATM records removed.

    Meeko processes ATOM and HETATM identically, so modified residues stored
    as HETATM (e.g. CSD, PTR) that share a residue number with a protein ATOM
    residue cause ``each residue key must have exactly 1 resname`` errors.
    Docking receptors should only contain protein atoms anyway.
    """
    st = gemmi.read_pdb(str(receptor_pdb))
    filter_residues(st, lambda res: res.het_flag != "H")
    return write_pdb(st, dest)


def _merge_sdfs_into_pdb(pdb: Path, sdfs: list[Path], output: Path) -> Path:
    """Append cofactor SDF atoms to a protein PDB as HETATM records.

    Each SDF is read with RDKit (hydrogens preserved) and its chains are merged
    into the protein structure via gemmi so meeko assigns AutoDock atom types to
    protein and cofactor(s) together in one pass.
    """
    st = gemmi.read_pdb(str(pdb))
    for sdf in sdfs:
        supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
        mol = next(iter(supplier), None)
        if mol is None:
            raise ValueError(f"Could not read cofactor SDF: {sdf}")
        pdb_block = Chem.MolToPDBBlock(mol) or ""
        append_chains(st, gemmi.read_pdb_string(pdb_block))
    return write_pdb(st, output)


def convert_receptor_pdb_to_pdbqt(
    receptor_pdb: Path,
    output_path: Path | None = None,
    *,
    cofactor_sdfs: list[Path] | None = None,
) -> Path:
    """Convert a receptor PDB to rigid receptor PDBQT using the Meeko Python API.

    The polymer is built directly from the PDB string with Gasteiger charges
    (the approach validated against non-standard-residue receptors such as
    3OLL, which the previous set_template/CYX-retry path could not convert).
    Receptor hydrogen addition, protonation-state decisions, and structural
    cleanup are expected to happen upstream (e.g. in stack_protein_prep).
    """
    receptor_pdb = require_file(Path(receptor_pdb), "Receptor PDB")

    if receptor_pdb.suffix.lower() != ".pdb":
        raise ValueError(f"Expected a .pdb receptor input, got: {receptor_pdb}")

    if output_path is None:
        output_path = receptor_pdb.with_suffix(".pdbqt")

    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Strip HETATM before Meeko so modified residues (CSD, PTR, …) stored as
    # HETATM with the same residue number as a protein ATOM don't cause errors.
    stripped = _strip_hetatm(receptor_pdb, output_path.parent / f"{receptor_pdb.stem}_protein_only.pdb")

    # Merge cofactors after stripping so meeko assigns AutoDock atom types to
    # both protein and cofactor atoms in one pass.
    if cofactor_sdfs:
        pdb_for_meeko = _merge_sdfs_into_pdb(
            stripped,
            cofactor_sdfs,
            output_path.parent / f"{receptor_pdb.stem}_with_cofactor.pdb",
        )
        LOGGER.info("Merged %d cofactor(s) into receptor PDB for meeko.", len(cofactor_sdfs))
    else:
        pdb_for_meeko = stripped

    pdb_string = pdb_for_meeko.read_text(encoding="utf-8")

    LOGGER.info("Preparing receptor with Meeko Python API: %s → %s", receptor_pdb.name, output_path.name)

    templates = ResidueChemTemplates.create_from_defaults()
    mk_prep = MoleculePreparation.from_config({"charge_model": "gasteiger"})
    try:
        polymer = Polymer.from_pdb_string(pdb_string, templates, mk_prep)
    except PolymerCreationError as exc:
        raise RuntimeError(f"Meeko could not build polymer from PDB: {exc}") from exc

    pdbqt_string, _flex = PDBQTWriterLegacy.write_from_polymer(polymer)

    if not pdbqt_string.strip():
        raise RuntimeError(
            f"Meeko produced an empty PDBQT for receptor: {receptor_pdb}\n"
            "Check that the PDB contains valid protein ATOM records."
        )

    output_path.write_text(pdbqt_string, encoding="utf-8")
    LOGGER.info("Meeko receptor PDBQT written: %s", output_path.name)
    return output_path


def merge_pdb_structures(base_pdb: Path, extra_pdb: Path, output_pdb: Path) -> Path:
    """Merge two PDB files by appending all chains from extra_pdb into base_pdb.

    Chain names that conflict with base_pdb are automatically renamed. The
    merged structure is written to output_pdb.
    """
    base_st = gemmi.read_pdb(str(base_pdb))
    append_chains(base_st, gemmi.read_pdb(str(extra_pdb)))
    return write_pdb(base_st, output_pdb)
