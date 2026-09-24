"""Membrane protein module.

Estimate membrane geometry parameters such as ``mthick`` and ``mctrdz``
for PB calculations directly from lipid phosphate atoms.
"""

from __future__ import annotations

import logging
import string
from dataclasses import dataclass
from typing import TYPE_CHECKING

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff

from gbsa_pipeline.mmbsa import PBParams

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path
    from typing import Any

logger = logging.getLogger(__name__)
__all__ = [
    "DEFAULT_LIPID_RESNAMES",
    "MembraneGeometry",
    "estimate_membrane_geometry",
    "extract_protein_ligand_system",
    "extract_receptor_pdb",
    "lipid_headgroup_restraint_atoms",
]


# Common phospholipid residue names used in PDB files.
DEFAULT_LIPID_RESNAMES: frozenset[str] = frozenset(
    {
        "DPP",
        "DPPC",
        "POP",
        "POPC",
        "POPE",
        "POPG",
        "POPS",
        "POPI",
        "DOP",
        "DOPC",
        "DOPE",
        "DMP",
        "DMPC",
        "DLPC",
    }
)

# A bilayer has exactly two leaflets.
_N_LEAFLETS = 2

# Index of the z-coordinate in a [x, y, z] positions array.
_Z_AXIS = 2


@dataclass(frozen=True)
class MembraneGeometry:
    """Bilayer geometry measured from lipid phosphate atoms.

    ``mctrdz`` is an absolute z-coordinate in the coordinate frame of the
    structure. ``mthick`` is the phosphate-to-phosphate bilayer thickness.
    """

    mctrdz: float
    mthick: float
    n_phosphates: int

    def pb_params(self, **overrides: Any) -> PBParams:
        """Build membrane-ready PBParams from this geometry."""
        kwargs: dict[str, Any] = {
            "memopt": 1,
            "mctrdz": self.mctrdz,
            "mthick": self.mthick,
            "eneopt": 1,
        }
        kwargs.update(overrides)

        return PBParams(**kwargs)


def estimate_membrane_geometry(
    universe: mda.Universe,
    lipid_resnames: Sequence[str] = tuple(DEFAULT_LIPID_RESNAMES),
    cutoff: float | None = None,
) -> MembraneGeometry:
    """Measure bilayer parameters from lipid phosphate atoms.

    Phosphate atoms are selected directly through MDAnalysis's own selection
    language, combining lipid residue names with a "P*" atom-name wildcard
    (real force-field topologies number phosphate atoms, e.g. "P8", "P31"),
    then grouped into two leaflets using LeafletFinder, a distance-based
    graph clustering.

    Assumes the bilayer normal is (approximately) the z-axis of ``universe``'s
    coordinate frame -- the convention used by essentially all membrane
    simulation builders, and required by gmx_MMPBSA's own implicit-membrane
    PB solver. A ValueError is raised if the two leaflets aren't primarily
    separated along z.
    """
    resnames = " ".join(sorted(set(lipid_resnames)))
    selection = f"resname {resnames} and name P*"
    phosphates = universe.select_atoms(selection)

    if len(phosphates) == 0:
        raise ValueError(
            f"No phosphate atoms belonging to {sorted(set(lipid_resnames))} were found. "
            "Check the lipid residue names and pass lipid_resnames explicitly."
        )
    if cutoff is None:
        try:
            cutoff, _ = optimize_cutoff(universe, selection, dmin=10.0, dmax=30.0, step=0.5)
        except Exception as exc:
            raise ValueError(
                "Could not auto-detect a LeafletFinder cutoff that splits "
                f"{sorted(set(lipid_resnames))} phosphates into two balanced leaflets "
                "between 10-30 Angstrom. Pass an explicit cutoff."
            ) from exc

    finder = LeafletFinder(universe, phosphates, cutoff=cutoff)
    groups = finder.groups()

    if len(groups) != _N_LEAFLETS:
        raise ValueError(
            "Phosphate atoms did not split into two leaflets. "
            "The structure may not contain a symmetric bilayer, or the "
            "lipid residue names may be incorrect."
        )

    upper, lower = groups
    separation = upper.centroid() - lower.centroid()

    # A bilayer whose true center sits at/near the periodic boundary has its two
    # leaflets near OPPOSITE box edges in unwrapped coordinates (e.g. one leaflet
    # near z=15, the other near z=box_z-15). The naive z-separation between them
    # then approaches the full box height instead of the true membrane
    # thickness -- confirmed on this project's own 2rh1/POPC testdata, where it
    # reported mthick=124 (box_z=161) instead of the true ~36.7. Wrapping the
    # z-component into (-box_z/2, box_z/2] recovers the true, shorter
    # periodic-image separation regardless of where the bilayer sits relative
    # to the box edges; it is a no-op when the bilayer does not straddle the
    # boundary (separation already within that range). Universes without box
    # information (e.g. a synthetic/merged system in a unit test) have no
    # periodic image to wrap against, so the raw separation is used as-is.
    box_z = float(universe.dimensions[_Z_AXIS]) if universe.dimensions is not None else None
    if box_z is not None:
        separation[_Z_AXIS] -= box_z * round(separation[_Z_AXIS] / box_z)

    lateral = float(np.linalg.norm(separation[:_Z_AXIS]))
    normal_component = abs(float(separation[_Z_AXIS]))

    if lateral > normal_component:
        raise ValueError("Rotate model so lipid layer along z axis")

    mthick = normal_component
    mctrdz = float(lower.centroid()[_Z_AXIS] + separation[_Z_AXIS] / 2)
    if box_z is not None:
        mctrdz %= box_z

    return MembraneGeometry(
        mctrdz=mctrdz,
        mthick=mthick,
        n_phosphates=len(phosphates),
    )


def lipid_headgroup_restraint_atoms(
    system: Any,
    lipid_resnames: Sequence[str],
) -> list[int]:
    """Absolute atom indices of lipid phosphate atoms, for use as a BSS restraint list.

    BSS's builtin "backbone"/"heavy"/"all" restraint keywords have no concept
    of a membrane -- "heavy" would restrain every non-hydrogen lipid tail atom
    too, freezing the whole bilayer instead of letting it relax around a fixed
    protein and headgroups. Restraining only phosphate atoms (the same "P*"
    name-prefix convention used by :func:`estimate_membrane_geometry`) mirrors
    CHARMM-GUI's standard equilibration protocol, which restrains lipid
    headgroups -- not the full lipid -- alongside the protein backbone during
    early NVT/NPT equilibration, then releases them before production.

    BSS accepts a restraint keyword *or* an explicit atom-index list for
    ``BSS.Protocol.Equilibration``, not both at once. To restrain protein
    backbone and lipid headgroups together, resolve "backbone" via
    ``system.getRestraintAtoms("backbone")`` first and pass the union of that
    with this function's result as the explicit list.

    ``system`` is a ``BSS._SireWrappers.System``. Atom indices are counted by
    accumulating each molecule's atom count in system order, matching the
    "absolute index" convention ``getRestraintAtoms`` itself returns.
    """
    resnames = set(lipid_resnames)
    indices: list[int] = []
    atom_offset = 0
    for mol in system.getMolecules():
        residues = mol.getResidues()
        if {res.name() for res in residues} & resnames:
            indices.extend(
                atom_offset + atom.index()
                for res in residues
                if res.name() in resnames
                for atom in res.getAtoms()
                if atom.name().startswith("P")
            )
        atom_offset += mol.nAtoms()

    if not indices:
        raise ValueError(
            f"No lipid phosphate atoms belonging to {sorted(resnames)} were found for "
            "restraints. Check the lipid residue names."
        )

    return indices


def extract_protein_ligand_system(
    system: Any,
    n_solute_molecules: int,
    n_protein_molecules: int,
) -> Any:
    """Strip lipids, water, and ions from a production membrane system, keeping only protein + ligand.

    A full bilayer patch (hundreds of lipids) makes the gmx_MMPBSA "complex"
    large enough to overflow AmberTools 24's 32-bit sander PB solver. Membrane
    geometry (``mctrdz``/``mthick``) is computed separately, from the original
    unstripped structure via :func:`estimate_membrane_geometry`, so reducing
    the complex here loses no membrane context -- this mirrors the official
    gmx_MMPBSA ``Protein_membrane`` example, whose Receptor is likewise
    protein-only despite a much larger full solvated structure.

    ``n_solute_molecules`` is the pre-built ``[membrane]`` system's molecule
    count (protein + lipids, before the ligand was merged in);
    ``n_protein_molecules`` is how many of those are protein. Molecules are
    identified by position, not name or number, matching the convention used
    throughout ``_stage_mmbsa``.

    Returns a new ``BSS._SireWrappers.System``; ``system`` itself is left
    untouched (``removeMolecules`` mutates in place, so this works on a copy).
    """
    reduced = system.copy()
    molecules = list(reduced.getMolecules())
    to_remove = molecules[n_protein_molecules:n_solute_molecules] + molecules[n_solute_molecules + 1 :]
    reduced.removeMolecules(to_remove)
    return reduced


def extract_receptor_pdb(
    gro_file: Path,
    output_pdb: Path,
) -> Path:
    """Extract a protein-only receptor PDB from a pre-built [membrane] system.

    Docking tools treat the receptor as rigid with no membrane representation
    -- lipids/water only matter for box placement and downstream MD, not for
    the docking score itself, so they are stripped here.

    GRO files carry no chain IDs, and a system built by excising a fusion
    protein (e.g. T4-lysozyme from a GPCR's ICL3) can leave the receptor as
    multiple polypeptide chains that are each independently renumbered from
    residue 1. Writing them all under one blank chain ID makes residue
    numbers collide across chains, which breaks downstream PDB parsers (e.g.
    Meeko's Polymer, which requires each chain:resid to be unique). A resid
    decrease is therefore treated as a new-chain boundary and each detected
    chain is given its own letter.
    """
    universe = mda.Universe(str(gro_file))
    protein = universe.select_atoms("protein")

    if protein.n_atoms == 0:
        raise ValueError(f"No protein atoms found in {gro_file}.")

    chain_ids = []
    current_chain = 0
    prev_resid = None
    for resid in protein.resids:
        if prev_resid is not None and resid < prev_resid:
            current_chain += 1
        chain_ids.append(string.ascii_uppercase[current_chain])
        prev_resid = resid

    universe.add_TopologyAttr("chainIDs")
    protein.chainIDs = chain_ids

    protein.write(str(output_pdb))
    return output_pdb
