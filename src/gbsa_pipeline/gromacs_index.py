"""Module to generate GROMACS index files from residue-name-based atom selection.

This module separates *selecting* which atoms belong to the Receptor/Ligand
groups from *writing* them out in the standard GROMACS ``.ndx`` format. The
output file is intended to be passed to ``gmx_MMPBSA`` via the ``-ci`` flag or
to standard GROMACS tools such as ``gmx trjconv`` and ``gmx make_ndx``.

Selection is driven entirely by residue name, read directly from a coordinate
file (``.gro``/``.pdb``) via MDAnalysis -- a bare ``.gro`` already carries
resnames, so no ``.top``/``.tpr`` file is needed. Receptor = everything that
is not solvent/ions and not the ligand; Ligand = the ligand's residue name.
No molecule-position or -ordering assumption is made anywhere: the module's
previous position/number-based implementation assumed water/ions always sit
after the ligand in molecule order, which does not hold in general.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.selection import ProteinSelection

from gbsa_pipeline._constants import ION_RESIDUE_NAMES, WATER_RESIDUE_NAMES

if TYPE_CHECKING:
    from collections.abc import Sequence
    from io import TextIOWrapper
    from pathlib import Path

    import sire.system

# Residue names GMXMMPBSA.make_top.cleantop() strips from the "-cp" topology
# before matching it against our Receptor/Ligand index. Hardcoded here
# because it's a local variable inside cleantop(), not a public constant, in
# the pinned gmx_MMPBSA==1.5.0.3:
# https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/blob/v1.5.0.3/GMXMMPBSA/make_top.py#L645-L704
# Notably missing "HOH" -- cleantop() will NOT strip crystallographic water
# named "HOH" (common in CHARMM-GUI-style prebuilt systems).
#
# gmx_MMPBSA>=1.7.0 exposes this as a real, importable constant
# (GMXMMPBSA.make_top.solvent_ion_residues), but upgrading to it requires
# downgrading ambertools (gmx_MMPBSA==1.7.0 needs ambertools<24; this project
# pins >=24.8,<25 to work around a real AmberTools 24 sander bug used
# elsewhere in this codebase). I do not know what would spit errors. I put
# it in the comment because I don't know what would be the proper way to
# communicate it over GitHub.
_CLEANTOP_STRIPPED_RESNAMES: frozenset[str] = frozenset(
    {
        "NA",
        "CL",
        "SOL",
        "SOD",
        "Na+",
        "CLA",
        "Cl-",
        "POT",
        "K+",
        "TIP3P",
        "TIP3",
        "TP3",
        "TIPS3P",
        "TIP3o",
        "TIP4P",
        "TIP4PEW",
        "T4E",
        "TIP4PD",
        "TIP5P",
        "SPC",
        "SPC/E",
        "SPCE",
        "WAT",
        "OPC",
    }
)

# Broader, human-recognizable solvent/ion names -- a superset of what
# cleantop() actually strips (e.g. it also recognizes "HOH", which cleantop()
# does not). Built from the package's shared water/ion name constants (see
# gbsa_pipeline._constants) rather than a fresh list, so this stays in sync
# with the names other stages already recognize as solvent.
_LOOKS_LIKE_SOLVENT_RESNAMES: frozenset[str] = WATER_RESIDUE_NAMES | ION_RESIDUE_NAMES


def identify_ligand_resname(system: sire.system.System) -> str:
    """Identify the ligand's residue name in a merged protein(+lipid)+ligand system.

    Identifies the ligand *by composition*, not by molecule position/index:
    every molecule is classified as protein (a standard amino-acid residue
    set, via MDAnalysis's own maintained ``ProteinSelection.prot_res`` table)
    or solvent/ions (``_LOOKS_LIKE_SOLVENT_RESNAMES``); whatever residue name
    is left over is the ligand. A fixed molecule index is not a reliable way
    to find the ligand -- GROMACS round-trips only guarantee append-only
    ordering, and the exact position varies between a [system] run, a
    [membrane] run, and however many lipids a given membrane patch has.

    Call this on a protein+ligand system (a [system] run's production system,
    or a [membrane] run's already lipid-reduced complex from
    :func:`~gbsa_pipeline.membrane.extract_protein_ligand_system`) -- lipids
    are not classified as protein or solvent here, so a system that still
    contains them would be misidentified.

    Raises ValueError if zero or more than one non-protein, non-solvent
    residue name is found.
    """
    candidates: set[str] = set()
    for mol in system:
        resnames = {res.name().value() for res in mol.residues()}
        if resnames <= ProteinSelection.prot_res or resnames <= _LOOKS_LIKE_SOLVENT_RESNAMES:
            continue
        candidates.update(resnames)

    if not candidates:
        raise ValueError("Could not identify a ligand: every molecule looks like protein or solvent/ions.")
    if len(candidates) > 1:
        raise ValueError(
            f"Ambiguous ligand identification: found multiple non-protein, non-solvent "
            f"residue names {sorted(candidates)}. Expected exactly one."
        )

    return next(iter(candidates))


def select_receptor_and_ligand_atoms(coord_file: Path, ligand_resname: str) -> tuple[list[int], list[int]]:
    """Select Receptor/Ligand atom indices for gmx_MMPBSA, from residue names alone.

    Works identically for both a [system] (soluble) run and a [membrane] run:
    Receptor is "every atom that is not solvent/ions and not the ligand" --
    protein for a soluble run, protein and lipids for a membrane run, since
    neither is a recognized solvent/ion name and gmx_MMPBSA's own topology
    cleaning (``GMXMMPBSA.make_top.cleantop``) never strips lipids either.
    No molecule-position or -ordering assumption is made: unlike this
    module's previous ``select_receptor_and_ligand_atoms_by_position``,
    solvent/ions are not assumed to sit after the ligand (or after anything
    else) in the file.

    ``coord_file`` is a ``.gro`` or ``.pdb`` -- resnames are already present
    in a bare coordinate file, so no ``.top``/``.tpr`` is read here.
    """
    if ligand_resname in _CLEANTOP_STRIPPED_RESNAMES:
        raise ValueError(
            f"Ligand resname {ligand_resname!r} is a name gmx_MMPBSA's own topology "
            "cleaning (GMXMMPBSA.make_top.cleantop) strips as solvent/ions -- it would be "
            "removed from the cleaned topology before the Receptor/Ligand index is ever "
            "applied. Rename the ligand's residue to something else."
        )

    u = mda.Universe(str(coord_file))
    resnames = u.atoms.resnames

    stray = (set(resnames) - _CLEANTOP_STRIPPED_RESNAMES) & _LOOKS_LIKE_SOLVENT_RESNAMES
    if stray:
        raise ValueError(
            f"Residue(s) {sorted(stray)} look like solvent/ions but are not names "
            "gmx_MMPBSA's own topology cleaning (GMXMMPBSA.make_top.cleantop) recognizes -- "
            "they will survive into the cleaned topology and the atom counts will no longer "
            f"match the Receptor/Ligand index. Rename them to a recognized name (one of "
            f"{sorted(_CLEANTOP_STRIPPED_RESNAMES)}) before the MMPBSA stage."
        )

    # A boolean mask rather than a selection string: some solvent/ion names
    # (e.g. "Na+", "SPC/E") are not guaranteed to survive every selection
    # grammar, so matching plain resname arrays sidesteps that entirely.
    solvent_mask = np.isin(resnames, sorted(_CLEANTOP_STRIPPED_RESNAMES))
    ligand_mask = resnames == ligand_resname

    ligand = u.atoms[ligand_mask]
    if ligand.n_atoms == 0:
        raise RuntimeError(f"No atoms with resname {ligand_resname!r} found in {coord_file}.")

    receptor = u.atoms[~solvent_mask & ~ligand_mask]
    if receptor.n_atoms == 0:
        raise RuntimeError("Protein/lipid atoms not found in system.")

    # MDAnalysis atom.ix is 0-based; GROMACS/gmx_MMPBSA index files are 1-based.
    receptor_atoms = sorted(int(i) + 1 for i in receptor.ix)
    ligand_atoms = sorted(int(i) + 1 for i in ligand.ix)

    return receptor_atoms, ligand_atoms


def write_index(
    receptor_atoms: Sequence[int],
    ligand_atoms: Sequence[int],
    index_file: Path,
) -> None:
    """Write a GROMACS index file with Receptor and Ligand atom groups.

    The groups are written as ``[ Receptor ]`` and ``[ Ligand ]`` sections,
    which are the names expected by gmx_MMPBSA when the ``-cg`` flag is used
    with group numbers 0 and 1. Raises ``RuntimeError`` if either group is
    empty, to fail early rather than silently produce an incomplete index
    file. See the GROMACS index file format documentation at
    https://manual.gromacs.org/documentation/current/reference-manual/file-formats.html#ndx
    for the format specification.
    """
    if not receptor_atoms:
        raise RuntimeError("Protein/lipid atoms not found in system.")

    if not ligand_atoms:
        raise RuntimeError("Ligand atoms not found in system.")

    with index_file.open("w") as f:
        f.write("[ Receptor ]\n")
        _write_group(f, receptor_atoms)

        f.write("\n[ Ligand ]\n")
        _write_group(f, ligand_atoms)


def _write_group(f: TextIOWrapper, atoms: Sequence[int], per_line: int = 15) -> None:
    """Write a single index group body, 15 atom indices per line."""
    for i in range(0, len(atoms), per_line):
        f.write(" ".join(map(str, atoms[i : i + per_line])) + "\n")
