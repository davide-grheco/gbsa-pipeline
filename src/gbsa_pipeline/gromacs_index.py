"""Module to generate GROMACS index files for a Sire molecular system.

This module separates *selecting* which atoms belong to the Receptor/Ligand
groups (one function per system-type convention -- soluble vs. membrane) from
*writing* them out in the standard GROMACS ``.ndx`` format. The output file is
intended to be passed to ``gmx_MMPBSA`` via the ``-ci`` flag or to standard
GROMACS tools such as ``gmx trjconv`` and ``gmx make_ndx``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import sire

if TYPE_CHECKING:
    from collections.abc import Sequence
    from io import TextIOWrapper
    from pathlib import Path

    import sire.system


def select_receptor_and_ligand_atoms_by_number(
    system: sire.system.System,
    protein: sire.mol.Molecule,
    ligand: sire.mol.Molecule,
) -> tuple[list[int], list[int]]:
    """Select Receptor/Ligand atom indices by matching molecule number.

    Used for a [system] (soluble) run: protein and ligand are each a single,
    specific molecule identified by number. Molecules are matched by their
    ``number()`` identifier rather than by Python object identity, so this
    works correctly with sire systems where iteration may yield new wrapper
    objects around the same underlying C++ molecule.
    """
    protein_num = protein.number()
    ligand_num = ligand.number()

    receptor_atoms: list[int] = []
    ligand_atoms: list[int] = []
    atom_counter = 1  # GROMACS uses 1-based indexing

    for mol in system:
        natoms = len(mol.atoms())
        start = atom_counter
        end = atom_counter + natoms
        num = mol.number()
        if num == protein_num:
            receptor_atoms.extend(range(start, end))
        elif num == ligand_num:
            ligand_atoms.extend(range(start, end))
        atom_counter = end

    return receptor_atoms, ligand_atoms


def select_receptor_and_ligand_atoms_by_position(
    system: sire.system.System,
    n_solute_molecules: int,
    ligand: sire.mol.Molecule,
) -> tuple[list[int], list[int]]:
    """Select Receptor/Ligand atom indices for a [membrane] run.

    Receptor = every molecule before the ligand (protein + lipids); Ligand =
    the molecule at position ``n_solute_molecules``. Water/ions (appended
    later by solvation) are excluded from both. Molecules are identified by
    position, not number, since GROMACS round-trips only ever append new
    molecules, never reorder existing ones.

    Lipids must stay in Receptor: gmx_MMPBSA's own topology cleaning
    (``GMXMMPBSA.make_top.cleantop``) strips only water/ions from the ``-cp``
    topology, never lipids, then requires the Receptor+Ligand selection to
    cover that cleaned topology exactly.
    """
    ligand_num = ligand.number()

    receptor_atoms: list[int] = []
    ligand_atoms: list[int] = []
    atom_counter = 1  # GROMACS uses 1-based indexing

    for position, mol in enumerate(system):
        natoms = len(mol.atoms())
        start = atom_counter
        end = atom_counter + natoms
        if mol.number() == ligand_num:
            ligand_atoms.extend(range(start, end))
        elif position < n_solute_molecules:
            receptor_atoms.extend(range(start, end))
        atom_counter = end

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
