"""Shared gemmi structure-traversal and PDB-editing helpers."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

import gemmi

from gbsa_pipeline._constants import WATER_RESIDUE_NAMES

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator
    from pathlib import Path


def iter_chains(structure: gemmi.Structure) -> Iterator[gemmi.Chain]:
    """Yield every chain in a structure, across all its models."""
    for model in structure:
        yield from model


def iter_residues(model: gemmi.Model) -> Iterator[gemmi.Residue]:
    """Yield every residue in a model, across all its chains."""
    for chain in model:
        yield from chain


def residue_is_water(residue: gemmi.Residue) -> bool:
    """Return True when the residue name is a known water residue name."""
    return residue.name.upper() in WATER_RESIDUE_NAMES


def filter_residues(structure: gemmi.Structure, keep: Callable[[gemmi.Residue], bool]) -> gemmi.Structure:
    """Delete every residue for which ``keep`` is False, across all models and chains.

    Mutates ``structure`` in place and returns it.
    """
    for chain in iter_chains(structure):
        # Iterate backwards so deleting a residue does not shift the indices
        # still to be visited.
        for i in reversed(range(len(chain))):
            if not keep(chain[i]):
                del chain[i]
    return structure


def append_chains(base: gemmi.Structure, donor: gemmi.Structure) -> None:
    """Append all chains of ``donor``'s first model into ``base``'s first model.

    Chain names that conflict with ``base`` are renamed automatically.
    """
    for chain in donor[0]:
        base[0].add_chain(chain, unique_name=True)


def write_pdb(structure: gemmi.Structure, output_pdb: Path) -> Path:
    """Write ``structure`` to ``output_pdb``, creating parent directories."""
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    structure.write_pdb(str(output_pdb))
    return output_pdb


def write_crystal_waters_pdb(protein_pdb: Path, output_pdb: Path) -> Path | None:
    """Write the crystallographic waters of ``protein_pdb`` to ``output_pdb``.

    Returns ``None`` when the input contains no water residues.
    """
    waters = filter_residues(gemmi.read_pdb(str(protein_pdb)), residue_is_water)

    if not any(len(chain) > 0 for chain in iter_chains(waters)):
        with contextlib.suppress(FileNotFoundError):
            output_pdb.unlink()
        return None

    return write_pdb(waters, output_pdb)
