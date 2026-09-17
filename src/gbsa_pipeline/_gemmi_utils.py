"""Shared gemmi structure-traversal helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Iterator

    import gemmi


def _iter_residues(model: gemmi.Model) -> Iterator[gemmi.Residue]:
    """Yield every residue in a model, across all its chains."""
    for chain in model:
        yield from chain
