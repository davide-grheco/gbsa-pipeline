"""Format-agnostic spatial helpers: distance queries and clash detection."""

from __future__ import annotations

from collections.abc import Hashable
from typing import TypeVar

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

_Coords = tuple[float, float, float]
_ResKey = tuple[int, str]  # canonical key for GRO residues; kept for _gro_io.py annotations

_K = TypeVar("_K", bound=Hashable)


def contact_pairs(
    a_coords: np.ndarray,
    b_coords: np.ndarray,
    cutoff: float,
) -> list[tuple[int, int, float]]:
    """Return ``(i, j, distance)`` for every pair of atoms within ``cutoff``.

    Both arrays are ``(N, 3)`` in the same unit system.  The boundary is
    inclusive (pairs at exactly ``cutoff`` are included).
    """
    if a_coords.size == 0 or b_coords.size == 0:
        return []
    dists = cdist(a_coords, b_coords)
    return [(int(i), int(j), float(dists[i, j])) for i, j in np.argwhere(dists <= cutoff)]


def _find_clashing_residues(
    candidate_entries: list[tuple[_K, _Coords]],
    reference_coords: list[_Coords],
    cutoff: float,
) -> set[_K]:
    """Return residue keys whose atoms come within ``cutoff`` of any reference coordinate.

    ``candidate_entries`` is a list of ``(residue_key, coords)`` pairs (e.g. water
    molecules).  ``reference_coords`` is the set of coordinates to check against
    (e.g. solute heavy atoms).  Both are in the same unit system; ``cutoff`` uses
    the same units.
    """
    if cutoff <= 0:
        raise ValueError("cutoff must be positive.")
    if not reference_coords or not candidate_entries:
        return set()
    tree = cKDTree(np.asarray(reference_coords))
    keys = [key for key, _ in candidate_entries]
    cand_arr = np.asarray([coords for _, coords in candidate_entries])
    hits = tree.query_ball_point(cand_arr, cutoff)
    return {keys[i] for i, neighbors in enumerate(hits) if neighbors}
