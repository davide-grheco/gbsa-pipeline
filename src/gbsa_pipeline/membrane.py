"""Membrane protein module.

Estimate membrane geometry parameters such as ``mthick`` and ``mctrdz``
for PB calculations directly from lipid phosphate atoms.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import gemmi
import MDAnalysis as mda  # noqa: N813 -- `mda` is the standard alias used throughout MDAnalysis's own docs
import numpy as np
from MDAnalysis.analysis.leaflet import LeafletFinder

if TYPE_CHECKING:
    from collections.abc import Sequence
    from pathlib import Path
    from typing import Any


from gbsa_pipeline.mmbsa import PBParams

logger = logging.getLogger(__name__)


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

# Minimum meaningful number of lipids per leaflet.
_MIN_PHOSPHATES_PER_LEAFLET = 5

# A bilayer has exactly two leaflets.
_N_LEAFLETS = 2


def _is_phosphate_atom(atom: gemmi.Atom) -> bool:
    """Whether an atom is a phosphate atom."""
    return atom.element.name == "P"


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
    structure: Path,
    lipid_resnames: Sequence[str] = tuple(DEFAULT_LIPID_RESNAMES),
    cutoff: float = 15.0,
) -> MembraneGeometry:
    """Measure bilayer parameters from lipid phosphate atoms.

    Phosphate atoms are found by gemmi (robust against e.g. Pt/P confusion)
    and grouped into two leaflets using MDAnalysis's LeafletFinder, a
    distance-based graph clustering.
    """
    resnames = frozenset(lipid_resnames)
    struct = gemmi.read_structure(str(structure))

    coords = [
        [atom.pos.x, atom.pos.y, atom.pos.z]
        for model in struct
        for chain in model
        for residue in chain
        if residue.name.strip() in resnames
        for atom in residue
        if _is_phosphate_atom(atom)
    ]

    if not coords:
        raise ValueError(
            f"No phosphate atoms belonging to {sorted(resnames)} were found "
            f"in {structure}. Check the lipid residue names and pass "
            "lipid_resnames explicitly."
        )

    positions = np.array(coords, dtype=np.float32)
    universe = mda.Universe.empty(len(positions), trajectory=True)
    universe.atoms.positions = positions

    finder = LeafletFinder(universe, universe.atoms, cutoff=cutoff)
    groups = finder.groups()

    if len(groups) != _N_LEAFLETS or min(len(group) for group in groups) < _MIN_PHOSPHATES_PER_LEAFLET:
        raise ValueError(
            "Phosphate atoms did not split into two comparable leaflets. "
            "The structure may not contain a symmetric bilayer, or the "
            "lipid residue names may be incorrect."
        )

    upper, lower = groups
    mthick = abs(float(upper.positions[:, 2].mean()) - float(lower.positions[:, 2].mean()))

    return MembraneGeometry(
        mctrdz=float(positions[:, 2].mean()),
        mthick=mthick,
        n_phosphates=len(positions),
    )
