"""Membrane protein module.

Estimate membrane geometry parameters such as ``mthick`` and ``mctrdz``
for PB calculations directly from lipid phosphate atoms.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.leaflet import LeafletFinder

from gbsa_pipeline._gemmi_utils import _iter_residues
from gbsa_pipeline.mmbsa import PBParams

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    import gemmi

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

# index opf the Z coordinate in a [x,y,z] position arraz
_Z_AXIS = 2


def _is_phosphate_atom(atom: gemmi.Atom) -> bool:
    """Whether an atom is a phosphorus atom."""
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
    structure: gemmi.Structure,
    lipid_resnames: Sequence[str] = tuple(DEFAULT_LIPID_RESNAMES),
    cutoff: float = 15.0,
) -> MembraneGeometry:
    """Measure bilayer parameters from lipid phosphate atoms.

    Phosphorus atoms are identified by comparing each atom's element symbol
    to "P" and grouped into two leaflets using MDAnalysis's LeafletFinder, a
    distance-based graph clustering.
    """
    resnames = frozenset(lipid_resnames)

    coords = [
        [atom.pos.x, atom.pos.y, atom.pos.z]
        for model in structure
        for residue in _iter_residues(model)
        if residue.name.strip() in resnames
        for atom in residue
        if _is_phosphate_atom(atom)
    ]

    if not coords:
        raise ValueError(
            f"No phosphate atoms belonging to {sorted(resnames)} were found "
            f"in structure '{structure.name}'. Check the lipid residue names "
            "and pass lipid_resnames explicitly."
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
    mthick = abs(float(upper.positions[:, _Z_AXIS].mean()) - float(lower.positions[:, _Z_AXIS].mean()))

    return MembraneGeometry(
        mctrdz=float(positions[:, _Z_AXIS].mean()),
        mthick=mthick,
        n_phosphates=len(positions),
    )
