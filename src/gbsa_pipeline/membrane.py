"""Membrane protein module.

Estimate membrane geometry parameters such as ``mthick`` and ``mctrdz``
for PB calculations directly from lipid phosphate atoms.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from MDAnalysis.analysis.leaflet import LeafletFinder

from gbsa_pipeline.mmbsa import PBParams
from gbsa_pipeline.parametrization import parametrize_ligand_only

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    import MDAnalysis as mda

logger = logging.getLogger(__name__)
__all__ = [
    "DEFAULT_LIPID_RESNAMES",
    "MembraneGeometry",
    "estimate_membrane_geometry",
    "parametrize_ligand_only",
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
    cutoff: float = 15.0,
) -> MembraneGeometry:
    """Measure bilayer parameters from lipid phosphate atoms.

    ``universe`` must already be parsed (e.g. via ``MDAnalysis.Universe(path)``)
    so that file-parsing errors are handled by the caller.

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
    phosphates = universe.select_atoms(f"resname {resnames} and name P*")

    if len(phosphates) == 0:
        raise ValueError(
            f"No phosphate atoms belonging to {sorted(set(lipid_resnames))} were found. "
            "Check the lipid residue names and pass lipid_resnames explicitly."
        )

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
    lateral = float(np.linalg.norm(separation[:_Z_AXIS]))
    normal_component = abs(float(separation[_Z_AXIS]))

    if lateral > normal_component:
        raise ValueError("Rotate model so lipid layer along z axis")

    mthick = normal_component

    return MembraneGeometry(
        mctrdz=float(phosphates.positions[:, _Z_AXIS].mean()),
        mthick=mthick,
        n_phosphates=len(phosphates),
    )
