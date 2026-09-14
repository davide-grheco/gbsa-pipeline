"""Membrane protein module.

Estimate membrane geometry parameters such as ``mthick`` and ``mctrdz``
for PB calculations directly from lipid phosphate atoms.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

# from pathlib import Path
# from typing import Any
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from collections.abc import Sequence
    from typing import Any

    from patlib import Path

import gemmi

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
) -> MembraneGeometry:
    """Measure bilayer parameters from lipid phosphate atoms.

    The phosphate atoms are separated into upper and lower leaflets using
    their overall mean z-coordinate.
    """
    resnames = frozenset(lipid_resnames)
    struct = gemmi.read_structure(str(structure))

    phosphate_z: list[float] = [
        atom.pos.z
        for model in struct
        for chain in model
        for residue in chain
        if residue.name.strip() in resnames
        for atom in residue
        if atom.name.strip().upper().startswith("P")
    ]

    if not phosphate_z:
        raise ValueError(
            f"No phosphate atoms belonging to {sorted(resnames)} were found "
            f"in {structure}. Check the lipid residue names and pass "
            "lipid_resnames explicitly."
        )

    mean_z = sum(phosphate_z) / len(phosphate_z)

    upper = [z for z in phosphate_z if z >= mean_z]
    lower = [z for z in phosphate_z if z < mean_z]

    if len(upper) < _MIN_PHOSPHATES_PER_LEAFLET or len(lower) < _MIN_PHOSPHATES_PER_LEAFLET:
        raise ValueError(
            "Phosphate atoms did not split into two comparable leaflets. "
            "The structure may not contain a symmetric bilayer, or the "
            "lipid residue names may be incorrect."
        )

    upper_mean = sum(upper) / len(upper)
    lower_mean = sum(lower) / len(lower)
    mthick = abs(upper_mean - lower_mean)

    return MembraneGeometry(
        mctrdz=mean_z,
        mthick=mthick,
        n_phosphates=len(phosphate_z),
    )
