"""Package-wide constants shared across modules."""

from __future__ import annotations

WATER_RESIDUE_NAMES: frozenset[str] = frozenset({"HOH", "WAT", "TIP3", "TIP3P", "SOL"})

# Monatomic ion residue names. "CA" here is the calcium ion residue — distinct
# from the CA (alpha-carbon) atom name, which lives in a different namespace.
ION_RESIDUE_NAMES: frozenset[str] = frozenset({"NA", "CL", "K", "MG", "CA", "ZN"})

# Bulk solvent + ions: residues that should never be position-restrained.
SOLVENT_RESIDUE_NAMES: frozenset[str] = WATER_RESIDUE_NAMES | ION_RESIDUE_NAMES
