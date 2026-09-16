"""Package-wide constants shared across modules."""

from __future__ import annotations

WATER_RESIDUE_NAMES: frozenset[str] = frozenset({"HOH", "WAT", "TIP3", "TIP3P", "SOL"})

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
