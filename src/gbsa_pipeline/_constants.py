"""Package-wide constants shared across modules."""

from __future__ import annotations

WATER_RESIDUE_NAMES: frozenset[str] = frozenset({"HOH", "WAT", "TIP3", "TIP3P", "SOL"})

# Monatomic ion residue names. "CA" here is the calcium ion residue — distinct
# from the CA (alpha-carbon) atom name, which lives in a different namespace.
ION_RESIDUE_NAMES: frozenset[str] = frozenset({"NA", "CL", "K", "MG", "CA", "ZN"})

# Bulk solvent + ions: residues that should never be position-restrained.
SOLVENT_RESIDUE_NAMES: frozenset[str] = WATER_RESIDUE_NAMES | ION_RESIDUE_NAMES

# Bulk monovalent counter-ion names across common naming conventions (GROMACS,
# CHARMM-GUI, AMBER). Deliberately excludes divalent metals (Mg2+, Ca2+,
# Zn2+, ...): those are commonly structural/catalytic, not bulk solvent ions,
# and gmx_MMPBSA's own topology cleaning (GMXMMPBSA.make_top.cleantop) does
# not strip them either. Named separately from ION_RESIDUE_NAMES above (a
# different, broader concept used for restraint exclusion) since the two
# lists serve different, incompatible purposes -- merging them would make
# gromacs_index.py's Receptor/Ligand atom-count matching wrong.
CLEANTOP_ION_RESIDUE_NAMES: frozenset[str] = frozenset({"NA", "CL", "SOD", "Na+", "CLA", "Cl-", "POT", "K+"})
