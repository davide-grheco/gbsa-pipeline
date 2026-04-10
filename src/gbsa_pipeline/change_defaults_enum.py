# /home/grheco/repositorios/gbsa-pipeline/src/gbsa_pipeline/change_defaults_enum.py

"""Classes for available options in GROMACS parameters."""

from __future__ import annotations

from enum import StrEnum


class Integrator(StrEnum):
    """Integration or minimization algorithm."""

    # MD / dynamics
    LEAP_FROG = "md"
    VELOCITY_VERLET = "md-vv"
    VELOCITY_VERLET_AVEK = "md-vv-avek"
    STOCHASTIC_DYNAMICS = "sd"
    BROWNIAN_DYNAMICS = "bd"

    # energy minimization / special modes
    STEEPEST_DESCENT = "steep"
    CONJUGATE_GRADIENT = "cg"
    LBFGS = "l-bfgs"
    NORMAL_MODE = "nm"
    TEST_PARTICLE_INSERTION = "tpi"

    # backward-compatible aliases
    LEAP_FROG_STOCHASTIC = "sd"
    LANGEVIN = "bd"


class CommMode(StrEnum):
    """Center-of-mass motion removal mode."""

    LINEAR = "Linear"
    ANGULAR = "Angular"
    LINEAR_ACC_CORRECTION = "Linear-acceleration-correction"
    NONE = "None"


class NghCutoffScheme(StrEnum):
    """Neighbor-list and cutoff handling scheme."""

    VERLET = "Verlet"
    GROUP = "group"  # legacy token; no longer supported by current GROMACS


class CoulombType(StrEnum):
    """Electrostatics method for long-range Coulomb interactions."""

    CUT_OFF = "Cut-off"
    EWALD = "Ewald"
    PME = "PME"
    PM3AD = "PM3-AD"
    REACTION_FIELD = "Reaction-field"


class CoulombModifier(StrEnum):
    """Modifier applied to short-range Coulomb interactions near cutoff."""

    POTENTIAL_SHIFT = "Potential-shift"
    NONE = "None"


class VDWType(StrEnum):
    """Method for computing van der Waals interactions."""

    CUT_OFF = "Cut-off"
    PME = "PME"
    SHIFT = "Shift"
    SWITCH = "Switch"
    USER = "User"


class DispCorr(StrEnum):
    """Long-range dispersion correction mode."""

    NO = "no"
    ENERGY_PRESSURE = "EnerPres"
    ENERGY = "Energy"


class VDWModifier(StrEnum):
    """Modifier applied to Lennard-Jones interactions near cutoff."""

    POTENTIAL_SHIFT = "Potential-shift"
    NONE = "None"
    FORCE_SWITCH = "Force-switch"
    POTENTIAL_SWITCH = "Potential-switch"

    # backward-compatible alias for old typo/spelling
    POTENTIAL_SHIFT_UNDERSCORE = "Potential-shift"


class LJPMECombination(StrEnum):
    """Combination rule for Lennard-Jones PME interactions."""

    GEOMETRIC = "Geometric"
    LORENTZ_BERTHELOT = "Lorentz-Berthelot"


class EnsembleTempSetting(StrEnum):
    """Temperature control strategy for ensemble definition."""

    AUTO = "auto"
    CONSTANT = "constant"
    VARIABLE = "variable"
    NOT_AVAILABLE = "not-available"


class Thermostat(StrEnum):
    """Thermostat algorithm for temperature coupling."""

    NO = "no"
    BERENDSEN = "berendsen"
    NOSE_HOOVER = "nose-hoover"
    ANDERSEN = "andersen"
    ANDERSEN_MASSIVE = "andersen-massive"
    V_RESCALE = "v-rescale"

    # backward-compatible alias
    VRESCALE = "v-rescale"


class Barostat(StrEnum):
    """Barostat algorithm for pressure coupling."""

    NO = "no"
    BERENDSEN = "Berendsen"
    CRESCALE = "C-rescale"
    PARRINELLO_RAHMAN = "Parrinello-Rahman"
    MARTYNA_TUCKERMAN_TK = "MTTK"


class PCoupleType(StrEnum):
    """Pressure coupling geometry mode."""

    ISOTROPIC = "isotropic"
    SEMIISOTROPIC = "semiisotropic"
    ANISOTROPIC = "anisotropic"
    SURFACE_TENSION = "surface-tension"


class VelocityGeneration(StrEnum):
    """Velocity generation flag at simulation start."""

    NO = "no"
    YES = "yes"


class Constraints(StrEnum):
    """Bond/angle constraint type applied during simulation."""

    NONE = "none"
    H_BONDS = "h-bonds"
    ALL_BONDS = "all-bonds"
    H_ANGLES = "h-angles"
    ALL_ANGLES = "all-angles"

    # backward-compatible aliases
    HYDROGENS_BONDS = "h-bonds"
    HANGLES = "h-angles"
    AANGLES = "all-angles"


class ConstraintsAlgorithms(StrEnum):
    """Constraint solver algorithm."""

    LINCS = "LINCS"
    SHAKE = "SHAKE"
