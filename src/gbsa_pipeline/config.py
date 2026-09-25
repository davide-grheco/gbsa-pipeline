"""Top-level RunConfig model for driving the pipeline from a TOML file."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Self

import tomllib
from pydantic import BaseModel, ConfigDict, Field, model_validator

from gbsa_pipeline.mdp import GromacsParams
from gbsa_pipeline.membrane import DEFAULT_LIPID_RESNAMES
from gbsa_pipeline.parametrization import ParametrizationConfig, ParametrizationInput
from gbsa_pipeline.solvation_box import BoxShape, SolvationParams


class SystemConfig(BaseModel):
    """[system] section — input files, charge settings, and system type.

    Either ``protein`` (a bare protein PDB to parametrize from scratch) or
    both ``gro_file``/``top_file`` (an already-built system, e.g. a protein
    embedded in a lipid bilayer) must be set, but not both. ``membrane=True``
    marks a pre-built system as a lipid bilayer, so the pipeline branches
    into membrane-aware ligand-merge and solvation stages; see
    :class:`MembraneConfig` for the bilayer-specific parameters that apply
    in that case.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    protein: Path | None = None
    gro_file: Path | None = None
    top_file: Path | None = None
    ligand: Path | None = None
    net_charge: int | None = None
    membrane: bool = False

    @model_validator(mode="after")
    def _validate_source(self) -> Self:
        """Exactly one of protein or (gro_file + top_file) must be set."""
        has_protein = self.protein is not None
        has_prebuilt = self.gro_file is not None or self.top_file is not None

        if has_protein == has_prebuilt:
            raise ValueError("Exactly one of protein or gro_file+top_file must be set.")
        if has_prebuilt and (self.gro_file is None or self.top_file is None):
            raise ValueError("gro_file and top_file must both be set together.")
        if self.membrane and self.protein is not None:
            raise ValueError("membrane=True requires gro_file/top_file, not a bare protein.")
        return self


class SolvationConfig(SolvationParams):
    """[solvation] section — solvent box settings with pipeline defaults."""

    shape: BoxShape = BoxShape.TRUNCATED_OCTAHEDRON
    padding: float | None = Field(default=None, ge=0.0)
    ion_concentration: float | None = Field(default=0.15, ge=0.0)


class MembraneConfig(BaseModel):
    """[membrane] section — lipid bilayer geometry.

    Only meaningful when ``[system]`` has ``membrane = true``.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    lipid_resnames: frozenset[str] = frozenset(DEFAULT_LIPID_RESNAMES)
    z_padding_nm: float = Field(default=1.5, ge=0.0)  # only used when solvation.solvate is True


class MinimizationConfig(BaseModel):
    """[minimization] section — energy minimization settings."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    nsteps: int = 10_000
    emtol: float = 10.0


class EquilibrationConfig(BaseModel):
    """[equilibration] section — NVT heating settings."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    simulation_time_ps: float = 50.0


class NptConfig(BaseModel):
    """[npt_equilibration] section — NPT equilibration time."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    simulation_time_ps: float = 100.0


class RunConfig(BaseModel):
    """Top-level configuration for a complete GBSA pipeline run.

    Load from a TOML file with :meth:`from_toml`. Each section maps to a
    nested model. The ``[md]`` section accepts any field of
    :class:`~gbsa_pipeline.mdp.GromacsParams`.

    Stages (in order):
    1. Parametrize  2. Solvate (BSS)  3. SD minimization  4. CG minimization
    5. NVT restrained heating  6. NPT restrained  7. NPT unrestrained
    8. Production MD

    Example:
    -------
    ```toml
    [system]
    protein = "protein.pdb"
    ligand  = "ligand.sdf"

    [solvation]
    water_model = "tip3p"
    padding = 1.0

    [equilibration]
    simulation_time_ps = 50.0

    [npt_equilibration]
    simulation_time_ps = 100.0

    [md]
    nsteps = 250000
    dt = 0.002
    tcoupl = "v-rescale"
    ref_t = 300.0
    ```

    A pre-built membrane system instead of a bare protein PDB:

    ```toml
    [system]
    gro_file = "system.gro"
    top_file = "system.top"
    ligand   = "ligand.sdf"
    membrane = true

    [membrane]
    z_padding_nm = 1.5

    [solvation]
    solvate = true
    ```
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    system: SystemConfig
    membrane: MembraneConfig | None = None
    forcefield: ParametrizationConfig = Field(default_factory=ParametrizationConfig)
    solvation: SolvationConfig = Field(default_factory=SolvationConfig)
    minimization: MinimizationConfig = Field(default_factory=MinimizationConfig)
    equilibration: EquilibrationConfig = Field(default_factory=EquilibrationConfig)
    npt_equilibration: NptConfig = Field(default_factory=NptConfig)
    md: GromacsParams = Field(default_factory=GromacsParams)

    @classmethod
    def from_toml(cls, path: Path) -> RunConfig:
        """Load and validate a :class:`RunConfig` from a TOML file.

        Parameters
        ----------
        path:
            Path to the ``.toml`` configuration file.

        Returns:
        -------
        RunConfig
            Validated configuration object.
        """
        with open(path, "rb") as f:
            data: dict[str, Any] = tomllib.load(f)
        return cls.model_validate(data)

    def to_parametrization_input(self, work_dir: Path) -> ParametrizationInput:
        """Build a :class:`~gbsa_pipeline.parametrization.ParametrizationInput` from this config.

        Parameters
        ----------
        work_dir:
            Directory where parametrization output files will be written.

        Returns:
        -------
        ParametrizationInput
            Ready to pass to :func:`~gbsa_pipeline.parametrization.parametrize`.

        Raises:
        ------
        ValueError
            If ``system.protein`` (pre-built systems use a different code
            path) or ``system.ligand`` is unset.
        """
        if self.system.protein is None:
            raise ValueError(
                "to_parametrization_input requires system.protein to be set; "
                "pre-built systems (gro_file/top_file) use a different code path."
            )
        if self.system.ligand is None:
            raise ValueError(
                "system.ligand must be set to run the parametrization stage. "
                "Provide a ligand SDF path in the [system] section of your config."
            )
        return ParametrizationInput(
            protein_pdb=self.system.protein,
            ligand_sdf=self.system.ligand,
            config=self.forcefield,
            net_charge=self.system.net_charge,
            work_dir=work_dir,
        )
