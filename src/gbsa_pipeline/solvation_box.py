"""Validated solvation-box parameters, shared result types, and BSS solvation helper."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import TYPE_CHECKING, Any, Self

import BioSimSpace as BSS
from pydantic import Field, field_validator, model_validator

from gbsa_pipeline._paths import require_file
from gbsa_pipeline._pydantic import StrictModel

if TYPE_CHECKING:
    from pathlib import Path


@dataclass(frozen=True)
class SolvatedComplex:
    """Solvated protein-ligand complex produced by a solvation helper.

    Carries the paths to the GROMACS files written for inspection, checkpointing,
    and direct loading by downstream MD stages.
    """

    gro_file: Path
    top_file: Path

    def load_bss(self) -> Any:
        """Load this complex as a BioSimSpace System for MD stages."""
        gro_file = require_file(self.gro_file, "SolvatedComplex coordinate file")
        top_file = require_file(self.top_file, "SolvatedComplex topology file")
        import BioSimSpace as BSS  # noqa: PLC0415

        return BSS.IO.readMolecules([str(gro_file), str(top_file)])


class WaterModel(StrEnum):
    """Supported water models for solvation.

    The enum values are the user-facing strings accepted by the solvation
    configuration layer. They are intentionally lower-case so Pydantic can parse
    simple config-file values such as ``"tip3p"`` directly. Execution helpers
    should consume this enum, not raw strings, so water-model lookup tables stay
    typed and mypy can validate them. Only water models currently handled by the
    solvation helpers should be listed here.
    """

    TIP3P = "tip3p"
    TIP4P = "tip4p"
    SPC = "spc"
    SPCE = "spce"
    TIP5P = "tip5p"

    @property
    def gmx_water_name(self) -> str:
        """String accepted by gmx solvate / Modeller.addSolvent(model=...)."""
        _names = {
            WaterModel.TIP3P: "tip3p",
            WaterModel.TIP4P: "tip4pew",
            WaterModel.SPC: "spce",
            WaterModel.SPCE: "spce",
            WaterModel.TIP5P: "tip5p",
        }
        return _names[self]

    @property
    def openmm_xml(self) -> str:
        """OpenMM force-field XML file for this water model."""
        _xml = {
            WaterModel.TIP3P: "amber14/tip3p.xml",
            WaterModel.TIP4P: "amber14/tip4pew.xml",
            WaterModel.SPC: "amber14/spce.xml",
            WaterModel.SPCE: "amber14/spce.xml",
            WaterModel.TIP5P: "tip5p.xml",
        }
        return _xml[self]


class BoxShape(StrEnum):
    """Supported solvent-box shapes."""

    CUBIC = "cubic"
    TRUNCATED_OCTAHEDRON = "truncated_octahedron"

    @property
    def gmx_bt(self) -> str:
        """Value for gmx editconf -bt."""
        if self is BoxShape.TRUNCATED_OCTAHEDRON:
            return "octahedron"
        return self.value


class SolvationParams(StrictModel):
    """Validated parameters for solvent-box construction.

    This model is the input boundary for user-facing solvation settings. It
    accepts simple strings for water model and box shape, but stores them as
    typed enum values after validation. This keeps execution modules such as
    ``solvation_openmm`` free from local string coercion while still allowing
    config-style inputs. ``box_size`` may be ``None`` when padding-based box
    construction is used.
    """

    water_model: WaterModel = WaterModel.TIP3P
    shape: BoxShape = BoxShape.CUBIC
    padding: float | None = Field(default=1.0, ge=0.0)
    box_size: float | None = Field(default=8.0, gt=0.0)
    neutralize: bool = True
    ion_concentration: float | None = Field(default=None, ge=0.0)

    @field_validator("water_model", "shape", mode="before")
    @classmethod
    def _normalise_enum_input(cls, value: object) -> object:
        """Normalize simple string input before enum parsing.

        Pydantic performs the actual enum validation after this method returns.
        This validator only trims whitespace and lower-cases user-provided
        strings so config files and CLI-style inputs are slightly more forgiving.
        Existing enum values pass through unchanged. Unsupported values still
        fail through the normal Pydantic enum validation error.
        """
        if isinstance(value, str):
            return value.strip().lower()
        return value

    @model_validator(mode="after")
    def _validate_box_definition(self) -> Self:
        """Validate that either padding or explicit box size is available.

        OpenMM and BioSimSpace can construct a solvent box from a padding
        distance or from an explicit box size. A missing ``box_size`` is valid
        when ``padding`` is present. If both values are missing, downstream
        solvation cannot define the simulation box and should fail before an
        external tool is called. This keeps the failure at the parameter model
        boundary.
        """
        if self.padding is None and self.box_size is None:
            raise ValueError("Either padding or box_size must be set.")
        return self


def _run_bss_solvent(
    bss: Any,
    system: Any,
    params: SolvationParams,
    work_dir: Path | str | None,
    *,
    shell: Any = None,
    box: Any = None,
    angles: Any = None,
) -> Any:
    """Call the BSS solvent function for ``params.water_model``.

    ``ion_concentration=None`` means no added salt and ``work_dir=None`` lets
    BSS pick a directory — resolved explicitly instead of relying on BSS defaults.
    """
    solvent = getattr(bss.Solvent, params.water_model.value)
    return solvent(
        molecule=system,
        is_neutral=params.neutralize,
        ion_conc=params.ion_concentration if params.ion_concentration is not None else 0,
        work_dir=str(work_dir) if work_dir is not None else None,
        shell=shell,
        box=box,
        angles=angles if angles is not None else [90 * bss.Units.Angle.degree] * 3,
    )


def run_solvation(
    system: Any,
    params: SolvationParams,
    work_dir: Path | str | None = None,
) -> Any:
    """Solvate a molecular system with BioSimSpace.

    This helper preserves the older BioSimSpace-based solvation entry point used
    by existing tests and callers. Padding-based solvation is mapped to
    BioSimSpace's ``shell`` argument, while explicit ``box_size`` values are
    mapped to BioSimSpace box vectors. This avoids constructing boxes that are
    accidentally too small for the input system. ``ion_conc`` is passed as a
    plain molar float because BioSimSpace validates that value directly.
    """
    import BioSimSpace as BSS  # noqa: PLC0415

    if params.padding is not None:
        shell = params.padding * BSS.Units.Length.nanometer
        return _run_bss_solvent(BSS, system, params, work_dir, shell=shell)

    if params.box_size is not None:
        box, angles = _make_bss_box(BSS, params.shape, params.box_size)
        return _run_bss_solvent(BSS, system, params, work_dir, box=box, angles=angles)

    raise AssertionError("unreachable: SolvationParams guarantees padding or box_size")


def _make_bss_box(bss: Any, shape: BoxShape, size_nm: float) -> tuple[Any, Any]:
    """Create a BioSimSpace box from validated solvation parameters.

    BioSimSpace box constructors return the box vectors and angles expected by
    the solvent helpers. This local helper keeps the shape mapping in the
    BioSimSpace compatibility layer instead of leaking BioSimSpace naming into
    the parameter model. ``size_nm`` is interpreted as a nanometer box size to
    match the OpenMM solvation helper. Unsupported shapes fail explicitly even
    though the enum currently exposes only implemented values.
    """
    size = size_nm * bss.Units.Length.nanometer

    if shape is BoxShape.CUBIC:
        return bss.Box.cubic(size)

    if shape is BoxShape.TRUNCATED_OCTAHEDRON:
        return bss.Box.truncatedOctahedron(size)

    raise ValueError(f"Unsupported solvation box shape: {shape!s}")


def solvate_membrane(
    system: Any,
    params: SolvationParams,
    z_padding_nm: float,
    work_dir: Path | None = None,
) -> Any:
    """Solvate a pre-built membrane-system with BioSimSpace.

    Unlike run_solvation (isotropic padding), this preserves x&y from the inputs system's own box extending only the z-vector.
    """
    dimensions = system._sire_object.property("space").dimensions()  # Å
    x, y, z = (dimension.value() / 10 for dimension in dimensions)  # nm
    new_box = [
        x * BSS.Units.Length.nanometer,
        y * BSS.Units.Length.nanometer,
        (z + 2 * z_padding_nm) * BSS.Units.Length.nanometer,
    ]

    system.setBox(new_box, angles=[90 * BSS.Units.Angle.degree] * 3)

    return _run_bss_solvent(BSS, system, params, work_dir)
