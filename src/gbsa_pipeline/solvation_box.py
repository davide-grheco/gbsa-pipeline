"""Validated solvation-box parameters and BioSimSpace solvation helper."""

from __future__ import annotations

from enum import StrEnum
from typing import TYPE_CHECKING, Any, Self

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

if TYPE_CHECKING:
    from pathlib import Path


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


class BoxShape(StrEnum):
    """Supported solvent-box shapes.

    The enum values are the user-facing strings accepted by the solvation
    parameter model. Execution code should consume this enum directly instead
    of repeating string normalization in multiple places. This keeps validation
    at the configuration boundary and makes invalid values fail before solvation
    starts. Only shapes implemented by the downstream solvation helpers should
    be exposed here.
    """

    CUBIC = "cubic"
    TRUNCATED_OCTAHEDRON = "truncated_octahedron"


class SolvationParams(BaseModel):
    """Validated parameters for solvent-box construction.

    This model is the input boundary for user-facing solvation settings. It
    accepts simple strings for water model and box shape, but stores them as
    typed enum values after validation. This keeps execution modules such as
    ``solvation_openmm`` free from local string coercion while still allowing
    config-style inputs. ``box_size`` may be ``None`` when padding-based box
    construction is used.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    water_model: WaterModel = WaterModel.TIP3P
    shape: BoxShape = BoxShape.CUBIC
    padding: float | None = Field(default=1.0, ge=0.0)
    box_size: float | None = Field(default=8.0, gt=0.0)
    neutralize: bool = True
    ion_concentration: float | None = Field(default=None, ge=0.0)

    def __init__(
        self,
        *,
        water_model: WaterModel | str = WaterModel.TIP3P,
        shape: BoxShape | str = BoxShape.CUBIC,
        padding: float | None = 1.0,
        box_size: float | None = 8.0,
        neutralize: bool = True,
        ion_concentration: float | None = None,
        **data: Any,
    ) -> None:
        """Create validated solvation parameters from enum or string input.

        The explicit constructor keeps existing mypy-checked call sites working
        when they pass strings such as ``"tip3p"`` or ``"cubic"``. Pydantic
        still performs the real validation and stores enum values on the model.
        Extra keyword arguments are passed through so Pydantic can report
        forbidden fields through the normal model validation path. Downstream
        code can therefore use ``params.water_model`` and ``params.shape`` as
        typed enums.
        """
        super().__init__(
            water_model=water_model,
            shape=shape,
            padding=padding,
            box_size=box_size,
            neutralize=neutralize,
            ion_concentration=ion_concentration,
            **data,
        )

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

    solvent = _get_bss_solvent_function(BSS, params.water_model)

    kwargs: dict[str, Any] = {
        "molecule": system,
        "is_neutral": params.neutralize,
    }

    if params.padding is not None:
        kwargs["shell"] = params.padding * BSS.Units.Length.nanometer
    else:
        if params.box_size is None:
            raise ValueError("BioSimSpace run_solvation requires params.box_size when padding is None.")

        box, angles = _make_bss_box(BSS, params.shape, params.box_size)
        kwargs["box"] = box
        kwargs["angles"] = angles

    if params.ion_concentration is not None:
        kwargs["ion_conc"] = params.ion_concentration

    if work_dir is not None:
        kwargs["work_dir"] = str(work_dir)

    return solvent(**kwargs)


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


def _get_bss_solvent_function(bss: Any, water_model: WaterModel) -> Any:
    """Return the BioSimSpace solvent function for a validated water model.

    The parameter model stores water models as ``WaterModel`` values, while
    BioSimSpace exposes solvent builders as functions on ``BSS.Solvent``. This
    helper keeps that API mapping local to the BioSimSpace compatibility path.
    It avoids string-indexed dictionaries so mypy can still track enum use in
    the OpenMM solvation implementation. Unsupported models fail explicitly.
    """
    if water_model is WaterModel.TIP3P:
        return bss.Solvent.tip3p

    if water_model is WaterModel.TIP4P:
        return bss.Solvent.tip4p

    if water_model is WaterModel.SPC:
        return bss.Solvent.spc

    if water_model is WaterModel.SPCE:
        return bss.Solvent.spce

    if water_model is WaterModel.TIP5P:
        return bss.Solvent.tip5p

    raise ValueError(f"Unsupported water model: {water_model!s}")
