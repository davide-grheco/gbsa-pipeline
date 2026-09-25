"""Tests for the shared BSS solvent-call helper."""

from __future__ import annotations

from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

from gbsa_pipeline.solvation_box import SolvationParams, _run_bss_solvent

if TYPE_CHECKING:
    from pathlib import Path


def _fake_bss(captured: dict[str, Any]) -> SimpleNamespace:
    return SimpleNamespace(
        Solvent=SimpleNamespace(tip3p=lambda **kwargs: captured.update(kwargs)),
        Units=SimpleNamespace(Angle=SimpleNamespace(degree=1)),
    )


def test_run_bss_solvent_resolves_unset_optionals_to_bss_defaults() -> None:
    captured: dict[str, Any] = {}
    system = object()
    params = SolvationParams(neutralize=False, ion_concentration=None)

    _run_bss_solvent(_fake_bss(captured), system, params, None)

    assert captured == {
        "molecule": system,
        "is_neutral": False,
        "ion_conc": 0,
        "work_dir": None,
        "shell": None,
        "box": None,
        "angles": [90, 90, 90],
    }


def test_run_bss_solvent_passes_explicit_arguments_through(tmp_path: Path) -> None:
    captured: dict[str, Any] = {}
    params = SolvationParams(neutralize=True, ion_concentration=0.15)

    _run_bss_solvent(_fake_bss(captured), object(), params, tmp_path, box="box", angles="angles")

    assert captured["ion_conc"] == 0.15
    assert captured["work_dir"] == str(tmp_path)
    assert captured["box"] == "box"
    assert captured["angles"] == "angles"
