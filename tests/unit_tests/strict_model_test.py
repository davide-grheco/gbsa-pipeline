"""Tests for the shared StrictModel base."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from gbsa_pipeline._pydantic import StrictModel


class _ExampleModel(StrictModel):
    value: int = 1


def test_strict_model_rejects_unknown_fields() -> None:
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        _ExampleModel(value=2, unknown=3)  # type: ignore[call-arg]


def test_strict_model_is_frozen() -> None:
    model = _ExampleModel()
    with pytest.raises(ValidationError, match="frozen"):
        model.value = 5


def test_strict_model_validates_defaults() -> None:
    class _BadDefault(StrictModel):
        value: int = "not an int"  # type: ignore[assignment]

    with pytest.raises(ValidationError):
        _BadDefault()
