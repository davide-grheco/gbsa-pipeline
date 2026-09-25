"""Shared Pydantic base model for the pipeline's input-model convention."""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict


class StrictModel(BaseModel):
    """Frozen, extra-forbidding base model with validated defaults.

    Unknown fields are rejected, instances are immutable after validation, and
    defaults go through the same validation as user-supplied values.
    """

    model_config = ConfigDict(frozen=True, extra="forbid", validate_default=True)
