"""GBSA pipeline package.

Optimize GBSA calculations for speed
"""

from __future__ import annotations

from gbsa_pipeline._internal.cli import get_parser, main

__all__: list[str] = ["get_parser", "main"]
