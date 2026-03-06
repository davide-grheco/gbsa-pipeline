"""Command-line interface for the GBSA pipeline."""

from __future__ import annotations

import argparse
import faulthandler
import logging
import sys
from pathlib import Path

from gbsa_pipeline.config import RunConfig
from gbsa_pipeline.pipeline import run_pipeline


def main(argv: list[str] | None = None) -> None:
    """Entry point for the ``gbsa-pipeline`` command.

    Parameters
    ----------
    argv:
        Argument list (defaults to ``sys.argv[1:]`` when ``None``).
    """
    parser = argparse.ArgumentParser(
        prog="gbsa-pipeline",
        description="Run a GBSA MD simulation from a TOML config file.",
    )
    parser.add_argument(
        "config",
        type=Path,
        help="Path to the TOML configuration file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory (default: <config_dir>/gbsa_output).",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )

    args = parser.parse_args(argv)

    faulthandler.enable(file=sys.stderr)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        stream=sys.stderr,
    )

    config = RunConfig.from_toml(args.config)
    output_dir: Path = args.output_dir or args.config.parent / "gbsa_output"
    run_pipeline(config, output_dir)
