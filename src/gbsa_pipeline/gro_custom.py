"""Preparing custom GROMACS protocol for GBSA optimization."""

import logging
from pathlib import Path

import BioSimSpace as BSS


class GromacsCustom(BSS.Protocol.Custom):
    """Thin wrapper for a custom GROMACS protocol."""

    def __init__(self, config: str | Path):
        """Constructin an object."""
        super().__init__(config)


def run_gro_custom(parameters: str | Path, system: BSS._SireWrappers.System) -> BSS._SireWrappers.System:
    """Creates custom GROMACS protocol and runs simulation."""
    custom_protocol = GromacsCustom(parameters)
    logging.info("Created Protocol")
    proces = BSS.Process.Gromacs(system, protocol=custom_protocol)
    logging.info("Starting proces")
    proces.start()
    proces.wait()
    logging.info("Proces finished")
    customized = proces.getSystem(block=True)

    return customized
