"""Preparing custom GROMACS protocol for GBSA optimization."""

import BioSimSpace as BSS


class GROMACS_custom(BSS.Protocol.Custom):
    """Thin wrapper for a custom GROMACS protocol."""

    def __init__(self, config):
        super().__init__(config)
