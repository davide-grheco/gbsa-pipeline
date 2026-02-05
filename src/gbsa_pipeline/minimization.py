# Preparing minimization protocol and running

from __future__ import annotations

import BioSimSpace as BSS


def preparing_protocol(nsteps: int):
    """Create a minimisation protocol."""
    return BSS.Protocol.Minimisation(steps=nsteps)


def setting_process(protocol, solvated: BSS._SireWrappers.System):
    """Run minimisation and return the minimised system."""
    min_process = BSS.MD.run(
        solvated,
        protocol,
        engine="auto",
        gpu_support=False,
        auto_start=True,
        name="min",
        property_map={},
    )
    return min_process.getSystem(block=True)


def run_minimization(nsteps: int, solvated: BSS._SireWrappers.System):
    """Convenience wrapper: build protocol + run minimisation."""
    protocol = preparing_protocol(nsteps)
    minimised_system = setting_process(protocol, solvated)
    return minimised_system
