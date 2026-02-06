# Preparing minimization protocol and running

from __future__ import annotations

import BioSimSpace as BSS


def preparing_protocol(nsteps: int):
    """Create a minimisation protocol."""
    return BSS.Protocol.Minimisation(steps=nsteps)


def setting_minimization(protocol, solvated: BSS._SireWrappers.System):
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

    protocol = preparing_protocol(nsteps)
    min_process = BSS.Process.Gromacs(protocol, solvated)

    min_process.start()
    min_process.wait()
    minimized = min_process.getSystem()

    return minimized
