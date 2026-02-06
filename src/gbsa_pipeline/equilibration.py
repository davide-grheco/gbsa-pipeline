# Preparing minimization protocol and running

from __future__ import annotations

import BioSimSpace as BSS


heating_protocol = BSS.Protocol.Equilibration(
    temperature_start=0 * BSS.Units.Temperature.kelvin,
    temperature_end=300 * BSS.Units.Temperature.kelvin,
    restraint="backbone",
)


def setting_heating(protocol, minimized: BSS._SireWrappers.System):
    """Run heating and return the equlibrated system."""
    heat_process = BSS.MD.run(
        minimized,
        heating_protocol,
        engine="Gromacs",
        gpu_support=False,
        auto_start=True,
        name="NVT",
        property_map={},
    )
    return heat_process.getSystem(block=True)


def run_heating(nsteps: int, minimized: BSS._SireWrappers.System, proces):

    protocol = heating_protocol()
    heating_process = BSS.Process.Gromacs(protocol, minimized)

    heating_process.start()
    heating_process.wait()
    equilibrated = heating_process.getSystem()

    return equilibrated

    getOutput(name="NVT.txt", block="AUTO", file_link=False)
