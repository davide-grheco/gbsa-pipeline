from typing import  Callable

import BioSimSpace as BSS



def water_function(water_model: str) -> Callable:
    water_params = {
        "tip3p": BSS.Solvent.tip3p,
        "tip4p": BSS.Solvent.tip4p,
        "tip5p": BSS.Solvent.tip5p,
        "spc":   BSS.Solvent.spc,
        "spce":  BSS.Solvent.spce,
    }

    water_model = water_model.lower()

    if water_model not in water_params:

        raise ValueError(f"Water model is not supported. ")

    return water_model



def box_parameters(box_size: float):
    box, angles = BSS.Box.truncatedOctahedron(
        box_size * BSS.Units.Length.nanometer
    )
    return box, angles


def solvate_system(system: BSS._SireWrappers.System, solvent: Callable, box, angles) -> BSS._SireWrappers.System:
    return solvent(system, box=box, angles=angles)


def run_solvation(system: BSS._SireWrappers.System, water_model: str, box_size: float) -> BSS._SireWrappers.System:
    solvent = water_function(water_model)
    box, angles = box_parameters(box_size)
    solvated = solvate_system(system, solvent, box, angles)
    return solvated


