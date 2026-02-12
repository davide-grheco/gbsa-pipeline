"""Module defining functions to add water box and ions to a olute system."""

from typing import Callable

import BioSimSpace as BSS


def water_function(water_model: str) -> Callable:
    """Checking which water model was selected."""
    water_params = {
        "tip3p": BSS.Solvent.tip3p,
        "tip4p": BSS.Solvent.tip4p,
        "tip5p": BSS.Solvent.tip5p,
        "spc": BSS.Solvent.spc,
        "spce": BSS.Solvent.spce,
    }

    water_model = water_model.lower().strip()

    if water_model not in water_params:
        raise ValueError("Water model is not supported. ")

    return water_params[water_model]


def box_parameters(box_size: float) -> tuple[BSS.Box, BSS.Types.Angle]:
    """Defining box shape and angles."""
    box, angles = BSS.Box.truncatedOctahedron(box_size * BSS.Units.Length.nanometer)
    return box, angles


def solvate_system(
    system: BSS._SireWrappers.System,
    solvent: Callable,
    box: BSS.Box,
    angles: BSS.Types.Angle,
) -> BSS._SireWrappers.System:
    """Adding water molecules at defined shape of box and angles."""
    return solvent(system, box=box, angles=angles)


def run_solvation(system: BSS._SireWrappers.System, water_model: str, box_size: float) -> BSS._SireWrappers.System:
    """Prepares the full solvation procedure for a solute."""
    solvent = water_function(water_model)
    box, angles = box_parameters(box_size)

    solvated = solvate_system(system, solvent, box, angles)
    return solvated
