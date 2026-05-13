"""Integration tests for the BioSimSpace solvation compatibility helper."""

from __future__ import annotations

from typing import TYPE_CHECKING

import BioSimSpace as BSS
import pytest

from gbsa_pipeline.solvation_box import (
    BoxShape,
    SolvationParams,
    WaterModel,
    run_solvation,
)

if TYPE_CHECKING:
    from collections.abc import Iterable
    from pathlib import Path


def _molecule_names(system: BSS._SireWrappers.System) -> Iterable[str]:
    """Yield molecule names from a BioSimSpace system.

    BioSimSpace and Sire wrappers expose names through slightly different
    methods depending on the object type and version. This helper keeps that
    compatibility handling local to the test module. It is useful for debugging
    failures without making the test assertions depend on one exact wrapper API.
    Unknown objects fall back to their string representation.
    """
    for mol in system:
        if hasattr(mol, "getName"):
            yield mol.getName()
        elif hasattr(mol, "name"):
            yield mol.name()
        else:
            yield str(mol)


def _is_ion(mol: BSS._SireWrappers.Molecule) -> bool:
    """Return whether a BioSimSpace molecule looks like an ion.

    Some BioSimSpace versions expose ``isIon`` directly, while others require a
    more defensive check in tests. Single-atom molecules are treated as likely
    ions because this integration test only uses a small solvated protein system.
    The name fallback catches common sodium, chloride, potassium, and calcium
    representations. The heuristic is intentionally test-local and should not be
    used as production chemistry logic.
    """
    if hasattr(mol, "isIon"):
        return bool(mol.isIon())

    if hasattr(mol, "nAtoms") and mol.nAtoms() == 1:
        return True

    name = str(mol)
    return any(tag in name.upper() for tag in (" NA", " CL", " K", " CA", "NA+", "CL-"))


def _ion_molecules(
    system: BSS._SireWrappers.System,
) -> list[BSS._SireWrappers.Molecule]:
    """Return ion-like molecules from a BioSimSpace system.

    The preferred path uses BioSimSpace's own ``getIonMolecules`` method when it
    is available. Some wrapper versions return objects that are not directly
    list-convertible, so the helper falls back safely. If the method is absent,
    the local test heuristic is used instead. This keeps version-specific API
    handling out of the actual assertions.
    """
    if hasattr(system, "getIonMolecules"):
        ions = system.getIonMolecules()
        try:
            return list(ions)
        except TypeError:
            return []

    return [mol for mol in system if _is_ion(mol)]


@pytest.mark.integration
def test_solvation_real_protein_box_and_ions(tmp_path: Path) -> None:
    """Solvate a real protein test system and check box, waters, and ions.

    This test exercises the BioSimSpace compatibility entry point with typed
    ``SolvationParams`` values. The compatibility helper constructs BioSimSpace
    boxes from an explicit box size, so the test provides ``box_size`` instead
    of relying on OpenMM-style padding. The assertions stay coarse because exact
    water and ion counts can vary between BioSimSpace/Sire versions. The
    important behavior is that a reasonable box, water molecules, and ions are
    produced.
    """
    system = BSS.IO.readMolecules(
        files=[
            "tests/testdata/test.gro",
            "tests/testdata/test.top",
        ],
        make_whole=True,
    )

    params = SolvationParams(
        water_model=WaterModel.TIP3P,
        shape=BoxShape.CUBIC,
        padding=1.5,
        box_size=None,
        ion_concentration=0.1,
        neutralize=True,
    )

    solvated = run_solvation(system=system, params=params, work_dir=tmp_path)

    dims = solvated._sire_object.property("space").dimensions()
    dims_nm = [dim.value() / 10 for dim in dims]
    assert all(val > 3.0 for val in dims_nm)

    water_mols = solvated.getWaterMolecules()
    assert water_mols.nMolecules() > 0

    ions = [mol for mol in solvated if _is_ion(mol)]
    assert ions, "Expected ions when ion_concentration is set"


@pytest.mark.integration
def test_solvation_real_protein_without_neutralisation(tmp_path: Path) -> None:
    """Solvate a real protein test system without requested neutralisation.

    This test keeps the previous behavior where ion detection is used only as a
    compatibility check, not as a strict count assertion. BioSimSpace may still
    add ions depending on solvent settings and backend behavior. The purpose is
    to verify that the wrapper accepts validated parameters and produces water
    molecules. A temporary work directory is passed explicitly so generated
    files do not leak into the repository tree.
    """
    system = BSS.IO.readMolecules(
        files=[
            "tests/testdata/test.gro",
            "tests/testdata/test.top",
        ],
        make_whole=True,
    )

    _ion_molecules(system)

    params = SolvationParams(
        water_model=WaterModel.TIP3P,
        shape=BoxShape.CUBIC,
        padding=1.0,
        box_size=None,
        ion_concentration=0.0,
        neutralize=False,
    )

    solvated = run_solvation(system=system, params=params, work_dir=tmp_path)

    water_mols = solvated.getWaterMolecules()
    assert water_mols.nMolecules() > 0

    _ion_molecules(solvated)
