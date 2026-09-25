"""Tests for SDF loading in mol_utils and its tleap consumer."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from gbsa_pipeline.mol_utils import load_first_sdf_molecule
from gbsa_pipeline.tleap import sdf_formal_charge

if TYPE_CHECKING:
    from pathlib import Path

_METHANE_RECORD = """\
methane
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""

_AMMONIUM_RECORD = """\
ammonium
     RDKit          3D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  CHG  1   1   1
M  END
$$$$
"""

_GARBAGE_RECORD = """\
garbage


not a counts line
$$$$
"""


def test_load_first_sdf_molecule_skips_unparsable_records(tmp_path: Path) -> None:
    sdf = tmp_path / "ligand.sdf"
    sdf.write_text(_GARBAGE_RECORD + _METHANE_RECORD, encoding="utf-8")

    molecule = load_first_sdf_molecule(sdf)

    assert molecule.GetNumAtoms() == 1
    assert molecule.GetAtomWithIdx(0).GetSymbol() == "C"


def test_load_first_sdf_molecule_no_valid_molecule_raises(tmp_path: Path) -> None:
    sdf = tmp_path / "ligand.sdf"
    sdf.write_text(_GARBAGE_RECORD, encoding="utf-8")

    with pytest.raises(ValueError, match="Could not read any molecule"):
        load_first_sdf_molecule(sdf)


def test_sdf_formal_charge_sums_atom_charges(tmp_path: Path) -> None:
    sdf = tmp_path / "ligand.sdf"
    sdf.write_text(_AMMONIUM_RECORD, encoding="utf-8")

    assert sdf_formal_charge(sdf) == 1
