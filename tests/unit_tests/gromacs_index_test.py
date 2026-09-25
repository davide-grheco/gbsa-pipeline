"""Unit tests for gromacs_index: resname-based atom selection and index-file writing."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from gbsa_pipeline.gromacs_index import (
    identify_ligand_resname,
    select_receptor_and_ligand_atoms,
    write_index,
)

if TYPE_CHECKING:
    from pathlib import Path

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_index(path: Path) -> str:
    return path.read_text()


def _gro_line(resid: int, resname: str, atomname: str, atomnum: int) -> str:
    """One fixed-width GROMACS .gro atom line (resid, resname, atomname, atomnum, xyz)."""
    x = 0.1 * atomnum
    return f"{resid:5d}{resname:<5s}{atomname:>5s}{atomnum:5d}{x:8.3f}{0.0:8.3f}{0.0:8.3f}"


def _write_gro(tmp_path: Path, atoms: list[tuple[int, str, str]]) -> Path:
    """Write a minimal, MDAnalysis-readable .gro file.

    ``atoms`` is a list of ``(resid, resname, atomname)``, one per atom, in
    file order.
    """
    lines = ["Test", f"{len(atoms):5d}"]
    for i, (resid, resname, atomname) in enumerate(atoms, start=1):
        lines.append(_gro_line(resid, resname, atomname, i))
    lines.append("   5.00000   5.00000   5.00000")

    gro_file = tmp_path / "test.gro"
    gro_file.write_text("\n".join(lines) + "\n")
    return gro_file


class _FakeResName:
    def __init__(self, name: str) -> None:
        self._name = name

    def value(self) -> str:
        return self._name


class _FakeResidue:
    def __init__(self, name: str) -> None:
        self._name = name

    def name(self) -> _FakeResName:
        return _FakeResName(self._name)


class _FakeMol:
    def __init__(self, resnames: list[str]) -> None:
        self._resnames = resnames

    def residues(self) -> list[_FakeResidue]:
        return [_FakeResidue(n) for n in self._resnames]


class _FakeSystem:
    def __init__(self, molecules: list[_FakeMol]) -> None:
        self._mols = molecules

    def __iter__(self):
        return iter(self._mols)


# ---------------------------------------------------------------------------
# identify_ligand_resname
# ---------------------------------------------------------------------------


def test_identify_ligand_resname_soluble_system() -> None:
    """Protein + ligand + water + ions -- ligand is the only non-protein, non-solvent residue."""
    system = _FakeSystem(
        [
            _FakeMol(["ALA", "GLY", "LEU"]),
            _FakeMol(["MOL"]),
            _FakeMol(["SOL"]),
            _FakeMol(["NA"]),
        ]
    )

    assert identify_ligand_resname(system) == "MOL"


def test_identify_ligand_resname_membrane_reduced_system() -> None:
    """Protein + ligand only (lipids/water already stripped) -- still finds the ligand."""
    system = _FakeSystem([_FakeMol(["ASP", "VAL"]), _FakeMol(["LIG"])])

    assert identify_ligand_resname(system) == "LIG"


def test_identify_ligand_resname_does_not_depend_on_position() -> None:
    """Ligand first, protein second -- still correctly identified (no positional assumption)."""
    system = _FakeSystem([_FakeMol(["LIG"]), _FakeMol(["ASP", "VAL"]), _FakeMol(["SOL"])])

    assert identify_ligand_resname(system) == "LIG"


def test_identify_ligand_resname_raises_when_no_candidate() -> None:
    """Everything looks like protein or solvent -- no ligand-like residue found."""
    system = _FakeSystem([_FakeMol(["ALA", "GLY"]), _FakeMol(["SOL"]), _FakeMol(["NA"])])

    with pytest.raises(ValueError, match="Could not identify"):
        identify_ligand_resname(system)


def test_identify_ligand_resname_raises_when_ambiguous() -> None:
    """Two distinct non-protein, non-solvent residue names -- can't tell which is the ligand."""
    system = _FakeSystem([_FakeMol(["ALA"]), _FakeMol(["LIG"]), _FakeMol(["COFACTOR"])])

    with pytest.raises(ValueError, match="Ambiguous"):
        identify_ligand_resname(system)


# ---------------------------------------------------------------------------
# select_receptor_and_ligand_atoms
# ---------------------------------------------------------------------------


def test_select_receptor_and_ligand_atoms_soluble_convention(tmp_path: Path) -> None:
    """[system] (soluble) convention: protein + ligand + water + ions."""
    gro_file = _write_gro(
        tmp_path,
        [
            (1, "ALA", "CA"),
            (1, "ALA", "CB"),
            (2, "GLY", "CA"),
            (3, "LIG", "C1"),
            (3, "LIG", "C2"),
            (4, "SOL", "OW"),
            (5, "NA", "NA"),
        ],
    )

    receptor, ligand = select_receptor_and_ligand_atoms(gro_file, "LIG")

    assert receptor == [1, 2, 3]
    assert ligand == [4, 5]


def test_select_receptor_and_ligand_atoms_membrane_convention_lipids_in_receptor(tmp_path: Path) -> None:
    """[membrane] convention: lipids must join Receptor, not be dropped or excluded."""
    gro_file = _write_gro(
        tmp_path,
        [
            (1, "ALA", "CA"),
            (1, "ALA", "CB"),
            (2, "GLY", "CA"),
            (3, "POP", "P8"),
            (3, "POP", "C1"),
            (4, "POP", "P8"),
            (4, "POP", "C1"),
            (5, "LIG", "C1"),
            (5, "LIG", "C2"),
        ],
    )

    receptor, ligand = select_receptor_and_ligand_atoms(gro_file, "LIG")

    assert receptor == [1, 2, 3, 4, 5, 6, 7]
    assert ligand == [8, 9]


def test_select_receptor_and_ligand_atoms_does_not_assume_solvent_comes_last(tmp_path: Path) -> None:
    """Water placed BEFORE the ligand in the file must still be excluded correctly.

    The previous position-based implementation assumed solvent/ions always
    sit after the ligand; this must hold regardless of file order.
    """
    gro_file = _write_gro(
        tmp_path,
        [
            (1, "SOL", "OW"),
            (2, "NA", "NA"),
            (3, "ALA", "CA"),
            (4, "GLY", "CA"),
            (5, "LIG", "C1"),
            (6, "SOL", "OW"),
        ],
    )

    receptor, ligand = select_receptor_and_ligand_atoms(gro_file, "LIG")

    assert receptor == [3, 4]
    assert ligand == [5]


def test_select_receptor_and_ligand_atoms_raises_on_hoh_water(tmp_path: Path) -> None:
    """A prebuilt system naming crystallographic water "HOH" is rejected with a clear message."""
    gro_file = _write_gro(
        tmp_path,
        [
            (1, "ALA", "CA"),
            (1, "ALA", "CB"),
            (2, "GLY", "CA"),
            (3, "LIG", "C1"),
            (3, "LIG", "C2"),
            (4, "HOH", "O"),
        ],
    )

    with pytest.raises(ValueError, match="HOH"):
        select_receptor_and_ligand_atoms(gro_file, "LIG")


def test_select_receptor_and_ligand_atoms_raises_when_ligand_resname_is_solvent_name(tmp_path: Path) -> None:
    """A ligand accidentally named "SOL" collides with a name cleantop() strips."""
    gro_file = _write_gro(tmp_path, [(1, "ALA", "CA"), (2, "SOL", "OW")])

    with pytest.raises(ValueError, match="SOL"):
        select_receptor_and_ligand_atoms(gro_file, "SOL")


def test_select_receptor_and_ligand_atoms_raises_when_ligand_absent(tmp_path: Path) -> None:
    gro_file = _write_gro(tmp_path, [(1, "ALA", "CA"), (2, "SOL", "OW")])

    with pytest.raises(RuntimeError, match="LIG"):
        select_receptor_and_ligand_atoms(gro_file, "LIG")


def test_select_receptor_and_ligand_atoms_raises_when_receptor_empty(tmp_path: Path) -> None:
    """Everything besides the ligand is recognized solvent -- no receptor atoms remain."""
    gro_file = _write_gro(tmp_path, [(1, "LIG", "C1"), (2, "SOL", "OW")])

    with pytest.raises(RuntimeError, match="Protein"):
        select_receptor_and_ligand_atoms(gro_file, "LIG")


# ---------------------------------------------------------------------------
# write_index
# ---------------------------------------------------------------------------


def test_write_index_writes_receptor_and_ligand_groups(tmp_path: Path) -> None:
    out = tmp_path / "test.ndx"
    write_index([1, 2, 3], [4, 5], out)

    content = _read_index(out)
    assert "[ Receptor ]" in content
    assert "[ Ligand ]" in content
    assert "1 2 3" in content
    assert "4 5" in content


def test_write_index_line_wrapping(tmp_path: Path) -> None:
    """16 receptor atoms - first line has 15 atoms, second has 1."""
    out = tmp_path / "test.ndx"
    write_index(list(range(1, 17)), [17], out)

    content = _read_index(out)
    lines = [ln for ln in content.splitlines() if ln and not ln.startswith("[")]
    first_line_nums = lines[0].split()
    assert len(first_line_nums) == 15
    second_line_nums = lines[1].split()
    assert len(second_line_nums) == 1
    assert second_line_nums[0] == "16"


def test_write_index_raises_when_receptor_empty(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="Protein"):
        write_index([], [1, 2], tmp_path / "test.ndx")


def test_write_index_raises_when_ligand_empty(tmp_path: Path) -> None:
    with pytest.raises(RuntimeError, match="Ligand"):
        write_index([1, 2, 3], [], tmp_path / "test.ndx")
