"""GROMACS GRO/TOP I/O helpers: parsing, clash removal, topology patching."""

from __future__ import annotations

import logging
from collections import Counter
from typing import TYPE_CHECKING, NamedTuple

import MDAnalysis as mda

if TYPE_CHECKING:
    from pathlib import Path

LOGGER = logging.getLogger(__name__)

_GRO_MIN_LINE_LEN = 44  # coordinates end at column 44


class _GROAtom(NamedTuple):
    atom_idx: int  # 1-based atom index
    res_num: int
    res_name: str
    atom_name: str
    x: float  # nm
    y: float  # nm
    z: float  # nm


# ---------------------------------------------------------------------------
# Line-level parsers
# ---------------------------------------------------------------------------


def _parse_gro_atom_line(line: str) -> _GROAtom:
    """Parse one GRO atom line into a :class:`_GROAtom`.

    Raises ``ValueError`` on malformed input so callers that process
    intermediate files fail explicitly rather than silently producing
    wrong geometry.
    """
    try:
        return _GROAtom(
            atom_idx=int(line[15:20]),
            res_num=int(line[0:5]),
            res_name=line[5:10].strip(),
            atom_name=line[10:15].strip(),
            x=float(line[20:28]),
            y=float(line[28:36]),
            z=float(line[36:44]),
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Could not parse GRO atom line: {line!r}") from exc


# ---------------------------------------------------------------------------
# File-level parsers
# ---------------------------------------------------------------------------


def _parse_gro(gro_path: Path) -> list[_GROAtom]:
    """Read a GROMACS GRO file and return one :class:`_GROAtom` per atom.

    Short lines (< 44 characters) are skipped silently — they indicate
    truncated or velocity-only trailing records that do not carry coordinates.
    """
    with gro_path.open(encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()
    n_atoms = int(lines[1])
    return [_parse_gro_atom_line(line) for line in lines[2 : 2 + n_atoms] if len(line) >= _GRO_MIN_LINE_LEN]


# ---------------------------------------------------------------------------
# Clash detection and solvent cleanup
# ---------------------------------------------------------------------------


def _write_cleaned_gro(
    input_gro: Path,
    output_gro: Path,
    cutoff_nm: float,
    water_resnames: set[str],
) -> dict[str, int]:
    """Write a GRO file with clashing solvent waters removed.

    Removes whole water residues with any atom within ``cutoff_nm`` of a
    non-water atom (minimum-image aware). The MDAnalysis writer renumbers atom
    serials and updates the atom count. Returns a ``{resname: count}`` dict of
    removed molecules for topology patching.
    """
    universe = mda.Universe(str(input_gro))
    water = "resname " + " ".join(sorted(water_resnames))
    clashing = universe.select_atoms(f"byres (({water}) and around {cutoff_nm * 10.0} (not ({water})))")

    (universe.atoms - clashing).write(str(output_gro))

    return dict(Counter(residue.resname for residue in clashing.residues))


# ---------------------------------------------------------------------------
# Topology patching
# ---------------------------------------------------------------------------


def _update_topology_water_counts(
    input_top: Path,
    output_top: Path,
    removed_counts: dict[str, int],
) -> None:
    """Write a topology with [ molecules ] water counts reduced to match a cleaned GRO.

    Raises ``ValueError`` if waters were removed but no matching entry can be
    found in the ``[ molecules ]`` section.
    """
    if not removed_counts:
        output_top.write_text(input_top.read_text(encoding="utf-8", errors="replace"), encoding="utf-8")
        return

    lines = input_top.read_text(encoding="utf-8", errors="replace").splitlines()
    in_molecules = False
    remaining = dict(removed_counts)
    output_lines: list[str] = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            in_molecules = stripped.strip("[]").strip().lower() == "molecules"
            output_lines.append(line)
            continue

        if not in_molecules or not stripped or stripped.startswith((";", "#")):
            output_lines.append(line)
            continue

        body, *comment_parts = line.split(";", 1)
        comment = ";" + comment_parts[0] if comment_parts else ""
        fields = body.split()
        if len(fields) < 2 or fields[0] not in remaining:  # noqa: PLR2004
            output_lines.append(line)
            continue

        molecule_name = fields[0]
        try:
            old_count = int(fields[1])
        except ValueError as exc:
            raise ValueError(f"Could not parse molecule count in topology line: {line!r}") from exc

        new_count = old_count - remaining[molecule_name]
        if new_count < 0:
            raise ValueError(
                f"Cannot remove {remaining[molecule_name]} {molecule_name} molecules: topology count is {old_count}."
            )

        prefix = line[: line.find(molecule_name)] if molecule_name in line else ""
        output_lines.append(f"{prefix}{molecule_name:<16} {new_count}{(' ' + comment) if comment else ''}".rstrip())
        del remaining[molecule_name]

    if remaining:
        raise ValueError(
            "Removed solvent waters but could not update matching topology molecule counts: "
            + ", ".join(f"{name}={count}" for name, count in sorted(remaining.items()))
        )

    output_top.write_text("\n".join(output_lines) + "\n", encoding="utf-8")
