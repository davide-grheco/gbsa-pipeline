"""Diagnostic helpers for GROMACS MD stage failures.

Three entry points are exposed:

``check_posre_consistency``
    Reads a GROMACS GRO file and a ``posre_*.itp`` position-restraint file,
    maps every restrained atom index back to its residue/atom name, and warns
    when restrained atoms are not the expected backbone atoms (N, CA, C, O).

``analyze_crash_frames``
    Searches a stage work directory for GROMACS crash-dump PDB files
    (``step*b.pdb`` / ``step*c.pdb``), identifies exploded or NaN atoms,
    measures pre-crash protein-water contacts, and writes a plain-text report
    alongside the failed-stage log.

``find_extreme_atoms``
    Scans a GRO file and returns atoms whose absolute coordinate exceeds a
    threshold.  Used as a quick sanity check between stages.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING, NamedTuple

if TYPE_CHECKING:
    from pathlib import Path

import MDAnalysis as mda
import numpy as np

from gbsa_pipeline._constants import SOLVENT_RESIDUE_NAMES
from gbsa_pipeline._gro_io import _parse_gro
from gbsa_pipeline._spatial import contact_pairs

logger = logging.getLogger(__name__)

# Atom names considered "backbone" for position-restraint validation.
_BACKBONE_ATOM_NAMES: frozenset[str] = frozenset({"N", "CA", "C", "O"})


# ---------------------------------------------------------------------------
# Small data types
# ---------------------------------------------------------------------------


class PosreCheckResult(NamedTuple):
    """Result of a position-restraint consistency check."""

    ok: bool
    n_restrained: int
    unexpected: list[tuple[int, str, str]]  # (atom_idx, res_name, atom_name)
    first_twenty: list[tuple[int, str, str]]  # (atom_idx, res_name, atom_name)


# ---------------------------------------------------------------------------
# PDB parser (GROMACS crash dumps)
# ---------------------------------------------------------------------------


def _parse_crash_pdb(pdb_path: Path) -> list[tuple[int, str, str, float, float, float]]:
    """Parse a GROMACS step*.pdb crash dump.

    Returns list of (atom_index, res_name, atom_name, x, y, z) in Angstrom.
    NaN/Inf coordinates are preserved so callers can detect them.
    """
    records: list[tuple[int, str, str, float, float, float]] = []
    with pdb_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                atom_idx = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:21].strip()
                x_str = line[30:38].strip()
                y_str = line[38:46].strip()
                z_str = line[46:54].strip()
                x = float(x_str) if x_str not in ("", "nan", "NaN") else math.nan
                y = float(y_str) if y_str not in ("", "nan", "NaN") else math.nan
                z = float(z_str) if z_str not in ("", "nan", "NaN") else math.nan
                records.append((atom_idx, res_name, atom_name, x, y, z))
            except (ValueError, IndexError):
                continue
    return records


# ---------------------------------------------------------------------------
# Position-restraint validator
# ---------------------------------------------------------------------------


def _parse_posre_indices(posre_path: Path) -> list[int]:
    """Return the 1-based restrained atom indices from a posre ITP file.

    Restraint entries are lines whose first field is a positive integer
    followed by at least one more field (function type, force constants).
    """
    indices: list[int] = []
    with posre_path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            fields = line.split()
            if len(fields) >= 2 and fields[0].isdigit() and int(fields[0]) > 0:  # noqa: PLR2004
                indices.append(int(fields[0]))
    return indices


def check_posre_consistency(
    gro_path: Path,
    posre_path: Path,
    expected_atom_names: frozenset[str] = _BACKBONE_ATOM_NAMES,
) -> PosreCheckResult:
    """Validate that posre ITP indices point to expected backbone atoms.

    Reads ``posre_path`` for the list of restrained 1-based atom indices and
    maps each index to the corresponding atom in ``gro_path``.  Returns a
    ``PosreCheckResult`` with:

    - ``ok`` — True when every restrained atom is in ``expected_atom_names``
      and no solvent/ligand residue is restrained.
    - ``n_restrained`` — total number of restrained atoms.
    - ``unexpected`` — list of (index, res_name, atom_name) for atoms that
      do not match the expected set.
    - ``first_twenty`` — first 20 restrained atoms for inspection.
    """
    if not gro_path.exists():
        logger.warning("check_posre_consistency: GRO not found: %s", gro_path)
        return PosreCheckResult(ok=False, n_restrained=0, unexpected=[], first_twenty=[])
    if not posre_path.exists():
        logger.warning("check_posre_consistency: posre ITP not found: %s", posre_path)
        return PosreCheckResult(ok=False, n_restrained=0, unexpected=[], first_twenty=[])

    # Map 1-based GRO atom index → (res_name, atom_name); unknown indices
    # resolve to ("?", "?"), which never matches the expected atom names.
    atoms = mda.Universe(str(gro_path)).atoms
    by_index = {
        atom_id: (res_name, atom_name)
        for atom_id, res_name, atom_name in zip(
            atoms.ids.tolist(), atoms.resnames.tolist(), atoms.names.tolist(), strict=True
        )
    }

    entries: list[tuple[int, str, str]] = []
    for idx in _parse_posre_indices(posre_path):
        res_name, atom_name = by_index.get(idx, ("?", "?"))
        entries.append((idx, res_name, atom_name))

    unexpected = [
        (idx, res_name, atom_name)
        for idx, res_name, atom_name in entries
        if res_name in SOLVENT_RESIDUE_NAMES or atom_name not in expected_atom_names
    ]

    ok = len(unexpected) == 0
    result = PosreCheckResult(
        ok=ok,
        n_restrained=len(entries),
        unexpected=unexpected,
        first_twenty=entries[:20],
    )

    if ok:
        logger.debug(
            "posre check PASSED: %d restrained atoms all in %s",
            result.n_restrained,
            sorted(expected_atom_names),
        )
    else:
        logger.warning(
            "posre check FAILED: %d of %d restrained atoms unexpected.\n  First 5 unexpected: %s",
            len(unexpected),
            result.n_restrained,
            result.unexpected[:5],
        )

    return result


# ---------------------------------------------------------------------------
# Extreme-atom scanner
# ---------------------------------------------------------------------------


def find_extreme_atoms(
    gro_path: Path,
    threshold_nm: float = 10.0,
) -> list[tuple[int, str, str, float, float, float]]:
    """Return atoms whose absolute coordinate exceeds ``threshold_nm``.

    Useful as a between-stage sanity check.  Returns list of
    (atom_index, res_name, atom_name, x, y, z) in nm.
    """
    if not gro_path.exists():
        return []
    atoms = _parse_gro(gro_path)
    return [
        (a.atom_idx, a.res_name, a.atom_name, a.x, a.y, a.z)
        for a in atoms
        if abs(a.x) > threshold_nm or abs(a.y) > threshold_nm or abs(a.z) > threshold_nm
    ]


# ---------------------------------------------------------------------------
# Crash-frame analyzer
# ---------------------------------------------------------------------------


def analyze_crash_frames(work_dir: Path) -> str:
    """Analyze GROMACS step*.pdb crash dumps and write crash_report.txt.

    Locates ``step*b.pdb`` (pre-crash) and ``step*c.pdb`` (post-constraint)
    files in ``work_dir``.  When found it:

    - identifies atoms with NaN or extreme (|coord| > 100 Å) coordinates
    - computes displacement between the b→c frames for each common atom
    - lists the 10 atoms with the largest displacement
    - lists contacts below 1.0 Å between non-solvent and solvent atoms in
      the pre-crash frame
    - writes a readable report to ``work_dir/crash_report.txt``

    Returns the report string (also written to disk).
    """
    step_pdbs_b = sorted(work_dir.glob("step*b.pdb"))
    step_pdbs_c = sorted(work_dir.glob("step*c.pdb"))

    if not step_pdbs_b and not step_pdbs_c:
        return ""

    lines: list[str] = [
        f"GROMACS crash-frame analysis: {work_dir}",
        "=" * 72,
    ]

    for b_path in step_pdbs_b:
        c_stem = b_path.stem[:-1] + "c"
        c_path = work_dir / (c_stem + ".pdb")

        lines.append(f"\nPre-crash frame:  {b_path.name}")
        b_records = _parse_crash_pdb(b_path)

        # --- NaN / extreme atoms in b frame --------------------------------
        nan_atoms = [r for r in b_records if math.isnan(r[3]) or math.isnan(r[4]) or math.isnan(r[5])]
        extreme_atoms = [
            r
            for r in b_records
            if not math.isnan(r[3]) and (abs(r[3]) > 100.0 or abs(r[4]) > 100.0 or abs(r[5]) > 100.0)  # noqa: PLR2004
        ]

        if nan_atoms:
            lines.append(f"  NaN atoms in pre-crash frame ({len(nan_atoms)}):")
            for idx, rname, aname, x, y, z in nan_atoms[:10]:
                lines.append(f"    atom {idx:6d}  {rname:5s} {aname:6s}  ({x}, {y}, {z})")

        if extreme_atoms:
            lines.append(f"  Extreme coordinates |coord| > 100 Å ({len(extreme_atoms)} atoms):")
            for idx, rname, aname, x, y, z in extreme_atoms[:10]:
                lines.append(f"    atom {idx:6d}  {rname:5s} {aname:6s}  ({x:.1f}, {y:.1f}, {z:.1f}) Å")

        if not nan_atoms and not extreme_atoms:
            lines.append("  No NaN or extreme atoms in pre-crash frame.")

        # --- Contacts below 1.0 Å between protein and water ---------------
        close_contacts = _find_close_contacts_pdb(b_records, threshold_ang=1.0)
        if close_contacts:
            lines.append(f"  Close contacts < 1.0 Å ({len(close_contacts)}):")
            for r1, r2, d in close_contacts[:15]:
                i1, rn1, an1, _x1, _y1, _z1 = r1
                i2, rn2, an2, _x2, _y2, _z2 = r2
                lines.append(f"    atom {i1:6d} {rn1:5s} {an1:6s} -- atom {i2:6d} {rn2:5s} {an2:6s}  d = {d:.3f} A")
        else:
            lines.append("  No close contacts < 1.0 Å between protein and solvent.")

        # --- Displacement b→c -------------------------------------------
        if c_path.exists():
            lines.append(f"\nPost-constraint frame: {c_path.name}")
            c_records = _parse_crash_pdb(c_path)

            nan_c = [r for r in c_records if math.isnan(r[3]) or math.isnan(r[4]) or math.isnan(r[5])]
            if nan_c:
                lines.append(f"  NaN atoms in post-constraint frame ({len(nan_c)}):")
                for idx, rname, aname, _x, _y, _z in nan_c[:15]:
                    lines.append(f"    atom {idx:6d}  {rname:5s} {aname:6s}")

            displacements = _compute_displacements(b_records, c_records)
            if displacements:
                displacements.sort(key=lambda t: -t[2])
                lines.append("  Top 10 largest b→c displacements:")
                for idx, (rname, aname), disp in displacements[:10]:
                    lines.append(f"    atom {idx:6d}  {rname:5s} {aname:6s}  Δ = {disp:.2f} Å")

    report = "\n".join(lines)
    report_path = work_dir / "crash_report.txt"
    try:
        report_path.write_text(report + "\n", encoding="utf-8")
        logger.info("Crash report written to %s", report_path)
    except OSError as exc:
        logger.warning("Could not write crash report: %s", exc)

    return report


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _find_close_contacts_pdb(
    records: list[tuple[int, str, str, float, float, float]],
    threshold_ang: float = 1.0,
) -> list[tuple[tuple, tuple, float]]:
    """Return pairs with distance below threshold (A), protein-solvent only."""
    solvent = [r for r in records if r[1] in SOLVENT_RESIDUE_NAMES and not any(math.isnan(v) for v in r[3:6])]
    protein = [r for r in records if r[1] not in SOLVENT_RESIDUE_NAMES and not any(math.isnan(v) for v in r[3:6])]

    if not protein or not solvent:
        return []

    protein_coords = np.array([[r[3], r[4], r[5]] for r in protein])
    solvent_coords = np.array([[r[3], r[4], r[5]] for r in solvent])
    return [
        (protein[i], solvent[j], dist) for i, j, dist in contact_pairs(protein_coords, solvent_coords, threshold_ang)
    ]


def _compute_displacements(
    b_records: list[tuple[int, str, str, float, float, float]],
    c_records: list[tuple[int, str, str, float, float, float]],
) -> list[tuple[int, tuple[str, str], float]]:
    """Compute per-atom displacement (Å) between two PDB frames."""
    c_by_idx: dict[int, tuple[float, float, float]] = {r[0]: (r[3], r[4], r[5]) for r in c_records}
    result = []
    for idx, rname, aname, bx, by, bz in b_records:
        if math.isnan(bx) or math.isnan(by) or math.isnan(bz):
            continue
        if idx not in c_by_idx:
            continue
        cx, cy, cz = c_by_idx[idx]
        if math.isnan(cx) or math.isnan(cy) or math.isnan(cz):
            disp = math.inf
        else:
            disp = float(np.linalg.norm([bx - cx, by - cy, bz - cz]))
        result.append((idx, (rname, aname), disp))
    return result
