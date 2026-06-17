"""Mol2 cap-stripping helpers for capped-dipeptide RESP parametrization."""

from __future__ import annotations

import contextlib
import logging
import warnings
from typing import TYPE_CHECKING

import gemmi
import networkx as nx

if TYPE_CHECKING:
    from pathlib import Path
import parmed as pmd

logger = logging.getLogger(__name__)

MIN_GREEK_ATOM_NAME_LENGTH = 2

# GAFF nitrogen and carbon atom type sets, used to detect ACE/NME cap atoms.
_GAFF_N_TYPES = {
    "n",
    "n1",
    "n2",
    "n3",
    "n4",
    "na",
    "nb",
    "nc",
    "nd",
    "nh",
    "no",
    "ns",
    "nt",
    "nu",
    "nv",
}
_GAFF_C_TYPES = {
    "c",
    "ca",
    "c1",
    "c2",
    "c3",
    "c4",
    "c5",
    "c6",
    "c7",
    "c8",
    "cc",
    "cd",
    "ce",
    "cf",
    "cg",
    "ch",
    "ci",
    "ck",
    "cm",
    "cn",
    "co",
    "cp",
    "cq",
    "cr",
    "cu",
    "cv",
    "cw",
    "cx",
    "cy",
    "cz",
}

# AMBER ff14SB backbone atom types for residue templates.
_AMBER_BACKBONE_TYPE = {
    "backbone_N": "N",
    "backbone_H": "H",
    "backbone_CA": "CX",
    "backbone_HA": "H1",
    "backbone_CB": "2C",
    "backbone_HB": "H1",
    "backbone_C": "C",
    "backbone_O": "O",
}


def _pdb_sidechain_names_by_depth(protein_pdb: Path, resname: str) -> dict[tuple[str, int], list[str]]:
    """Return a mapping (element, depth_from_CA) → [atom_names] for a residue in the PDB."""
    st = gemmi.read_pdb(str(protein_pdb))
    pdb_atoms: list[tuple[str, str]] = []
    for model in st:
        for chain in model:
            for residue in chain:
                if residue.name.upper() != resname:
                    continue
                for atom in residue:
                    pdb_atoms.append((atom.name, atom.element.name))

    backbone_names = {
        "N",
        "H",
        "CA",
        "HA",
        "C",
        "O",
        "HN",
        "1H",
        "2H",
        "3H",
        "H1",
        "H2",
        "H3",
        "OXT",
        "HB",
        "HB2",
        "HB3",
    }
    greek = {"A": 0, "B": 1, "G": 2, "D": 3, "E": 4, "Z": 5, "H": 6}
    depth: dict[str, int] = {"CA": 0}
    for aname, _ in pdb_atoms:
        if len(aname) >= MIN_GREEK_ATOM_NAME_LENGTH and aname[1:2].upper() in greek:
            depth[aname] = greek[aname[1:2].upper()]
        else:
            depth[aname] = 0

    result: dict[tuple[str, int], list[str]] = {}
    for aname, elem in pdb_atoms:
        if aname in backbone_names:
            continue
        key = (elem.upper(), depth.get(aname, 99))
        result.setdefault(key, []).append(aname)
    return result


def _gaff_type_to_element(gaff_type: str) -> str:
    """Map a GAFF atom type string to its element symbol."""
    t = gaff_type.lower()
    if t.startswith("cl"):
        return "Cl"
    if t.startswith("br"):
        return "Br"
    if t.startswith("s"):
        return "S"
    if t.startswith("o"):
        return "O"
    if t.startswith("n"):
        return "N"
    if t.startswith("p"):
        return "P"
    if t.startswith("h"):
        return "H"
    if t.startswith("c"):
        return "C"
    if t.startswith("f"):
        return "F"
    return t[0].upper()


def _strip_mol2_dipeptide_caps(
    mol2_path: Path,
    output_mol2: Path,
    protein_pdb: Path | None = None,
) -> Path:
    """Strip ACE/NME caps from a capped-dipeptide mol2 and rename backbone atoms.

    Raises ``ValueError`` when no ACE residue is found, so callers that process
    bare residue templates (e.g. MCPB.py CS1-4 files) can fall back gracefully.
    """
    structure = pmd.load_file(str(mol2_path), structure=True)

    res_names = {r.name.upper() for r in structure.residues}
    if "ACE" not in res_names:
        raise ValueError(f"No ACE residue found by ParmEd in {mol2_path}")

    cap_idx: set[int] = {a.idx for a in structure.atoms if a.residue.name.upper() in {"ACE", "NME"}}

    adj_pmd: dict[int, list[pmd.Atom]] = {a.idx: [] for a in structure.atoms}
    for bond in structure.bonds:
        adj_pmd[bond.atom1.idx].append(bond.atom2)
        adj_pmd[bond.atom2.idx].append(bond.atom1)

    backbone_n = None
    for atom in structure.atoms:
        if atom.idx in cap_idx or atom.type.lower() not in _GAFF_N_TYPES:
            continue
        if any(nb.idx in cap_idx for nb in adj_pmd[atom.idx]):
            backbone_n = atom
            break
    if backbone_n is None:
        raise ValueError(f"Could not identify backbone N in {mol2_path}")

    backbone_ca = next(
        (nb for nb in adj_pmd[backbone_n.idx] if nb.idx not in cap_idx and nb.type.lower() == "c3"),
        None,
    )
    if backbone_ca is None:
        raise ValueError(f"Could not identify backbone CA in {mol2_path}")

    backbone_c = backbone_o = None
    for nb in adj_pmd[backbone_ca.idx]:
        if nb is backbone_n or nb.idx in cap_idx:
            continue
        if nb.type.lower() in _GAFF_C_TYPES and nb.type.lower() != "c3":
            o_nbs = [x for x in adj_pmd[nb.idx] if x.type.lower() == "o" and x is not backbone_ca]
            if o_nbs:
                backbone_c = nb
                backbone_o = o_nbs[0]
                break
    if backbone_c is None or backbone_o is None:
        raise ValueError(f"Could not identify backbone C/O in {mol2_path}")

    backbone_ha = next(
        (nb for nb in adj_pmd[backbone_ca.idx] if nb.type.lower() == "h1" and nb.idx not in cap_idx),
        None,
    )
    backbone_h = next(
        (nb for nb in adj_pmd[backbone_n.idx] if nb.type.lower() in {"hn", "h"} and nb.idx not in cap_idx),
        None,
    )
    backbone_cb = next(
        (
            nb
            for nb in adj_pmd[backbone_ca.idx]
            if nb is not backbone_n
            and nb is not backbone_c
            and nb.idx not in cap_idx
            and nb.type.lower() == "c3"
            and nb is not backbone_ha
        ),
        None,
    )
    backbone_hb_atoms: list[pmd.Atom] = []
    if backbone_cb is not None:
        backbone_hb_atoms = [
            nb for nb in adj_pmd[backbone_cb.idx] if nb.type.lower() in {"h1", "hc", "hx"} and nb.idx not in cap_idx
        ]

    backbone_n.name = "N"
    backbone_n.type = _AMBER_BACKBONE_TYPE["backbone_N"]

    backbone_ca.name = "CA"
    backbone_ca.type = _AMBER_BACKBONE_TYPE["backbone_CA"]

    backbone_c.name = "C"
    backbone_c.type = _AMBER_BACKBONE_TYPE["backbone_C"]

    backbone_o.name = "O"
    backbone_o.type = _AMBER_BACKBONE_TYPE["backbone_O"]

    if backbone_ha:
        backbone_ha.name = "HA"
        backbone_ha.type = _AMBER_BACKBONE_TYPE["backbone_HA"]

    if backbone_h:
        backbone_h.name = "H"
        backbone_h.type = _AMBER_BACKBONE_TYPE["backbone_H"]

    if backbone_cb:
        backbone_cb.name = "CB"
        backbone_cb.type = _AMBER_BACKBONE_TYPE["backbone_CB"]

    for i, hb in enumerate(backbone_hb_atoms, start=2):
        hb.name = f"HB{i}"
        hb.type = _AMBER_BACKBONE_TYPE["backbone_HB"]

    resname = structure.residues[0].name.strip() if structure.residues else "UNK"
    pdb_names_by_depth: dict[tuple[str, int], list[str]] = {}
    if protein_pdb is not None:
        with contextlib.suppress(Exception):
            pdb_names_by_depth = _pdb_sidechain_names_by_depth(protein_pdb, resname)

    if pdb_names_by_depth:
        named_idx = {
            a.idx
            for a in [
                backbone_n,
                backbone_ca,
                backbone_c,
                backbone_o,
                backbone_ha,
                backbone_h,
                backbone_cb,
                *backbone_hb_atoms,
            ]
            if a is not None
        }
        graph_core = nx.Graph()
        graph_core.add_edges_from(
            (b.atom1.idx, b.atom2.idx)
            for b in structure.bonds
            if b.atom1.idx not in cap_idx and b.atom2.idx not in cap_idx
        )
        mol2_depth: dict[int, int] = dict(nx.single_source_shortest_path_length(graph_core, backbone_ca.idx))

        sc_by_elem_depth: dict[tuple[str, int], list[pmd.Atom]] = {}
        for atom in structure.atoms:
            if atom.idx in cap_idx or atom.idx in named_idx:
                continue
            elem = _gaff_type_to_element(atom.type)
            sc_by_elem_depth.setdefault((elem, mol2_depth.get(atom.idx, 99)), []).append(atom)

        pdb_names_used: set[str] = set()
        for (elem, depth), sc_atoms in sorted(sc_by_elem_depth.items()):
            available = [n for n in pdb_names_by_depth.get((elem, depth), []) if n not in pdb_names_used]
            for atom, pdb_name in zip(sc_atoms, available):
                atom.name = pdb_name
                pdb_names_used.add(pdb_name)

    structure.strip(":ACE,NME")  # noqa: B005
    output_mol2.parent.mkdir(parents=True, exist_ok=True)
    structure.save(str(output_mol2), format="mol2", overwrite=True)
    return output_mol2


def _strip_mol2_or_original(
    mol2: Path,
    work_dir: Path,
    protein_pdb: Path | None = None,
) -> Path:
    """Strip ACE/NME caps from mol2, falling back to the original path on any failure.

    Returns ``work_dir/<stem>_stripped.mol2`` on success, or the original ``mol2``
    path unchanged when stripping fails (e.g. the file is already a bare residue
    template produced by MCPB.py). A :class:`UserWarning` is emitted in the
    fallback case so callers are informed without raising.
    """
    stripped = work_dir / f"{mol2.stem}_stripped.mol2"
    try:
        _strip_mol2_dipeptide_caps(mol2, stripped, protein_pdb=protein_pdb)
    except Exception as exc:  # noqa: BLE001
        warnings.warn(f"Cap stripping skipped for {mol2.name}: {exc}", stacklevel=2)
        return mol2
    else:
        return stripped
