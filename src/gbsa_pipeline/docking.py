# /home/grheco/repositorios/gbsa-pipeline/src/gbsa_pipeline/docking.py

"""Docking helpers for preparing ligands and running AutoDock Vina.

This module is intended as a reusable docking utility layer for the pipeline.

Main responsibilities
---------------------
- normalize ligand input from common formats,
- prepare ligands for Vina through RDKit + Meeko,
- convert receptors from PDB to PDBQT with Open Babel,
- export docked PDBQT poses back to SDF through mk_export.py,
- optionally rebuild bond orders with RDKit AssignBondOrdersFromTemplate,
- generate docking boxes directly from ligand coordinates,
- run native-ligand redocking and derive a final query box from the
  redocked native pose.

Design principles
-----------------
1. Keep the workflow Vina-only.
2. Prefer explicit intermediate files.
3. Fail early with clear error messages.
4. Write full subprocess logs to disk while keeping terminal output compact.

Typical workflow
----------------
1. Read a native ligand from SDF / PDB / MOL / MOL2 / SMILES / RDKit Mol.
2. Normalize it to SDF if needed.
3. Add hydrogens, ensure 3D coordinates, and geometry-optimize the ligand.
4. Prepare the ligand to PDBQT using Meeko.
5. Redock the native ligand with Vina.
6. Export the best docked native pose back to SDF.
7. Build the final query docking box from that redocked native pose.
8. Prepare a query ligand and dock it into the native-derived box.

Notes on force fields
---------------------
Geometry optimization prefers MMFF94 when possible and falls back to UFF.
This is intended for typical small organic ligands. The choice is only for
pre-docking geometry cleanup; it is not part of the actual Vina scoring step.
"""

from __future__ import annotations

import logging
import re
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from subprocess import CompletedProcess, run
from typing import TYPE_CHECKING, Any, Protocol

from meeko import MoleculePreparation, PDBQTWriterLegacy
from pydantic import BaseModel, ConfigDict, Field, field_validator
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import (
    MMFFGetMoleculeProperties,
    MMFFOptimizeMolecule,
    UFFOptimizeMolecule,
)

if TYPE_CHECKING:
    from collections.abc import Mapping

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

if not LOGGER.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(message)s"))
    LOGGER.addHandler(handler)

LOGGER.propagate = False


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


class DockingBox(BaseModel):
    """Docking box definition in Angstrom."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    center: tuple[float, float, float]
    size: tuple[float, float, float]


class DockingRequest(BaseModel):
    """Normalized request object for a Vina docking run."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    receptor: Path
    ligands: list[Path]
    box: DockingBox
    seed: int | None = None
    workdir: Path | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)

    @field_validator("receptor")
    @classmethod
    def _check_receptor_exists(cls, path: Path) -> Path:
        """Validate receptor existence and supported suffix."""
        path = Path(path)

        if not path.exists():
            raise ValueError(f"Receptor file does not exist: {path}")

        if not path.is_file():
            raise ValueError(f"Receptor path is not a file: {path}")

        if path.suffix.lower() not in {".pdb", ".pdbqt"}:
            raise ValueError(f"Receptor must be .pdb or .pdbqt for vina docking: {path}")

        return path.resolve()

    @field_validator("ligands")
    @classmethod
    def _check_ligands(cls, ligands: list[Path]) -> list[Path]:
        """Validate that at least one ligand file exists."""
        if not ligands:
            raise ValueError("At least one ligand file required.")

        checked: list[Path] = []
        for ligand_entry in ligands:
            ligand_path = Path(ligand_entry)

            if not ligand_path.exists():
                raise ValueError(f"Ligand file missing: {ligand_path}")

            if not ligand_path.is_file():
                raise ValueError(f"Ligand path not file: {ligand_path}")

            checked.append(ligand_path.resolve())

        return checked


@dataclass(frozen=True)
class DockedPose:
    """Single docked pose plus compact metadata."""

    ligand: Path
    pose_path: Path
    score: float | None
    rank: int | None
    engine: str
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class DockingResult:
    """Collection of docking outputs for one request."""

    poses: list[DockedPose]
    engine: str
    parameters: Mapping[str, Any]
    raw_outputs: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class NativeRedockBoxResult:
    """Result bundle for native redocking + query-box derivation."""

    resolved_input_sdf: Path
    before_meeko_sdf: Path
    native_pdbqt: Path
    initial_redock_box: DockingBox
    docking_result: DockingResult
    best_pose_pdbqt: Path
    best_pose_export_sdf: Path
    query_box: DockingBox
    best_score: float | None


class DockingEngine(Protocol):
    """Protocol implemented by docking backends."""

    name: str

    def dock(self, request: DockingRequest) -> DockingResult:
        """Run docking for the provided request."""


# ---------------------------------------------------------------------------
# Small logging / text helpers
# ---------------------------------------------------------------------------


def _box(title: str, char: str = "=") -> str:
    """Return a compact ASCII title box."""
    line = char * max(len(title), 12)
    return f"{line}\n{title}\n{line}"


def _section(title: str) -> None:
    """Log a section header."""
    LOGGER.info("")
    LOGGER.info(_box(title))


def _step(prefix: str, message: str) -> None:
    """Log one normal progress line."""
    LOGGER.info("[%s] %s", prefix, message)


def _warn(prefix: str, message: str) -> None:
    """Log one warning line."""
    LOGGER.warning("[%s] warning: %s", prefix, message)


def _error(prefix: str, message: str) -> None:
    """Log one error line."""
    LOGGER.error("[%s] error: %s", prefix, message)


def _summarize_stderr(stderr: str, max_lines: int = 4) -> str:
    """Return a short one-line summary of stderr."""
    lines = [line.strip() for line in stderr.splitlines() if line.strip()]

    if not lines:
        return "no stderr output"

    preview = lines[:max_lines]
    text = " | ".join(preview)

    if len(lines) > max_lines:
        text += " | ..."

    return text


def _write_process_log(
    log_path: Path,
    process: CompletedProcess[str],
    *,
    command: list[str],
    title: str,
) -> None:
    """Write complete subprocess stdout/stderr to a plain-text log file."""
    log_path = Path(log_path).resolve()
    log_path.parent.mkdir(parents=True, exist_ok=True)

    text = (
        f"{title}\n"
        f"{'=' * len(title)}\n\n"
        f"Command:\n{' '.join(command)}\n\n"
        f"Return code:\n{process.returncode}\n\n"
        f"STDOUT\n"
        f"------\n"
        f"{process.stdout or ''}\n\n"
        f"STDERR\n"
        f"------\n"
        f"{process.stderr or ''}\n"
    )
    log_path.write_text(text, encoding="utf-8")


# ---------------------------------------------------------------------------
# Meeko / Vina parsing helpers
# ---------------------------------------------------------------------------


def _extract_pdbqt_string_from_meeko_result(result: Any) -> str:
    """Extract the PDBQT string from Meeko's version-dependent return value."""
    if isinstance(result, str):
        return result

    if isinstance(result, tuple):
        if not result:
            raise ValueError("Meeko returned an empty tuple from write_string().")

        pdbqt_string = result[0]

        if not isinstance(pdbqt_string, str):
            raise TypeError(
                f"Expected first element of Meeko write_string() result to be a str, got {type(pdbqt_string).__name__}."
            )

        return pdbqt_string

    raise TypeError(f"Unexpected return type from Meeko write_string(): {type(result).__name__}")


def _parse_vina_best_score(stdout: str) -> float | None:
    """Parse the best affinity from the first row of Vina's output table."""
    for line in stdout.splitlines():
        stripped = line.strip()
        match = re.match(r"^1\s+(-?\d+(?:\.\d+)?)\s+", stripped)
        if match:
            return float(match.group(1))
    return None


def _build_compact_pose_metadata(
    *,
    returncode: int,
    log_file: Path,
    output_exists: bool,
    receptor_used: Path,
) -> dict[str, Any]:
    """Build compact metadata for a single docked pose."""
    return {
        "returncode": returncode,
        "log_file": str(log_file),
        "output_exists": output_exists,
        "receptor_used": str(receptor_used),
    }


# ---------------------------------------------------------------------------
# Ligand input normalization
# ---------------------------------------------------------------------------


def _looks_like_path_string(text: str) -> bool:
    """Heuristically distinguish path-like strings from SMILES strings."""
    lower = text.lower()
    known_suffixes = (".sdf", ".pdb", ".mol", ".mol2", ".pdbqt")
    return lower.endswith(known_suffixes) or "/" in text or "\\" in text


def read_sdf_molecules(path: Path) -> list[Chem.Mol]:
    """Read all valid molecules from an SDF file."""
    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"SDF file not found: {path}")
    if not path.is_file():
        raise ValueError(f"SDF path is not a file: {path}")

    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    mols = [mol for mol in supplier if mol is not None]

    if not mols:
        raise ValueError(f"Could not read any valid molecules from SDF: {path}")

    return mols


def read_first_sdf_molecule(path: Path) -> Chem.Mol:
    """Read the first valid molecule from an SDF file."""
    return read_sdf_molecules(path)[0]


def read_pdb_ligand_molecule(
    path: Path,
    *,
    sanitize: bool = True,
    remove_hs: bool = False,
) -> Chem.Mol:
    """Read a ligand-like molecule from a PDB file using RDKit."""
    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"PDB ligand file not found: {path}")
    if not path.is_file():
        raise ValueError(f"PDB ligand path is not a file: {path}")

    mol = Chem.MolFromPDBFile(
        str(path),
        sanitize=sanitize,
        removeHs=remove_hs,
    )
    if mol is None:
        raise ValueError(f"Could not read ligand molecule from PDB: {path}")

    return mol


def read_mol_ligand_molecule(path: Path) -> Chem.Mol:
    """Read a MOL ligand file using RDKit."""
    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"MOL ligand file not found: {path}")
    if not path.is_file():
        raise ValueError(f"MOL ligand path is not a file: {path}")

    mol = Chem.MolFromMolFile(str(path), sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read ligand molecule from MOL: {path}")

    return mol


def read_mol2_ligand_molecule(path: Path) -> Chem.Mol:
    """Read a MOL2 ligand file using RDKit."""
    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"MOL2 ligand file not found: {path}")
    if not path.is_file():
        raise ValueError(f"MOL2 ligand path is not a file: {path}")

    mol = Chem.MolFromMol2File(str(path), sanitize=True, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read ligand molecule from MOL2: {path}")

    return mol


def write_sdf_molecules(molecules: list[Chem.Mol], path: Path) -> Path:
    """Write one or more RDKit molecules to an SDF file."""
    path = Path(path).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    writer = Chem.SDWriter(str(path))
    try:
        for mol in molecules:
            writer.write(mol)
    finally:
        writer.close()

    return path


def write_single_sdf(mol: Chem.Mol, path: Path) -> Path:
    """Write a single RDKit molecule to SDF."""
    return write_sdf_molecules([mol], path)


def _load_ligand_file_as_mol(path: Path) -> Chem.Mol:
    """Load a ligand file into RDKit based on its suffix."""
    suffix = path.suffix.lower()

    if suffix == ".sdf":
        return read_first_sdf_molecule(path)
    if suffix == ".pdb":
        return read_pdb_ligand_molecule(path)
    if suffix == ".mol":
        return read_mol_ligand_molecule(path)
    if suffix == ".mol2":
        return read_mol2_ligand_molecule(path)

    raise ValueError(f"Unsupported ligand file type. Expected one of: .sdf, .pdb, .mol, .mol2. Got: {path}")


def _load_ligand_input_as_mol(
    ligand: str | Path | Chem.Mol,
    *,
    name: str | None = None,
) -> tuple[Chem.Mol, Path | None]:
    """Load ligand input into an RDKit molecule."""
    if isinstance(ligand, Chem.Mol):
        mol = Chem.Mol(ligand)
        source_path = None

    elif isinstance(ligand, Path):
        source_path = ligand.expanduser().resolve()
        mol = _load_ligand_file_as_mol(source_path)

    elif isinstance(ligand, str):
        if _looks_like_path_string(ligand):
            source_path = Path(ligand).expanduser().resolve()
            if not source_path.exists():
                raise FileNotFoundError(f"Ligand file not found: {source_path}")
            mol = _load_ligand_file_as_mol(source_path)
        else:
            mol = Chem.MolFromSmiles(ligand)
            if mol is None:
                raise ValueError("Failed to parse SMILES.")
            source_path = None
    else:
        raise TypeError(f"Unsupported ligand input type: {type(ligand).__name__}")

    if name is not None:
        mol.SetProp("_Name", name)
    elif source_path is not None and not mol.HasProp("_Name"):
        mol.SetProp("_Name", source_path.stem)
    elif not mol.HasProp("_Name"):
        mol.SetProp("_Name", "LIG")

    return mol, source_path


def ensure_ligand_sdf_for_meeko(
    ligand: str | Path | Chem.Mol,
    output_sdf: Path,
    *,
    name: str | None = None,
) -> tuple[Chem.Mol, Path]:
    """Ensure that a ligand exists as an SDF on disk for Meeko workflows."""
    mol, source_path = _load_ligand_input_as_mol(ligand, name=name)

    if source_path is not None and source_path.suffix.lower() == ".sdf":
        return mol, source_path

    output_sdf = Path(output_sdf).expanduser().resolve()
    write_single_sdf(mol, output_sdf)
    return mol, output_sdf


# ---------------------------------------------------------------------------
# Ligand cleanup and geometry preparation
# ---------------------------------------------------------------------------


def add_hydrogens_and_3d_if_needed(
    mol: Chem.Mol,
    *,
    random_seed: int = 20260410,
) -> Chem.Mol:
    """Return a copy of the molecule with explicit hydrogens and 3D coordinates."""
    if mol is None:
        raise ValueError("Ligand molecule is None.")

    work = Chem.Mol(mol)
    work = Chem.AddHs(work, addCoords=True)

    needs_3d = work.GetNumConformers() == 0
    if not needs_3d:
        try:
            conf = work.GetConformer()
            needs_3d = not conf.Is3D()
        except ValueError:
            needs_3d = True

    if needs_3d:
        embed_status = EmbedMolecule(work, randomSeed=random_seed)
        if embed_status != 0:
            raise RuntimeError("RDKit failed to embed 3D coordinates for ligand.")

    return work


def ligand_geoopt(
    mol: Chem.Mol,
    *,
    random_seed: int = 20260410,
    mmff_variant: str = "MMFF94",
    max_iters: int = 2000,
) -> Chem.Mol:
    """Return a geometry-cleaned ligand for docking."""
    if mol is None:
        raise ValueError("Ligand molecule is None.")

    work = add_hydrogens_and_3d_if_needed(mol, random_seed=random_seed)

    mmff_props = MMFFGetMoleculeProperties(work, mmffVariant=mmff_variant)
    if mmff_props is not None:
        opt_status = MMFFOptimizeMolecule(
            work,
            mmffVariant=mmff_variant,
            maxIters=max_iters,
        )
        if opt_status not in (0, 1):
            _warn("rdkit", f"MMFF optimization returned non-standard status {opt_status}")
        _step("rdkit", f"ligand geometry optimized with {mmff_variant}")
        return work

    _warn("rdkit", "MMFF unavailable; falling back to UFF")

    uff_status = UFFOptimizeMolecule(work, maxIters=max_iters)
    if uff_status not in (0, 1):
        _warn("rdkit", f"UFF optimization returned non-standard status {uff_status}")

    _step("rdkit", "ligand geometry optimized with UFF")
    return work


# ---------------------------------------------------------------------------
# Box generation
# ---------------------------------------------------------------------------


def build_box_from_mol(
    mol: Chem.Mol,
    *,
    padding: float = 4.0,
    min_size: float = 0.0,
    use_heavy_atoms_only: bool = True,
) -> DockingBox:
    """Build a Vina search box from ligand coordinates."""
    if mol is None:
        raise ValueError("Ligand molecule is None.")

    box_mol = Chem.Mol(mol)
    if use_heavy_atoms_only:
        stripped = Chem.RemoveHs(Chem.Mol(mol))
        if stripped.GetNumAtoms() > 0:
            box_mol = stripped

    if box_mol.GetNumConformers() == 0:
        raise ValueError("Ligand molecule has no 3D conformer for box generation.")

    conf = box_mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []

    for atom_index in range(box_mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(atom_index)
        xs.append(float(pos.x))
        ys.append(float(pos.y))
        zs.append(float(pos.z))

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    center = (
        (min_x + max_x) / 2.0,
        (min_y + max_y) / 2.0,
        (min_z + max_z) / 2.0,
    )
    size = (
        max(max_x - min_x + 2.0 * padding, min_size),
        max(max_y - min_y + 2.0 * padding, min_size),
        max(max_z - min_z + 2.0 * padding, min_size),
    )

    return DockingBox(center=center, size=size)


def build_box_from_sdf(
    sdf_path: Path,
    *,
    padding: float = 4.0,
    min_size: float = 0.0,
    use_heavy_atoms_only: bool = True,
) -> DockingBox:
    """Build a Vina box from the first molecule in an SDF file."""
    mol = read_first_sdf_molecule(sdf_path)
    return build_box_from_mol(
        mol,
        padding=padding,
        min_size=min_size,
        use_heavy_atoms_only=use_heavy_atoms_only,
    )


# ---------------------------------------------------------------------------
# PDBQT -> SDF export and chemistry reconstruction
# ---------------------------------------------------------------------------


def _detect_mk_export_output_flag(mk_export_binary: str) -> str:
    """Detect which output flag the installed mk_export.py supports."""
    help_process: CompletedProcess[str] = run(  # noqa: S603
        [mk_export_binary, "-h"],
        capture_output=True,
        text=True,
        check=False,
    )
    help_text = (help_process.stdout or "") + "\n" + (help_process.stderr or "")

    if "--write_sdf" in help_text:
        return "-s"
    if re.search(r"(^|\s)-o(\s|,|$)", help_text):
        return "-o"

    raise RuntimeError(
        "Could not determine mk_export.py output flag from help text.\n"
        f"STDOUT:\n{help_process.stdout}\n"
        f"STDERR:\n{help_process.stderr}"
    )


def _copy_heavy_atom_coordinates_if_possible(
    target_mol: Chem.Mol,
    source_mol: Chem.Mol,
) -> Chem.Mol:
    """Copy heavy-atom coordinates from source to target when possible."""
    if source_mol.GetNumConformers() == 0:
        return target_mol
    if target_mol.GetNumConformers() > 0:
        return target_mol

    source_no_h = Chem.RemoveHs(Chem.Mol(source_mol))
    if target_mol.GetNumAtoms() != source_no_h.GetNumAtoms():
        return target_mol

    source_conf = source_no_h.GetConformer()
    new_conf = Chem.Conformer(target_mol.GetNumAtoms())

    for atom_idx in range(target_mol.GetNumAtoms()):
        pos = source_conf.GetAtomPosition(atom_idx)
        new_conf.SetAtomPosition(atom_idx, pos)

    target_mol.AddConformer(new_conf, assignId=True)
    return target_mol


def assign_bond_orders_from_template_mol(
    template_mol: Chem.Mol,
    target_mol: Chem.Mol,
    *,
    add_hydrogens: bool = False,
) -> Chem.Mol:
    """Rebuild target bond orders from a chemistry template."""
    if template_mol is None:
        raise ValueError("template_mol is None.")
    if target_mol is None:
        raise ValueError("target_mol is None.")

    template_no_h = Chem.RemoveHs(Chem.Mol(template_mol))
    target_no_h = Chem.RemoveHs(Chem.Mol(target_mol))

    try:
        rebuilt = AllChem.AssignBondOrdersFromTemplate(template_no_h, target_no_h)
    except Exception as exc:
        raise RuntimeError(
            "RDKit AssignBondOrdersFromTemplate failed. Template and exported molecule likely do not match."
        ) from exc

    rebuilt = _copy_heavy_atom_coordinates_if_possible(rebuilt, target_mol)

    if add_hydrogens:
        rebuilt = Chem.AddHs(rebuilt, addCoords=True)

    return rebuilt


def export_pdbqt_to_sdf(
    pdbqt_path: Path,
    output_sdf: Path,
    *,
    mk_export_binary: str = "mk_export.py",
    template_mol: Chem.Mol | None = None,
    template_bond_orders: bool = False,
    add_hydrogens_after_template: bool = True,
) -> Path:
    """Export a PDBQT back to SDF using Meeko mk_export.py."""
    if shutil.which(mk_export_binary) is None:
        raise RuntimeError(f"Meeko export executable not found in PATH: {mk_export_binary}")

    pdbqt_path = Path(pdbqt_path).expanduser().resolve()
    output_sdf = Path(output_sdf).expanduser().resolve()
    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    if not pdbqt_path.exists():
        raise FileNotFoundError(f"PDBQT file not found: {pdbqt_path}")
    if not pdbqt_path.is_file():
        raise ValueError(f"PDBQT path is not a file: {pdbqt_path}")

    output_flag = _detect_mk_export_output_flag(mk_export_binary)

    raw_output_sdf = output_sdf
    if template_bond_orders:
        if template_mol is None:
            raise ValueError("template_mol is required when template_bond_orders=True.")
        raw_output_sdf = output_sdf.with_name(f"{output_sdf.stem}_raw{output_sdf.suffix}")

    process: CompletedProcess[str] = run(  # noqa: S603
        [mk_export_binary, str(pdbqt_path), output_flag, str(raw_output_sdf)],
        capture_output=True,
        text=True,
        check=False,
    )

    if process.returncode != 0:
        raise RuntimeError(
            "mk_export.py failed.\n"
            f"PDBQT: {pdbqt_path}\n"
            f"Expected SDF: {raw_output_sdf}\n"
            f"STDOUT:\n{process.stdout}\n"
            f"STDERR:\n{process.stderr}"
        )

    if not raw_output_sdf.exists():
        produced_sdfs = sorted(raw_output_sdf.parent.glob("*.sdf"))
        raise RuntimeError(
            "mk_export.py returned success but expected SDF was not created.\n"
            f"Expected: {raw_output_sdf}\n"
            f"Produced SDF files: {[str(p) for p in produced_sdfs]}\n"
            f"STDOUT:\n{process.stdout}\n"
            f"STDERR:\n{process.stderr}"
        )

    if not template_bond_orders:
        return raw_output_sdf

    raw_molecules = read_sdf_molecules(raw_output_sdf)
    rebuilt_molecules: list[Chem.Mol] = []

    for mol in raw_molecules:
        rebuilt = assign_bond_orders_from_template_mol(
            template_mol=template_mol,
            target_mol=mol,
            add_hydrogens=add_hydrogens_after_template,
        )
        if mol.HasProp("_Name"):
            rebuilt.SetProp("_Name", mol.GetProp("_Name"))
        rebuilt_molecules.append(rebuilt)

    write_sdf_molecules(rebuilt_molecules, output_sdf)
    return output_sdf


# ---------------------------------------------------------------------------
# Receptor preparation
# ---------------------------------------------------------------------------


def convert_receptor_pdb_to_pdbqt(
    receptor_pdb: Path,
    output_path: Path | None = None,
    *,
    obabel_binary: str = "obabel",
    preserve_atom_names: bool = True,
    preserve_atom_indices: bool = True,
    preserve_hydrogens: bool = True,
) -> Path:
    """Convert a receptor PDB to rigid receptor PDBQT using Open Babel."""
    receptor_pdb = Path(receptor_pdb).resolve()

    if not receptor_pdb.exists():
        raise FileNotFoundError(f"Receptor PDB not found: {receptor_pdb}")

    if receptor_pdb.suffix.lower() != ".pdb":
        raise ValueError(f"Expected a .pdb receptor input, got: {receptor_pdb}")

    if shutil.which(obabel_binary) is None:
        raise RuntimeError(f"Open Babel executable not found in PATH: {obabel_binary}")

    if output_path is None:
        output_path = receptor_pdb.with_suffix(".pdbqt")

    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    log_path = output_path.with_suffix(".obabel.log")

    cmd = [
        obabel_binary,
        str(receptor_pdb),
        "-O",
        str(output_path),
        "-xr",
    ]

    if preserve_atom_names:
        cmd.append("-xn")

    if preserve_atom_indices:
        cmd.append("-xp")

    if preserve_hydrogens:
        cmd.append("-xh")

    _step("obabel", f"converting receptor: {receptor_pdb.name} -> {output_path.name}")

    process: CompletedProcess[str] = run(  # noqa: S603
        cmd,
        capture_output=True,
        text=True,
        check=False,
    )

    _write_process_log(
        log_path,
        process,
        command=cmd,
        title=f"Open Babel receptor conversion log for {receptor_pdb.name}",
    )

    if process.returncode != 0:
        raise RuntimeError(
            "Open Babel receptor conversion failed.\n"
            f"Receptor: {receptor_pdb}\n"
            f"Log: {log_path}\n"
            f"stderr summary: {_summarize_stderr(process.stderr)}"
        )

    if not output_path.exists():
        raise RuntimeError(
            "Open Babel reported success but receptor PDBQT output is missing.\n"
            f"Expected: {output_path}\n"
            f"Log: {log_path}"
        )

    if process.stderr.strip():
        _warn("obabel", f"finished with warnings; full details in {log_path.name}")
    else:
        _step("obabel", f"finished successfully; log written to {log_path.name}")

    return output_path


# ---------------------------------------------------------------------------
# Ligand preparation for docking
# ---------------------------------------------------------------------------


def prepare_ligand_with_meeko(
    ligand: str | Path | Chem.Mol,
    output_path: Path,
    name: str | None = None,
    *,
    rigid_macrocycles: bool = False,
    do_geometry_optimization: bool = True,
    mmff_variant: str = "MMFF94",
    max_geoopt_iters: int = 2000,
    random_seed: int = 20260410,
    input_sdf_path: Path | None = None,
) -> Path:
    """Prepare a ligand as PDBQT using RDKit and Meeko."""
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    _step("meeko", f"preparing ligand -> {output_path.name}")

    if input_sdf_path is None:
        input_sdf_path = output_path.with_name(f"{output_path.stem}_input_for_meeko.sdf")

    mol, ligand_sdf_path = ensure_ligand_sdf_for_meeko(
        ligand,
        input_sdf_path,
        name=name,
    )

    _step("meeko", f"ligand SDF for Meeko: {ligand_sdf_path}")

    if do_geometry_optimization:
        mol = ligand_geoopt(
            mol,
            random_seed=random_seed,
            mmff_variant=mmff_variant,
            max_iters=max_geoopt_iters,
        )
    else:
        mol = add_hydrogens_and_3d_if_needed(
            mol,
            random_seed=random_seed,
        )

    prep = MoleculePreparation(
        merge_these_atom_types=(),
        rigid_macrocycles=rigid_macrocycles,
    )
    mol_setups = prep.prepare(mol)

    if not mol_setups:
        raise RuntimeError("Meeko produced no molecule setups.")

    meeko_result = PDBQTWriterLegacy.write_string(mol_setups[0])
    pdbqt_string = _extract_pdbqt_string_from_meeko_result(meeko_result)

    if not pdbqt_string.strip():
        raise RuntimeError("Generated ligand PDBQT string is empty.")

    output_path.write_text(pdbqt_string, encoding="utf-8")
    _step("meeko", f"ligand PDBQT written: {output_path.name}")

    return output_path


def prepare_ligand_for_vina(
    ligand_input: str | Path | Chem.Mol,
    *,
    name: str,
    workdir: Path,
    rigid_macrocycles: bool = False,
    random_seed: int = 20260410,
    mmff_variant: str = "MMFF94",
    max_geoopt_iters: int = 2000,
) -> dict[str, Any]:
    """Prepare one ligand for a standard Vina workflow."""
    workdir = Path(workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    input_for_meeko_sdf = workdir / f"{name}_input_for_meeko.sdf"
    input_mol, resolved_input_sdf = ensure_ligand_sdf_for_meeko(
        ligand_input,
        input_for_meeko_sdf,
        name=name,
    )

    before_meeko_mol = ligand_geoopt(
        input_mol,
        random_seed=random_seed,
        mmff_variant=mmff_variant,
        max_iters=max_geoopt_iters,
    )
    before_meeko_sdf = workdir / f"{name}_before_meeko.sdf"
    write_single_sdf(before_meeko_mol, before_meeko_sdf)

    after_meeko_pdbqt = workdir / f"{name}_after_meeko.pdbqt"
    prepare_ligand_with_meeko(
        before_meeko_mol,
        after_meeko_pdbqt,
        name=name,
        rigid_macrocycles=rigid_macrocycles,
        do_geometry_optimization=False,
        mmff_variant=mmff_variant,
        max_geoopt_iters=max_geoopt_iters,
        random_seed=random_seed,
        input_sdf_path=before_meeko_sdf,
    )

    return {
        "input_mol": input_mol,
        "resolved_input_sdf": resolved_input_sdf,
        "before_meeko_mol": before_meeko_mol,
        "before_meeko_sdf": before_meeko_sdf,
        "after_meeko_pdbqt": after_meeko_pdbqt,
    }


# ---------------------------------------------------------------------------
# Vina engine
# ---------------------------------------------------------------------------


class VinaEngine:
    """Minimal AutoDock Vina wrapper with compact terminal output."""

    name = "vina"

    def __init__(
        self,
        binary: str = "vina",
        obabel_binary: str = "obabel",
    ):
        """Initialize Vina and Open Babel executable names."""
        if shutil.which(binary) is None:
            _warn("vina", f"binary not found in PATH: {binary}")

        if shutil.which(obabel_binary) is None:
            _warn("obabel", f"binary not found in PATH: {obabel_binary}")

        self.binary = binary
        self.obabel_binary = obabel_binary

    def _build_command(
        self,
        *,
        receptor: Path,
        ligand: Path,
        output: Path,
        box: DockingBox,
        seed: int | None = None,
        num_modes: int | None = None,
        exhaustiveness: int | None = None,
        energy_range: float | None = None,
        extra_flags: Mapping[str, Any] | None = None,
    ) -> list[str]:
        """Build the Vina command line for a single ligand."""
        cmd: list[str] = [
            self.binary,
            "--receptor",
            str(Path(receptor)),
            "--ligand",
            str(Path(ligand)),
            "--center_x",
            str(box.center[0]),
            "--center_y",
            str(box.center[1]),
            "--center_z",
            str(box.center[2]),
            "--size_x",
            str(box.size[0]),
            "--size_y",
            str(box.size[1]),
            "--size_z",
            str(box.size[2]),
            "--out",
            str(Path(output)),
        ]

        if seed is not None:
            cmd += ["--seed", str(seed)]

        if num_modes is not None:
            cmd += ["--num_modes", str(num_modes)]

        if exhaustiveness is not None:
            cmd += ["--exhaustiveness", str(exhaustiveness)]

        if energy_range is not None:
            cmd += ["--energy_range", str(energy_range)]

        if extra_flags:
            for key, value in extra_flags.items():
                cmd.append(str(key))
                if value is not None:
                    cmd.append(str(value))

        return cmd

    def _prepare_receptor_for_docking(
        self,
        receptor: Path,
        workdir: Path,
    ) -> Path:
        """Return the receptor PDBQT to use for docking."""
        receptor = Path(receptor).resolve()

        if receptor.suffix.lower() == ".pdbqt":
            _step("vina", f"using receptor PDBQT directly: {receptor.name}")
            return receptor

        if receptor.suffix.lower() == ".pdb":
            return convert_receptor_pdb_to_pdbqt(
                receptor,
                output_path=workdir / f"{receptor.stem}.pdbqt",
                obabel_binary=self.obabel_binary,
            )

        raise ValueError(f"Unsupported receptor file type for docking: {receptor}")

    def dock(self, request: DockingRequest) -> DockingResult:
        """Run docking for all ligands in the request."""
        _section("DOCKING")

        workdir = request.workdir or Path.cwd()
        workdir = Path(workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)

        _step("vina", f"workdir: {workdir}")

        receptor_for_docking = self._prepare_receptor_for_docking(
            request.receptor,
            workdir,
        )

        poses: list[DockedPose] = []
        raw_outputs: dict[str, Any] = {}

        for ligand_entry in request.ligands:
            ligand_path = Path(ligand_entry).resolve()

            out_file = workdir / f"{ligand_path.stem}_vina_out.pdbqt"
            log_file = workdir / f"{ligand_path.stem}_vina.log"

            cmd = self._build_command(
                receptor=receptor_for_docking,
                ligand=ligand_path,
                output=out_file,
                box=request.box,
                seed=request.seed,
                num_modes=request.parameters.get("num_modes"),
                exhaustiveness=request.parameters.get("exhaustiveness"),
                energy_range=request.parameters.get("energy_range"),
                extra_flags=request.parameters.get("extra_flags"),
            )

            _step("vina", f"docking ligand: {ligand_path.name}")

            process: CompletedProcess[str] = run(  # noqa: S603
                cmd,
                cwd=workdir,
                capture_output=True,
                text=True,
                check=False,
            )

            _write_process_log(
                log_file,
                process,
                command=cmd,
                title=f"Vina docking log for {ligand_path.name}",
            )

            best_score = _parse_vina_best_score(process.stdout)

            compact_metadata = _build_compact_pose_metadata(
                returncode=process.returncode,
                log_file=log_file,
                output_exists=out_file.exists(),
                receptor_used=receptor_for_docking,
            )

            if process.returncode == 0 and out_file.exists():
                if best_score is not None:
                    _step(
                        "vina",
                        f"success: {out_file.name} written; best score {best_score:.3f} kcal/mol; log: {log_file.name}",
                    )
                else:
                    _step(
                        "vina",
                        f"success: {out_file.name} written; score not parsed; log: {log_file.name}",
                    )
            elif process.returncode == 0 and not out_file.exists():
                _warn(
                    "vina",
                    f"process returned 0 but output missing; log: {log_file.name}",
                )
            else:
                _error(
                    "vina",
                    f"failed: return code {process.returncode}; "
                    f"stderr summary: {_summarize_stderr(process.stderr)}; "
                    f"log: {log_file.name}",
                )

            poses.append(
                DockedPose(
                    ligand=ligand_path,
                    pose_path=out_file,
                    score=best_score,
                    rank=1 if best_score is not None else None,
                    engine=self.name,
                    metadata=compact_metadata,
                )
            )

            raw_outputs[str(ligand_path)] = compact_metadata

        return DockingResult(
            poses=poses,
            engine=self.name,
            parameters=request.parameters,
            raw_outputs=raw_outputs,
        )


# ---------------------------------------------------------------------------
# Higher-level native-redocking workflow
# ---------------------------------------------------------------------------


def redock_native_and_build_box(
    *,
    receptor: Path,
    native_ligand: str | Path | Chem.Mol,
    workdir: Path,
    native_name: str = "native_ligand",
    vina_binary: str = "vina",
    obabel_binary: str = "obabel",
    mk_export_binary: str = "mk_export.py",
    template_bond_orders: bool = True,
    rigid_macrocycles: bool = False,
    seed: int = 42,
    mmff_variant: str = "MMFF94",
    max_geoopt_iters: int = 2000,
    initial_box: DockingBox | None = None,
    initial_box_padding: float = 5.0,
    initial_box_min_size: float = 20.0,
    query_box_padding: float = 4.0,
    query_box_min_size: float = 0.0,
    native_num_modes: int = 1,
    native_exhaustiveness: int = 16,
    native_energy_range: float = 3.0,
) -> NativeRedockBoxResult:
    """Prepare a native ligand, redock it, and build a final query box."""
    workdir = Path(workdir).expanduser().resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    native_bundle = prepare_ligand_for_vina(
        native_ligand,
        name=native_name,
        workdir=workdir,
        rigid_macrocycles=rigid_macrocycles,
        random_seed=seed,
        mmff_variant=mmff_variant,
        max_geoopt_iters=max_geoopt_iters,
    )

    if initial_box is None:
        initial_box = build_box_from_sdf(
            native_bundle["before_meeko_sdf"],
            padding=initial_box_padding,
            min_size=initial_box_min_size,
        )

    request = DockingRequest(
        receptor=receptor,
        ligands=[native_bundle["after_meeko_pdbqt"]],
        box=initial_box,
        seed=seed,
        workdir=workdir,
        parameters={
            "num_modes": native_num_modes,
            "exhaustiveness": native_exhaustiveness,
            "energy_range": native_energy_range,
        },
    )

    engine = VinaEngine(binary=vina_binary, obabel_binary=obabel_binary)
    result = engine.dock(request)

    successful_poses = [
        pose for pose in result.poses if pose.pose_path.exists() and pose.metadata.get("returncode") == 0
    ]
    if not successful_poses:
        raise RuntimeError("Native redocking produced no successful output pose.")

    best_pose = min(
        successful_poses,
        key=lambda pose: float("inf") if pose.score is None else pose.score,
    )

    best_pose_export_sdf = workdir / f"{native_name}_after_docking_export.sdf"
    export_pdbqt_to_sdf(
        best_pose.pose_path,
        best_pose_export_sdf,
        mk_export_binary=mk_export_binary,
        template_mol=native_bundle["before_meeko_mol"],
        template_bond_orders=template_bond_orders,
        add_hydrogens_after_template=True,
    )

    query_box = build_box_from_sdf(
        best_pose_export_sdf,
        padding=query_box_padding,
        min_size=query_box_min_size,
    )

    return NativeRedockBoxResult(
        resolved_input_sdf=native_bundle["resolved_input_sdf"],
        before_meeko_sdf=native_bundle["before_meeko_sdf"],
        native_pdbqt=native_bundle["after_meeko_pdbqt"],
        initial_redock_box=initial_box,
        docking_result=result,
        best_pose_pdbqt=best_pose.pose_path,
        best_pose_export_sdf=best_pose_export_sdf,
        query_box=query_box,
        best_score=best_pose.score,
    )
