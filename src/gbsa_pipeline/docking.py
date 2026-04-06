# /home/grheco/repositorios/gbsa-pipeline/src/gbsa_pipeline/docking.py

"""Docking helpers for preparing ligands and running AutoDock Vina.

This module provides:
- lightweight request/result models for docking
- ligand preparation to PDBQT using RDKit + Meeko
- optional receptor conversion from PDB to receptor PDBQT via Open Babel
- a small Vina engine wrapper that writes full subprocess output to log files
  while keeping terminal output compact

Design notes
------------
- terminal messages are intentionally short and human-readable
- full subprocess stdout/stderr always goes to log files
- this module writes its own log messages to stdout, not stderr
- failures are represented in returned metadata and/or raised exceptions,
  depending on where they happen

Important:
---------
This module is a reusable library layer. It is not a full CLI command.
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
from rdkit.Chem.rdDistGeom import EmbedMolecule
from rdkit.Chem.rdForceFieldHelpers import UFFOptimizeMolecule

if TYPE_CHECKING:
    from collections.abc import Mapping

LOGGER = logging.getLogger(__name__)
LOGGER.setLevel(logging.INFO)

if not LOGGER.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(message)s"))
    LOGGER.addHandler(handler)

LOGGER.propagate = False


class DockingBox(BaseModel):
    """Docking-box center and size in Angstrom."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    center: tuple[float, float, float]
    size: tuple[float, float, float]


class DockingRequest(BaseModel):
    """Normalized docking request for a receptor, one or more ligands, and a box."""

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
    """Single docked pose plus compact metadata about the docking run."""

    ligand: Path
    pose_path: Path
    score: float | None
    rank: int | None
    engine: str
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class DockingResult:
    """Collection of produced poses and compact run metadata."""

    poses: list[DockedPose]
    engine: str
    parameters: Mapping[str, Any]
    raw_outputs: Mapping[str, Any] = field(default_factory=dict)


class DockingEngine(Protocol):
    """Protocol implemented by docking backends."""

    name: str

    def dock(self, request: DockingRequest) -> DockingResult:
        """Run docking for the provided request."""


def _box(title: str, char: str = "=") -> str:
    line = char * max(len(title), 12)
    return f"{line}\n{title}\n{line}"


def _section(title: str) -> None:
    LOGGER.info("")
    LOGGER.info(_box(title))


def _step(prefix: str, message: str) -> None:
    LOGGER.info("[%s] %s", prefix, message)


def _warn(prefix: str, message: str) -> None:
    LOGGER.warning("[%s] warning: %s", prefix, message)


def _error(prefix: str, message: str) -> None:
    LOGGER.error("[%s] error: %s", prefix, message)


def _summarize_stderr(stderr: str, max_lines: int = 4) -> str:
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
    log_path.write_text(text)


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
    """Parse the best affinity from the first row of Vina's result table."""
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
    """Build compact per-pose metadata stored in the returned result."""
    return {
        "returncode": returncode,
        "log_file": str(log_file),
        "output_exists": output_exists,
        "receptor_used": str(receptor_used),
    }


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


def prepare_ligand_with_meeko(
    ligand: str | Chem.Mol,
    output_path: Path,
    name: str | None = None,
) -> Path:
    """Prepare a ligand as ligand PDBQT using RDKit and Meeko."""
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    _step("meeko", f"preparing ligand -> {output_path.name}")

    if isinstance(ligand, str):
        mol = Chem.MolFromSmiles(ligand)

        if mol is None:
            raise ValueError("Failed to parse SMILES.")

        mol.SetProp("_Name", name or "LIG")
        mol = Chem.AddHs(mol)

        embed_status = EmbedMolecule(mol)
        if embed_status != 0:
            raise RuntimeError(f"RDKit failed to embed 3D coordinates for ligand: {ligand}")

        optimize_status = UFFOptimizeMolecule(mol)
        if optimize_status not in (0, 1):
            _warn(
                "meeko",
                f"UFF optimization returned non-standard status {optimize_status}",
            )

    else:
        mol = ligand

        if mol is None:
            raise ValueError("Ligand molecule is None.")

        if name is not None:
            mol.SetProp("_Name", name)
        elif not mol.HasProp("_Name"):
            mol.SetProp("_Name", "LIG")

    prep = MoleculePreparation()
    mol_setups = prep.prepare(mol)

    if not mol_setups:
        raise RuntimeError("Meeko produced no molecule setups.")

    meeko_result = PDBQTWriterLegacy.write_string(mol_setups[0])
    pdbqt_string = _extract_pdbqt_string_from_meeko_result(meeko_result)

    if not pdbqt_string.strip():
        raise RuntimeError("Generated ligand PDBQT string is empty.")

    output_path.write_text(pdbqt_string)
    _step("meeko", f"ligand PDBQT written: {output_path.name}")

    return output_path


class VinaEngine:
    """Minimal AutoDock Vina wrapper with compact console output and log files."""

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
