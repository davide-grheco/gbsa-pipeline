# /home/grheco/repositorios/gbsa-pipeline/src/gbsa_pipeline/unigbsa.py

"""Helpers for preparing and running gmx_MMPBSA / Uni-GBSA style calculations."""

from __future__ import annotations

import json
import subprocess
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Any

PathLike = str | Path


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _fmt(v: Any) -> str:
    """Format a Python value for an MMPBSA namelist."""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (int, float)):
        return str(v)
    if v is None:
        raise ValueError("None should be filtered before formatting")
    s = str(v)
    if s == "" or any(c.isspace() for c in s) or "," in s:
        return f'"{s}"'
    return s


def _as_kv(dc: Any) -> dict[str, Any]:
    """Convert a configuration dataclass to a flat key-value mapping."""
    d: dict[str, Any] = {}
    extra: dict[str, Any] = {}

    for f in fields(dc):
        val = getattr(dc, f.name)
        if f.name == "extra":
            extra = val or {}
            continue
        if val is not None:
            d[f.name] = val

    d.update(extra)
    return d


def _render_namelist(name: str, kv: Mapping[str, Any]) -> str:
    """Render one Fortran-style namelist block."""
    lines = [f"&{name}"]
    for k, v in kv.items():
        lines.append(f"  {k:<20} = {_fmt(v)}")
    lines.append("/")
    return "\n".join(lines)


def _to_path(value: PathLike) -> Path:
    """Convert a path-like value to ``Path``."""
    return value if isinstance(value, Path) else Path(value)


def _jsonify(value: Any) -> Any:
    """Convert nested objects into JSON-serializable structures."""
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value) and not isinstance(value, type):
        return {k: _jsonify(v) for k, v in asdict(value).items()}
    if isinstance(value, Mapping):
        return {str(k): _jsonify(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonify(v) for v in value]
    return value


# -----------------------------------------------------------------------------
# Config dataclasses
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class GeneralParams:
    """General control parameters for an MMPBSA input file."""

    sys_name: str = ""
    startframe: int = 1
    endframe: int = 9_999_999
    interval: int = 1
    temperature: float = 298.15
    keep_files: int = 2
    netcdf: int = 0
    solvated_trajectory: int = 1
    verbose: int = 1
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class GBParams:
    """Generalized-Born parameter block for an MMPBSA input file."""

    igb: int = 5
    intdiel: float = 1.0
    extdiel: float = 78.5
    saltcon: float = 0.0
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PBParams:
    """Poisson-Boltzmann parameter block for an MMPBSA input file."""

    ipb: int = 2
    indi: float = 1.0
    exdi: float = 80.0
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class MMPBSAConfig:
    """Complete MMPBSA namelist configuration."""

    general: GeneralParams = field(default_factory=GeneralParams)
    gb: GBParams | None = field(default_factory=GBParams)
    pb: PBParams | None = field(default_factory=PBParams)

    def to_text(self) -> str:
        """Render the MMPBSA configuration as a namelist text file."""
        parts = [_render_namelist("general", _as_kv(self.general))]

        if self.gb is not None:
            parts.append(_render_namelist("gb", _as_kv(self.gb)))

        if self.pb is not None:
            parts.append(_render_namelist("pb", _as_kv(self.pb)))

        return "\n\n".join(parts) + "\n"

    def write(self, path: PathLike) -> Path:
        """Write the MMPBSA configuration to disk and return the output path."""
        path = _to_path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.to_text(), encoding="utf-8")
        return path


# -----------------------------------------------------------------------------
# Request / Result
# -----------------------------------------------------------------------------


@dataclass(frozen=True)
class GmxMmpbsaRequest:
    """Request object describing one gmx_MMPBSA run."""

    complex_structure: PathLike
    trajectory: PathLike
    topology: PathLike
    index_file: PathLike
    receptor_group: int
    ligand_group: int
    output_dir: PathLike

    config: MMPBSAConfig = field(default_factory=MMPBSAConfig)
    gmx_mmpbsa: str = "gmx_MMPBSA"
    extra_args: tuple[str, ...] = ()
    input_filename: str = "mmpbsa.in"


@dataclass(frozen=True)
class GmxMmpbsaResult:
    """Structured result returned from a gmx_MMPBSA run."""

    ok: bool
    returncode: int
    command: list[str]
    cwd: str

    stdout_log: str
    stderr_log: str
    result_json: str

    input_file: str | None = None
    complex_structure: str | None = None
    trajectory: str | None = None
    topology: str | None = None
    index_file: str | None = None
    receptor_group: int | None = None
    ligand_group: int | None = None
    final_results_mmpbsa: str | None = None
    final_results_decomp: str | None = None
    info_file: str | None = None

    def write_json(self, path: PathLike) -> None:
        """Write the run result as formatted JSON."""
        path = _to_path(path)
        path.write_text(json.dumps(_jsonify(self), indent=2), encoding="utf-8")


# -----------------------------------------------------------------------------
# Runner
# -----------------------------------------------------------------------------


def run_gmx_mmpbsa(request: GmxMmpbsaRequest) -> GmxMmpbsaResult:
    """Run gmx_MMPBSA from a prepared request object."""
    output_dir = _to_path(request.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    input_file = output_dir / request.input_filename
    stdout_log = output_dir / "gmx_mmpbsa.stdout.log"
    stderr_log = output_dir / "gmx_mmpbsa.stderr.log"
    result_json = output_dir / "result.json"
    final_results_mmpbsa = output_dir / "FINAL_RESULTS_MMPBSA.dat"
    final_results_decomp = output_dir / "FINAL_DECOMP_MMPBSA.dat"
    info_file = output_dir / "gmx_MMPBSA_info.dat"

    request.config.write(input_file)

    cmd = [
        request.gmx_mmpbsa,
        "-O",
        "-i",
        str(input_file),
        "-cs",
        str(request.complex_structure),
        "-ct",
        str(request.trajectory),
        "-cp",
        str(request.topology),
        "-ci",
        str(request.index_file),
        "-cg",
        str(request.receptor_group),
        str(request.ligand_group),
        *request.extra_args,
    ]

    completed = subprocess.run(  # noqa: S603
        cmd,
        cwd=str(output_dir),
        text=True,
        capture_output=True,
        check=False,
    )

    stdout_log.write_text(completed.stdout or "", encoding="utf-8")
    stderr_log.write_text(completed.stderr or "", encoding="utf-8")

    result = GmxMmpbsaResult(
        ok=(completed.returncode == 0),
        returncode=completed.returncode,
        command=cmd,
        cwd=str(output_dir),
        stdout_log=str(stdout_log),
        stderr_log=str(stderr_log),
        result_json=str(result_json),
        input_file=str(input_file),
        complex_structure=str(_to_path(request.complex_structure)),
        trajectory=str(_to_path(request.trajectory)),
        topology=str(_to_path(request.topology)),
        index_file=str(_to_path(request.index_file)),
        receptor_group=request.receptor_group,
        ligand_group=request.ligand_group,
        final_results_mmpbsa=str(final_results_mmpbsa),
        final_results_decomp=str(final_results_decomp),
        info_file=str(info_file),
    )

    result.write_json(result_json)
    return result


# -----------------------------------------------------------------------------
# Legacy wrapper
# -----------------------------------------------------------------------------


def run_gmx_mmpbsa_from_gromacs(
    input_file: PathLike,
    complex_structure: PathLike,
    trajectory: PathLike,
    topology: PathLike,
    index_file: PathLike,
    receptor_group: int,
    ligand_group: int,
    output_dir: PathLike,
    gmx_mmpbsa: str = "gmx_MMPBSA",
    extra_args: Sequence[str] = (),
) -> subprocess.CompletedProcess[str]:
    """Run gmx_MMPBSA directly from GROMACS-style input files."""
    output_dir = _to_path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        gmx_mmpbsa,
        "-O",
        "-i",
        str(input_file),
        "-cs",
        str(complex_structure),
        "-ct",
        str(trajectory),
        "-cp",
        str(topology),
        "-ci",
        str(index_file),
        "-cg",
        str(receptor_group),
        str(ligand_group),
        *extra_args,
    ]

    return subprocess.run(  # noqa: S603
        cmd,
        cwd=str(output_dir),
        text=True,
        capture_output=True,
        check=False,
    )
