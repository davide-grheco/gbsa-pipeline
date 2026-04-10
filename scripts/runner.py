# /home/grheco/repositorios/gbsa-pipeline/scripts/runner.py

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
import traceback
from collections import defaultdict
from contextlib import contextmanager, redirect_stderr, redirect_stdout, suppress
from datetime import UTC, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, TextIO

if TYPE_CHECKING:
    from collections.abc import Iterator

from rdkit import Chem

import gbsa_pipeline.docking as docking_module
import gbsa_pipeline.md as md_module
from gbsa_pipeline.docking import (
    DockingBox,
    DockingRequest,
    VinaEngine,
    prepare_ligand_with_meeko,
)
from gbsa_pipeline.md import run_md_from_docking
from gbsa_pipeline.unigbsa import (
    GeneralParams,
    GmxMmpbsaRequest,
    MMPBSAConfig,
    run_gmx_mmpbsa,
)

STATE_FILENAME = "pipeline_state.json"
LOGS_DIRNAME = "logs"
SUMMARY_FILENAME = "run_summary.txt"
STEP_ORDER = {"docking": 0, "md": 1, "gbsa": 2}

MIN_GRO_TOTAL_LINES = 3
MIN_GRO_ATOM_LINE_LENGTH = 20

MD_BANNER_PENDING_NONE = 0
MD_BANNER_PENDING_TOP = 1
MD_BANNER_PENDING_BOTTOM = 2

PROTEIN_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "HID",
    "HIE",
    "HIP",
    "ASH",
    "GLH",
    "CYM",
    "CYX",
    "LYN",
    "ACE",
    "NME",
    "NMA",
    "NASP",
    "CASP",
    "NGLU",
    "CGLU",
    "NHID",
    "NHIE",
    "NHIP",
    "CHID",
    "CHIE",
    "CHIP",
    "NCYS",
    "CCYS",
    "NLYS",
    "CLYS",
    "NARG",
    "CARG",
}

SOLVENT_ION_RESNAMES = {
    "HOH",
    "WAT",
    "SOL",
    "TIP3",
    "TIP4",
    "TIP5",
    "SPC",
    "SPCE",
    "NA",
    "NA+",
    "CL",
    "CL-",
    "K",
    "K+",
    "CA",
    "CA2",
    "MG",
    "MG2",
    "ZN",
    "ZN2",
    "MN",
    "MN2",
    "FE",
    "FE2",
    "FE3",
    "CU",
    "CU1",
    "CU2",
}


def _require_text_stream(stream: TextIO | None, label: str) -> TextIO:
    """Return a guaranteed text stream or raise a clear runtime error."""
    if stream is None:
        raise RuntimeError(f"{label} is not available.")
    return stream


def _current_stdout() -> TextIO:
    """Return the currently active stdout as a guaranteed text stream."""
    return _require_text_stream(sys.stdout, "sys.stdout")


def _current_stderr() -> TextIO:
    """Return the currently active stderr as a guaranteed text stream."""
    return _require_text_stream(sys.stderr, "sys.stderr")


def _original_stdout() -> TextIO:
    """Return the original interpreter stdout as a guaranteed text stream."""
    return _require_text_stream(sys.__stdout__, "sys.__stdout__")


def _original_stderr() -> TextIO:
    """Return the original interpreter stderr as a guaranteed text stream."""
    return _require_text_stream(sys.__stderr__, "sys.__stderr__")


class TeeStream:
    """Write output to terminal and one or more log files."""

    def __init__(self, *streams: TextIO) -> None:
        """Store all target streams."""
        self.streams = streams

    def write(self, data: str) -> int:
        """Write data to all configured streams."""
        for stream in self.streams:
            stream.write(data)
        return len(data)

    def flush(self) -> None:
        """Flush all configured streams."""
        for stream in self.streams:
            stream.flush()

    def isatty(self) -> bool:
        """Return False so redirected output is treated as non-interactive."""
        return False


class StageLoggerBinding:
    """Temporarily bind imported module loggers to redirected stdout."""

    def __init__(self) -> None:
        """Initialize restore bookkeeping for temporary logger redirection."""
        self._restores: list[tuple[logging.StreamHandler[TextIO], TextIO]] = []
        self._added_handlers: list[tuple[logging.Logger, logging.StreamHandler[TextIO]]] = []

    def bind(self) -> None:
        """Redirect imported module stream handlers to current stdout."""
        self._bind_existing_stream_handlers(docking_module.LOGGER)
        self._ensure_simple_stdout_handler(md_module.LOGGER)
        self._bind_existing_stream_handlers(md_module.LOGGER)

    def restore(self) -> None:
        """Restore original logger streams and remove temporary handlers."""
        for handler, old_stream in self._restores:
            handler.setStream(old_stream)
        for logger, handler in self._added_handlers:
            with suppress(Exception):
                logger.removeHandler(handler)

    def _bind_existing_stream_handlers(self, logger: logging.Logger) -> None:
        for handler in logger.handlers:
            if isinstance(handler, logging.StreamHandler):
                typed_handler = handler
                self._restores.append((typed_handler, typed_handler.stream))
                typed_handler.setStream(_current_stdout())

    def _ensure_simple_stdout_handler(self, logger: logging.Logger) -> None:
        if logger.handlers:
            return

        handler: logging.StreamHandler[TextIO] = logging.StreamHandler(_current_stdout())
        handler.setLevel(logging.INFO)
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.propagate = False
        self._added_handlers.append((logger, handler))


def _now() -> str:
    """Return a timezone-aware current timestamp string."""
    return datetime.now(UTC).strftime("%Y-%m-%d %H:%M:%S %Z")


def _elapsed_text(seconds: float) -> str:
    """Format elapsed seconds as MM:SS.xx or HH:MM:SS.xx."""
    minutes, sec = divmod(seconds, 60.0)
    hours, minutes = divmod(minutes, 60.0)
    if hours >= 1:
        return f"{int(hours):02d}:{int(minutes):02d}:{sec:05.2f}"
    return f"{int(minutes):02d}:{sec:05.2f}"


def _box_lines(
    title: str,
    lines: list[str] | None = None,
    *,
    heavy: bool = False,
    center_title: bool = False,
    style: str = "standard",
) -> str:
    """Build a framed text box."""
    lines = lines or []
    width = max([len(title), *(len(line) for line in lines), 40])

    if style == "wide_header":
        tl, tr, bl, br, hz, vt = "╒", "╕", "╘", "╛", "═", "│"
        sep_left, sep_mid, sep_right = "╞", "═", "╡"
    elif heavy:
        tl, tr, bl, br, hz, vt = "╔", "╗", "╚", "╝", "═", "║"
        sep_left, sep_mid, sep_right = "╟", "─", "╢"
    else:
        tl, tr, bl, br, hz, vt = "┌", "┐", "└", "┘", "─", "│"
        sep_left, sep_mid, sep_right = "├", "─", "┤"

    title_text = title.center(width) if center_title else title.ljust(width)

    text: list[str] = [
        f"{tl}{hz * (width + 2)}{tr}",
        f"{vt} {title_text} {vt}",
    ]
    if lines:
        text.append(f"{sep_left}{sep_mid * (width + 2)}{sep_right}")
        text.extend(f"{vt} {line.ljust(width)} {vt}" for line in lines)
    text.append(f"{bl}{hz * (width + 2)}{br}")
    return "\n".join(text)


def _print_box(
    title: str,
    lines: list[str] | None = None,
    *,
    heavy: bool = False,
    center_title: bool = False,
    style: str = "standard",
) -> None:
    """Print a framed text box."""
    print(_box_lines(title, lines, heavy=heavy, center_title=center_title, style=style))


def _print_spacer(lines: int = 1) -> None:
    """Print one or more blank lines."""
    for _ in range(lines):
        print()


def _install_md_box_style() -> None:
    """Patch md module banner logging to use the same boxed style as the runner."""
    original_info: Callable[..., Any] = md_module.LOGGER.info
    state = {"pending": MD_BANNER_PENDING_NONE}

    def styled_info(msg: object, *args: Any, **kwargs: Any) -> None:
        text_msg = str(msg)

        if text_msg == "====================================================":
            if state["pending"] == MD_BANNER_PENDING_NONE:
                state["pending"] = MD_BANNER_PENDING_TOP
                return
            if state["pending"] == MD_BANNER_PENDING_BOTTOM:
                _print_spacer(1)
                print(_box_lines("MD", center_title=True, style="wide_header"))
                _print_spacer(1)
                state["pending"] = MD_BANNER_PENDING_NONE
                return

        if state["pending"] == MD_BANNER_PENDING_TOP and text_msg == "MD":
            state["pending"] = MD_BANNER_PENDING_BOTTOM
            return

        if state["pending"] != MD_BANNER_PENDING_NONE:
            original_info("====================================================")
            original_info("MD")
            state["pending"] = MD_BANNER_PENDING_NONE

        original_info(msg, *args, **kwargs)

    md_module.LOGGER.info = styled_info  # type: ignore[method-assign]


def _install_docking_box_style() -> None:
    """Patch docking module banner rendering to match runner style."""

    def _styled_docking_box(title: str, char: str = "=") -> str:
        width = max(len(title), 18) + 12
        return "\n".join(
            [
                f"╒{'═' * width}╕",
                f"│{title.center(width)}│",
                f"╘{'═' * width}╛",
            ]
        )

    docking_module._box = _styled_docking_box


def _print_stage_banner(title: str) -> None:
    """Print a centered wide header for a stage."""
    _print_spacer(1)
    _print_box(title.upper(), center_title=True, style="wide_header")
    _print_spacer(1)


def _write_json(path: Path, data: dict[str, Any]) -> None:
    """Write formatted JSON to disk."""
    path = path.expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _read_json(path: Path) -> dict[str, Any]:
    """Read JSON from disk."""
    path = path.expanduser().resolve()
    return json.loads(path.read_text(encoding="utf-8"))


def _ensure_state_defaults(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Fill missing state keys with defaults."""
    defaults = {
        "run_dir": str(run_dir),
        "state_file": str(run_dir / STATE_FILENAME),
        "logs_dir": str(run_dir / LOGS_DIRNAME),
        "created_at": state.get("created_at", _now()),
        "last_updated": _now(),
        "status": state.get("status", "initialized"),
        "failed_stage": state.get("failed_stage"),
        "selected_stages": state.get("selected_stages", []),
        "docking": state.get("docking"),
        "md": state.get("md"),
        "gbsa": state.get("gbsa"),
        "progress": state.get("progress", {}),
    }

    for stage in STEP_ORDER:
        defaults["progress"].setdefault(
            stage,
            {
                "status": "pending",
                "started_at": None,
                "finished_at": None,
                "elapsed_seconds": None,
                "elapsed_text": None,
                "log_file": None,
                "error": None,
            },
        )

    return defaults


def _load_state(run_dir: Path) -> dict[str, Any]:
    """Load state file if present, otherwise create default in-memory state."""
    state_path = run_dir / STATE_FILENAME
    if state_path.exists():
        return _ensure_state_defaults(run_dir, _read_json(state_path))
    return _ensure_state_defaults(run_dir, {})


def _save_state(run_dir: Path, state: dict[str, Any]) -> Path:
    """Persist current state to disk."""
    state["last_updated"] = _now()
    state_path = run_dir / STATE_FILENAME
    _write_json(state_path, state)
    return state_path


def _stage_is_selected(stage: str, start: str, stop: str) -> bool:
    """Return whether a stage falls between selected start and stop stages."""
    return STEP_ORDER[start] <= STEP_ORDER[stage] <= STEP_ORDER[stop]


def _selected_stage_names(start: str, stop: str) -> list[str]:
    """Return ordered stage names included in the current run."""
    return [stage for stage in STEP_ORDER if _stage_is_selected(stage, start, stop)]


def _require_existing_file(path_str: str | Path | None, label: str) -> Path:
    """Validate that a required file path exists and points to a file."""
    if path_str is None:
        raise RuntimeError(f"Missing path for {label}.")

    path = Path(path_str).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")
    return path


def _load_first_sdf_molecule(path: Path) -> Chem.Mol:
    """Load the first valid molecule from an SDF file."""
    path = path.expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"SDF file not found: {path}")
    if not path.is_file():
        raise ValueError(f"SDF path is not a file: {path}")

    try:
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    except Exception as exc:
        raise RuntimeError(f"Could not open SDF with RDKit: {path}") from exc

    for mol in supplier:
        if mol is not None:
            return mol

    raise ValueError(f"Could not read a valid molecule from SDF: {path}")


def _normalize_ligand_for_docking(ligand_path: Path, docking_dir: Path) -> Path:
    """Convert an SDF ligand to PDBQT for docking or reuse an existing PDBQT."""
    ligand_path = ligand_path.expanduser().resolve()

    if ligand_path.suffix.lower() == ".pdbqt":
        return ligand_path

    if ligand_path.suffix.lower() == ".sdf":
        mol = _load_first_sdf_molecule(ligand_path)
        return prepare_ligand_with_meeko(
            mol,
            output_path=docking_dir / f"{ligand_path.stem}.pdbqt",
            name=ligand_path.stem,
        )

    raise ValueError(f"Unsupported ligand format for docking. Expected .sdf or .pdbqt, got: {ligand_path}")


def _auto_box_from_sdf(sdf_path: Path, padding: float, min_size: float) -> DockingBox:
    """Derive a simple docking box from ligand coordinates."""
    mol = _load_first_sdf_molecule(sdf_path)

    if mol.GetNumConformers() == 0:
        raise ValueError(f"SDF has no 3D conformer: {sdf_path}")

    conf = mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []

    for atom_index in range(mol.GetNumAtoms()):
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


def _default_run_dir(inpdb: Path, inpsdf: Path) -> Path:
    """Build the default run directory from input stems."""
    base = Path.cwd() / "pipeline_runs"
    name = f"{inpdb.stem}__{inpsdf.stem}"
    return (base / name).resolve()


def _serialize_docking_pose(pose: Any) -> dict[str, Any]:
    """Convert a DockedPose-like object to JSON-serializable form."""
    return {
        "ligand": str(pose.ligand),
        "pose_path": str(pose.pose_path),
        "score": pose.score,
        "rank": pose.rank,
        "engine": pose.engine,
        "metadata": dict(pose.metadata),
    }


def _load_saved_docking_manifest(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Load docking manifest from state or disk."""
    if state.get("docking"):
        return state["docking"]

    manifest_path = run_dir / "docking" / "result.json"
    if manifest_path.exists():
        manifest = _read_json(manifest_path)
        state["docking"] = manifest
        _save_state(run_dir, state)
        return manifest

    raise RuntimeError("No saved docking result found. Run docking first for this run directory.")


def _load_saved_md_manifest(run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Load MD manifest from state or disk."""
    if state.get("md"):
        return state["md"]

    manifest_path = run_dir / "md" / "md" / "result.json"
    if manifest_path.exists():
        manifest = _read_json(manifest_path)
        state["md"] = manifest
        _save_state(run_dir, state)
        return manifest

    raise RuntimeError("No saved MD result found. Run docking/MD first for this run directory.")


def _split_fixed_gro_fields(line: str) -> tuple[str, str, str, int]:
    """Parse a fixed-width GROMACS GRO atom line."""
    resnum = line[0:5].strip()
    resname = line[5:10].strip()
    atomname = line[10:15].strip()
    atomnum = int(line[15:20].strip())
    return resnum, resname, atomname, atomnum


def _format_ndx_block(name: str, atom_numbers: list[int], width: int = 15) -> str:
    """Render one GROMACS index group block."""
    lines = [f"[ {name} ]"]
    for start in range(0, len(atom_numbers), width):
        chunk = atom_numbers[start : start + width]
        lines.append(" ".join(str(x) for x in chunk))
    lines.append("")
    return "\n".join(lines)


def _write_auto_gbsa_index_from_gro(gro_file: Path, output_ndx: Path) -> tuple[Path, str]:
    """Build a minimal GBSA index file with RECEPTOR and LIGAND groups."""
    gro_file = gro_file.expanduser().resolve()
    if not gro_file.exists():
        raise FileNotFoundError(f"GRO file not found for index generation: {gro_file}")

    lines = gro_file.read_text(encoding="utf-8").splitlines()
    if len(lines) < MIN_GRO_TOTAL_LINES:
        raise ValueError(f"Invalid GRO file: {gro_file}")

    atom_lines = lines[2:-1]
    protein_atoms: list[int] = []
    ligand_candidates: dict[str, list[int]] = defaultdict(list)

    for line in atom_lines:
        if len(line) < MIN_GRO_ATOM_LINE_LENGTH:
            continue

        _, resname, _, atomnum = _split_fixed_gro_fields(line)

        if resname in PROTEIN_RESNAMES:
            protein_atoms.append(atomnum)
            continue

        if resname in SOLVENT_ION_RESNAMES:
            continue

        ligand_candidates[resname].append(atomnum)

    if not protein_atoms:
        raise RuntimeError("Could not identify any protein atoms while building GBSA index.")
    if not ligand_candidates:
        raise RuntimeError("Could not identify any ligand-like residue while building GBSA index.")

    ligand_resname, ligand_atoms = max(
        ligand_candidates.items(),
        key=lambda item: len(item[1]),
    )

    output_ndx = output_ndx.expanduser().resolve()
    output_ndx.parent.mkdir(parents=True, exist_ok=True)

    ndx_text = _format_ndx_block("RECEPTOR", protein_atoms) + _format_ndx_block("LIGAND", ligand_atoms)
    output_ndx.write_text(ndx_text, encoding="utf-8")

    return output_ndx, ligand_resname


def _read_ndx_group_numbers(index_file: Path) -> dict[str, int]:
    """Read 0-based group numbers from a GROMACS .ndx file."""
    index_file = index_file.expanduser().resolve()

    group_names: list[str] = []
    for line in index_file.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            group_names.append(stripped[1:-1].strip())

    return {name: i for i, name in enumerate(group_names)}


def _get_required_group_numbers(
    index_file: Path,
    receptor_name: str = "RECEPTOR",
    ligand_name: str = "LIGAND",
) -> tuple[int, int]:
    """Read required receptor and ligand group numbers from an index file."""
    groups = _read_ndx_group_numbers(index_file)

    if receptor_name not in groups:
        raise RuntimeError(f"Group {receptor_name!r} not found in index file: {index_file}")
    if ligand_name not in groups:
        raise RuntimeError(f"Group {ligand_name!r} not found in index file: {index_file}")

    return groups[receptor_name], groups[ligand_name]


def _summary_lines_from_state(state: dict[str, Any]) -> list[str]:
    """Build human-readable summary lines from pipeline state."""
    lines: list[str] = []
    for stage in STEP_ORDER:
        progress = state["progress"][stage]
        status = progress["status"]
        elapsed = progress["elapsed_text"] or "--"
        log_file = progress["log_file"] or "--"
        lines.append(f"{stage:<8} | {status:<7} | {elapsed:<10} | {log_file}")
    return lines


@contextmanager
def _stage_output_context(stage: str, run_dir: Path) -> Iterator[Path]:
    """Redirect stdout/stderr for one stage to terminal and log files."""
    logs_dir = run_dir / LOGS_DIRNAME
    logs_dir.mkdir(parents=True, exist_ok=True)

    stage_log_path = logs_dir / f"{stage}.log"
    pipeline_log_path = logs_dir / "pipeline.log"

    with (
        stage_log_path.open("a", encoding="utf-8") as stage_log,
        pipeline_log_path.open("a", encoding="utf-8") as pipeline_log,
    ):
        tee_stdout = TeeStream(_original_stdout(), stage_log, pipeline_log)
        tee_stderr = TeeStream(_original_stderr(), stage_log, pipeline_log)
        binding = StageLoggerBinding()

        with redirect_stdout(tee_stdout), redirect_stderr(tee_stderr):
            binding.bind()
            try:
                yield stage_log_path
            finally:
                binding.restore()


def _run_stage_wrapped(
    *,
    stage: str,
    args: argparse.Namespace,
    run_dir: Path,
    state: dict[str, Any],
    runner: Callable[[argparse.Namespace, Path, dict[str, Any]], dict[str, Any]],
) -> dict[str, Any]:
    """Run one pipeline stage with logging, timing, and state persistence."""
    start_wall = _now()
    start_perf = time.perf_counter()
    log_path = run_dir / LOGS_DIRNAME / f"{stage}.log"

    state["progress"][stage].update(
        {
            "status": "running",
            "started_at": start_wall,
            "finished_at": None,
            "elapsed_seconds": None,
            "elapsed_text": None,
            "log_file": str(log_path),
            "error": None,
        }
    )
    _save_state(run_dir, state)

    with _stage_output_context(stage, run_dir) as actual_log_path:
        _print_spacer(1)
        _print_box(
            f"STEP {stage.upper()} STARTED",
            [
                f"timestamp : {start_wall}",
                f"log file  : {actual_log_path}",
            ],
            heavy=False,
            center_title=True,
            style="wide_header",
        )

        try:
            result = runner(args, run_dir, state)
        except Exception as exc:
            traceback.print_exc()
            elapsed = time.perf_counter() - start_perf
            finished = _now()
            state["progress"][stage].update(
                {
                    "status": "failed",
                    "finished_at": finished,
                    "elapsed_seconds": round(elapsed, 3),
                    "elapsed_text": _elapsed_text(elapsed),
                    "error": f"{type(exc).__name__}: {exc}",
                    "log_file": str(actual_log_path),
                }
            )
            state["status"] = "failed"
            state["failed_stage"] = stage
            _save_state(run_dir, state)

            _print_spacer(1)
            _print_box(
                f"STEP {stage.upper()} FAILED",
                [
                    f"timestamp : {finished}",
                    f"elapsed   : {_elapsed_text(elapsed)}",
                    f"error     : {type(exc).__name__}: {exc}",
                    f"log file  : {actual_log_path}",
                ],
                heavy=False,
            )
            raise
        else:
            elapsed = time.perf_counter() - start_perf
            finished = _now()
            state["progress"][stage].update(
                {
                    "status": "done",
                    "finished_at": finished,
                    "elapsed_seconds": round(elapsed, 3),
                    "elapsed_text": _elapsed_text(elapsed),
                    "error": None,
                    "log_file": str(actual_log_path),
                }
            )
            state["status"] = "running"
            state["failed_stage"] = None
            _save_state(run_dir, state)

            _print_spacer(1)
            _print_box(
                f"STEP {stage.upper()} FINISHED",
                [
                    f"timestamp : {finished}",
                    f"elapsed   : {_elapsed_text(elapsed)}",
                    f"log file  : {actual_log_path}",
                ],
                heavy=False,
            )
            return result


def run_docking_stage(args: argparse.Namespace, run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Run docking stage and persist a compact manifest."""
    docking_dir = run_dir / "docking"
    docking_dir.mkdir(parents=True, exist_ok=True)

    protein_pdb = _require_existing_file(args.inpdb, "inpdb")
    ligand_input = _require_existing_file(args.inpsdf, "inpsdf")
    ligand_for_docking = _normalize_ligand_for_docking(ligand_input, docking_dir)
    auto_box = _auto_box_from_sdf(
        ligand_input,
        padding=args.box_padding,
        min_size=args.box_min_size,
    )

    _print_spacer(1)
    _print_box(
        "Docking inputs",
        [
            f"protein       : {protein_pdb}",
            f"ligand input  : {ligand_input}",
            f"ligand pdbqt  : {ligand_for_docking}",
            f"auto center   : {tuple(round(x, 3) for x in auto_box.center)}",
            f"auto size     : {tuple(round(x, 3) for x in auto_box.size)}",
        ],
        heavy=False,
    )

    request = DockingRequest(
        receptor=protein_pdb,
        ligands=[ligand_for_docking],
        box=auto_box,
        seed=args.seed,
        workdir=docking_dir,
        parameters={
            "num_modes": args.num_modes,
            "exhaustiveness": args.exhaustiveness,
            "energy_range": args.energy_range,
        },
    )

    engine = VinaEngine(binary=args.vina_binary, obabel_binary=args.obabel_binary)
    result = engine.dock(request)

    successful_poses = [
        pose for pose in result.poses if pose.pose_path.exists() and pose.metadata.get("returncode") == 0
    ]
    if not successful_poses:
        raise RuntimeError("Docking finished without a successful output pose.")

    selected_pose = min(
        successful_poses,
        key=lambda pose: float("inf") if pose.score is None else pose.score,
    )

    manifest_path = docking_dir / "result.json"
    manifest = {
        "status": "success",
        "manifest_path": str(manifest_path),
        "work_dir": str(docking_dir),
        "protein_pdb": str(protein_pdb),
        "ligand_input": str(ligand_input),
        "ligand_for_docking": str(ligand_for_docking),
        "engine": result.engine,
        "box": {
            "center": list(auto_box.center),
            "size": list(auto_box.size),
            "padding": args.box_padding,
            "min_size": args.box_min_size,
        },
        "parameters": dict(result.parameters),
        "selected_pose": str(selected_pose.pose_path),
        "selected_score": selected_pose.score,
        "poses": [_serialize_docking_pose(pose) for pose in result.poses],
    }

    _write_json(manifest_path, manifest)
    state["docking"] = manifest
    _save_state(run_dir, state)

    _print_spacer(1)
    _print_box(
        "Docking outputs",
        [
            f"manifest      : {manifest_path}",
            f"selected pose : {selected_pose.pose_path}",
            f"score         : {selected_pose.score}",
        ],
        heavy=False,
    )
    return manifest


def run_md_stage(args: argparse.Namespace, run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Run MD stage from saved docking outputs."""
    docking_manifest = _load_saved_docking_manifest(run_dir, state)
    docked_pose = _require_existing_file(
        docking_manifest.get("selected_pose"),
        "selected docked pose",
    )
    protein_pdb = _require_existing_file(args.inpdb, "inpdb")

    _print_spacer(1)
    _print_box(
        "MD inputs",
        [
            f"protein pdb   : {protein_pdb}",
            f"docked pose   : {docked_pose}",
            f"work dir      : {run_dir / 'md'}",
        ],
        heavy=False,
    )

    md_work_dir = run_dir / "md"
    result = run_md_from_docking(
        protein_pdb=protein_pdb,
        docked_ligand_pose=docked_pose,
        work_dir=md_work_dir,
        ligand_net_charge=args.ligand_net_charge,
        meeko_export_binary="mk_export.py",
        solvate=True,
        save_manifest=True,
        save_parametrized_pdb=True,
        save_solvated_pdb=True,
        save_md_final_pdb=True,
        save_md_native_outputs=True,
        md_output_stem="md_final",
    )

    if result.manifest_path is not None and result.manifest_path.exists():
        manifest = _read_json(result.manifest_path)
    else:
        manifest = result.to_dict()

    state["md"] = manifest
    _save_state(run_dir, state)

    _print_spacer(1)
    _print_box(
        "MD outputs",
        [
            f"manifest      : {result.manifest_path}",
            f"final pdb     : {result.final_md_pdb}",
            f"trajectory    : {result.native_md_outputs.trajectory_xtc}",
            f"run dir       : {result.native_md_outputs.run_directory}",
        ],
        heavy=False,
    )
    return manifest


def run_gbsa_stage(args: argparse.Namespace, run_dir: Path, state: dict[str, Any]) -> dict[str, Any]:
    """Run GBSA stage from saved MD outputs."""
    _print_stage_banner("GBSA")
    md_manifest = _load_saved_md_manifest(run_dir, state)

    native_outputs = md_manifest["md"]["native_outputs"]
    solvation = md_manifest["solvation"]
    parametrization = md_manifest["parametrization"]

    complex_structure = (
        native_outputs.get("portable_run_input_tpr")
        or native_outputs.get("final_gro")
        or solvation.get("gro_file")
        or parametrization.get("gro_file")
    )
    trajectory = native_outputs.get("trajectory_xtc") or native_outputs.get("trajectory_trr")
    topology = solvation.get("top_file") or parametrization.get("top_file")
    gro_for_index = native_outputs.get("final_gro") or solvation.get("gro_file") or parametrization.get("gro_file")

    complex_structure_path = _require_existing_file(complex_structure, "GBSA complex structure")
    trajectory_path = _require_existing_file(trajectory, "GBSA trajectory")
    topology_path = _require_existing_file(topology, "GBSA topology")
    gro_for_index_path = _require_existing_file(gro_for_index, "GRO file for GBSA index generation")

    output_dir = run_dir / "gbsa"
    output_dir.mkdir(parents=True, exist_ok=True)

    index_file_path, ligand_resname = _write_auto_gbsa_index_from_gro(
        gro_file=gro_for_index_path,
        output_ndx=output_dir / "gbsa_auto_index.ndx",
    )
    receptor_group, ligand_group = _get_required_group_numbers(index_file_path)

    _print_spacer(1)
    _print_box(
        "GBSA inputs",
        [
            f"complex       : {complex_structure_path}",
            f"trajectory    : {trajectory_path}",
            f"topology      : {topology_path}",
            f"index file    : {index_file_path}",
            f"ligand resname: {ligand_resname}",
            f"groups        : RECEPTOR={receptor_group}, LIGAND={ligand_group}",
        ],
        heavy=False,
    )

    config = MMPBSAConfig(
        general=GeneralParams(
            startframe=args.startframe,
            endframe=args.endframe,
            interval=args.interval,
        )
    )

    request = GmxMmpbsaRequest(
        complex_structure=complex_structure_path,
        trajectory=trajectory_path,
        topology=topology_path,
        index_file=index_file_path,
        receptor_group=receptor_group,
        ligand_group=ligand_group,
        output_dir=output_dir,
        config=config,
        gmx_mmpbsa=args.gmx_mmpbsa,
        extra_args=tuple(args.gbsa_extra_args),
    )
    result = run_gmx_mmpbsa(request)

    manifest_path = output_dir / "result.json"
    manifest = _read_json(manifest_path)
    state["gbsa"] = manifest
    _save_state(run_dir, state)

    _print_spacer(1)
    _print_box(
        "GBSA outputs",
        [
            f"manifest      : {manifest_path}",
            f"stdout log    : {result.stdout_log}",
            f"stderr log    : {result.stderr_log}",
            f"return code   : {result.returncode}",
            f"ok            : {result.ok}",
        ],
        heavy=False,
    )

    if not result.ok:
        raise RuntimeError("gmx_MMPBSA failed. See logs in the gbsa directory.")

    return manifest


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Run docking -> md -> gbsa with minimal CLI input, resumable JSON state, and end-user logging.",
    )

    parser.add_argument("--inpdb", required=True, help="Protein PDB path.")
    parser.add_argument("--inpsdf", required=True, help="Ligand SDF path.")
    parser.add_argument(
        "--step",
        choices=["docking", "md", "gbsa"],
        default="docking",
        help="Stage to start from. Default: docking",
    )
    parser.add_argument(
        "--stop-after",
        choices=["docking", "md", "gbsa"],
        default="gbsa",
        help="Last stage to run. Default: gbsa",
    )
    parser.add_argument(
        "--run-dir",
        default=None,
        help="Optional manual run directory. If omitted, one is created automatically.",
    )

    parser.add_argument("--box-padding", type=float, default=5.0)
    parser.add_argument("--box-min-size", type=float, default=20.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--num-modes", type=int, default=9)
    parser.add_argument("--exhaustiveness", type=int, default=16)
    parser.add_argument("--energy-range", type=float, default=3.0)
    parser.add_argument("--vina-binary", default="vina")
    parser.add_argument("--obabel-binary", default="obabel")

    parser.add_argument("--ligand-net-charge", type=int, default=None)
    parser.add_argument("--gmx-mmpbsa", default="gmx_MMPBSA")
    parser.add_argument("--startframe", type=int, default=1)
    parser.add_argument("--endframe", type=int, default=9_999_999)
    parser.add_argument("--interval", type=int, default=1)
    parser.add_argument("--gbsa-extra-args", nargs="*", default=[])
    return parser


def _write_summary_file(run_dir: Path, state: dict[str, Any]) -> Path:
    """Write final summary text file."""
    summary_path = run_dir / LOGS_DIRNAME / SUMMARY_FILENAME
    lines = [
        _box_lines(
            "RUN SUMMARY",
            [
                f"run dir      : {run_dir}",
                f"state file   : {run_dir / STATE_FILENAME}",
                f"status       : {state['status']}",
                f"failed stage : {state.get('failed_stage') or '--'}",
            ],
            heavy=True,
            center_title=True,
            style="wide_header",
        ),
        "",
        _box_lines(
            "Stage overview",
            _summary_lines_from_state(state),
            heavy=False,
        ),
    ]
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return summary_path


def main() -> None:
    """Parse CLI arguments and run selected pipeline stages."""
    parser = build_parser()
    args = parser.parse_args()

    if STEP_ORDER[args.step] > STEP_ORDER[args.stop_after]:
        parser.error("--step cannot be later than --stop-after.")

    inpdb = Path(args.inpdb).expanduser().resolve()
    inpsdf = Path(args.inpsdf).expanduser().resolve()
    run_dir = Path(args.run_dir).expanduser().resolve() if args.run_dir else _default_run_dir(inpdb, inpsdf)
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / LOGS_DIRNAME).mkdir(parents=True, exist_ok=True)

    state = _load_state(run_dir)
    selected_stages = _selected_stage_names(args.step, args.stop_after)
    state["selected_stages"] = selected_stages
    state["status"] = "running"
    _save_state(run_dir, state)

    _install_docking_box_style()
    _install_md_box_style()

    _print_spacer(1)
    _print_box(
        "DOCKING → MD → GBSA PIPELINE",
        [
            f"timestamp    : {_now()}",
            f"protein pdb  : {inpdb}",
            f"ligand sdf   : {inpsdf}",
            f"run dir      : {run_dir}",
            f"start step   : {args.step}",
            f"stop after   : {args.stop_after}",
            f"logs dir     : {run_dir / LOGS_DIRNAME}",
        ],
        heavy=True,
        center_title=True,
        style="wide_header",
    )

    overall_start = time.perf_counter()

    try:
        if _stage_is_selected("docking", args.step, args.stop_after):
            _run_stage_wrapped(
                stage="docking",
                args=args,
                run_dir=run_dir,
                state=state,
                runner=run_docking_stage,
            )

        if _stage_is_selected("md", args.step, args.stop_after):
            _run_stage_wrapped(
                stage="md",
                args=args,
                run_dir=run_dir,
                state=state,
                runner=run_md_stage,
            )

        if _stage_is_selected("gbsa", args.step, args.stop_after):
            _run_stage_wrapped(
                stage="gbsa",
                args=args,
                run_dir=run_dir,
                state=state,
                runner=run_gbsa_stage,
            )

    except Exception:
        state["status"] = "failed"
        _save_state(run_dir, state)
        summary_path = _write_summary_file(run_dir, state)
        _print_spacer(2)
        _print_box(
            "RUN FAILED",
            [
                f"elapsed      : {_elapsed_text(time.perf_counter() - overall_start)}",
                f"state file   : {run_dir / STATE_FILENAME}",
                f"summary file : {summary_path}",
                f"pipeline log : {run_dir / LOGS_DIRNAME / 'pipeline.log'}",
            ],
            heavy=True,
            center_title=True,
            style="wide_header",
        )
        raise
    else:
        state["status"] = "success"
        _save_state(run_dir, state)
        summary_path = _write_summary_file(run_dir, state)
        _print_spacer(2)
        _print_box(
            "RUN FINISHED",
            [
                f"elapsed      : {_elapsed_text(time.perf_counter() - overall_start)}",
                f"state file   : {run_dir / STATE_FILENAME}",
                f"summary file : {summary_path}",
                f"pipeline log : {run_dir / LOGS_DIRNAME / 'pipeline.log'}",
            ],
            heavy=True,
            center_title=True,
            style="wide_header",
        )
        _print_spacer(1)
        _print_box(
            "Stage overview",
            _summary_lines_from_state(state),
            heavy=False,
            center_title=True,
            style="wide_header",
        )


if __name__ == "__main__":
    main()
