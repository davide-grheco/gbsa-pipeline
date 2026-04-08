# /home/grheco/repositorios/gbsa-pipeline/src/gbsa_pipeline/md.py

"""Combined MD stage for gbsa-pipeline.

Overview
--------
This module provides one importable entry point for the MD stage so external
workflow tools only need to resolve file paths and parameters, then call this
module.

Current workflow
----------------
1. clean the protein PDB for parametrization
   - remove crystallographic waters (fresh solvation is done later)
   - map common AMBER-style terminal residue names back to standard names
2. prepare a ligand SDF suitable for parametrization
   - if ligand_pose is already an SDF, copy it
   - if ligand_pose is a PDBQT or DLG/DLG.GZ, export it to SDF with Meeko's
     Python API
3. parameterize protein + ligand via gbsa_pipeline.parametrization
4. solvate with OpenMM/ParmEd via gbsa_pipeline.solvation_openmm
5. load the solvated GROMACS system into BioSimSpace
6. run a first short GROMACS MD process
7. persist both:
   - user-friendly structural snapshots (PDB)
   - native GROMACS run artefacts (trajectory, topology/run files, logs, energy)
8. save a machine-readable manifest

Important design choices
------------------------
- Docking-derived PDBQT is treated as a pose container, not as a chemistry source.
- No Meeko CLI tool is used. Ligand export uses Meeko's Python API directly.
- Input crystal waters are removed before parametrization because this workflow
  re-solvates the system explicitly.
- Solvation is enabled by default so the first MD run has a neutralized/solvated box.
- For downstream GBSA, the module now preserves native MD trajectory artefacts
  instead of exposing only the last-frame PDB snapshot.

What this module does NOT do
----------------------------
- no signac logic
- no grubicy logic
- no parent/child job resolution
- no CLI parsing
"""

from __future__ import annotations

import contextlib
import gzip
import json
import logging
import shutil
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import BioSimSpace as BSS
from meeko import PDBQTMolecule, RDKitMolCreate

from gbsa_pipeline.change_defaults import GromacsCustom, GromacsParams
from gbsa_pipeline.parametrization import (
    ParametrisedComplex,
    ParametrizationConfig,
    ParametrizationInput,
    parametrize,
)
from gbsa_pipeline.solvation_box import SolvationParams
from gbsa_pipeline.solvation_openmm import SolvatedComplex, solvate_openmm

LOGGER = logging.getLogger(__name__)

PDB_RESNAME_END_COL = 21

_AMBER_TO_STANDARD_RESNAME: dict[str, str] = {
    "NASP": "ASP",
    "CASP": "ASP",
    "NGLU": "GLU",
    "CGLU": "GLU",
    "NHID": "HIS",
    "NHIE": "HIS",
    "NHIP": "HIS",
    "CHID": "HIS",
    "CHIE": "HIS",
    "CHIP": "HIS",
    "NCYS": "CYS",
    "CCYS": "CYS",
    "CYM": "CYS",
    "CYX": "CYS",
    "NLYS": "LYS",
    "CLYS": "LYS",
    "NARG": "ARG",
    "CARG": "ARG",
}


@dataclass(frozen=True)
class MDInput:
    """Inputs for the combined docking-pose -> parametrization -> solvation -> MD stage."""

    protein_pdb: Path
    ligand_pose: Path
    work_dir: Path

    parametrization_config: ParametrizationConfig = field(default_factory=ParametrizationConfig)
    gromacs_params: GromacsParams = field(default_factory=GromacsParams)
    solvation_params: SolvationParams = field(default_factory=SolvationParams)

    ligand_net_charge: int | None = None

    # kept for API compatibility with older callers; ignored intentionally
    meeko_export_binary: str = "mk_export.py"

    solvate: bool = True
    save_manifest: bool = True
    save_parametrized_pdb: bool = True
    save_solvated_pdb: bool = True
    save_md_final_pdb: bool = True
    save_md_native_outputs: bool = True
    md_output_stem: str = "md_final"


@dataclass(frozen=True)
class NativeMDOutputs:
    """Stable copy of native GROMACS run outputs.

    These are the files that downstream analysis/GBSA workflows typically need.
    Any field may be None if the corresponding artefact was not produced or
    could not be located in the BioSimSpace/GROMACS work directory.
    """

    run_directory: Path | None
    trajectory_xtc: Path | None = None
    trajectory_trr: Path | None = None
    final_gro: Path | None = None
    portable_run_input_tpr: Path | None = None
    checkpoint_cpt: Path | None = None
    energy_edr: Path | None = None
    md_log: Path | None = None
    index_ndx: Path | None = None
    mdp: Path | None = None

    def to_dict(self) -> dict[str, str | None]:
        """Return a JSON-serializable mapping of native MD output paths."""
        return {
            "run_directory": str(self.run_directory) if self.run_directory else None,
            "trajectory_xtc": str(self.trajectory_xtc) if self.trajectory_xtc else None,
            "trajectory_trr": str(self.trajectory_trr) if self.trajectory_trr else None,
            "final_gro": str(self.final_gro) if self.final_gro else None,
            "portable_run_input_tpr": str(self.portable_run_input_tpr) if self.portable_run_input_tpr else None,
            "checkpoint_cpt": str(self.checkpoint_cpt) if self.checkpoint_cpt else None,
            "energy_edr": str(self.energy_edr) if self.energy_edr else None,
            "md_log": str(self.md_log) if self.md_log else None,
            "index_ndx": str(self.index_ndx) if self.index_ndx else None,
            "mdp": str(self.mdp) if self.mdp else None,
        }


@dataclass(frozen=True)
class MDResult:
    """Structured result for the combined MD stage."""

    work_dir: Path
    protein_dir: Path
    ligand_dir: Path
    parametrization_dir: Path
    solvation_dir: Path
    md_dir: Path

    prepared_protein_pdb: Path
    normalized_ligand_sdf: Path
    parametrized_complex: ParametrisedComplex
    solvated_complex: SolvatedComplex | None

    parametrized_pdb: Path | None
    solvated_pdb: Path | None
    final_md_pdb: Path | None
    native_md_outputs: NativeMDOutputs
    manifest_path: Path | None

    def to_dict(self) -> dict[str, Any]:
        """Convert the successful MD result into a JSON-serializable dictionary."""
        pickle_path = self.parametrization_dir / "complex.pickle"
        return {
            "status": "success",
            "work_dir": str(self.work_dir),
            "protein": {
                "prepared_pdb": str(self.prepared_protein_pdb),
            },
            "ligand": {
                "normalized_sdf": str(self.normalized_ligand_sdf),
            },
            "parametrization": {
                "directory": str(self.parametrization_dir),
                "gro_file": str(self.parametrized_complex.gro_file),
                "top_file": str(self.parametrized_complex.top_file),
                "pickle_file": str(pickle_path) if pickle_path.exists() else None,
                "parametrized_pdb": str(self.parametrized_pdb) if self.parametrized_pdb is not None else None,
                "config": _safe_model_dump(self.parametrized_complex.config),
            },
            "solvation": {
                "enabled": self.solvated_complex is not None,
                "directory": str(self.solvation_dir),
                "gro_file": str(self.solvated_complex.gro_file) if self.solvated_complex is not None else None,
                "top_file": str(self.solvated_complex.top_file) if self.solvated_complex is not None else None,
                "solvated_pdb": str(self.solvated_pdb) if self.solvated_pdb is not None else None,
            },
            "md": {
                "directory": str(self.md_dir),
                "final_pdb": str(self.final_md_pdb) if self.final_md_pdb else None,
                "native_outputs": self.native_md_outputs.to_dict(),
            },
            "manifest_path": str(self.manifest_path) if self.manifest_path else None,
        }


def _safe_model_dump(obj: Any) -> dict[str, Any]:
    """Best-effort model dump for pydantic/dataclass-like config objects."""
    if hasattr(obj, "model_dump"):
        return obj.model_dump(mode="json")
    if hasattr(obj, "__dict__"):
        return dict(obj.__dict__)
    return {"value": str(obj)}


def _write_json(path: Path, data: dict[str, Any]) -> None:
    """Write JSON data to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def _write_text_log(path: Path, title: str, body: str) -> Path:
    """Write a plain-text debug log."""
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    text = f"{title}\n{'=' * len(title)}\n\n{body}"
    path.write_text(text, encoding="utf-8")
    return path


def _validate_existing_file(path: Path, label: str) -> Path:
    """Validate that a path exists and is a file."""
    path = Path(path).expanduser().resolve()

    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")

    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")

    return path


def _read_text_maybe_gzip(path: Path) -> str:
    """Read a plain-text or .gz-compressed file."""
    path = Path(path).resolve()
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return handle.read()
    return path.read_text(encoding="utf-8")


def _detect_is_dlg(path: Path) -> bool:
    """Return True if the file should be treated as an AutoDock-GPU DLG result."""
    name = Path(path).name.lower()
    return name.endswith((".dlg", ".dlg.gz"))


def _prepare_protein_pdb_for_parametrization(
    protein_pdb: Path,
    protein_dir: Path,
) -> tuple[Path, Path]:
    """Write a temporary protein PDB suitable for the OpenMM parametrization path.

    Current behavior
    ----------------
    - removes waters (HOH/WAT), because this workflow re-solvates explicitly
    - maps common AMBER-style terminal residue names back to standard names
    """
    protein_pdb = _validate_existing_file(protein_pdb, "protein_pdb")
    protein_dir.mkdir(parents=True, exist_ok=True)

    output_pdb = protein_dir / "protein_for_parametrization.pdb"
    log_path = protein_dir / "protein_preparation.log"

    n_removed_waters = 0
    n_mapped_names = 0
    log_lines: list[str] = [
        f"Input protein PDB: {protein_pdb}",
        f"Output PDB: {output_pdb}",
    ]

    fixed_lines: list[str] = []

    for line in protein_pdb.read_text(encoding="utf-8").splitlines():
        updated_line = line

        if updated_line.startswith(("ATOM  ", "HETATM")) and len(updated_line) >= PDB_RESNAME_END_COL:
            resname = updated_line[17:21].strip()

            if resname in {"HOH", "WAT"}:
                n_removed_waters += 1
                continue

            mapped = _AMBER_TO_STANDARD_RESNAME.get(resname, resname)
            if mapped != resname:
                n_mapped_names += 1
                log_lines.append(f"Mapped residue name {resname} -> {mapped}")
                updated_line = f"{updated_line[:17]}{mapped:>3s} {updated_line[21:]}"

        fixed_lines.append(updated_line)

    output_pdb.write_text("\n".join(fixed_lines) + "\n", encoding="utf-8")

    log_lines.append(f"Removed water records: {n_removed_waters}")
    log_lines.append(f"Mapped residue-name records: {n_mapped_names}")

    _write_text_log(
        log_path,
        "Protein preparation log",
        "\n".join(log_lines),
    )

    return output_pdb, log_path


def _export_pose_pdbqt_to_sdf_with_meeko(
    *,
    ligand_pose: Path,
    ligand_dir: Path,
) -> tuple[Path, Path]:
    """Export a docked PDBQT pose to SDF using Meeko's Python API.

    This follows the same core logic as Meeko's own export implementation.
    """
    ligand_pose = _validate_existing_file(ligand_pose, "ligand_pose")
    ligand_dir.mkdir(parents=True, exist_ok=True)

    output_sdf = ligand_dir / "ligand_for_parametrization.sdf"
    log_path = ligand_dir / "ligand_pose_to_sdf.meeko.log"

    is_dlg = _detect_is_dlg(ligand_pose)
    input_string = _read_text_maybe_gzip(ligand_pose)

    LOGGER.info(
        "[ligand] exporting %s -> %s with Meeko Python API",
        ligand_pose.name,
        output_sdf.name,
    )

    pdbqt_mol = PDBQTMolecule(
        input_string,
        is_dlg=is_dlg,
        skip_typing=True,
    )

    output_string, failures = RDKitMolCreate.write_sd_string(
        pdbqt_mol,
        only_cluster_leads=False,
    )

    log_lines = [
        f"Input file: {ligand_pose}",
        f"Output file: {output_sdf}",
        f"is_dlg: {is_dlg}",
        f"Number of failures: {len(failures)}",
    ]

    for idx in failures:
        warnings.warn(
            f"molecule {idx} not converted to RDKit/SD File",
            stacklevel=2,
        )
        log_lines.append(f"Failed molecule index: {idx}")

    mol_index_annotations = pdbqt_mol._atom_annotations.get("mol_index", [])
    if len(failures) == len(mol_index_annotations):
        msg = "\nCould not convert to RDKit. Maybe meeko was not used for preparing\n"
        msg += "the input PDBQT for docking, and the SMILES string is missing?\n"
        msg += "Except for standard protein sidechains, all ligands and flexible residues\n"
        msg += "require a REMARK SMILES line in the PDBQT, which is added automatically by Meeko."

        log_lines.append("")
        log_lines.append("Conversion failed completely.")
        log_lines.append(msg)

        _write_text_log(
            log_path,
            f"Meeko export log for {ligand_pose.name}",
            "\n".join(log_lines),
        )
        raise RuntimeError(msg)

    if not output_string.strip():
        log_lines.append("")
        log_lines.append("Output SDF string was empty.")
        _write_text_log(
            log_path,
            f"Meeko export log for {ligand_pose.name}",
            "\n".join(log_lines),
        )
        raise RuntimeError(f"Meeko conversion returned an empty SDF string for ligand pose {ligand_pose}")

    output_sdf.write_text(output_string, encoding="utf-8")

    log_lines.append("")
    log_lines.append("Conversion succeeded.")
    log_lines.append(f"Wrote SDF: {output_sdf}")

    _write_text_log(
        log_path,
        f"Meeko export log for {ligand_pose.name}",
        "\n".join(log_lines),
    )

    return output_sdf, log_path


def _prepare_ligand_for_parametrization(
    *,
    ligand_pose: Path,
    ligand_dir: Path,
) -> tuple[Path, list[Path]]:
    """Produce the ligand SDF used for parametrization.

    Rules
    -----
    - .sdf: copy directly
    - .pdbqt: export with Meeko Python API
    - .dlg / .dlg.gz: export with Meeko Python API
    - otherwise: fail
    """
    ligand_pose = _validate_existing_file(ligand_pose, "ligand_pose")
    ligand_dir.mkdir(parents=True, exist_ok=True)

    ligand_name = ligand_pose.name.lower()

    if ligand_pose.suffix.lower() == ".sdf":
        output_sdf = ligand_dir / "ligand_for_parametrization.sdf"
        shutil.copy2(ligand_pose, output_sdf)
        return output_sdf, []

    if ligand_pose.suffix.lower() == ".pdbqt" or ligand_name.endswith((".dlg", ".dlg.gz")):
        output_sdf, log_path = _export_pose_pdbqt_to_sdf_with_meeko(
            ligand_pose=ligand_pose,
            ligand_dir=ligand_dir,
        )
        return output_sdf, [log_path]

    raise ValueError(
        f"Unsupported ligand_pose type for MD parametrization: {ligand_pose}. Expected .sdf, .pdbqt, .dlg, or .dlg.gz."
    )


def _load_bss_system_from_gromacs(gro_file: Path, top_file: Path) -> Any:
    """Load a BioSimSpace system from a GROMACS .gro/.top pair."""
    return BSS.IO.readMolecules(
        files=[str(gro_file), str(top_file)],
        make_whole=True,
    )


def _save_system_as_pdb(system: Any, output_path: Path) -> Path:
    """Save a BioSimSpace system as PDB."""
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    BSS.IO.saveMolecules(str(output_path), system, fileformat="PDB")
    return output_path


def _copy_if_exists(src: Path | None, dst: Path) -> Path | None:
    """Copy a file if it exists; otherwise return None."""
    if src is None:
        return None

    src = Path(src).expanduser().resolve()
    if not src.exists() or not src.is_file():
        return None

    dst = Path(dst).expanduser().resolve()
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def _find_first_existing_file(directory: Path, patterns: list[str]) -> Path | None:
    """Return the first matching file found for the provided glob patterns."""
    directory = Path(directory).expanduser().resolve()
    if not directory.exists():
        return None

    for pattern in patterns:
        matches = sorted(path for path in directory.glob(pattern) if path.is_file())
        if matches:
            return matches[0]

    return None


def _get_process_work_dir(proc: Any) -> Path | None:
    """Best-effort access to the BioSimSpace process work directory.

    Different BioSimSpace versions expose this slightly differently, so this
    helper is intentionally defensive.
    """
    candidates: list[Any] = []

    if hasattr(proc, "workDir"):
        with contextlib.suppress(Exception):
            candidates.append(proc.workDir())

    if hasattr(proc, "_work_dir"):
        with contextlib.suppress(Exception):
            candidates.append(proc._work_dir)

    for candidate in candidates:
        if candidate is None:
            continue

        try:
            path = Path(str(candidate)).expanduser().resolve()
        except (OSError, RuntimeError, TypeError, ValueError):
            continue

        if path.exists() and path.is_dir():
            return path

    return None


def _collect_native_md_outputs(
    *,
    proc: Any,
    md_dir: Path,
) -> NativeMDOutputs:
    """Collect native GROMACS artefacts from the BioSimSpace process work directory.

    Why this exists
    ---------------
    For GBSA or later trajectory analysis, a final-frame PDB is not enough.
    We therefore locate the actual GROMACS outputs generated by the process and
    copy stable copies into `work_dir/md/`.

    Notes:
    -----
    - The exact filenames can vary between BioSimSpace/GROMACS versions.
    - The function therefore searches by extension/pattern instead of assuming
      one fixed basename.
    - Missing files are tolerated and returned as None.
    """
    md_dir = Path(md_dir).expanduser().resolve()
    md_dir.mkdir(parents=True, exist_ok=True)

    run_dir = _get_process_work_dir(proc)
    if run_dir is None:
        LOGGER.warning("Could not determine BioSimSpace/GROMACS work directory. Native MD outputs will be unavailable.")
        return NativeMDOutputs(run_directory=None)

    archive_dir = md_dir / "native_outputs"
    archive_dir.mkdir(parents=True, exist_ok=True)

    xtc_src = _find_first_existing_file(run_dir, ["*.xtc"])
    trr_src = _find_first_existing_file(run_dir, ["*.trr"])
    gro_src = _find_first_existing_file(run_dir, ["*.gro"])
    tpr_src = _find_first_existing_file(run_dir, ["*.tpr"])
    cpt_src = _find_first_existing_file(run_dir, ["*.cpt"])
    edr_src = _find_first_existing_file(run_dir, ["*.edr"])
    log_src = _find_first_existing_file(run_dir, ["*.log"])
    ndx_src = _find_first_existing_file(run_dir, ["*.ndx"])
    mdp_src = _find_first_existing_file(run_dir, ["*.mdp"])

    outputs = NativeMDOutputs(
        run_directory=run_dir,
        trajectory_xtc=_copy_if_exists(xtc_src, archive_dir / "trajectory.xtc"),
        trajectory_trr=_copy_if_exists(trr_src, archive_dir / "trajectory.trr"),
        final_gro=_copy_if_exists(gro_src, archive_dir / "md_final.gro"),
        portable_run_input_tpr=_copy_if_exists(tpr_src, archive_dir / "md.tpr"),
        checkpoint_cpt=_copy_if_exists(cpt_src, archive_dir / "md.cpt"),
        energy_edr=_copy_if_exists(edr_src, archive_dir / "md.edr"),
        md_log=_copy_if_exists(log_src, archive_dir / "md.log"),
        index_ndx=_copy_if_exists(ndx_src, archive_dir / "index.ndx"),
        mdp=_copy_if_exists(mdp_src, archive_dir / "md.mdp"),
    )

    LOGGER.info("Native GROMACS work dir : %s", run_dir)
    if outputs.trajectory_xtc:
        LOGGER.info("Trajectory XTC         : %s", outputs.trajectory_xtc)
    if outputs.trajectory_trr:
        LOGGER.info("Trajectory TRR         : %s", outputs.trajectory_trr)
    if outputs.portable_run_input_tpr:
        LOGGER.info("Run input TPR          : %s", outputs.portable_run_input_tpr)
    if outputs.energy_edr:
        LOGGER.info("Energy EDR             : %s", outputs.energy_edr)
    if outputs.md_log:
        LOGGER.info("MD log                 : %s", outputs.md_log)

    return outputs


def _build_failure_payload(
    inp: MDInput,
    exc: Exception,
    manifest_path: Path,
) -> dict[str, Any]:
    """Build a machine-readable failure manifest."""
    return {
        "status": "failed",
        "inputs": {
            "protein_pdb": str(inp.protein_pdb),
            "ligand_pose": str(inp.ligand_pose),
            "work_dir": str(inp.work_dir),
            "ligand_net_charge": inp.ligand_net_charge,
            "solvate": inp.solvate,
            "md_output_stem": inp.md_output_stem,
            "save_md_final_pdb": inp.save_md_final_pdb,
            "save_md_native_outputs": inp.save_md_native_outputs,
        },
        "parametrization_config": _safe_model_dump(inp.parametrization_config),
        "solvation_params": _safe_model_dump(inp.solvation_params),
        "gromacs_params": _safe_model_dump(inp.gromacs_params),
        "error": {
            "type": type(exc).__name__,
            "message": str(exc),
        },
        "manifest_path": str(manifest_path),
    }


def run_md(inp: MDInput) -> MDResult:
    """Run the combined MD stage.

    This is the main importable entry point for workflow wrappers.
    """
    protein_pdb = _validate_existing_file(inp.protein_pdb, "protein_pdb")
    ligand_pose = _validate_existing_file(inp.ligand_pose, "ligand_pose")
    work_dir = Path(inp.work_dir).expanduser().resolve()

    protein_dir = work_dir / "protein"
    ligand_dir = work_dir / "ligand"
    parametrization_dir = work_dir / "parametrization"
    solvation_dir = work_dir / "solvation"
    md_dir = work_dir / "md"
    manifest_path = md_dir / "result.json"

    protein_dir.mkdir(parents=True, exist_ok=True)
    ligand_dir.mkdir(parents=True, exist_ok=True)
    parametrization_dir.mkdir(parents=True, exist_ok=True)
    solvation_dir.mkdir(parents=True, exist_ok=True)
    md_dir.mkdir(parents=True, exist_ok=True)

    LOGGER.info("====================================================")
    LOGGER.info("MD")
    LOGGER.info("====================================================")
    LOGGER.info("Protein PDB : %s", protein_pdb)
    LOGGER.info("Ligand pose : %s", ligand_pose)
    LOGGER.info("Work dir    : %s", work_dir)

    try:
        LOGGER.info("[1/5] Preparing protein and ligand for parametrization")
        prepared_protein_pdb, protein_log = _prepare_protein_pdb_for_parametrization(
            protein_pdb,
            protein_dir,
        )
        ligand_sdf, ligand_logs = _prepare_ligand_for_parametrization(
            ligand_pose=ligand_pose,
            ligand_dir=ligand_dir,
        )

        LOGGER.info("Protein PDB : %s", prepared_protein_pdb)
        LOGGER.info("Protein log : %s", protein_log)
        LOGGER.info("Ligand SDF  : %s", ligand_sdf)
        for log_path in ligand_logs:
            LOGGER.info("Ligand log  : %s", log_path)

        LOGGER.info("[2/5] Parametrizing system")
        complex_ = parametrize(
            ParametrizationInput(
                protein_pdb=prepared_protein_pdb,
                ligand_sdf=ligand_sdf,
                config=inp.parametrization_config,
                net_charge=inp.ligand_net_charge,
                work_dir=parametrization_dir,
            )
        )
        LOGGER.info("GRO         : %s", complex_.gro_file)
        LOGGER.info("TOP         : %s", complex_.top_file)

        LOGGER.info("[3/5] Preparing PDB snapshots")
        parametrized_pdb: Path | None = None
        unsolvated_bss_system = _load_bss_system_from_gromacs(
            gro_file=complex_.gro_file,
            top_file=complex_.top_file,
        )
        if inp.save_parametrized_pdb:
            parametrized_pdb = _save_system_as_pdb(
                unsolvated_bss_system,
                parametrization_dir / "complex_parametrized.pdb",
            )
            LOGGER.info("Param PDB   : %s", parametrized_pdb)

        solvated_complex: SolvatedComplex | None = None
        solvated_pdb: Path | None = None

        LOGGER.info("[4/5] Solvation")
        if inp.solvate:
            solvated_complex = solvate_openmm(
                parametrized=complex_,
                params=inp.solvation_params,
                output_gro=solvation_dir / "solvated.gro",
                output_top=solvation_dir / "solvated.top",
            )
            LOGGER.info("Solv GRO    : %s", solvated_complex.gro_file)
            LOGGER.info("Solv TOP    : %s", solvated_complex.top_file)

            bss_system = solvated_complex.load_bss()
            if inp.save_solvated_pdb:
                solvated_pdb = _save_system_as_pdb(
                    bss_system,
                    solvation_dir / "solvated.pdb",
                )
                LOGGER.info("Solv PDB    : %s", solvated_pdb)
        else:
            LOGGER.info("Solvation   : skipped")
            bss_system = unsolvated_bss_system

        LOGGER.info("[5/5] Running first GROMACS MD")
        protocol = GromacsCustom(params=inp.gromacs_params)
        proc = BSS.Process.Gromacs(bss_system, protocol=protocol)
        proc.start()
        proc.wait()
        LOGGER.info("GROMACS run finished")

        final_md_pdb: Path | None = None
        final_system = proc.getSystem(block=True)

        if inp.save_md_final_pdb:
            final_md_pdb = _save_system_as_pdb(
                final_system,
                md_dir / f"{inp.md_output_stem}.pdb",
            )
            LOGGER.info("Final PDB   : %s", final_md_pdb)
        else:
            LOGGER.info("Final PDB   : skipped")

        if inp.save_md_native_outputs:
            native_md_outputs = _collect_native_md_outputs(
                proc=proc,
                md_dir=md_dir,
            )
        else:
            LOGGER.info("Native MD outputs : skipped")
            native_md_outputs = NativeMDOutputs(run_directory=None)

        result = MDResult(
            work_dir=work_dir,
            protein_dir=protein_dir,
            ligand_dir=ligand_dir,
            parametrization_dir=parametrization_dir,
            solvation_dir=solvation_dir,
            md_dir=md_dir,
            prepared_protein_pdb=prepared_protein_pdb,
            normalized_ligand_sdf=ligand_sdf,
            parametrized_complex=complex_,
            solvated_complex=solvated_complex,
            parametrized_pdb=parametrized_pdb,
            solvated_pdb=solvated_pdb,
            final_md_pdb=final_md_pdb,
            native_md_outputs=native_md_outputs,
            manifest_path=manifest_path if inp.save_manifest else None,
        )

        if inp.save_manifest:
            _write_json(manifest_path, result.to_dict())
            LOGGER.info("Manifest    : %s", manifest_path)

    except Exception as exc:
        LOGGER.exception("MD stage failed")
        if inp.save_manifest:
            _write_json(manifest_path, _build_failure_payload(inp, exc, manifest_path))
        raise
    else:
        LOGGER.info("MD stage finished successfully")
        return result


def run_md_from_docking(
    *,
    protein_pdb: str | Path,
    docked_ligand_pose: str | Path,
    work_dir: str | Path,
    parametrization_config: ParametrizationConfig | None = None,
    gromacs_params: GromacsParams | None = None,
    solvation_params: SolvationParams | None = None,
    ligand_net_charge: int | None = None,
    meeko_export_binary: str = "mk_export.py",
    solvate: bool = True,
    save_manifest: bool = True,
    save_parametrized_pdb: bool = True,
    save_solvated_pdb: bool = True,
    save_md_final_pdb: bool = True,
    save_md_native_outputs: bool = True,
    md_output_stem: str = "md_final",
) -> MDResult:
    """Run the preferred docking-pose to MD workflow.

    - dock with Vina -> PDBQT or DLG
    - export pose with Meeko -> SDF
    - parametrize that SDF
    - solvate
    - run first MD
    - keep both final snapshot and native trajectory outputs.
    """
    inp = MDInput(
        protein_pdb=Path(protein_pdb),
        ligand_pose=Path(docked_ligand_pose),
        work_dir=Path(work_dir),
        parametrization_config=parametrization_config or ParametrizationConfig(),
        gromacs_params=gromacs_params or GromacsParams(),
        solvation_params=solvation_params or SolvationParams(),
        ligand_net_charge=ligand_net_charge,
        meeko_export_binary=meeko_export_binary,
        solvate=solvate,
        save_manifest=save_manifest,
        save_parametrized_pdb=save_parametrized_pdb,
        save_solvated_pdb=save_solvated_pdb,
        save_md_final_pdb=save_md_final_pdb,
        save_md_native_outputs=save_md_native_outputs,
        md_output_stem=md_output_stem,
    )
    return run_md(inp)


__all__ = [
    "MDInput",
    "MDResult",
    "NativeMDOutputs",
    "run_md",
    "run_md_from_docking",
]
