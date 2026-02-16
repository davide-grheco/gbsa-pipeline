from __future__ import annotations

import subprocess
from collections.abc import Sequence
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Optional


def _fmt(v: Any) -> str:
    """Format values in a way gmx_MMPBSA input likes."""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, (int, float)):
        return str(v)
    if v is None:
        raise ValueError("None should be filtered before formatting")
    s = str(v)
    # Quote strings that could confuse parsing
    if s == "" or any(c.isspace() for c in s) or "," in s:
        return f'"{s}"'
    return s


def _as_kv(dc) -> dict[str, Any]:
    """Dataclass -> dict excluding None + merging extra."""
    d: dict[str, Any] = {}
    extra = {}
    for f in fields(dc):
        val = getattr(dc, f.name)
        if f.name == "extra":
            extra = val or {}
            continue
        if val is not None:
            d[f.name] = val
    # extra overrides explicit fields
    d.update(extra)
    return d


def _render_namelist(name: str, kv: dict[str, Any]) -> str:
    lines = [f"&{name}"]
    for k, v in kv.items():
        lines.append(f"  {k:<20} = {_fmt(v)}")
    lines.append("/")
    return "\n".join(lines)


# --------- Dataclasses reflecting your defaults ---------


@dataclass(frozen=True)
class GeneralParams:
    sys_name: str = ""
    startframe: int = 1
    endframe: int = 9_999_999
    interval: int = 1

    # Only relevant when gmx_MMPBSA builds AMBER topology itself (no -cp)
    forcefields: str | None = None
    ions_parameters: int | None = None
    PBRadii: int | None = None

    temperature: float = 298.15

    qh_entropy: int = 0
    interaction_entropy: int = 0
    ie_segment: int = 25
    c2_entropy: int = 0

    assign_chainID: int = 0
    exp_ki: float = 0.0

    full_traj: int = 0
    gmx_path: str = ""

    keep_files: int = 2
    netcdf: int = 0
    solvated_trajectory: int = 1
    verbose: int = 1

    # Any additional keys supported by your gmx_MMPBSA version
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class GBParams:
    igb: int = 5
    intdiel: float = 1.0
    extdiel: float = 78.5
    saltcon: float = 0.0

    surften: float = 0.0072
    surfoff: float = 0.0
    molsurf: int = 0
    msoffset: float = 0.0
    probe: float = 1.4

    # QM/MM options
    ifqnt: int = 0
    qm_theory: str = ""
    qm_residues: str = ""
    qmcharge_com: int = 0
    qmcharge_lig: int = 0
    qmcharge_rec: int = 0
    qmcut: float = 9999.0
    scfconv: float = 1e-08
    peptide_corr: int = 0
    writepdb: int = 1
    verbosity: int = 0

    alpb: int = 0
    arad_method: int = 1

    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PBParams:
    ipb: int = 2
    inp: int = 2
    sander_apbs: int = 0

    indi: float = 1.0
    exdi: float = 80.0
    emem: float = 4.0

    smoothopt: int = 1
    istrng: float = 0.0
    radiopt: int = 1
    prbrad: float = 1.4
    iprob: float = 2.0
    sasopt: int = 0
    arcres: float = 0.25

    memopt: int = 0
    mprob: float = 2.7
    mthick: float = 40.0
    mctrdz: float = 0.0
    poretype: int = 1

    npbopt: int = 0
    solvopt: int = 1
    accept: float = 0.001
    linit: int = 1000

    fillratio: float = 4.0
    scale: float = 2.0
    nbuffer: float = 0.0
    nfocus: int = 2
    fscale: int = 8
    npbgrid: int = 1
    bcopt: int = 5

    eneopt: int = 2
    frcopt: int = 0
    scalec: int = 0

    cutfd: float = 5.0
    cutnb: float = 0.0
    nsnba: int = 1

    decompopt: int = 2
    use_rmin: int = 1

    sprob: float = 0.557
    vprob: float = 1.3
    rhow_effect: float = 1.129
    use_sav: int = 1
    cavity_surften: float = 0.0378
    cavity_offset: float = -0.5692

    maxsph: int = 400
    maxarcdot: int = 1500
    npbverb: int = 0

    extra: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class MMPBSAConfig:
    general: GeneralParams = field(default_factory=GeneralParams)
    gb: Optional[GBParams] = field(default_factory=GBParams)
    pb: Optional[PBParams] = field(default_factory=PBParams)

    # If later you want to support rism/decomposition/nmode/etc without hardcoding:
    other_namelists: dict[str, dict[str, Any]] = field(default_factory=dict)

    def to_text(self) -> str:
        parts = []
        parts.append(_render_namelist("general", _as_kv(self.general)))
        if self.gb is not None:
            parts.append(_render_namelist("gb", _as_kv(self.gb)))
        if self.pb is not None:
            parts.append(_render_namelist("pb", _as_kv(self.pb)))
        for name, kv in self.other_namelists.items():
            parts.append(_render_namelist(name, kv))
        return "\n\n".join(parts).rstrip() + "\n"

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.write_text(self.to_text())
        return path


def run_gmx_mmpbsa_from_gromacs(
    input_file: str | Path,  # mmpbsa.in
    complex_structure: str | Path,  # complex.gro (or complex.pdb)
    trajectory: str | Path,  # traj.xtc
    topology: str | Path,  # topol.top
    index_file: str | Path,  # index.ndx
    receptor_group: int,  # group id in index_file (starts at 0)
    ligand_group: int,  # group id in index_file (starts at 0)
    output_dir: str | Path,
    gmx_mmpbsa: str = "gmx_MMPBSA",
    extra_args: Sequence[str] = (),
) -> subprocess.CompletedProcess[str]:
    """Run gmx_MMPBSA using GROMACS inputs (gro/top/xtc/ndx).

    This uses:
      -cs  complex_structure
      -ct  trajectory
      -cp  topology
      -ci  index_file
      -cg  receptor_group ligand_group
      -i   input_file
    """
    output_dir = Path(output_dir)
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
    ]
    cmd.extend(extra_args)

    return subprocess.run(
        cmd,
        cwd=str(output_dir),
        text=True,
        capture_output=True,
        check=False,
    )
