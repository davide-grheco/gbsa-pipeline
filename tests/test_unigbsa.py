# /home/grheco/repositorios/gbsa-pipeline/tests/test_unigbsa.py

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

from gbsa_pipeline.unigbsa import (
    GBParams,
    GmxMmpbsaRequest,
    MMPBSAConfig,
    PBParams,
    run_gmx_mmpbsa,
    run_gmx_mmpbsa_from_gromacs,
)

if TYPE_CHECKING:
    import pytest


def test_mmpbsa_config_write_creates_expected_input_file(tmp_path: Path) -> None:
    """MMPBSAConfig.write should render a valid gmx_MMPBSA input file."""
    config = MMPBSAConfig(
        gb=GBParams(igb=8, saltcon=0.15),
        pb=PBParams(ipb=1, indi=2.0),
    )

    input_path = tmp_path / "mmpbsa.in"
    written = config.write(input_path)

    assert written == input_path
    assert input_path.exists()

    text = input_path.read_text(encoding="utf-8")

    assert "&general" in text
    assert "&gb" in text
    assert "&pb" in text

    assert "igb" in text
    assert "8" in text

    assert "saltcon" in text
    assert "0.15" in text

    assert "ipb" in text
    assert "indi" in text
    assert "2.0" in text


def test_run_gmx_mmpbsa_writes_logs_and_result_json(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """run_gmx_mmpbsa should write input/log/result artifacts and return a serializable result."""
    complex_structure = tmp_path / "complex.gro"
    trajectory = tmp_path / "traj.xtc"
    topology = tmp_path / "topol.top"
    index_file = tmp_path / "index.ndx"

    complex_structure.write_text("dummy complex\n", encoding="utf-8")
    trajectory.write_text("dummy trajectory\n", encoding="utf-8")
    topology.write_text("dummy topology\n", encoding="utf-8")
    index_file.write_text("dummy index\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run(
        cmd: list[str],
        cwd: str,
        text: bool,
        capture_output: bool,
        check: bool,
    ) -> subprocess.CompletedProcess[str]:
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        captured["text"] = text
        captured["capture_output"] = capture_output
        captured["check"] = check

        return subprocess.CompletedProcess(
            args=cmd,
            returncode=0,
            stdout="mock stdout\n",
            stderr="mock stderr\n",
        )

    monkeypatch.setattr("gbsa_pipeline.unigbsa.subprocess.run", fake_run)

    output_dir = tmp_path / "gbsa"

    request = GmxMmpbsaRequest(
        complex_structure=complex_structure,
        trajectory=trajectory,
        topology=topology,
        index_file=index_file,
        receptor_group=1,
        ligand_group=2,
        output_dir=output_dir,
    )

    result = run_gmx_mmpbsa(request)

    input_file = output_dir / "mmpbsa.in"
    stdout_log = output_dir / "gmx_mmpbsa.stdout.log"
    stderr_log = output_dir / "gmx_mmpbsa.stderr.log"
    result_json = output_dir / "result.json"

    assert result.ok is True
    assert result.returncode == 0

    assert input_file.exists()
    assert stdout_log.exists()
    assert stderr_log.exists()
    assert result_json.exists()

    assert stdout_log.read_text(encoding="utf-8") == "mock stdout\n"
    assert stderr_log.read_text(encoding="utf-8") == "mock stderr\n"

    saved = json.loads(result_json.read_text(encoding="utf-8"))

    assert saved["ok"] is True
    assert saved["returncode"] == 0
    assert saved["receptor_group"] == 1
    assert saved["ligand_group"] == 2

    assert saved["input_file"] == str(input_file)
    assert saved["stdout_log"] == str(stdout_log)
    assert saved["stderr_log"] == str(stderr_log)
    assert saved["result_json"] == str(result_json)

    assert saved["final_results_mmpbsa"] == str(output_dir / "FINAL_RESULTS_MMPBSA.dat")
    assert saved["final_results_decomp"] == str(output_dir / "FINAL_DECOMP_MMPBSA.dat")
    assert saved["info_file"] == str(output_dir / "gmx_MMPBSA_info.dat")

    assert captured["cwd"] == str(output_dir)
    assert captured["text"] is True
    assert captured["capture_output"] is True
    assert captured["check"] is False

    cmd = captured["cmd"]
    assert isinstance(cmd, list)
    assert cmd[0] == "gmx_MMPBSA"
    assert "-O" in cmd
    assert "-i" in cmd
    assert "-cs" in cmd
    assert "-ct" in cmd
    assert "-cp" in cmd
    assert "-ci" in cmd
    assert "-cg" in cmd

    cg_index = cmd.index("-cg")
    assert cmd[cg_index + 1] == "1"
    assert cmd[cg_index + 2] == "2"


def test_run_gmx_mmpbsa_passes_extra_args(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Extra CLI arguments should be appended to the gmx_MMPBSA command."""
    complex_structure = tmp_path / "complex.gro"
    trajectory = tmp_path / "traj.xtc"
    topology = tmp_path / "topol.top"
    index_file = tmp_path / "index.ndx"

    complex_structure.write_text("x\n", encoding="utf-8")
    trajectory.write_text("x\n", encoding="utf-8")
    topology.write_text("x\n", encoding="utf-8")
    index_file.write_text("x\n", encoding="utf-8")

    recorded_cmd: list[str] = []

    def fake_run(
        cmd: list[str],
        cwd: str,
        text: bool,
        capture_output: bool,
        check: bool,
    ) -> subprocess.CompletedProcess[str]:
        recorded_cmd[:] = cmd
        return subprocess.CompletedProcess(
            args=cmd,
            returncode=0,
            stdout="",
            stderr="",
        )

    monkeypatch.setattr("gbsa_pipeline.unigbsa.subprocess.run", fake_run)

    request = GmxMmpbsaRequest(
        complex_structure=complex_structure,
        trajectory=trajectory,
        topology=topology,
        index_file=index_file,
        receptor_group=3,
        ligand_group=7,
        output_dir=tmp_path / "gbsa",
        extra_args=("-nogui", "-eo", "energy.csv"),
    )

    run_gmx_mmpbsa(request)

    assert recorded_cmd[-3:] == ["-nogui", "-eo", "energy.csv"]


def test_run_gmx_mmpbsa_from_gromacs_backward_compatible(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The legacy wrapper should still call subprocess.run and return CompletedProcess."""
    input_file = tmp_path / "mmpbsa.in"
    complex_structure = tmp_path / "complex.gro"
    trajectory = tmp_path / "traj.xtc"
    topology = tmp_path / "topol.top"
    index_file = tmp_path / "index.ndx"
    output_dir = tmp_path / "legacy_run"

    for path in [input_file, complex_structure, trajectory, topology, index_file]:
        path.write_text("x\n", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_run(
        cmd: list[str],
        cwd: str,
        text: bool,
        capture_output: bool,
        check: bool,
    ) -> subprocess.CompletedProcess[str]:
        captured["cmd"] = cmd
        captured["cwd"] = cwd
        return subprocess.CompletedProcess(
            args=cmd,
            returncode=5,
            stdout="legacy stdout",
            stderr="legacy stderr",
        )

    monkeypatch.setattr("gbsa_pipeline.unigbsa.subprocess.run", fake_run)

    completed = run_gmx_mmpbsa_from_gromacs(
        input_file=input_file,
        complex_structure=complex_structure,
        trajectory=trajectory,
        topology=topology,
        index_file=index_file,
        receptor_group=4,
        ligand_group=9,
        output_dir=output_dir,
        extra_args=("-nogui",),
    )

    assert isinstance(completed, subprocess.CompletedProcess)
    assert completed.returncode == 5
    assert completed.stdout == "legacy stdout"
    assert completed.stderr == "legacy stderr"

    assert output_dir.exists()
    assert captured["cwd"] == str(output_dir)

    cmd = captured["cmd"]
    assert isinstance(cmd, list)
    assert cmd[0] == "gmx_MMPBSA"
    assert cmd[-1] == "-nogui"


def test_run_gmx_mmpbsa_records_failure_in_result_json(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Non-zero exit codes should be reflected in the returned result and result.json."""
    complex_structure = tmp_path / "complex.gro"
    trajectory = tmp_path / "traj.xtc"
    topology = tmp_path / "topol.top"
    index_file = tmp_path / "index.ndx"

    for path in [complex_structure, trajectory, topology, index_file]:
        path.write_text("x\n", encoding="utf-8")

    def fake_run(
        cmd: list[str],
        cwd: str,
        text: bool,
        capture_output: bool,
        check: bool,
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args=cmd,
            returncode=17,
            stdout="",
            stderr="something failed\n",
        )

    monkeypatch.setattr("gbsa_pipeline.unigbsa.subprocess.run", fake_run)

    result = run_gmx_mmpbsa(
        GmxMmpbsaRequest(
            complex_structure=complex_structure,
            trajectory=trajectory,
            topology=topology,
            index_file=index_file,
            receptor_group=0,
            ligand_group=1,
            output_dir=tmp_path / "gbsa_fail",
        )
    )

    assert result.ok is False
    assert result.returncode == 17

    saved = json.loads(Path(result.result_json).read_text(encoding="utf-8"))
    assert saved["ok"] is False
    assert saved["returncode"] == 17
