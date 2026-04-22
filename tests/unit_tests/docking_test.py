"""Tests for the lightweight docking adapter layer."""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest

from gbsa_pipeline.docking import (
    DockingBox,
    DockingRequest,
    VinaEngine,
    prepare_ligand_with_meeko,
)

TESTDATA = Path(__file__).parents[1] / "testdata"
DOCKING_TESTDATA = TESTDATA / "docking"


def test_meeko_smiles_to_pdbqt(tmp_path: Path) -> None:
    output = tmp_path / "ligand.pdbqt"

    path = prepare_ligand_with_meeko("CCO", output, name="ETH")

    assert path == output
    assert output.exists()

    content = output.read_text(encoding="utf-8")
    assert "ROOT" in content
    assert "ATOM" in content


def test_vina_build_command(tmp_path: Path) -> None:
    engine = VinaEngine(binary="vina")
    box = DockingBox(center=(0.0, 1.0, 2.0), size=(10.0, 10.0, 10.0))

    receptor = tmp_path / "receptor.pdbqt"
    ligand = tmp_path / "ligand.pdbqt"
    output = tmp_path / "out.pdbqt"

    receptor.write_text("", encoding="utf-8")
    ligand.write_text("", encoding="utf-8")

    cmd = engine._build_command(
        receptor=receptor,
        ligand=ligand,
        output=output,
        box=box,
        seed=42,
        num_modes=5,
        exhaustiveness=3,
        energy_range=4.5,
        extra_flags={"--cpu": 2},
    )

    assert cmd[:2] == ["vina", "--receptor"]
    assert "--seed" in cmd
    assert "42" in cmd
    assert "--cpu" in cmd
    assert "2" in cmd
    assert "--num_modes" in cmd
    assert "5" in cmd
    assert "--exhaustiveness" in cmd
    assert "3" in cmd
    assert "--energy_range" in cmd
    assert "4.5" in cmd


@pytest.mark.integration
def test_vina_binary_smoke(tmp_path: Path) -> None:
    if shutil.which("vina") is None:
        pytest.skip("vina not available in PATH")

    receptor = TESTDATA / "protein1.pdbqt"
    if not receptor.exists():
        pytest.skip(f"missing receptor test file: {receptor}")

    engine = VinaEngine(binary="vina")
    box = DockingBox(center=(7.0, 20.0, 14.0), size=(10.0, 10.0, 10.0))

    ligand = tmp_path / "ligand.pdbqt"
    prepare_ligand_with_meeko("CCO", ligand, name="ETH")

    request = DockingRequest(
        receptor=receptor,
        ligands=[ligand],
        box=box,
        workdir=tmp_path,
    )

    result = engine.dock(request=request)

    assert result.engine == "vina"
    assert len(result.poses) == 1
    assert result.poses[0].pose_path.exists()
    assert result.poses[0].metadata["returncode"] == 0
    assert result.poses[0].score is not None
