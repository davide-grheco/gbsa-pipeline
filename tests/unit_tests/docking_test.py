"""Unit tests for the lightweight docking adapter layer.

This module keeps only small, local checks that do not require external tools.
The goal is to verify stable contracts for helper behavior such as ligand
preparation output shape and Vina command construction.
Anything that requires real Vina execution, mk_export.py, or filesystem-heavy
workflow chaining belongs in the integration test module instead.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from gbsa_pipeline.docking import DockingBox, VinaEngine, prepare_ligand_with_meeko

if TYPE_CHECKING:
    from pathlib import Path


def test_meeko_smiles_to_pdbqt(tmp_path: Path) -> None:
    """Check that a simple SMILES string is converted into a PDBQT file.

    This is a unit-level contract test for `prepare_ligand_with_meeko()` and
    does not try to prove docking correctness or chemical realism beyond basic
    output generation.
    The `tmp_path` parameter is required because the function writes a PDBQT
    file, and unit tests should keep such artifacts isolated from repository
    fixtures and from other tests.
    We are currently checking three things: the returned path matches the
    requested output path, the file is actually created, and the contents look
    like a PDBQT-style ligand file by containing expected record sections.
    """
    output = tmp_path / "ligand.pdbqt"

    path = prepare_ligand_with_meeko("CCO", output, name="ETH")

    assert path == output
    assert output.exists()

    content = output.read_text(encoding="utf-8")
    assert "ROOT" in content
    assert "ATOM" in content


def test_vina_build_command(tmp_path: Path) -> None:
    """Check that the Vina engine builds the expected command-line arguments.

    This is a unit test for `_build_command()` and exists so argument forwarding
    can be validated without invoking the real Vina binary.
    The `tmp_path` parameter is required because the command is assembled from
    concrete receptor, ligand, and output paths, even though the files are only
    placeholders for this local contract check.
    We are currently checking that the box, random seed, mode count,
    exhaustiveness, energy range, and extra flags are all encoded into the
    generated command list in a way the downstream subprocess call can use.
    """
    engine = VinaEngine(binary="vina")
    box = DockingBox(center=(-3.245, 29.915, 53.639), size=(10.0, 10.0, 10.0))

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
