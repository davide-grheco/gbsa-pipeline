# /home/grheco/repositorios/gbsa-pipeline/tests/test_md.py

from __future__ import annotations

import json
from pathlib import Path  # noqa: TC003
from types import SimpleNamespace
from typing import TYPE_CHECKING, Any

import gbsa_pipeline.md as md_module

if TYPE_CHECKING:
    import pytest


class DummyParametrizedComplex:
    """Minimal parametrized-complex stub used by MD tests."""

    def __init__(self, gro_file: Path, top_file: Path, config: Any) -> None:
        """Store the parametrization outputs expected by the MD module."""
        self.gro_file = gro_file
        self.top_file = top_file
        self.config = config


class DummySolvatedComplex:
    """Minimal solvated-complex stub used by MD tests."""

    def __init__(self, gro_file: Path, top_file: Path) -> None:
        """Store the solvated GROMACS outputs expected by the MD module."""
        self.gro_file = gro_file
        self.top_file = top_file

    def load_bss(self) -> object:
        """Return a BioSimSpace-like solvated system placeholder."""
        return {"kind": "solvated_bss_system"}


class DummyProcess:
    """Minimal BioSimSpace process stub used by MD tests."""

    def __init__(self, system: object, protocol: object) -> None:
        """Store process inputs and lifecycle flags."""
        self.system = system
        self.protocol = protocol
        self.started = False
        self.waited = False

    def start(self) -> None:
        """Mark the process as started."""
        self.started = True

    def wait(self) -> None:
        """Mark the process as waited on."""
        self.waited = True

    def getSystem(self, block: bool = True) -> object:  # noqa: N802, FBT002
        """Return a BioSimSpace-like final system placeholder."""
        return {"kind": "final_bss_system", "block": block}


def test_run_md_from_docking_writes_manifest_and_native_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test the happy path of ``gbsa_pipeline.md.run_md_from_docking()``.

    This test does not run real Meeko/OpenMM/BioSimSpace/GROMACS work.
    Instead, it mocks the heavy steps and verifies that the MD module:

    - accepts the new API arguments
    - returns a structured result
    - writes the manifest
    - preserves native MD outputs needed for downstream GBSA
    """
    protein_pdb = tmp_path / "protein_input.pdb"
    ligand_pose = tmp_path / "ligand_pose.pdbqt"
    work_dir = tmp_path / "md_work"

    protein_pdb.write_text(
        ("ATOM      1  N   ALA A   1      11.104  13.207  12.011  1.00 20.00           N\nEND\n"),
        encoding="utf-8",
    )
    ligand_pose.write_text(
        "REMARK  test ligand pose\n",
        encoding="utf-8",
    )

    prepared_protein_pdb = work_dir / "protein" / "protein_for_parametrization.pdb"
    prepared_protein_log = work_dir / "protein" / "protein_preparation.log"
    normalized_ligand_sdf = work_dir / "ligand" / "ligand_for_parametrization.sdf"
    ligand_log = work_dir / "ligand" / "ligand_pose_to_sdf.meeko.log"

    parametrized_gro = work_dir / "parametrization" / "complex.gro"
    parametrized_top = work_dir / "parametrization" / "complex.top"

    final_md_pdb = work_dir / "md" / "md_final.pdb"
    manifest_path = work_dir / "md" / "result.json"

    native_outputs = md_module.NativeMDOutputs(
        run_directory=work_dir / "md" / "native_run_dir",
        trajectory_xtc=work_dir / "md" / "native_outputs" / "trajectory.xtc",
        trajectory_trr=work_dir / "md" / "native_outputs" / "trajectory.trr",
        final_gro=work_dir / "md" / "native_outputs" / "md_final.gro",
        portable_run_input_tpr=work_dir / "md" / "native_outputs" / "md.tpr",
        checkpoint_cpt=work_dir / "md" / "native_outputs" / "md.cpt",
        energy_edr=work_dir / "md" / "native_outputs" / "md.edr",
        md_log=work_dir / "md" / "native_outputs" / "md.log",
        index_ndx=work_dir / "md" / "native_outputs" / "index.ndx",
        mdp=work_dir / "md" / "native_outputs" / "md.mdp",
    )

    def fake_prepare_protein(
        protein_pdb_arg: Path,
        protein_dir_arg: Path,
    ) -> tuple[Path, Path]:
        assert protein_pdb_arg == protein_pdb.resolve()
        assert protein_dir_arg == (work_dir / "protein").resolve()
        prepared_protein_pdb.parent.mkdir(parents=True, exist_ok=True)
        prepared_protein_pdb.write_text("ATOM\nEND\n", encoding="utf-8")
        prepared_protein_log.write_text("protein prep ok\n", encoding="utf-8")
        return prepared_protein_pdb, prepared_protein_log

    def fake_prepare_ligand(
        *,
        ligand_pose: Path,
        ligand_dir: Path,
    ) -> tuple[Path, list[Path]]:
        assert ligand_pose == ligand_pose.resolve()
        assert ligand_dir == (work_dir / "ligand").resolve()
        normalized_ligand_sdf.parent.mkdir(parents=True, exist_ok=True)
        normalized_ligand_sdf.write_text("$$$$\n", encoding="utf-8")
        ligand_log.parent.mkdir(parents=True, exist_ok=True)
        ligand_log.write_text("ligand prep ok\n", encoding="utf-8")
        return normalized_ligand_sdf, [ligand_log]

    def fake_parametrize(inp: Any) -> DummyParametrizedComplex:
        parametrized_gro.parent.mkdir(parents=True, exist_ok=True)
        parametrized_gro.write_text("test gro\n", encoding="utf-8")
        parametrized_top.parent.mkdir(parents=True, exist_ok=True)
        parametrized_top.write_text("test top\n", encoding="utf-8")
        return DummyParametrizedComplex(
            gro_file=parametrized_gro,
            top_file=parametrized_top,
            config=inp.config,
        )

    def fake_load_bss_system_from_gromacs(gro_file: Path, top_file: Path) -> object:
        assert gro_file == parametrized_gro
        assert top_file == parametrized_top
        return {"kind": "unsolvated_bss_system"}

    def fake_save_system_as_pdb(system: object, output_path: Path) -> Path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(f"PDB for {system}\n", encoding="utf-8")
        return output_path

    def fake_solvate_openmm(
        *,
        parametrized: object,
        params: object,
        output_gro: Path,
        output_top: Path,
    ) -> DummySolvatedComplex:
        del parametrized, params
        output_gro.parent.mkdir(parents=True, exist_ok=True)
        output_gro.write_text("solvated gro\n", encoding="utf-8")
        output_top.parent.mkdir(parents=True, exist_ok=True)
        output_top.write_text("solvated top\n", encoding="utf-8")
        return DummySolvatedComplex(
            gro_file=output_gro,
            top_file=output_top,
        )

    def fake_collect_native_md_outputs(
        *,
        proc: object,
        md_dir: Path,
    ) -> md_module.NativeMDOutputs:
        del proc
        assert md_dir == (work_dir / "md").resolve()

        for path in [
            native_outputs.trajectory_xtc,
            native_outputs.trajectory_trr,
            native_outputs.final_gro,
            native_outputs.portable_run_input_tpr,
            native_outputs.checkpoint_cpt,
            native_outputs.energy_edr,
            native_outputs.md_log,
            native_outputs.index_ndx,
            native_outputs.mdp,
        ]:
            assert path is not None
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("native output\n", encoding="utf-8")

        assert native_outputs.run_directory is not None
        native_outputs.run_directory.mkdir(parents=True, exist_ok=True)
        return native_outputs

    monkeypatch.setattr(
        md_module,
        "_prepare_protein_pdb_for_parametrization",
        fake_prepare_protein,
    )
    monkeypatch.setattr(
        md_module,
        "_prepare_ligand_for_parametrization",
        fake_prepare_ligand,
    )
    monkeypatch.setattr(md_module, "parametrize", fake_parametrize)
    monkeypatch.setattr(
        md_module,
        "_load_bss_system_from_gromacs",
        fake_load_bss_system_from_gromacs,
    )
    monkeypatch.setattr(md_module, "_save_system_as_pdb", fake_save_system_as_pdb)
    monkeypatch.setattr(md_module, "solvate_openmm", fake_solvate_openmm)
    monkeypatch.setattr(
        md_module,
        "_collect_native_md_outputs",
        fake_collect_native_md_outputs,
    )

    fake_bss = SimpleNamespace(
        Process=SimpleNamespace(Gromacs=DummyProcess),
    )
    monkeypatch.setattr(md_module, "BSS", fake_bss)

    result = md_module.run_md_from_docking(
        protein_pdb=protein_pdb,
        docked_ligand_pose=ligand_pose,
        work_dir=work_dir,
        solvate=True,
        save_manifest=True,
        save_parametrized_pdb=True,
        save_solvated_pdb=True,
        save_md_final_pdb=True,
        save_md_native_outputs=True,
        md_output_stem="md_final",
    )

    assert result.prepared_protein_pdb == prepared_protein_pdb
    assert result.normalized_ligand_sdf == normalized_ligand_sdf
    assert result.parametrized_pdb == work_dir / "parametrization" / "complex_parametrized.pdb"
    assert result.solvated_pdb == work_dir / "solvation" / "solvated.pdb"
    assert result.final_md_pdb == final_md_pdb
    assert result.native_md_outputs.trajectory_xtc == native_outputs.trajectory_xtc
    assert result.native_md_outputs.portable_run_input_tpr == native_outputs.portable_run_input_tpr
    assert result.manifest_path == manifest_path

    assert manifest_path.exists()

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "success"
    assert manifest["protein"]["prepared_pdb"] == str(prepared_protein_pdb)
    assert manifest["ligand"]["normalized_sdf"] == str(normalized_ligand_sdf)
    assert manifest["md"]["final_pdb"] == str(final_md_pdb)
    assert manifest["md"]["native_outputs"]["trajectory_xtc"] == str(native_outputs.trajectory_xtc)
    assert manifest["md"]["native_outputs"]["portable_run_input_tpr"] == str(native_outputs.portable_run_input_tpr)
