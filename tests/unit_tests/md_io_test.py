"""Unit tests for BioSimSpace/GROMACS IO boundary helpers.

These tests cover the small file-based boundary added for external workflow
orchestration. They do not run GROMACS and do not require real molecular input
files, because the behavior under test is local path validation and delegation
to ``BioSimSpace.IO``. The BioSimSpace IO calls are monkeypatched so the tests
remain fast and focused. Full integration coverage for whether the written
``.gro`` and ``.top`` files are scientifically usable belongs to the existing
integration tests.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gbsa_pipeline import md_io


def test_load_bss_system_from_gromacs_validates_and_delegates_to_bss(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Load a system from existing GROMACS files through BioSimSpace IO.

    The helper should resolve both paths, validate that they exist as files, and
    then pass the matching ``.gro``/``.top`` pair to ``BSS.IO.readMolecules``.
    The test uses placeholder text files because BioSimSpace itself is not under
    test here. Returning a sentinel object verifies that the helper returns the
    BioSimSpace-loaded system unchanged. This keeps the test focused on the
    workflow boundary rather than on molecular file parsing.
    """
    gro_file = tmp_path / "system.gro"
    top_file = tmp_path / "system.top"
    gro_file.write_text("dummy gro\n", encoding="utf-8")
    top_file.write_text("dummy top\n", encoding="utf-8")

    expected_system = object()
    calls: list[list[str]] = []

    def fake_read_molecules(paths: list[str]) -> object:
        calls.append(paths)
        return expected_system

    monkeypatch.setattr(md_io.BSS.IO, "readMolecules", fake_read_molecules)

    loaded = md_io.load_bss_system_from_gromacs(gro_file, top_file)

    assert loaded is expected_system
    assert calls == [[str(gro_file.resolve()), str(top_file.resolve())]]


def test_load_bss_system_from_gromacs_reports_missing_coordinate_file(
    tmp_path: Path,
) -> None:
    """Fail early when the GROMACS coordinate file is missing.

    Workflow engines decide whether a parent action is complete from declared
    output files, but a library helper should still report missing files
    explicitly when called directly. This test verifies that the coordinate file
    is checked before BioSimSpace is invoked. The error message includes the
    problematic path so failed workflow actions are easier to diagnose. The
    topology file exists here to isolate the missing-coordinate case.
    """
    gro_file = tmp_path / "missing.gro"
    top_file = tmp_path / "system.top"
    top_file.write_text("dummy top\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="GROMACS coordinate file not found"):
        md_io.load_bss_system_from_gromacs(gro_file, top_file)


def test_load_bss_system_from_gromacs_reports_missing_topology_file(
    tmp_path: Path,
) -> None:
    """Fail early when the GROMACS topology file is missing.

    The ``.gro`` and ``.top`` files form a pair for BioSimSpace loading, so a
    missing topology should be reported before the lower-level loader runs. This
    keeps the workflow boundary deterministic and gives callers a clear product
    path to inspect. The coordinate file exists here to isolate the
    missing-topology case. No BioSimSpace monkeypatch is required because the
    function should fail before reaching BioSimSpace.
    """
    gro_file = tmp_path / "system.gro"
    top_file = tmp_path / "missing.top"
    gro_file.write_text("dummy gro\n", encoding="utf-8")

    with pytest.raises(FileNotFoundError, match="GROMACS topology file not found"):
        md_io.load_bss_system_from_gromacs(gro_file, top_file)


def test_save_bss_system_to_gromacs_creates_directory_and_delegates_to_bss(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Save a BioSimSpace system as matching GROMACS coordinate/topology files.

    The helper should create the output directory, construct ``.gro`` and
    ``.top`` paths from a shared prefix, and call ``BSS.IO.saveMolecules`` with
    the expected BioSimSpace file formats. The monkeypatched save function
    writes small placeholder files so the test can also verify the returned
    paths. The molecular system itself is a sentinel object because this unit
    test only covers IO delegation and path handling.
    """
    system = object()
    output_prefix = tmp_path / "nested" / "system"

    calls: list[tuple[Path, object, str]] = []

    def fake_save_molecules(
        path: str,
        saved_system: object,
        *,
        fileformat: str,
    ) -> None:
        output_path = Path(path)
        calls.append((output_path, saved_system, fileformat))
        output_path.write_text(f"saved as {fileformat}\n", encoding="utf-8")

    monkeypatch.setattr(md_io.BSS.IO, "saveMolecules", fake_save_molecules)

    gro_file, top_file = md_io.save_bss_system_to_gromacs(system, output_prefix)

    expected_gro = output_prefix.resolve().with_suffix(".gro")
    expected_top = output_prefix.resolve().with_suffix(".top")

    assert gro_file == expected_gro
    assert top_file == expected_top
    assert gro_file.exists()
    assert top_file.exists()
    assert calls == [
        (expected_gro, system, "gro87"),
        (expected_top, system, "grotop"),
    ]
