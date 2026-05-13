"""Integration test for the minimization helper."""

from __future__ import annotations

from pathlib import Path

import BioSimSpace as BSS
import pytest

from gbsa_pipeline.minimization import run_minimization


def _minimization_testdata_dir() -> Path:
    """Return the shared minimization testdata directory.

    The minimization integration test uses a prebuilt solvated GROMACS system as
    input. Keeping the path construction in one helper avoids duplicating
    repository-relative paths in the test body. The helper resolves paths from
    the location of this test file, so the test can be run from the repository
    root or through pytest collection. Individual file existence checks stay in
    the test body for clearer failure messages.
    """
    return Path(__file__).resolve().parents[1] / "testdata" / "minimization"


@pytest.mark.integration
def test_run_minimization(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Run the minimization helper on a prebuilt solvated system.

    This test exercises the public ``run_minimization`` entry point rather than
    constructing ``GromacsParams`` directly. That keeps the integration test
    focused on the production path that prepares and runs the GROMACS
    minimization protocol. The working directory is redirected to ``tmp_path`` so
    BioSimSpace/GROMACS artefacts do not leak into the repository tree. The
    assertion is intentionally coarse because this is a smoke-level integration
    test for successful execution.
    """
    testdata = _minimization_testdata_dir()
    gro_file = testdata / "example.gro"
    top_file = testdata / "example.top"

    assert gro_file.exists(), f"Missing testdata file: {gro_file}"
    assert top_file.exists(), f"Missing testdata file: {top_file}"

    monkeypatch.chdir(tmp_path)

    system = BSS.IO.readMolecules(
        [
            str(gro_file),
            str(top_file),
        ],
        make_whole=True,
    )

    minimized = run_minimization(
        nsteps=500,
        system=system,
        ignore_warnings=True,
    )

    assert minimized is not None
