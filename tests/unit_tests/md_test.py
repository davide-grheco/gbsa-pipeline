"""Unit tests for the small BioSimSpace MD helpers.

These tests focus on the local behavior of the gbsa-pipeline MD helper
functions and do not run external GROMACS processes. BioSimSpace process
creation is mocked because the unit test should verify that this module wires
the protocol, system, and optional work directory correctly. Real process
execution belongs in a later integration test once the complete MD workflow is
available. This keeps the first minimization helper test small, fast, and
aligned with the public function that is being reviewed.
"""

from pathlib import Path
from unittest.mock import Mock

import pytest

from gbsa_pipeline import md


def test_run_minimization_builds_bss_process_and_returns_minimized_system(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that the minimization helper delegates to BioSimSpace correctly.

    The helper should create a BioSimSpace ``Minimisation`` protocol and pass it
    together with the input system to ``BSS.Process.Gromacs``. The optional
    ``work_dir`` should be converted to a string because BioSimSpace process
    constructors commonly expect file-system paths in that form. The test uses
    mocks instead of a real Sire system because the unit-level contract is the
    local wiring, not molecular dynamics execution. The returned object should
    be exactly the system returned by ``getSystem(block=True)``.
    """
    input_system = Mock(name="input_system")
    minimized_system = Mock(name="minimized_system")
    minimization_protocol = Mock(name="minimization_protocol")
    process = Mock(name="gromacs_process")
    process.getSystem.return_value = minimized_system

    minimisation_factory = Mock(return_value=minimization_protocol)
    gromacs_factory = Mock(return_value=process)

    monkeypatch.setattr(md.BSS.Protocol, "Minimisation", minimisation_factory)
    monkeypatch.setattr(md.BSS.Process, "Gromacs", gromacs_factory)

    result = md.run_minimization(
        system=input_system,
        work_dir=tmp_path,
    )

    minimisation_factory.assert_called_once_with()
    gromacs_factory.assert_called_once_with(
        protocol=minimization_protocol,
        system=input_system,
        work_dir=str(tmp_path),
    )
    process.start.assert_called_once_with()
    process.wait.assert_called_once_with()
    process.getSystem.assert_called_once_with(block=True)
    assert result is minimized_system


def test_run_minimization_omits_work_dir_when_not_provided(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test that the minimization helper does not invent a work directory.

    The ``work_dir`` argument is optional because callers may prefer to let
    BioSimSpace create and manage its own process directory. When no path is
    provided, the helper should not pass ``work_dir=None`` explicitly because
    that can differ from omitting the keyword for third-party APIs. The test
    checks only this small branching behavior and reuses mocks to avoid any
    dependency on external simulation tools. The returned minimized system is
    still taken from ``getSystem(block=True)`` to keep the behavior identical to
    the explicit-work-directory case.
    """
    input_system = Mock(name="input_system")
    minimized_system = Mock(name="minimized_system")
    minimization_protocol = Mock(name="minimization_protocol")
    process = Mock(name="gromacs_process")
    process.getSystem.return_value = minimized_system

    minimisation_factory = Mock(return_value=minimization_protocol)
    gromacs_factory = Mock(return_value=process)

    monkeypatch.setattr(md.BSS.Protocol, "Minimisation", minimisation_factory)
    monkeypatch.setattr(md.BSS.Process, "Gromacs", gromacs_factory)

    result = md.run_minimization(system=input_system)

    minimisation_factory.assert_called_once_with()
    gromacs_factory.assert_called_once_with(
        protocol=minimization_protocol,
        system=input_system,
    )
    process.start.assert_called_once_with()
    process.wait.assert_called_once_with()
    process.getSystem.assert_called_once_with(block=True)
    assert result is minimized_system
