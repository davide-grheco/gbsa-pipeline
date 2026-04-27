"""Unit tests for the small BioSimSpace MD helpers.

These tests focus on the local behavior of the gbsa-pipeline MD helper
functions and do not run external GROMACS processes. BioSimSpace process
creation is mocked because the unit tests should verify that this module wires
the protocol, system, and optional work directory correctly. Real process
execution belongs in later integration tests once the complete MD workflow is
available. This keeps the first MD helper tests small, fast, and aligned with
the public functions that are being reviewed.
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


def test_run_heating_builds_bss_process_and_returns_equilibrated_system(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that the heating helper delegates to BioSimSpace correctly.

    The helper should create a BioSimSpace ``Equilibration`` protocol with the
    requested runtime, a 0 K start temperature, a 300 K end temperature, and the
    current conservative backbone restraint default. The optional ``work_dir``
    should be converted to a string and passed only as a process keyword,
    matching the minimization helper behavior. The test mocks BioSimSpace
    because unit coverage should verify local protocol wiring, not run GROMACS.
    The returned object should be exactly the system returned by
    ``getSystem(block=True)``.
    """
    minimized_system = Mock(name="minimized_system")
    equilibrated_system = Mock(name="equilibrated_system")
    simulation_time = Mock(name="simulation_time")
    heating_protocol = Mock(name="heating_protocol")
    process = Mock(name="gromacs_process")
    process.getSystem.return_value = equilibrated_system

    equilibration_factory = Mock(return_value=heating_protocol)
    gromacs_factory = Mock(return_value=process)

    monkeypatch.setattr(md.BSS.Protocol, "Equilibration", equilibration_factory)
    monkeypatch.setattr(md.BSS.Process, "Gromacs", gromacs_factory)

    result = md.run_heating(
        simulation_time=simulation_time,
        minimized=minimized_system,
        work_dir=tmp_path,
    )

    equilibration_factory.assert_called_once_with(
        runtime=simulation_time,
        temperature_start=0 * md.BSS.Units.Temperature.kelvin,
        temperature_end=300 * md.BSS.Units.Temperature.kelvin,
        restraint="backbone",
    )
    gromacs_factory.assert_called_once_with(
        protocol=heating_protocol,
        system=minimized_system,
        work_dir=str(tmp_path),
    )
    process.start.assert_called_once_with()
    process.wait.assert_called_once_with()
    process.getSystem.assert_called_once_with(block=True)
    assert result is equilibrated_system


def test_run_heating_omits_work_dir_when_not_provided(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test that the heating helper does not invent a work directory.

    The ``work_dir`` argument is optional for the same reason as in the
    minimization helper: BioSimSpace can create its own process directory when
    the caller does not need a specific location. The helper should omit the
    keyword entirely instead of passing ``work_dir=None``. This test keeps the
    assertion local to process construction and return behavior. It does not
    validate physical heating behavior, which belongs in an integration test
    with a real parametrized system and a working GROMACS installation.
    """
    minimized_system = Mock(name="minimized_system")
    equilibrated_system = Mock(name="equilibrated_system")
    simulation_time = Mock(name="simulation_time")
    heating_protocol = Mock(name="heating_protocol")
    process = Mock(name="gromacs_process")
    process.getSystem.return_value = equilibrated_system

    equilibration_factory = Mock(return_value=heating_protocol)
    gromacs_factory = Mock(return_value=process)

    monkeypatch.setattr(md.BSS.Protocol, "Equilibration", equilibration_factory)
    monkeypatch.setattr(md.BSS.Process, "Gromacs", gromacs_factory)

    result = md.run_heating(
        simulation_time=simulation_time,
        minimized=minimized_system,
    )

    equilibration_factory.assert_called_once_with(
        runtime=simulation_time,
        temperature_start=0 * md.BSS.Units.Temperature.kelvin,
        temperature_end=300 * md.BSS.Units.Temperature.kelvin,
        restraint="backbone",
    )
    gromacs_factory.assert_called_once_with(
        protocol=heating_protocol,
        system=minimized_system,
    )
    process.start.assert_called_once_with()
    process.wait.assert_called_once_with()
    process.getSystem.assert_called_once_with(block=True)
    assert result is equilibrated_system


def test_run_npt_equilibration_builds_bss_process_and_returns_equilibrated_system(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that the NPT equilibration helper delegates to BioSimSpace correctly.

    The helper should create a BioSimSpace ``Equilibration`` protocol with the
    requested runtime, a 300 K target temperature, a 1 atm pressure target, and
    the current conservative backbone restraint default. The optional
    ``work_dir`` should be converted to a string and passed only as a process
    keyword, matching the other MD helper behavior. The test mocks BioSimSpace
    because unit coverage should verify local protocol wiring, not pressure
    equilibration with a real GROMACS process. The returned object should be
    exactly the system returned by ``getSystem(block=True)``.
    """
    heated_system = Mock(name="heated_system")
    equilibrated_system = Mock(name="equilibrated_system")
    simulation_time = Mock(name="simulation_time")
    equilibration_protocol = Mock(name="equilibration_protocol")
    process = Mock(name="gromacs_process")
    process.getSystem.return_value = equilibrated_system

    equilibration_factory = Mock(return_value=equilibration_protocol)
    gromacs_factory = Mock(return_value=process)

    monkeypatch.setattr(md.BSS.Protocol, "Equilibration", equilibration_factory)
    monkeypatch.setattr(md.BSS.Process, "Gromacs", gromacs_factory)

    result = md.run_npt_equilibration(
        simulation_time=simulation_time,
        heated=heated_system,
        work_dir=tmp_path,
    )

    equilibration_factory.assert_called_once_with(
        runtime=simulation_time,
        temperature=300 * md.BSS.Units.Temperature.kelvin,
        pressure=1 * md.BSS.Units.Pressure.atm,
        restraint="backbone",
    )
    gromacs_factory.assert_called_once_with(
        protocol=equilibration_protocol,
        system=heated_system,
        work_dir=str(tmp_path),
    )
    process.start.assert_called_once_with()
    process.wait.assert_called_once_with()
    process.getSystem.assert_called_once_with(block=True)
    assert result is equilibrated_system


def test_run_npt_equilibration_omits_work_dir_when_not_provided(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Test that the NPT equilibration helper does not invent a work directory.

    The ``work_dir`` argument remains optional so the caller can either choose a
    predictable process directory or let BioSimSpace create one. When no path is
    provided, the helper should omit the keyword entirely instead of passing
    ``work_dir=None``. This test checks only local process construction and
    return behavior, keeping real NPT stability questions out of unit coverage.
    A later integration test can cover the complete minimization, heating, and
    pressure-equilibration sequence with a real parametrized test system.
    """
    heated_system = Mock(name="heated_system")
    equilibrated_system = Mock(name="equilibrated_system")
    simulation_time = Mock(name="simulation_time")
    equilibration_protocol = Mock(name="equilibration_protocol")
    process = Mock(name="gromacs_process")
    process.getSystem.return_value = equilibrated_system

    equilibration_factory = Mock(return_value=equilibration_protocol)
    gromacs_factory = Mock(return_value=process)

    monkeypatch.setattr(md.BSS.Protocol, "Equilibration", equilibration_factory)
    monkeypatch.setattr(md.BSS.Process, "Gromacs", gromacs_factory)

    result = md.run_npt_equilibration(
        simulation_time=simulation_time,
        heated=heated_system,
    )

    equilibration_factory.assert_called_once_with(
        runtime=simulation_time,
        temperature=300 * md.BSS.Units.Temperature.kelvin,
        pressure=1 * md.BSS.Units.Pressure.atm,
        restraint="backbone",
    )
    gromacs_factory.assert_called_once_with(
        protocol=equilibration_protocol,
        system=heated_system,
    )
    process.start.assert_called_once_with()
    process.wait.assert_called_once_with()
    process.getSystem.assert_called_once_with(block=True)
    assert result is equilibrated_system
