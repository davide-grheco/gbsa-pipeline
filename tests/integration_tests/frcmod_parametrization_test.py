"""Tests for frcmod_parametrization and AMBER/GROMACS round-trips.

Unit tests cover input validation.
Functional tests exercise build_amber_ff_xml against MCPB.py testdata.
Integration tests validate OpenMM and ParmEd round-trips for the Zn-finger system.
"""

from __future__ import annotations

import os
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import parmed as pmd
import pytest
from openff.toolkit import Molecule
from openmm import LangevinIntegrator, Platform, unit
from openmm.app import ForceField, NoCutoff, PDBFile, Simulation
from openmm.app.modeller import Modeller
from openmmforcefields.generators import GAFFTemplateGenerator
from parmed.modeller import ResidueTemplateContainer

from gbsa_pipeline.frcmod_parametrization import (
    AmberFFInput,
    AmberInput,
    build_amber_ff_xml,
    load_amber_complex,
)
from gbsa_pipeline.parametrization import (
    ParametrizationConfig,
    ParametrizationInput,
    parametrize,
)
from gbsa_pipeline.parametrization_enum import ChargeMethod, LigandFF, ProteinFF

# Expected values read from the prmtop via ParmEd
_EXPECTED_ATOMS = 12814
_EXPECTED_RESIDUES = 826
_ZN_TYPE = "M1"
_ZN_CHARGE = 0.4314
_ZN_EPSILON_KCAL = 0.01492
_ZN_RMIN_ANG = 1.395

# ---------------------------------------------------------------------------
# Test data paths
# ---------------------------------------------------------------------------

_TESTDATA = Path("tests/testdata").resolve()

_FRCMOD = _TESTDATA / "Model4_mcpbpy.frcmod"
_MCPB_PDB = _TESTDATA / "Model4_mcpbpy.pdb"
_COMPLEX3_SDF = _TESTDATA / "complex3.sdf"

_MOL2S = (
    _TESTDATA / "CM1.mol2",
    _TESTDATA / "CM2.mol2",
    _TESTDATA / "HD1.mol2",
    _TESTDATA / "HD2.mol2",
    _TESTDATA / "ZN1.mol2",
)

DRY_PRMTOP = _TESTDATA / "Model4_dry.prmtop"
DRY_INPCRD = _TESTDATA / "Model4_dry.inpcrd"

_PROTEIN_XMLS = "amber14-all.xml"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _active_env_prefix() -> Path:
    """Return the active Pixi/Conda prefix from CONDA_PREFIX."""
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix is None:
        pytest.skip("CONDA_PREFIX is not set; cannot locate active Pixi/Conda environment.")
    assert prefix is not None
    return Path(prefix).resolve()


def _ambertools_parm_dir() -> Path:
    """Return the active AmberTools parm directory."""
    parm_dir = _active_env_prefix() / "dat" / "leap" / "parm"
    if not parm_dir.exists():
        pytest.skip(f"AmberTools parm directory not found: {parm_dir}")
    return parm_dir


def _build_test_xml(tmp_path: Path) -> Path:
    """Build combined MCPB XML for the Zn-finger test system."""
    return build_amber_ff_xml(
        AmberFFInput(
            frcmod_files=(_FRCMOD,),
            residue_mol2s=_MOL2S,
            protein_ff=ProteinFF.FF14SB,
            ligand_ff=LigandFF.GAFF,
            output_xml=tmp_path / "combined.xml",
        )
    )


def _load_mcpb_pdb() -> PDBFile:
    """Load the MCPB test PDB."""
    return PDBFile(str(_MCPB_PDB))


def _build_openmm_forcefield(xml_path: Path) -> ForceField:
    """Build an OpenMM ForceField from the standard protein XML plus MCPB XML."""
    return ForceField(_PROTEIN_XMLS, str(xml_path))


def _build_openmm_system_from_mcpb_pdb(
    xml_path: Path,
) -> tuple[PDBFile, ForceField, Any]:
    """Create an OpenMM System from the MCPB test PDB and combined XML."""
    pdb = _load_mcpb_pdb()
    forcefield = _build_openmm_forcefield(xml_path)
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )
    return pdb, forcefield, system


def _run_cpu_minimization(
    topology: Any,
    system: Any,
    positions: Any,
    max_iterations: int = 100,
) -> tuple[Simulation, Any, Any]:
    """Run a short CPU energy minimization and return simulation and energies."""
    integrator = LangevinIntegrator(
        300 * unit.kelvin,
        1.0 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    simulation = Simulation(
        topology,
        system,
        integrator,
        Platform.getPlatformByName("CPU"),
    )
    simulation.context.setPositions(positions)

    initial_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy(maxIterations=max_iterations)
    final_energy = simulation.context.getState(getEnergy=True).getPotentialEnergy()

    return simulation, initial_energy, final_energy


def _build_parmed_openmm_parameterset() -> pmd.openmm.OpenMMParameterSet:
    """Build a ParmEd/OpenMM parameter set from active AmberTools + MCPB frcmod."""
    parm_dir = _ambertools_parm_dir()

    base_parm = pmd.amber.AmberParameterSet(
        str(parm_dir / "parm19.dat"),
        str(parm_dir / "frcmod.ff14SB"),
        str(parm_dir / "gaff.dat"),
        str(_FRCMOD),
    )
    ff = pmd.openmm.OpenMMParameterSet.from_parameterset(base_parm)

    for mol2_file in _MOL2S:
        mol2 = pmd.load_file(str(mol2_file))
        if isinstance(mol2, ResidueTemplateContainer):
            ff.residues.update(mol2.to_library())
        else:
            ff.residues[mol2.name] = mol2

    return ff


def _assert_zn_parameters(struct: pmd.Structure) -> None:
    """Assert that Zn parameters survived the conversion."""
    zn_atoms = [atom for atom in struct.atoms if atom.atomic_number == 30]
    assert len(zn_atoms) == 1, "expected exactly one Zn atom"

    zn = zn_atoms[0]
    assert abs(zn.charge - 0.431432) < 1e-4, f"Zn charge {zn.charge} != 0.431432"
    assert abs(zn.epsilon - 0.014917) < 1e-4, f"Zn epsilon {zn.epsilon} != 0.014917 kcal/mol"
    assert abs(zn.rmin - 1.3950) < 1e-3, f"Zn Rmin/2 {zn.rmin} != 1.395 Å"


def _assert_all_forcefield_types_present(struct: pmd.Structure) -> None:
    """Assert masses and bonded/nonbonded parameters are assigned."""
    for atom in struct.atoms:
        assert atom.mass is not None, f"Mass not assigned for atom {atom.name}"
        assert atom.charge is not None, f"Charge not assigned for atom {atom.name}"
        assert atom.epsilon is not None, f"Epsilon not assigned for atom {atom.name}"

    for bond in struct.bonds:
        assert bond.type is not None, f"Bond type missing for bond between {bond.atom1.name} and {bond.atom2.name}"

    for angle in struct.angles:
        assert angle.type is not None, (
            f"Angle type missing for angle between {angle.atom1.name}, {angle.atom2.name}, {angle.atom3.name}"
        )

    for dihedral in struct.dihedrals:
        assert dihedral.type is not None, (
            f"Torsion type missing for dihedral between {dihedral.atom1.name}, "
            f"{dihedral.atom2.name}, {dihedral.atom3.name}, {dihedral.atom4.name}"
        )

    for improper in struct.impropers:
        assert improper.type is not None, (
            f"Improper type missing for improper between {improper.atom1.name}, "
            f"{improper.atom2.name}, {improper.atom3.name}, {improper.atom4.name}"
        )


# ---------------------------------------------------------------------------
# AmberInput validation (unit tests)
# ---------------------------------------------------------------------------


def test_amber_input_missing_prmtop_raises(tmp_path: Path) -> None:
    """AmberInput raises ValueError when prmtop does not exist."""
    inpcrd = tmp_path / "x.inpcrd"
    inpcrd.touch()
    with pytest.raises(ValueError, match="path_not_file"):
        AmberInput(prmtop=Path("nonexistent.prmtop"), inpcrd=inpcrd)


def test_amber_input_missing_inpcrd_raises(tmp_path: Path) -> None:
    """AmberInput raises ValueError when inpcrd does not exist."""
    prmtop = tmp_path / "x.prmtop"
    prmtop.touch()
    with pytest.raises(ValueError, match="path_not_file"):
        AmberInput(prmtop=prmtop, inpcrd=Path("nonexistent.inpcrd"))


def test_amber_input_output_dir_defaults_to_none(tmp_path: Path) -> None:
    """AmberInput.output_dir defaults to None."""
    prmtop = tmp_path / "x.prmtop"
    inpcrd = tmp_path / "x.inpcrd"
    prmtop.touch()
    inpcrd.touch()
    assert AmberInput(prmtop=prmtop, inpcrd=inpcrd).output_dir is None


# ---------------------------------------------------------------------------
# AmberFFInput validation (unit tests)
# ---------------------------------------------------------------------------


def test_amber_ff_input_missing_frcmod_raises(tmp_path: Path) -> None:
    """AmberFFInput raises ValueError when a frcmod file does not exist."""
    with pytest.raises(ValueError, match="Files not found"):
        AmberFFInput(frcmod_files=(tmp_path / "nonexistent.frcmod",))


def test_amber_ff_input_missing_mol2_raises(tmp_path: Path) -> None:
    """AmberFFInput raises ValueError when a mol2 file does not exist."""
    with pytest.raises(ValueError, match="Files not found"):
        AmberFFInput(residue_mol2s=(tmp_path / "nonexistent.mol2",))


def test_amber_ff_input_defaults() -> None:
    """AmberFFInput defaults: empty tuples, FF14SB, GAFF2, no output_xml."""
    inp = AmberFFInput()
    assert inp.frcmod_files == ()
    assert inp.residue_mol2s == ()
    assert inp.protein_ff == ProteinFF.FF14SB
    assert inp.ligand_ff == LigandFF.GAFF2
    assert inp.output_xml is None


# ---------------------------------------------------------------------------
# build_amber_ff_xml — functional tests
# ---------------------------------------------------------------------------


def test_build_amber_ff_xml_creates_file(tmp_path: Path) -> None:
    """build_amber_ff_xml writes a non-empty XML file."""
    xml = _build_test_xml(tmp_path)
    assert xml.exists()
    assert xml.stat().st_size > 0


def test_build_amber_ff_xml_no_missing_atomtypes(tmp_path: Path) -> None:
    """All atom types referenced by residue templates are defined in the XML."""
    xml = _build_test_xml(tmp_path)

    tree = ET.parse(xml)  # noqa: S314
    defined = {entry.get("name") for entry in tree.findall(".//AtomTypes/Type")} - {None}
    used = {
        atom_type for atom in tree.findall(".//Residues/Residue/Atom") if (atom_type := atom.get("type")) is not None
    }

    missing_types = used - defined
    assert not missing_types, f"Missing atom types: {sorted(missing_types)}"


def test_build_amber_ff_xml_loadable_by_openmm(tmp_path: Path) -> None:
    """The generated XML can be loaded by OpenMM's ForceField without errors."""
    xml = _build_test_xml(tmp_path)
    _build_openmm_forcefield(xml)


def test_build_amber_ff_xml_output_xml_none_creates_temp() -> None:
    """When output_xml is None, a file in a temp directory is returned."""
    xml = build_amber_ff_xml(
        AmberFFInput(
            frcmod_files=(_FRCMOD,),
            residue_mol2s=_MOL2S,
            protein_ff=ProteinFF.FF14SB,
            ligand_ff=LigandFF.GAFF,
        )
    )
    assert xml.exists()
    assert xml.suffix == ".xml"


# ---------------------------------------------------------------------------
# MCPB system tests: build_amber_ff_xml → ForceField → createSystem
# ---------------------------------------------------------------------------


def test_build_amber_ff_xml_create_system(tmp_path: Path) -> None:
    """build_amber_ff_xml XML + amber14-all.xml can createSystem on MCPB PDB."""
    xml = _build_test_xml(tmp_path)
    pdb, _forcefield, system = _build_openmm_system_from_mcpb_pdb(xml)

    assert system.getNumParticles() == pdb.topology.getNumAtoms()
    assert system.getNumForces() > 0


@pytest.mark.integration
def test_build_amber_ff_xml_energy_minimization(tmp_path: Path) -> None:
    """System built from the generated XML can be energy-minimized."""
    xml = _build_test_xml(tmp_path)
    pdb, _forcefield, system = _build_openmm_system_from_mcpb_pdb(xml)

    _simulation, initial_energy, final_energy = _run_cpu_minimization(
        pdb.topology,
        system,
        pdb.positions,
    )

    assert final_energy < initial_energy


@pytest.mark.integration
def test_parametrize_mcpb_with_ligand(tmp_path: Path) -> None:
    """Frcmod → XML → parametrization with a ligand produces valid GROMACS files."""
    xml = build_amber_ff_xml(
        AmberFFInput(
            frcmod_files=(_FRCMOD,),
            residue_mol2s=_MOL2S,
            protein_ff=ProteinFF.FF14SB,
            ligand_ff=LigandFF.GAFF,
            output_xml=tmp_path / "mcpb.xml",
        )
    )

    result = parametrize(
        ParametrizationInput(
            protein_pdb=_MCPB_PDB,
            ligand_sdf=_COMPLEX3_SDF,
            config=ParametrizationConfig(
                protein_ff=ProteinFF.FF14SB,
                ligand_ff=LigandFF.GAFF,
                extra_ff_files=(xml,),
            ),
            work_dir=tmp_path,
        )
    )

    assert result.gro_file.exists()
    assert result.gro_file.stat().st_size > 0
    assert result.top_file.exists()
    assert result.top_file.stat().st_size > 0

    struct = pmd.load_file(str(result.top_file), xyz=str(result.gro_file))

    residue_names = {residue.name for residue in struct.residues}
    assert {"CM1", "CM2", "HD1", "HD2", "ZN1"} <= residue_names

    _assert_zn_parameters(struct)
    _assert_all_forcefield_types_present(struct)


# ---------------------------------------------------------------------------
# Reference OpenMM/ParmEd-only tests
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_load_frcmod_openmm(tmp_path: Path) -> None:
    """Pure OpenMM/ParmEd reference test for MCPB frcmod loading."""
    pdb = _load_mcpb_pdb()
    ff = _build_parmed_openmm_parameterset()

    xml_path = tmp_path / "combined.xml"
    ff.write(str(xml_path), write_unused=True)
    assert xml_path.exists()

    forcefield = ForceField(_PROTEIN_XMLS, str(xml_path))
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )

    assert system.getNumParticles() == pdb.topology.getNumAtoms()
    assert system.getNumForces() > 0

    _simulation, initial_energy, final_energy = _run_cpu_minimization(
        pdb.topology,
        system,
        pdb.positions,
    )

    assert final_energy < initial_energy
    assert not unit.is_quantity(final_energy) or (final_energy._value == final_energy._value)

    structure = pmd.openmm.load_topology(pdb.topology, system, pdb.positions)

    gro_file = tmp_path / "complex.gro"
    top_file = tmp_path / "complex.top"
    structure.save(str(top_file), format="gromacs")
    structure.save(str(gro_file))

    struct = pmd.load_file(str(top_file), xyz=str(gro_file))
    assert len(struct.atoms) == _EXPECTED_ATOMS
    assert len(struct.residues) == _EXPECTED_RESIDUES

    residue_names = {residue.name for residue in struct.residues}
    assert {"CM1", "CM2", "HD1", "HD2", "ZN1"} <= residue_names

    _assert_zn_parameters(struct)
    _assert_all_forcefield_types_present(struct)


@pytest.mark.integration
def test_parametrize_mcpb_with_ligand_openmm(tmp_path: Path) -> None:
    """Pure OpenMM/GAFF reference test for MCPB protein + ligand."""
    pdb = _load_mcpb_pdb()
    ff = _build_parmed_openmm_parameterset()

    xml_path = tmp_path / "combined.xml"
    ff.write(str(xml_path), write_unused=True)
    assert xml_path.exists()

    ligand = Molecule.from_file(str(_COMPLEX3_SDF))
    if not ligand.conformers:
        raise ValueError(
            f"Ligand SDF '{_COMPLEX3_SDF}' contains no 3-D conformers. "
            "Provide an SDF file with embedded 3-D coordinates."
        )

    kwargs: dict[str, Any] = {
        "partial_charge_method": ChargeMethod.AM1BCC.value,
        "normalize_partial_charges": False,
        "use_conformers": ligand.conformers,
    }
    ligand.assign_partial_charges(**kwargs)

    forcefield = ForceField(_PROTEIN_XMLS, str(xml_path))
    gaff = GAFFTemplateGenerator(molecules=[ligand], forcefield="gaff-1.81", cache=None)
    forcefield.registerTemplateGenerator(gaff.generator)

    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0].to_openmm())

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=None,
    )

    assert system.getNumParticles() == modeller.topology.getNumAtoms()
    assert system.getNumForces() > 0

    _simulation, initial_energy, final_energy = _run_cpu_minimization(
        modeller.topology,
        system,
        modeller.positions,
    )

    assert final_energy < initial_energy
    assert not unit.is_quantity(final_energy) or (final_energy._value == final_energy._value)


# ---------------------------------------------------------------------------
# Load ZN1_Blimp1 pre-built AMBER files
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_load_amber_complex_creates_gromacs_files(tmp_path: Path) -> None:
    """load_amber_complex writes non-empty .gro and .top files."""
    result = load_amber_complex(AmberInput(prmtop=DRY_PRMTOP, inpcrd=DRY_INPCRD, output_dir=tmp_path))
    assert result.gro_file.exists()
    assert result.gro_file.stat().st_size > 0
    assert result.top_file.exists()
    assert result.top_file.stat().st_size > 0


@pytest.mark.integration
def test_load_amber_complex_atom_and_residue_count(tmp_path: Path) -> None:
    """Exported topology has the expected number of atoms and residues."""
    result = load_amber_complex(AmberInput(prmtop=DRY_PRMTOP, inpcrd=DRY_INPCRD, output_dir=tmp_path))
    struct = pmd.load_file(str(result.top_file), xyz=str(result.gro_file))
    assert len(struct.atoms) == _EXPECTED_ATOMS
    assert len(struct.residues) == _EXPECTED_RESIDUES


@pytest.mark.integration
def test_load_amber_complex_zn_type(tmp_path: Path) -> None:
    """Exported topology contains exactly one Zn atom with type M1."""
    result = load_amber_complex(AmberInput(prmtop=DRY_PRMTOP, inpcrd=DRY_INPCRD, output_dir=tmp_path))
    struct = pmd.load_file(str(result.top_file), xyz=str(result.gro_file))
    zn_atoms = [atom for atom in struct.atoms if atom.atomic_number == 30]
    assert len(zn_atoms) == 1
    assert zn_atoms[0].type == _ZN_TYPE


@pytest.mark.integration
def test_load_amber_complex_zn_charge(tmp_path: Path) -> None:
    """Zn RESP charge survives the prmtop → GROMACS round-trip."""
    result = load_amber_complex(AmberInput(prmtop=DRY_PRMTOP, inpcrd=DRY_INPCRD, output_dir=tmp_path))
    struct = pmd.load_file(str(result.top_file), xyz=str(result.gro_file))
    zn = next(atom for atom in struct.atoms if atom.atomic_number == 30)
    assert abs(zn.charge - _ZN_CHARGE) < 1e-4


@pytest.mark.integration
def test_load_amber_complex_zn_lj_params(tmp_path: Path) -> None:
    """Zn LJ parameters survive the prmtop → GROMACS round-trip."""
    result = load_amber_complex(AmberInput(prmtop=DRY_PRMTOP, inpcrd=DRY_INPCRD, output_dir=tmp_path))
    struct = pmd.load_file(str(result.top_file), xyz=str(result.gro_file))
    zn = next(atom for atom in struct.atoms if atom.atomic_number == 30)
    assert abs(zn.epsilon - _ZN_EPSILON_KCAL) < 1e-4
    assert abs(zn.rmin - _ZN_RMIN_ANG) < 1e-3
