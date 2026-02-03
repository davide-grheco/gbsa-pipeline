
from os import PathLike
from rdkit import Chem
from rdkit.Chem import Mol
from molvs import Standardizer
import sire as sr


def load_ligand_sdf(sdf_path: PathLike) :
    """Load the first molecule from an SDF file.

    NOTE: RDKit reads a list of molecules from an SDF file.
    Only the first molecule is processed.
    """
    molecule_list = Chem.SDMolSupplier(str(sdf_path), remove_HS=false)


    if len(molecule_list) == 0 or molecule_list[0] is None:
        raise ValueError(f"No valid molecules found in {sdf_path}")

    mol = molecule_list[0]
    return mol


def ligand_standardizer(mol: Mol) -> Mol:

    """Standardize a ligand using MolVS."""
    s=Standardizer()
    mol_standard = s.standardize(mol)
    return mol_standard

def ligand_converter(mol: Mol):

    """Converting standardized ligand to sire format"""
    mol_standard = ligand_standardizer(mol)
    return sr.convert.to(mol_standard, "sire")