#!/usr/bin/env python3

#The script should check ligand files
# It should check what type of ligands (SDF, PDB, MOL) are given
# It should read the guest molecules and check their dimensionality


from rdkit import Chem
from pathlib import Path
import sys

def main(directory: str = ".", log_path: str = "run.log") -> int:
    d = Path(directory)
    log_file = Path(log_path)

    sdf_files = sorted(d.glob("*.sdf"))

###########################################################################################
#                       Printing error if the SDF diles are not there                     #
###########################################################################################
  
    if not sdf_files:
        msg = f"ERROR: No .sdf file found in directory: {d.resolve()}"
        print(msg, file=sys.stderr)  # 1) to screen (stderr)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with log_file.open("a", encoding="utf-8") as f:  # 2) to log file
            f.write(msg + "\n")
        return 1



###########################################################################################
#                       Reading the molecules from SDF file                               #
###########################################################################################

  
    name = sdf_files[0].name  # includes ".sdf"
    suppl = Chem.SDMolSupplier(name)               # Reading molecules
  
    return 0

if __name__ == "__main__":
    # optional args: directory and log file
    directory = sys.argv[1] if len(sys.argv) > 1 else "."
    log_path = sys.argv[2] if len(sys.argv) > 2 else "run.log"
    raise SystemExit(main(directory, log_path))

  n_records = 0
    n_fail = 0
    unique_names = set()

    for mol in suppl:
        if mol is None:
            n_fail += 1
            continue
        n_records += 1
        if mol.HasProp("_Name"):
            unique_names.add(mol.GetProp("_Name"))

    # In SDF: typically 1 record = 1 molecule entry (often 1 conformer per entry)
    n_molecules = n_records
    n_conformers = n_records  # if your conformers are stored as separate SDF records

    write_both(f"SDF: {sdf_path.name}", log_file)
    write_both(f"Molecules (records): {n_molecules}", log_file)
    write_both(f"Conformers (records): {n_conformers}", log_file)
    write_both(f"Parse failures: {n_fail}", log_file)
    write_both(f"Unique _Name entries: {len(unique_names)}", log_file)




