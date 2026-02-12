"""Performing first MD Run.

We start from SDF file, read ligand, standardize and hydrogen it;
We load protein;
We parametrize ligand and load the protein force field parameters;
We add solvent of a chosen water model and counter ions;
We minimize the system;
We perform equilibration with restraints by stepwise heating the system from 0 to 300 K.
"""

from __future__ import annotations

import logging

import BioSimSpace as BSS

from gbsa_pipeline.gro_custom import run_gro_custom


def main() -> None:
    """Reading the minimized system data."""
    test_system = BSS.IO.readMolecules(
        files=["tests/testdata/test.gro", "tests/testdata/test.top"],
        make_whole=True,
    )
    # Running customized MD RUN.
    customized = run_gro_custom("tests/testdata/paramaterer.config", system=test_system)

    BSS.IO.saveMolecules("customized.pdb", customized, fileformat="PDB")

    logging.info("Custom Run done")


if __name__ == "__main__":
    main()
