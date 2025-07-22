# --- Import Libraries ---
import os

import numpy as np
from matscipy.dislocation import get_elastic_constants, BCCEdge111Dislocation
from matscipy.calculators.eam import EAM

from ase.io import write

# --- Paths ---
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARENT_DIR = os.path.dirname(os.path.abspath(__file__))
FILE_PATH = os.path.abspath(__file__)

POTENTIAL_PATH = os.path.join('potentials', 'malerba.fs')
OUTPUT_DIR = os.path.join(PARENT_DIR, 'output')

# --- Parameters ---

DISLO_CYL_RADIUS = 150
DISLO_LEN = 50

def main():

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    eam_calc = EAM(POTENTIAL_PATH)
    alat, C11, C12, C44 = get_elastic_constants(calculator=eam_calc, symbol='Fe', verbose='False') # Get information from the potential file for file creation
    print(f"{alat:.3f} (Angstrom), {C11:.2f}, {C12:.2f}, {C44:.2f} (GPa)") # Print information

    Fe_edge = BCCEdge111Dislocation(alat, C11, C12, C44, symbol="Fe") # Create dislocation object

    edge_bulk, edge_dislo = Fe_edge.build_cylinder(radius=DISLO_CYL_RADIUS) # Create a single plane of cylinder around the dislocation

    edge_dislo_long = edge_dislo.repeat((1,1,DISLO_LEN)) # Replicate the cylinder along the dislocation axis (z)

    print(f"Number of atoms: {len(edge_dislo_long)}") # Find the number of atoms in the sim

    write(os.path.join(OUTPUT_DIR, 'dislo.lmp'), edge_dislo_long, format="lammps-data", specorder=['Fe']) # Write the file out to lammps input file

    return None

# --- Entrypoint ---

if __name__ == "__main__":

    main()