#!/usr/bin/env python3

"""
This script calculates:

1. Ensemble Deviation: Measures how different each ensemble structure (model) is from the crystal structure (xtal). This is quantified as the average backbone RMSD (N, CA, C, O atoms) between each ensemble member and the crystal structure. Higher deviation values indicate greater structural differences from the reference.

2. Ensemble Diversity: Measures how structurally varied the ensemble members are relative to each other. This is quantified as the average pairwise backbone RMSD (N, CA, C, O atoms) between all models in the ensemble. High diversity values suggest greater conformational variability within the ensemble.

**Important:** Alternate rotamers should be removed from the crystal structure before running this script to ensure consistency in backbone RMSD calculations.

Author: Niayesh Zarifi

Usage:
    python Calculate_Deviation_Diversity_Ens.py <crystal_structure.pdb> <ensemble_dir> <output_file.csv>

Arguments:
    crystal_structure.pdb : The reference crystal structure PDB file.
    ensemble_dir (str)    : Directory containing the ensemble PDB files.
    output_file.csv (str) : Name of the output CSV file (must end in .csv).

Outputs:
    - CSV file containing:
      * Backbone RMSD values of ensemble members vs. crystal structure (Deviation)
      * Pairwise backbone RMSD values among ensemble members (Diversity)
"""

import os
import sys
import RMSD
import csv
import itertools
import time

def deviation(ensemble, crystal_structure):
    """Calculate backbone RMSD of each ensemble member against the crystal structure."""
    RMSDs_wrt_xtal = []
    for structure in ensemble:
        rmsd = RMSD.main(structure, crystal_structure, backbone_only=True)
        RMSDs_wrt_xtal.append(rmsd)
    return RMSDs_wrt_xtal

def diversity(ensemble):
    """Calculate pairwise backbone RMSD among ensemble members."""
    all_iterations = list(itertools.combinations(ensemble, 2))
    pair_rmsds = [RMSD.main(pair[0], pair[1], backbone_only=True) for pair in all_iterations]
    return pair_rmsds

def write_rmsd(results, output_file):
    """Write RMSD results to a CSV file."""
    with open(output_file, 'w', newline='') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        for value in results:
            wr.writerow([value])

def average(lst):
    """Calculate the average of a list of values."""
    return sum(lst) / len(lst) if lst else 0.0

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python Calculate_Deviation_Diversity_Ens.py <crystal_structure.pdb> <ensemble_dir> <output_file.csv>")
        sys.exit(1)
    
    crystal_structure = sys.argv[1]
    ensemble_dir = sys.argv[2]
    output_file = sys.argv[3]
    
    if not output_file.endswith(".csv"):
        print("Error: Output file must have a .csv extension.")
        sys.exit(1)
    
    # Collect ensemble PDB files
    ensemble = [os.path.join(ensemble_dir, f) for f in os.listdir(ensemble_dir) if f.endswith(".pdb")]
    
    if not ensemble:
        print("Error: No PDB files found in the ensemble directory.")
        sys.exit(1)
    
    # Compute deviation and diversity
    RMSDs_wrt_xtal = deviation(ensemble, crystal_structure)
    pair_rmsds = diversity(ensemble)
    
    # Write results
    write_rmsd(RMSDs_wrt_xtal, f"deviation_{output_file}")
    write_rmsd(pair_rmsds, f"diversity_{output_file}")
    
    # Compute mean values
    mean_deviation = average(RMSDs_wrt_xtal)
    mean_diversity = average(pair_rmsds)
    
    print(f"Mean Backbone Deviation: {mean_deviation:.3f} Å")
    print(f"Mean Backbone Diversity: {mean_diversity:.3f} Å")
    
    # Save summary results
    with open(output_file, 'w', newline='') as myfile:
        myfile.write("Deviation,Diversity\n")
        myfile.write(f"{mean_deviation},{mean_diversity}\n")
