#!/usr/bin/env python3

"""
This script is used to assess multiple designed protein sequences derived from a
sequence design application such as ProteinMPNN, ESM-IF1, or similar tools.

It processes multiple PDB files in a given directory and computes:
- The average pLDDT score across all residues.
- The average pLDDT for the 5 residues with the lowest pLDDT scores.
- The backbone RMSD (BB RMSD) compared to a reference PDB (optional),
  using the updated RMSD.py script with the --backbone flag.

If no reference PDB is provided, each structure is compared to itself for RMSD.
The results are saved in a user-specified CSV file.

Author: Niayesh Zarifi

Usage:
    python Analyse_AlphaFold_Outputs.py <input_pdb_dir> <output_csv> [ref_pdb]

Arguments:
    input_pdb_dir (str): Directory containing the generated PDB files.
    output_csv (str): Name of the output CSV file (must end in .csv).
    ref_pdb (str, optional): Reference PDB file for RMSD comparison.
                              If not provided, each PDB is compared to itself.
"""

import os
import sys
import statistics
import subprocess

def calculate_backbone_rmsd(structure_pdb, reference_pdb):
    """Calculate backbone RMSD using the updated RMSD.py script."""
    try:
        result = subprocess.run(["python", "RMSD.py", structure_pdb, reference_pdb, "--backbone"],
                                capture_output=True, text=True, check=True)
        rmsd_value = float(result.stdout.strip().split()[-2])  # Extract numeric RMSD value
        return rmsd_value
    except Exception as e:
        print(f"Error calculating backbone RMSD: {e}")
        return None

def analyze_af_outputs(input_pdb_dir, output_csv, ref_pdb=None):
    """Analyzes multiple AlphaFold output PDBs, calculates pLDDT metrics, and backbone RMSD."""
    
    if not output_csv.endswith(".csv"):
        print("Error: Output file must have a .csv extension.")
        sys.exit(1)
    
    if not os.path.isdir(input_pdb_dir):
        print("Error: The input must be a directory containing PDB files.")
        sys.exit(1)
    
    with open(output_csv, "w") as f:
        f.write("pdbfile,average pLDDT,median pLDDT,min pLDDT,backbone RMSD\n")
    
    for pdbfile in os.listdir(input_pdb_dir):
        input_pdb_path = os.path.join(input_pdb_dir, pdbfile)
        
        if os.path.isfile(input_pdb_path) and pdbfile.endswith(".pdb"):
            plddts = []
            
            with open(input_pdb_path, 'r') as infile:
                for line in infile:
                    tokens = line.split()
                    if len(tokens) > 10 and tokens[2] == 'CA':
                        plddts.append(float(tokens[10]))
            
            if plddts:
                avg = round(statistics.mean(plddts), 2)
                median = round(statistics.median(plddts), 2)
                avg_min = round(sum(sorted(plddts)[:5]) / 5, 2)
            else:
                print(f"Warning: No valid pLDDT values found in {pdbfile}. Skipping...")
                continue
            
            reference_pdb = ref_pdb if ref_pdb else input_pdb_path
            bb_rmsd = calculate_backbone_rmsd(input_pdb_path, reference_pdb)
            
            if bb_rmsd is not None:
                with open(output_csv, "a") as f:
                    f.write(f"{pdbfile},{avg},{median},{avg_min},{round(bb_rmsd,2)}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python Analyse_AlphaFold_Outputs.py <input_pdb_dir> <output_csv> [ref_pdb]")
        sys.exit(1)
    
    input_pdb_dir = sys.argv[1]
    output_csv = sys.argv[2]
    ref_pdb = sys.argv[3] if len(sys.argv) > 3 else None
    
    analyze_af_outputs(input_pdb_dir, output_csv, ref_pdb)

