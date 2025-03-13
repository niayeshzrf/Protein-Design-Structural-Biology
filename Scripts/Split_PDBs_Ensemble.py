#!/usr/bin/env python3

"""
PDB Ensemble Splitter
=====================
Author: Niayesh Zarifi

This script processes a multi-template PDB file containing an ensemble of structures
and extracts each member into a separate PDB file.

Usage:
------
    python Split_PDBs_Ensemble.py <input_pdb>

Example:
--------
    python Split_PDBs_Ensemble.py ensemble.pdb

Outputs:
--------
    - Individual PDB files for each model in the ensemble (e.g., ensemble_1.pdb, ensemble_2.pdb, etc.)

"""

import sys
import os

def read_pdb_file(file_path):
    """Reads a PDB file and returns its contents as a list of lines."""
    with open(file_path, 'r') as pdb_file:
        return pdb_file.readlines()

def split_ensemble_to_models(pdb_lines):
    """Splits the input PDB lines into separate models and stores them in a dictionary."""
    models = {}
    current_model = None

    for line in pdb_lines:
        if line.startswith("MODEL"):
            current_model = int(line.split()[1])
            models[current_model] = [line]
        elif line.startswith("ENDMDL"):
            if current_model is not None:
                models[current_model].append(line)
                current_model = None
        elif current_model is not None:
            models[current_model].append(line)

    return models

def write_individual_models(models, input_filename):
    """Writes each extracted model to its own PDB file, using input filename as prefix."""
    base_name = os.path.splitext(os.path.basename(input_filename))[0]  # Remove .pdb extension

    for model_id, model_lines in models.items():
        output_filename = f"{base_name}_{model_id}.pdb"
        with open(output_filename, 'w') as output_file:
            output_file.writelines(model_lines)
        print(f"Saved: {output_filename}")

def main():
    """Main function to process the input PDB file and split its ensemble."""
    if len(sys.argv) != 2:
        print("Usage: python Split_PDBs_Ensemble.py <input_pdb>")
        sys.exit(1)

    input_pdb = sys.argv[1]
    pdb_lines = read_pdb_file(input_pdb)
    models = split_ensemble_to_models(pdb_lines)
    write_individual_models(models, input_pdb)

if __name__ == "__main__":
    main()
