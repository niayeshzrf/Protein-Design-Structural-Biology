#!/usr/bin/env python3

"""
PDB Cleaner and Processor
=========================
This script is designed to process and clean PDB files by applying various user-defined modifications. It allows for structure refinement, molecular modeling preparation, and customization of protein structures for computational studies.

Author: Niayesh Zarifi

Features:
---------
- **Solvent Removal (`--remove-solvent`)**
  - Removes all water molecules (HOH) from the structure.
  - Useful for molecular docking and dry simulations where explicit solvent is not needed.

- **Ligand Removal (`--remove-ligands`)**
  - Removes all `HETATM` records that are not water molecules.
  - Useful when isolating the protein structure for simulations or modelling.

- **Hydrogen Removal (`--remove-hydrogens`)**
  - Removes all hydrogen atoms from the structure.
  - Automatically renumbers atoms to maintain sequential order.
  - Required for force field parameterization in many computational tools.

- **Alternate Conformation Handling (`--keep-highest-occup`)**
  - Retains only the conformation with the highest occupancy for each atom.
  - Modifies residue names to remove alternate location indicators (e.g., AMET → MET).
  - Sets occupancy to 1.00 for all retained atoms to standardize the structure.

- **Chain Removal (`--remove-chains <chains>`)**
  - Deletes specific chains from the structure.
  - Accepts multiple chains, separated by commas (e.g., A,B,C).
  - Allows users to isolate a single chain or a subset of chains for analysis.

- **Residue Removal (`--remove-residues <list>`)**
  - Deletes residues with chain specification (e.g., A/10-20, B/5-8, B/123, 123).
  - Supports range-based removal and single residue deletions.
  - Can remove residues across multiple chains or globally.

- **Residue Renumbering (`--renumber`)**
  - Restarts residue numbering from 1 for each chain.
  - Ensures ligands and solvent continue numbering from the last residue of the associated chain.
  - Maintains correct sequence continuity across chains and molecule types.

Usage:
------
    python Process_PDB.py <input_pdb> <output_pdb> [options]

Options:
--------
    --remove-solvent       Remove solvent molecules (e.g., HOH)
    --remove-ligands       Remove all ligands (HETATM records, except HOH if not using --remove-solvent)
    --remove-hydrogens     Remove hydrogen atoms and renumber atoms
    --keep-highest-occup   Keep only the conformation with the highest occupancy, rename residues correctly, renumber atoms, and set occupancy to 1
    --remove-chains <chains> Remove specific chains (comma-separated, e.g., A,B,C)
    --remove-residues <list> Remove specific residues with chain specification (e.g., A/10-20,B/5-8,B/123,123)
    --renumber             Renumber residues from 1 per chain, ensuring ligands and solvent continue numbering from the last residue of their associated chain

Examples:
---------
1️⃣ **Remove Water, Ligands, and Hydrogen Atoms**
```bash
python Process_PDB.py input.pdb output.pdb --remove-solvent --remove-ligands --remove-hydrogens
```

2️⃣ **Remove Chains and Renumber Residues**
```bash
python Process_PDB.py input.pdb output.pdb --remove-chains A,B --renumber
```

3️⃣ **Remove Specific Residues**
```bash
python Process_PDB.py input.pdb output.pdb --remove-residues A/10-20,B/5-8,B/123,123
```

4️⃣ **Keep Only the Highest Occupancy Conformation**
```bash
python Process_PDB.py input.pdb output.pdb --keep-highest-occup
```

5️⃣ **Complete Refinement and Renumbering**
```bash
python Process_PDB.py input.pdb output.pdb --remove-solvent --remove-ligands --remove-hydrogens --remove-chains B --remove-residues A/50-100 --renumber
```

Use Cases:
----------
✅ **Protein-Only Refinement:**
   - Strip down a structure to only protein atoms for computational workflows.

✅ **MD Simulations & Force Fields:**
   - Remove water, hydrogens, and ligands to prepare PDB files for molecular dynamics simulations.

✅ **Structure-Based Drug Design:**
   - Keep only ligands or remove them depending on docking vs. receptor-based modeling.

✅ **Resequence PDB Files:**
   - Ensure continuous and logical residue numbering for chains and ligands.

This script provides an **efficient and customizable** approach to preparing **PDB structures** for various **computational and modeling applications**. 🚀
"""

import sys
import re

def parse_pdb(pdb_file):
    """Read a PDB file and return a list of lines."""
    with open(pdb_file, 'r') as file:
        return file.readlines()

def remove_ligands(pdb_lines):
    """Remove all HETATM records that are not water molecules (HOH)."""
    return [line for line in pdb_lines if not (line.startswith("HETATM") and "HOH" not in line[17:20])]

def remove_solvent(pdb_lines):
    """Remove solvent molecules (HOH)."""
    return [line for line in pdb_lines if not (line.startswith("HETATM") and "HOH" in line[17:20])]

def remove_hydrogens_and_renumber(pdb_lines):
    """Remove hydrogen atoms and renumber atom serial numbers."""
    new_pdb = []
    atom_number = 1
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            if line[76:78].strip() != "H":  # Keep only non-hydrogen atoms
                new_line = f"{line[:6]}{atom_number:5d}{line[11:]}"
                new_pdb.append(new_line)
                atom_number += 1
        else:
            new_pdb.append(line)
    return new_pdb

def keep_highest_occupancy(pdb_lines):
    """Keep only the highest occupancy conformation for each atom, rename residues correctly, renumber atoms, and set occupancy to 1."""
    atom_dict = {}
    atom_number = 1  # Start atom numbering from 1
    updated_pdb = []
    
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            alt_loc = line[16]  # Alternate location indicator
            chain_id = line[21]
            res_number = line[22:26].strip()
            occupancy = float(line[54:60])
            key = (atom_name, chain_id, res_number)
            
            if key not in atom_dict or occupancy > atom_dict[key][0]:
                # Remove alternate location indicator, update residue name, renumber atom, and set occupancy to 1
                updated_line = f"{line[:6]}{atom_number:5d}{line[11:16]} {residue_name}{line[20:54]}  1.00{line[60:]}"
                atom_dict[key] = (occupancy, updated_line)
                atom_number += 1  # Increment atom number sequentially
    
    return [atom_dict[key][1] for key in atom_dict]

def remove_non_atom_hetatm(pdb_lines):
    """Remove anything that is not an ATOM or HETATM record."""
    return [line for line in pdb_lines if line.startswith("ATOM") or line.startswith("HETATM")]

def remove_chains(pdb_lines, chains):
    """Remove specific chains from the PDB file."""
    chains_to_remove = set(chains.split(","))
    return [line for line in pdb_lines if line[21] not in chains_to_remove]

def remove_residues(pdb_lines, residues):
    """Remove specific residues with chain specification from the PDB file."""
    residues_to_remove = set()
    for entry in residues.split(","):
        if "/" in entry:
            chain, res_range = entry.split("/")
            if "-" in res_range:
                start, end = map(int, res_range.split("-"))
                residues_to_remove.update((chain, str(res)) for res in range(start, end + 1))
            else:
                residues_to_remove.add((chain, res_range))
        else:
            try:
                residues_to_remove.add((None, str(int(entry))))  # Ensure it's an integer
            except ValueError:
                print(f"Invalid residue specification: {entry}")
                continue

    return [line for line in pdb_lines if not ((line[21], line[22:26].strip()) in residues_to_remove or (None, line[22:26].strip()) in residues_to_remove or (line[21], line[22:26].strip()) in residues_to_remove)]

def renumber_pdb(pdb_lines):
    """Renumber residues sequentially from 1 per chain, ensuring ligands and solvent continue numbering from the last residue of their associated chain."""
    new_pdb = []
    residue_map = {}
    current_resnum = {}
    last_resnum = {}
    prev_chain = None
    prev_resnum = None
    
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain_id = line[21]
            resnum = int(line[22:26].strip())
            
            # If it's a new chain and has not been seen before, start from 1
            if chain_id not in current_resnum:
                current_resnum[chain_id] = 1
            
            # If it's a known chain, continue from the last residue
            elif chain_id in last_resnum:
                current_resnum[chain_id] = last_resnum[chain_id] + 1
            
            # Assign a new residue number if we haven't seen this resnum before
            if (chain_id, resnum) not in residue_map:
                residue_map[(chain_id, resnum)] = current_resnum[chain_id]
                current_resnum[chain_id] += 1
                
            last_resnum[chain_id] = residue_map[(chain_id, resnum)]  # Track last residue number
            prev_chain = chain_id
            prev_resnum = resnum
            
            # Modify the residue number in the line
            new_line = line[:22] + f"{residue_map[(chain_id, resnum)]:4d}" + line[26:]
            new_pdb.append(new_line)
        else:
            new_pdb.append(line)
    
    return new_pdb

def process_pdb(input_pdb, output_pdb, options):
    """Apply user-selected operations to clean the PDB file."""
    pdb_lines = parse_pdb(input_pdb)

    # Default cleanup (removes anything other than ATOM and HETATM)
    pdb_lines = remove_non_atom_hetatm(pdb_lines)

    if "--remove-solvent" in options:
        pdb_lines = remove_solvent(pdb_lines)
    if "--remove-ligands" in options:
        pdb_lines = remove_ligands(pdb_lines)
    if "--remove-hydrogens" in options:
        pdb_lines = remove_hydrogens_and_renumber(pdb_lines)
    if "--keep-highest-occup" in options:
        pdb_lines = keep_highest_occupancy(pdb_lines)
    if "--remove-chains" in options:
        chain_index = options.index("--remove-chains") + 1
        if chain_index < len(options):
            pdb_lines = remove_chains(pdb_lines, options[chain_index])
    if "--remove-residues" in options:
        res_index = options.index("--remove-residues") + 1
        if res_index < len(options):
            pdb_lines = remove_residues(pdb_lines, options[res_index])
    if "--renumber" in options:
        pdb_lines = renumber_pdb(pdb_lines)

    with open(output_pdb, 'w') as file:
        file.writelines(pdb_lines)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python Process_PDB.py <input_pdb> <output_pdb> [options]")
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    options = sys.argv[3:]

    process_pdb(input_pdb, output_pdb, options)

