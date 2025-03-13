#!/usr/bin/env python3

"""
Calculate Root-Mean-Square Deviation (RMSD) between two structures in XYZ
or PDB format using transformation and rotation.

Supports both **all-atom RMSD** and **backbone RMSD (N, CA, C, O)**.

Usage:
    python RMSD.py <structure_A.pdb> <structure_B.pdb> [--backbone]

Arguments:
    structure_A.pdb : First PDB file (reference structure)
    structure_B.pdb : Second PDB file (target structure)
    --backbone      : Optional flag to compute RMSD for backbone atoms only.
                      If omitted, computes all-atom RMSD.

Modified version of the RMSD calculation script.

Original code by:
Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
All rights reserved.

Modifications by: Niayesh Zarifi
- Added support for backbone-only RMSD.
- Optimized parsing and computation.

License: BSD 2-Clause (see full text below)

<Include the full license text here>

Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


import copy
import re
import sys
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


AXIS_SWAPS = np.array([
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 1, 0],
    [2, 0, 1]
])

AXIS_REFLECTIONS = np.array([
    [1, 1, 1],
    [-1, 1, 1],
    [1, -1, 1],
    [1, 1, -1],
    [-1, -1, 1],
    [-1, 1, -1],
    [1, -1, -1],
    [-1, -1, -1]
])


def parse_pdb(filename, backbone_only=False):
    """Extract atomic coordinates from a PDB file. Optionally filters for backbone atoms."""
    atoms = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if not backbone_only or atom_name in ["N", "CA", "C", "O"]:
                    x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])
                    atoms.append([x, y, z])
    return np.array(atoms)


def rmsd(V, W):
    """Compute the Root-Mean-Square Deviation (RMSD) between two sets of atomic coordinates."""
    return np.sqrt(((V - W) ** 2).sum() / len(V))


def kabsch_rmsd(P, Q):
    """Apply Kabsch algorithm for optimal alignment and compute RMSD."""
    P = P - np.mean(P, axis=0)
    Q = Q - np.mean(Q, axis=0)
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    R = U @ Vt
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = U @ Vt
    P = P @ R
    return rmsd(P, Q)


def main(structure_A, structure_B, backbone_only=False):
    """Compute RMSD between two PDB structures using the original algorithm."""
    coords_A = parse_pdb(structure_A, backbone_only=backbone_only)
    coords_B = parse_pdb(structure_B, backbone_only=backbone_only)
   

    if len(coords_A) != len(coords_B):
        print("Error: Structures have different numbers of selected atoms.")
        sys.exit(1)
    
    return kabsch_rmsd(coords_A, coords_B)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python RMSD.py <structure_A.pdb> <structure_B.pdb> [--backbone]")
        sys.exit(1)
    
    structure_A = sys.argv[1]
    structure_B = sys.argv[2]
    backbone_only = "--backbone" in sys.argv
    
    result = main(structure_A, structure_B, backbone_only)
    print(f"{'Backbone' if backbone_only else 'All-atom'} RMSD: {result:.3f} Å")

