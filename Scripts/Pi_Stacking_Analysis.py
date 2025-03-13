#!/usr/bin/env python3

"""
Pi-Stacking Distance and Angle Calculation
===========================================
Author: Niayesh Zarifi

This script analyzes multiple PDB files in a given directory, calculating pi-stacking
interactions by computing centroid distances and angles between aromatic ring systems.

Usage:
------
    python Pi_Stacking_Analysis.py <pdb_directory> <output_csv> <atom_list>

Arguments:
----------
    - pdb_directory: Path to the directory containing PDB files.
    - output_csv: Name of the output CSV file to store results.
    - atom_list: List of atom groups for pi-stacking calculations in the format:
      "[['A48CG','A48CD2','A48CE2','A48CZ','A48CE1','A48CD1'],
        ['A301C5','A301C4','A301C3A','A301C7A','A301C7','A301C6']]"

Example:
--------
    python Pi_Stacking_Analysis.py ./pdbs/ results.csv "[['A48CG','A48CD2','A48CE2','A48CZ','A48CE1','A48CD1'],['A301C5','A301C4','A301C3A','A301C7A','A301C7','A301C6']]"

Output:
-------
    - A CSV file containing centroid distances and angles between aromatic groups.
"""


import numpy as np
from math import pow
from math import sqrt
from math import acos
from os import listdir
import sys

# extract specific pdb coordinates  -------------------------------------------
def Extract_PDB_Coords(_inputfile, _atoms):
    """Extracts and returns a specific set of atomic coordinates from a protein
    data bank (pdb) file. Requires two arguments: location and name of the
    input file (string), and a list of the atoms to extract (list of strings)
    formatted: ChainResidueAtomname, i.e. A123HA for chain A, residue 123, atom
    HA."""   
    
    # create an empty list for the extraction of coordinates
    _listlength = len(_atoms)
    _crdlist = [[] for _it in range(_listlength)]
    # open and read the pdb file line-by-line
    _pdb = open(_inputfile, 'r')
    for _line in _pdb:
        _record = _line[0:6]
        _record = _record.replace(' ', '')
        # identify if the entry is for a set of atomic coordinates
        if _record == 'ATOM' or _record == 'HETATM':
            _chain = _line[21]
            _chain = _chain.replace(' ', '')
            _residue = _line[22:26]
            _residue = _residue.replace(' ', '')
            _atom = _line[12:16]
            _atom = _atom.replace(' ', '')
            _name = _chain+_residue+_atom
            # compare atom to list of desired atoms
            _it = 0
            _read = True
            while _read and _it < _listlength:
                if _name == _atoms[_it]:
                    _read = False
                else:
                    _it = _it + 1
            # if atom matches, read coordiantes and place in crd list
            if _read is False:
                _x = _line[30:38]
                _x = _x.replace(' ', '')
                _x = float(_x)
                _y = _line[38:46]
                _y = _y.replace(' ', '')
                _y = float(_y)
                _z = _line[46:54]
                _z = _z.replace(' ', '')
                _z = float(_z)
                _crdlist[_it] = [_x, _y, _z]
    # close the pdb file and return the coordinate list
    _pdb.close()
    return(_crdlist)
#------------------------------------------------------------------------------


# compute geometric information  ----------------------------------------------
def Compute_Geometry(_crds):
    """Compute geometry for a set of atomic coordinates. Specifically, the
    centroid for each coordinate group, the distance between centroids of each
    coordinate group, and the angle between coordinate groups. The function
    expects a list of coordinates (float), nested according to their grouping.
    For example: [[[xi1,yi1,zi1],[xi2,yi2,zi2]],[[xj1,yj1,zj1],[xj2,yj2,zj2]]],
    specifies two groups (i and j) of coordinates, each containing two atoms
    (1 and 2)."""
    
    # compute centroid for each coordinate grouping
    _cents = []
    for _group in _crds:
        _cent = [0.0, 0.0, 0.0]
        _n = 0
        for _crd in _group:
            _n = _n + 1
            _cent[0] = _cent[0] + _crd[0]
            _cent[1] = _cent[1] + _crd[1]
            _cent[2] = _cent[2] + _crd[2]
        _cent = [_cent[0] / _n, _cent[1] / _n, _cent[2] / _n]
        _cents.append(_cent)
    # compute the distance between each pair of centroids
    _dists = []
    for _i in range(0, len(_crds) - 1):
        for _j in range(_i + 1, len(_crds)):
            _dist = 0.0
            _dist = _dist + pow(_cents[_i][0] - _cents[_j][0], 2)
            _dist = _dist + pow(_cents[_i][1] - _cents[_j][1], 2)
            _dist = _dist + pow(_cents[_i][2] - _cents[_j][2], 2)
            _dist = sqrt(_dist)
            _dists.append(_dist)
    # compute the norm for each coordinate grouping by calculation of each
    # matrix's singular value decomposition
    _norms = []
    for _group in _crds:
        _matrix = np.array(_group, dtype=float)
        _matrix = np.transpose(_matrix)
        _U, _s, _V = np.linalg.svd(_matrix)
        _norms.append([_U[0][2], _U[1][2], _U[2][2]])
    # compute the angle between norms for pairs of groups
    _angles = []
    for _i in range(0, len(_crds) - 1):
        for _j in range(_i + 1, len(_crds)):
            _angle = _norms[_i][0] * _norms[_j][0]\
                   + _norms[_i][1] * _norms[_j][1]\
                   + _norms[_i][2] * _norms[_j][2]
            _length_i = pow(_norms[_i][0], 2)\
                      + pow(_norms[_i][1], 2)\
                      + pow(_norms[_i][2], 2)
            _length_i = sqrt(_length_i)
            _length_j = pow(_norms[_j][0], 2)\
                      + pow(_norms[_j][1], 2)\
                      + pow(_norms[_j][2], 2)
            _length_j = sqrt(_length_j)
            _length = _length_i * _length_j
            _angle = (acos(_angle / _length)) * 57.2957795
            _angles.append(_angle)
    # return centroid distances and angles for each group pair
    return(_dists, _angles)
#------------------------------------------------------------------------------

# helper function to read directory and return files  -------------------------
def Read_Dir(_directory, _ext='.pdb'):
    """Return a list of the files in the directory (string) provided they share
    the same extension (string)."""
    
    # list files in the directory, with the extension, and return list
    _filelist = []
    _files = listdir(_directory)
    for _file in _files:
        if _file[-len(_ext):] == _ext:
            _filelist.append(_file)
    return(_filelist)
#------------------------------------------------------------------------------

# helper function for pdb analysis  -------------------------------------------
def Analyze_Pi_Stacking(_directory, _pdblist, _atoms):
    """Read and analyzed the pdb list."""
    
    _geoms = []
    for _pdbs in _pdblist:
        _geom = []
        _crds = []
        for _group in _atoms:
            _crd = Extract_PDB_Coords(_directory+_pdbs, _group)
            _crds.append(_crd)

        _dist, _angle = Compute_Geometry(_crds)

        _geom.append([_dist, _angle])
        _geoms.append(_geom)
    return(_geoms)
#------------------------------------------------------------------------------

def Save_Results(_csvfile, _geometry, _pdbs):
    _wfile = open(_csvfile, 'w')
    _wfile.write('PDB,Distance,Angle')
    for _geom, _pdb in zip(_geometry, _pdbs):
        _wfile.write('\n'+_pdb)
        for _geo in _geom:
            _wfile.write(','+str(_geo[0][0])+','+str(_geo[1][0]))
    _wfile.close()
    return(None)

# RUN SCRIPT HERE  ------------------------------------------------------------

def save_output(results, output_csv):
    """Saves computed pi-stacking distances and angles to a CSV file."""
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["PDB File", "Group 1", "Group 2", "Distance (Å)", "Angle (°)"])
        for row in results:
            writer.writerow([row[0], ",".join(row[1]), ",".join(row[2]), row[3], row[4]])

def main():
    """Main function to process the PDB directory and compute pi-stacking interactions."""
    if len(sys.argv) != 4:
        print("Usage: python pi_angles_modified.py <pdb_directory> <output_csv> <atom_list>")
        sys.exit(1)
    
    pdb_directory = sys.argv[1]
    output_csv = sys.argv[2]
    ring_atoms = eval(sys.argv[3])
    
    pdb_files = Read_Dir(pdb_directory)
    results = Analyze_Pi_Stacking(pdb_directory, pdb_files, ring_atoms)
    Save_Results(output_csv, results , pdb_files)
    print(f"Analysis complete. Results saved to {output_csv}")

if __name__ == "__main__":
    main()


