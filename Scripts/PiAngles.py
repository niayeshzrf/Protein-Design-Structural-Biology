import numpy as np
from math import pow
from math import sqrt
from math import acos
from os import listdir

# helper function _isint  -----------------------------------------------------
def _isint(_val):
    """Returns True if _val can be converted to an integer, False if not."""
    try:
        int(_val)
        return True
    except ValueError:
        return False
#------------------------------------------------------------------------------
        

# extract specific pdb coordinates  -------------------------------------------
def _extractpdbcoordinates(_inputfile, _atoms):
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
def _computegeometry(_crds):
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


# read and extract msd scores  ------------------------------------------------
def _extractmsdscores(_inputfile, _strs=[]):
    """Reads a triad msd output file and returns sequence scores across
    ensemble states. Requires the input file (string). Key _strs sets the top
    (n, int) number of sequences to extract, default = [] (all)."""
    
    # open the msd file and read line-by-line
    _msd = open(_inputfile, 'r')
    _read = False
    _seqs = []
    for _line in _msd:
        # format _line
        _line = _line.replace('\n', '')
        _line = _line.split()
        # parse sequence information
        if _read and _isint(_line[0]):
            _seqno = int(_line[0])
            # build sequence string if it falls within the top _strs count
            if _seqno < _strs:
                _score = float(_line[2])
                _statescores = [float(_it) for _it in _line[5:-2]]
                _seq = list(_line[3])
                _name = ''
                for _s, _w in zip(_seq, _wt):
                    if _s == '-':
                        _name = _name + _w
                    else:
                        _name = _name + _s
                _seqs.append([_name, _score, _statescores])
            # terminate sequence count exceeded
            else:
                _read = False
        # terminate reading sequences
        else:
            _read = False
        # identify data location, set _read  to True
        if _read is False and len(_line) > 0\
           and _line[0] == 'Index' and _line[1] == 'Tags':
            _read = True
            _states = (_line[5:-2])
        # identify the wt sequence
        if _read is False and len(_line) > 0 and _line[0] == 'WT':
            _wt = _line[1:]
    # close multi_design.stdout and return states and scores
    _msd.close()
    return(_states, _seqs)
#------------------------------------------------------------------------------


# helper function to read directory and return files  -------------------------
def _readdirectory(_directory, _ext='.pdb'):
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


# write a .csv file for the multistate design output  -------------------------
def _writemsdcsv(_csvfile, _state, _sequences):
    """Write a .csv file with the specific formatting. Col-1: Sequence, Col-2:
    Fitness, and Col-3 to Col-n: Scores across all states."""
    
    # open, write, and close csv file
    _wfile = open(_csvfile, 'w')
    _wfile.write('Sequence,Fitness')
    for _state in _states:
        _wfile.write(','+_state)
    for _seq in _sequences:
        _wfile.write('\n')
        _wfile.write(_seq[0]+','+str(_seq[1]))
        for _scores in _seq[2]:
            _wfile.write(','+str(_scores))
    _wfile.close()
    return None
#------------------------------------------------------------------------------


# prepare a list of pdb files  ------------------------------------------------
def _preppdblist(_sequences, _states, _existingpdbs):
    """Generate the list of pdb files to perform geometry analysis against. The
    function expects the sequence list returned by _extractmsdscores, the state
    list returned by _extractmsd scores, and the list of existing pdb files
    found in the directory."""
    
    _pdbs = []
    _soln = 0
    for _seq in _sequences:
        _mpi = 0
        _pdb = []
        for _state in _states:
            _name = _state[0:-4]+'_designed_mpi'+str(_mpi)+'_soln'+str(_soln)+'.pdb'
            if any([_name == _it for _it in _existingpdbs]):
                _pdb.append(_name)
                _mpi = _mpi + 1
        _pdbs.append(_pdb)
        _soln = _soln + 1
    return(_pdbs)
#------------------------------------------------------------------------------


# helper function for pdb analysis  -------------------------------------------
def _analyzepdbs(_directory, _pdblist, _atoms):
    """Read and analyzed the pdb list."""
    
    _geoms = []
    for _pdbs in _pdblist:
        _geom = []
        _crds = []
        for _group in _atoms:
            _crd = _extractpdbcoordinates(_directory+_pdbs, _group)
            _crds.append(_crd)

        _dist, _angle = _computegeometry(_crds)

        _geom.append([_dist, _angle])
        _geoms.append(_geom)
    return(_geoms)
#------------------------------------------------------------------------------


def _writegeomcsv(_csvfile, _geometry, _pdbs):
    _wfile = open(_csvfile, 'w')
    _wfile.write('PDB,Distance,Angle')
    for _geom, _pdb in zip(_geometry, _pdbs):
        _wfile.write('\n'+_pdb)
        for _geo in _geom:
            print(_geo)
            _wfile.write(','+str(_geo[0][0])+','+str(_geo[1][0]))
    _wfile.close()
    return(None)

# RUN SCRIPT HERE  ------------------------------------------------------------
_directory = '/Volumes/Nia_HardDrive/ensemble_refinement_KEs/PiStacking_withAngle/KE70/Evolved6NT/'

_atoms=[['A48CG','A48CD2','A48CE2','A48CZ','A48CE1','A48CD1'],['A301C5','A301C4','A301C3A','A301C7A','A301C7','A301C6']]
_pdblist = _readdirectory(_directory)

_geometries = _analyzepdbs(_directory, _pdblist, _atoms)
_writegeomcsv(_directory+'multi_design_geometry.csv', _geometries, _pdblist)
