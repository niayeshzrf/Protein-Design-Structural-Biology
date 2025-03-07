#!/usr/bin/env python 

from modeller import * 
from modeller.automodel import * 
from modeller.scripts import complete_pdb
import sys 

atom_files = sys.argv[1]
code = sys.argv[2]
outpdb = sys.argv[3]

log.verbose()
env = environ() 

# directories for input atom files 
atom_files = sys.argv[1]
env.io.atom_files_directory = ['.', atom_files]
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')

# set up the environment 
aln = alignment(env) 
mdl = model(env, file=code)
aln.append_model(mdl, atom_files=code, align_codes=code)

# application
mdl = complete_pdb(env, code, transfer_res_num = True) 

# add the mdl with added atoms 
aln.append_model(mdl, align_codes=code+"-1")

# generate molecular topology for the mutant 
mdl.clear_topology()
mdl.generate_topology(aln[code+"-1"])

# transfer all the coordinates you from the template native structure 
# to the mutant 
mdl.transfer_xyz(aln) 

# Build the remaining unknown coordinates for the mutant 
mdl.build(initialize_xyz=False, build_method="INTERNAL_COORDINATES")

# required to do a tranfer of res number 
mdl2 = model(env, file=code) 

# transfers from model 2 to model 1 
mdl.res_num_from(mdl2, aln) 

# write the mutant to a file
mdl.write(file=outpdb)

