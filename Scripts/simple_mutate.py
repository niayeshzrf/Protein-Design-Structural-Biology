#!/usr/bin/env python3 

from modeller import * 
import sys

atom_files = sys.argv[1]
code = sys.argv[2]
outpdb= sys.argv[3]

env = environ() 
env.io.atom_files_directory = [atom_files]

# Read the topology library to produce a mutant with hydrogens 
env.libs.topology.read(file='$(LIB)/top_allh.lib')

# Read the CHARMM parameter library: 
env.libs.parameters.read(file='$(LIB)/par.lib') 

# Read the original PDB file and copy its sequence to the alignment array
aln = alignment(env) 
mdl = model(env, file=code)
aln.append_model(mdl, atom_files=code, align_codes=code) 


# Select the residues to be mutated 
# Modeller re number residue, make sure you consider the first residue 
# of your protein 

first = mdl.residues[0].num
mutation = 84 - int(first)

# position to be mutated 
sel = selection(mdl.residues[mutation])

# Mutate the selected residues into CYS residue 
sel.mutate(residue_type = 'CYS')

# Add the mutated sequence to the alignment arrays (it is now the second 
# sequence in the alignment): 
aln.append_model(mdl, align_codes=code+'-1')

# Generate molecular topology for the mutant: 
mdl.clear_topology() 
mdl.generate_topology(aln[code+'-1'])

# transfer all the coordinates you can from the template native strucure 
# to the mutant 
mdl.transfer_xyz(aln)

# Build the remaining unknown coordinates for the mutant: 
mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

# required to do a transfer of res number 

mdl2 = model(env, file=code) 
#aln.append_model(mdl2, atom_files=code, align_codes=code) 

# Transfers from model 2 to model 1 
mdl.res_num_from(mdl2, aln)

# Write the mutant to a file: 
mdl.write(file=outpdb)
