#!/usr/bin/env python 

import __main__
__main__.pymol_argv=["pymol", "-qc"]

import pymol 
from pymol import cmd, stored

pymol.finish_launching()

import sys 
pdb=sys.argv[1]

pymol.cmd.load("../Xray/refine_3/"+pdb+"_noTSA/"+pdb+"_RT_noTSA_refine_001.pdb", pdb+"_noTSA") 
pymol.cmd.load("../Xray/refine_3/"+pdb+"_TSA/"+pdb+"_RT_TSA_refine_001.pdb", pdb+"_TSA")

# PDB PREPARATION 
pdblist=[pdb+"_noTSA",pdb+"_TSA"]

for p in pdblist: 
    pymol.cmd.remove("h. and "+p) 
    pymol.cmd.remove("resid 3+4+302 and "+p)
    pymol.cmd.remove("resn HOH or resn 6NT")


pymol.cmd.align(pdb+"_noTSA", pdb+"_TSA")




pymol.stored.residues_pdb1=[]
pymol.stored.residues_pdb2=[]


pymol.cmd.iterate(pdb+"_noTSA and n. CA and (not alt B and not alt C)" , "pymol.stored.residues_pdb1.append([resi, resn])") 

pymol.cmd.iterate(pdb+"_TSA and n. CA and (not alt B and not alt C)" , "pymol.stored.residues_pdb2.append([resi, resn])") 


# RMSD Backbone 
rmsd_bb=pymol.cmd.rms_cur("n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_noTSA and chain A","n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_TSA and chain A")

print("Backbone = {:.3f}".format(rmsd_bb))

# RMSD SC
rmsd_sc=pymol.cmd.rms_cur("not n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_noTSA and chain A","not n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_TSA and chain A")

print("Side chains = {:.3f}".format(rmsd_sc))

# RMSD Backbone act site 
act="21+44+50+83+84+127+172+236+237+265+267"

rmsd_bb_act=pymol.cmd.rms_cur("n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_noTSA and chain A and resi "+act,"n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_TSA and chain A and resi "+act)

print("Backbone active site = {:.3f}".format(rmsd_bb_act))

# RMSD SC act site 

rmsd_sc_act=pymol.cmd.rms_cur("not n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_noTSA and chain A and resi "+act,"not n. N+CA+C+O and (not alt B and not alt C) and "+pdb+"_TSA and chain A and resi "+act)

print("Side chains active site = {:.3f}".format(rmsd_sc_act))
