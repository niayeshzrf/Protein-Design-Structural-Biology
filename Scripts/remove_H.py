#!/usr/bin/env python 

# this script aims for removing hydrogen atoms 
# usage: ./remove_H.py pdb_in.pdb pdb_out.pdb

def remove_H(x):
    new_pdb=[]
    for line in open(x): 
        if line[76:78].strip()=="H" : 
            continue 
        else:
            new_pdb+=line
    return(new_pdb) 


###########
# program #
###########

import sys 

pdb_in = sys.argv[1]
pdb_out = sys.argv[2]

filout = open(pdb_out, "w") 
pdb_noH=remove_H(pdb_in) 

for p in pdb_noH: 
    filout.write(p) 

filout.close
