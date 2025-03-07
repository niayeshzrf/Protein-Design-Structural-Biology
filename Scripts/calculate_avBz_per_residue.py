#!/usr/bin/env python

############
# function #
############

def calculate_avBz(x,i): # I need the ensemble 

    bz_dico={}
    bz=[]
    count=0
    for line in open(x): 
        
        # keep residue between 4 and 302
        
        if (int(line[22:26].strip()) == i and re.search(" CA | C | O | N ", line[12:16])) and line[17:20].strip() != "GLY" :   
            continue
        if int(line[22:26].strip()) == i and re.search("H", line[12:16]) :   
            continue
        elif int(line[22:26].strip()) == i:   
            count+=1
            bz+=[float(line[60:66].strip())]
#    print(bz)
    bz_sum=sum(bz)/len(bz)
    
    return(bz_sum)


###########
# program #
###########

import sys
import re

structure=sys.argv[1]
output_path=sys.argv[2]
protein_name=structure.split('/')[-1][:-4]
print(protein_name)

filout=open(output_path+"/"+protein_name+"_per_residue.dat", "w") 
for i in range(4, 302): 
#for i in [50,127]: 
    prot=calculate_avBz(sys.argv[1],i)
    filout.write("{:3f}\n".format(prot))
filout.close()
