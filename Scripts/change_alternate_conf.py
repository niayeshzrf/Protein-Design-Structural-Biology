#!/usr/bin/env python 

# this script aims to replace the major alternate conformation
# in a pdb file

import sys 
import re 

###########
# Program #
###########

def open_pdb(x): 
    pdb=[]
    seen=[]
    for line in open(x): 
        alternate_conf=line[16:17].strip()
        residue_id=line[17:20].strip()+line[22:26].strip()+line[21:22].strip()

        #if re.search("SO4|ACT|ANISOU", line[17:20].strip()):
        if re.search("PCA|XYP|XYS", line[17:20].strip()):
            continue
        if not re.search("^ATOM|^HETATM", line[0:6].strip()):
            continue
        if residue_id not in seen : 
            seen.append(residue_id)
            pdb.append([line])
        
        else:
            # here you want to create a list into the list of resid
            pdb[-1].append(line)

    return(pdb)

def change_alternate_conf(prot): # list_of_list
    pdb_line=prot
    new_pdb=[]
    for i in range(len(pdb_line)):
        alt=[]
        alt_loc=pdb_line[i][0][16:17]
        if alt_loc != " ": 
            occ=[]
            
            for line in pdb_line[i] : 
                occ.append(line[55:60].strip())
                alt.append(line[16:17].strip())
            high=occ[0]
            low=occ[0]
            av=occ[0]

            for e in range(len(occ)) :  
                if occ[e] > high : 
                    high=occ[e]
                elif occ[e] < low : 
                    low=occ[e]
                else:
                    av=occ[e]

            # test if there are alternate conformations C
            for line in pdb_line[i]: 
                if high==low: 
                    new_pdb+=line

                elif "C" in alt:
                    if float(high) == float(line[54:60].strip()): 
                        new_pdb+=[line[0:16]+"A"+line[17:]]
                    elif float(low) == float(line[54:60].strip()): 
                        new_pdb+=[line[0:16]+"C"+line[17:]]
                    elif float(av) == float(line[54:60].strip()): 
                        new_pdb+=[line[0:16]+"B"+line[17:]]

                elif "C" not in alt:   
                    if float(high) == float(line[54:60].strip()): 
                        new_pdb+=[line[0:16]+"A"+line[17:]]
                    elif float(low) == float(line[54:60].strip()): 
                        new_pdb+=[line[0:16]+"B"+line[17:]]
                 

        else: 

            new_pdb+=pdb_line[i]


    return(new_pdb)
    

############
# Function #
############


pdb_input=sys.argv[1]
pdb_output=sys.argv[2]
new_pdb=open(pdb_output, "w")

prot=open_pdb(pdb_input) 
new_prot=change_alternate_conf(prot)


for e in new_prot: 
    new_pdb.write(e)
new_pdb.close()
