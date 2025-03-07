##!/usr/bin/env python2

import __main__
__main__.pymol_argv = ['pymol', '-qc']

import pymol as py
from pymol import cmd, stored 
import os 


pdb_list = ["Seq12_15_B_h_pm_run3"]
dir_energy = "Energies_formatted"

#if os.path.isdir("sasa") == True : 
#    os.remove("sasa/")

#os.mkdir("sasa_ligand")

def open_list(x): 
    final_list = []
    for line in open(x): 
        final_list.append(line.split(",")[0])

    return(final_list) 

def open_pdb(pdb): 

    # reference 
    #cmd.load("../../../xtal/HG4_TSA_template.pdb", "HG4") 
#    cmd.select("act_HG4", "resi 21+44+50+84+127+172+236+265+267+237 and HG4") 
#    cmd.remove("resn NBX and HG4") 

    # set the options to calculate SASA
    cmd.set("dot_solvent", 1)
    cmd.set("dot_density", 3)



    # write the reference in the output file
    filout=open("sasa/"+pdb+".csv", "w")
 #   filout.write("{},{:.3f}\n".format("HG4", float(cmd.get_area("act_HG4"))))
    #filout.write("{},{:.3f}\n".format("HG4", float(cmd.get_area("resn NBX and HG4"))))

    # list of the pdb 
    list_energy = open_list(dir_energy+"/"+pdb+"_results.csv") 
    # dir of the energies 
    dir_pdb = "Best_pdb/"+pdb+"/"

    for p in list_energy: 
        # designs
        cmd.load(dir_pdb+p, p) 

        cmd.remove("resname PTL")
        cmd.remove("resname PLN")
        #cmd.select("act"+p, "resi 21+44+50+84+127+172+236++237+265+267 and "+p) 
        #cmd.remove("resn NBX and "+p) 

        # set the options to calculate SASA
        cmd.set("dot_solvent", 1)
        cmd.set("dot_density", 3)

        #filout.write("{},{:.3f}\n".format(p, float(cmd.get_area("act"+p))))
        filout.write("{},{:.3f}\n".format(p, float(cmd.get_area("resn NBX and "+p))))
    filout.close()

###########
# program #
###########

#if __name__ == "__main__": 
for pdb in pdb_list: 
    open_pdb(pdb)
