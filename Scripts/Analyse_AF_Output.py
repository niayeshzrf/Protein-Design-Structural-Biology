from pymol import cmd
import os
import shutil
import statistics
p_input = "./renamed_best10_RFdiff1st/"
af2_out = "./JulystartwRF_AtoT_out/"
af2_filt = "./AFold_filt/"


with open('metrics.csv','w') as infile:
    infile.write('pdbfile, average plddt, median plddt, minimum pdltts, backbone RMSD\n')

for pdbfile in os.listdir(p_input):
    for i in range(1,5):
        i=str(i)
        if os.path.exists(af2_out+pdbfile[:-4]+"_"+i+"_Af2.pdb"):
            plddts = []
            cmd.load(p_input+pdbfile,pdbfile[:-4])
            cmd.load(af2_out+pdbfile[:-4]+"_"+i+"_Af2.pdb",pdbfile[:-4]+"_"+i+"_Af2")
            bb_rmsd = cmd.align(pdbfile[:-4]+"_"+i+"_Af2 and name c+n+o+ca",pdbfile[:-4]+" and name c+n+o+ca",cycles=0)[0]
            if bb_rmsd < 3:
                with open(af2_out+pdbfile[:-4]+"_"+i+"_Af2.pdb",'r') as infile:
                    lines = infile.readlines()
                    for line in lines:
                        if len(line.split()) > 10 and line.split()[2] == 'CA':
                            plddts.append(float(line.split()[10]))
                avg = round(statistics.mean(plddts),2)
                median = round(statistics.median(plddts),2)
                avg_min = round(sum(sorted(plddts)[0:5])/5,2)
                print(avg, median, avg_min)
                if avg_min>45 and avg>80 and median>80:
                    shutil.copy(af2_out+pdbfile[:-4]+"_"+i+"_Af2.pdb",af2_filt+pdbfile[:-4]+"_"+i+"_Af2.pdb")
                    with open('metrics.csv','a') as infile:
                        infile.write(pdbfile[:-4]+"_"+i + ',' +  str(avg) + ',' + str(median) + ',' +str(avg_min) + ',' + str(round(bb_rmsd,2)) + '\n')

