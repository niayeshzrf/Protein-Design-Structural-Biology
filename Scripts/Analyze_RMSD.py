#!/usr/bin/env python

import os
import RMSD
import csv
import itertools
import time
start_time = time.time()
directory = "."

templates = []
for filename in os.listdir(directory):
    if "model" in filename:
        templates.append(filename)
    if "xtal" in filename:
        reference = filename

def deviation(templates, reference):
    RMSDs_wrt_xtal = []
    for structure in templates:
        rmsd = RMSD.main(structure, reference)
        RMSDs_wrt_xtal.append(rmsd)
    return RMSDs_wrt_xtal
        
def write_rmsd(list_res, Outputfilename):
    with open(Outputfilename+'.csv','wb') as myfile:
        wr = csv.writer(myfile, quoting = csv.QUOTE_ALL)
        for value in list_res:
            wr.writerow([value])

def average(lst):
    return sum(lst)/len(lst)

def diversity(templates):
    all_iterations = []
    pair_rmsds = []
    for subset in itertools.combinations(templates, 2):
        all_iterations.append(subset)

    for pair in all_iterations:
        pair_rmsd = RMSD.main(pair[0], pair[1])
        pair_rmsds.append(pair_rmsd)
    return pair_rmsds
#Execute------------------------------------------------------
RMSDs_wrt_xtal = deviation(templates, reference)
pair_rmsds = diversity(templates)
write_rmsd(pair_rmsds, "diversity")
write_rmsd(RMSDs_wrt_xtal, "deviation")
mean_deviation = average(RMSDs_wrt_xtal)
mean_diversity = average(pair_rmsds)
print "Diversity:%f" %mean_diversity
print "Deviation:%f" %mean_deviation

with open("Results.csv",'wb') as myfile:
    myfile.write("Deviation,Diversity")
    myfile.write(str(mean_deviation + ',' + mean_diversity))



print time.time() - start_time
