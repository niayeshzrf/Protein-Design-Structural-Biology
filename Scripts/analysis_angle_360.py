#!/usr/bin/env python 

import sys 
import xlsxwriter
import os 

# write one table in one excel worksheet
# format of a table: structure chi1 chi2 d.chi1 d.chi2 
# you open each file and from each file 
# you should have a list of your chi angles and the protein 
# in each list: [structure, chi1, chi2]

# usage: ./analysis_angles.py angles_directory"

def chi_list(directory, x): 
    # x should contain a list of the file 
    line = []
    chi = []

    for pdb in x: 
        # this loop will iterate over the pdb files in the file_list 
        fo = open(directory+"/"+pdb, "r")
        line += [fo.readline().split()]
        fo.close()
    return(line) 

def reference(directory, ref): 
    # this is aimed to write the reference 
    line = []
    new_chi = []
    fo = open(directory+"/"+ref, "r")
    line += fo.readline().split()

    return(line)

def delta_chi(chi_list, reference): 
    delta_chi1 = 0.0
    delta_chi2 = 0.0
    delta_list = []
    # find a solution to deal with negative values
    for chi in chi_list : 
        chi1 = 0.0
        chi2 = 0.0 
        reference1 = 0.0
        reference2 = 0.0
        if len(chi) >= 3:   
            if float(chi[1]) < 0.0 : 
                chi1 = 360 + float(chi[1])
            elif float(chi[2]) < 0.0 : 
                chi2 = 360 + float(chi[2])
            elif float(reference[1]) < 0.0:
                reference1 = 360 + float(reference[1])
            elif float(reference[2]) < 0.0:
                reference2 = 360 + float(reference[2])
            else : 
                chi1 = chi[1]
                chi2 = chi[2]
                reference1 = reference[1]
                reference2 = reference[2]

            delta_chi1 = float(chi1) - float(reference1)
            delta_chi2 = float(chi2) - float(reference2)
            chi += [delta_chi1]  
            chi += [delta_chi2]

            # append new chi angles in a list 
            delta_list += [chi]

        # else in case the residue only have one chi angle     
        else: 
            delta_chi1 = float(chi[1]) - float(reference[1])
            chi += [delta_chi1]  
       
            # append new chi angles in a list 
            delta_list += [chi]
    return(delta_list)

# write in an excel workbook and in a woksheet
workbook = xlsxwriter.Workbook("../../HG-series/chi_angles/chi_angles_results.xlsx")

def write_in_excel(directory, residue_dir, x):

    # it needs the work directory, the residue directory 
    # and the file list from the directory 
    file_list = os.listdir(directory+"/"+residue_dir)

    # call the other functions
    p = chi_list(directory+"/"+residue_dir, file_list) 
    r = reference(directory+"/"+residue_dir, "HG4_xtal.dat")
    d = delta_chi(p, r) 

    # Here we want to create a worksheet for each residue  
    x = workbook.add_worksheet(residue_dir)

    if len(d[0]) > 3:

        row = 1
        col = 0

        # write the headers 
        x.write("A1", "Template") 
        x.write("B1", "chi1") 
        x.write("C1", "chi2") 
        x.write("D1", "delta_chi1") 
        x.write("E1", "delta_chi2") 

        # Iterate over the data in d 
        for structure, chi1, chi2, delta_chi1, delta_chi2 in d: 
            x.write(row, col, structure) 
            x.write(row, col + 1, float(chi1))
            x.write(row, col + 2, float(chi2)) 
            x.write(row, col + 3, delta_chi1)
            x.write(row, col + 4, delta_chi2) 
            row += 1

    else: 

        row = 1
        col = 0

        # write the headers 
        x.write("A1", "Template") 
        x.write("B1", "chi1") 
        x.write("C1", "delta_chi1") 

        # Iterate over the data in d 
        for structure, chi1, delta_chi1, in d: 
            x.write(row, col, structure) 
            x.write(row, col + 1, float(chi1))
            x.write(row, col + 2, delta_chi1)
            row += 1


# input for the functions directory containing all the chi angles   
directory = sys.argv[1]

# call the funtion 
list_residue = ["44", "50", "84", "127", "172", "236", "237", 
                "265", "267"]

for l in list_residue: 
    write_in_excel(directory, "residue_"+l, "worksheet")

# close the excel sheet 
workbook.close()


