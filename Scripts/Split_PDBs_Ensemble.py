#! /usr/bin/env python

#This script extracts an ensemble model into separate PDB files

import os
import sys
from stat import *
import logging
import argparse
from argparse import RawTextHelpFormatter

def generate_output_files_prefix(file_name):
	main_part_of_name, file_extension = os.path.splitext(file_name)
	if '.' in file_name:
		return main_part_of_name
	else:
		return file_name

def extract_models(ensemble_model):
	generated_output_files_prefix = generate_output_files_prefix(ensemble_model)
	
	the_multi_file_stream = open(ensemble_model, "r")
	model_number = 1
	new_file_text = ""
	for line in the_multi_file_stream:
		line = line.strip() 
		if line == "ENDMDL":
			output_file = open(generated_output_files_prefix + "_model_" + str(
			    model_number) + ".pdb", "w")
			output_file.write(new_file_text.rstrip('\r\n'))
			output_file.close()
			model_number += 1
			new_file_text = ""
		elif not line.startswith("MODEL"):
			new_file_text += line + '\n'
	
	the_multi_file_stream.close()
	model_number -= 1
	return (model_number, generated_output_files_prefix)


def THE_MAIN_FUNCTION(File_or_Directory):
	ensemble_model = File_or_Directory
	number_of_models = 0

	sys.stderr.write("Reading in your file...")
	number_of_models, output_files_prefix = (
		extract_models(ensemble_model))
	sys.stderr.write("\nConcluded. \n")
        if number_of_models < 2:
        	sys.stderr.write("Sorry. Only one model was recognized in the file. \
            		Please examine your file for the presence of multiple 'MODEL' and 'ENDMDL' indicators.")
        	sys.stderr.write("A file named '"+ output_files_prefix+"_model_1.pdb' was made in the process. ")
    	else:
        	sys.stderr.write("File split into "+ str(number_of_models)+" models. ")
        	sys.stderr.write("\nFiles with names '"+
            		output_files_prefix+"_model_1.pdb', '"+
            		output_files_prefix+"_model_2.pdb', etc., "+
            		"\nhave been created in same directory as the input file.\n\n")

parser = argparse.ArgumentParser(prog='split_PDBs_Ensemble.py', description="This script extracts an ensemble model into separate PDB files."
,formatter_class=RawTextHelpFormatter)

parser.add_argument("File_or_Directory", help="name of the ensemble model required")

if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

if os.path.isfile(args.File_or_Directory):
    	THE_MAIN_FUNCTION(args.File_or_Directory)
elif os.path.isdir(args.File_or_Directory):
    	for f in os.listdir(args.File_or_Directory):
        	pathname = os.path.join(args.File_or_Directory, f)
        	mode = os.stat(pathname).st_mode #need to import stat module
        	if S_ISREG(mode) and pathname[-4:].upper() == ".PDB":
            	# It's a PDB file, call the function
            		THE_MAIN_FUNCTION(pathname)
else:
	sys.stderr.write("SORRY. " + args.File_or_Directory + " IS NOT RECOGNIZED AS A PDB FILE.\n\n")
    	parser.print_help()
    	sys.exit(1)

