# script by Broom to measure deviation and diversity for some structures

import subprocess
import os
import sys

# function to generate the VMD TCL script
def GenerateVMDScript(Commands):

	VMDScript = ""

	# write the alignment function and do the alignment
	VMDScript += "proc AlignBackbone {numberOfFrames} {\n"
	VMDScript += "\tset reference [atomselect 0 \"name C CA N O\"]\n"
	VMDScript += "\tfor {set i 0} {$i < $numberOfFrames} {incr i} {\n"
	VMDScript += "\t\tset compare [atomselect 1 \"name C CA N O\" frame $i]\n"
	VMDScript += "\t\tset compareAll [atomselect 1 \"name C CA N O\" frame $i]\n"
	VMDScript += "\t\tset transMatrix [measure fit $compare $reference]\n"
	VMDScript += "\t\t$compareAll move $transMatrix\n"
	VMDScript += "\t}\n"
	VMDScript += "}\n"

	# write the main loop
	VMDScript += "set numberOfFrames [molinfo 1 get numframes]\n"
	VMDScript += "AlignBackbone $numberOfFrames\n"
	VMDScript += "set deviation 0.0\n"
	VMDScript += "set diversity 0.0\n"
	VMDScript += "set diversityCount 0\n"
	VMDScript += "for {set i 0} {$i < $numberOfFrames} {incr i} {\n"
	VMDScript += "\tset reference [atomselect 0 \"name C CA N O\"]\n"
	VMDScript += "\tset compareOne [atomselect 1 \"name C CA N O\" frame $i]\n"
	VMDScript += "\tset singleDeviation [measure rmsd $reference $compareOne]\n"
	VMDScript += "\tset deviation [expr $deviation + $singleDeviation]\n"
	VMDScript += "\tfor {set j [expr $i + 1]} {$j < $numberOfFrames} {incr j} {\n"
	VMDScript += "\t\tset compareTwo [atomselect 1 \"name C CA N O\" frame $j]\n"
	VMDScript += "\t\tset diversityCount [expr $diversityCount + 1]\n"
	VMDScript += "\t\tset singleDiversity [measure rmsd $compareOne $compareTwo]\n"
	VMDScript += "\t\tset diversity [expr $diversity + $singleDiversity]\n"
	VMDScript += "\t}\n"
	VMDScript += "}\n"

	# compute the final values
	VMDScript += "set deviation [expr $deviation / $numberOfFrames]\n"
	VMDScript += "set diversity [expr $diversity / $diversityCount]\n"

	# write the output
	VMDScript += "set OutFile [open RMSD_Output.dat w];\n"
	VMDScript += "puts $OutFile $numberOfFrames\n"
	VMDScript += "puts $OutFile $deviation\n"
	VMDScript += "puts $OutFile $diversity\n"
	VMDScript += "close $OutFile\n"

	return VMDScript


# get the command-line arguments
Commands = {}
for arg in sys.argv:
	arg = arg.split("=")
	if len(arg) == 2:
		Commands[arg[0]] = arg[1]

# if command-line arguments were not sensible report a brief help description
if ("ref" not in Commands) or ("dir" not in Commands):
	print "You must provide a reference and ensemble directory:"
	print "python GetDeviationAndDiversity.py ref=reference.pdb dir=/home/chica/ensemble/"

# generate the VMD TCL script
VMDScript = GenerateVMDScript(Commands)
VMDFile = open("GetDeviationAndDiversity.tcl", "w")
VMDFile.write(VMDScript)
VMDFile.close()

# call VMD
callLine = "vmd -dispdev text -f " + Commands["ref"] + " -f " + Commands["dir"].rstrip("/") + "/*.pdb -eofexit < GetDeviationAndDiversity.tcl"
subprocess.check_output(callLine, shell=True)

# load the metrics
outputFile = open("RMSD_Output.dat", "r")
ensembleSize = int(outputFile.readline())
deviation = float(outputFile.readline())
diversity = float(outputFile.readline())	
outputFile.close()

# report the metrics
print "\nEnsemble Size:", ensembleSize
print "Deviation:", deviation, "angstroms"
print "Diversity:", diversity, "angstroms\n"

# clean up
os.remove("GetDeviationAndDiversity.tcl")
os.remove("RMSD_Output.dat")
