#parent_dir, pdb_filename, and mtz_filename are hard-coded and will need to be edited below
#this script sets up all the required files, etc. and produces a shell script
#run the shell script to submit jobs to the cluster

import string
import re
import os

Tx_list = ["0.5", "1", "1.5", "2.0"]
ptls_list = ["1.0", "0.9", "0.8", "0.7", "0.6"]
tbath_offset_list = ["2.5", "5.0", "7.5", "10.0"]


parent_dir = os.getcwd()+"/"
pdb_filename = "PDBFile_Name.updated.pdb"
mtz_filename = "PDBFile_Name.mtz"


directories = []
qsub_lines = []

for Tx in Tx_list:
  for ptls in ptls_list:
    for tbath_offset in tbath_offset_list:

      dir_name = "ENSEMBLE_ref_tx_%s_ptls_%s_tbath_%s" % (Tx, ptls, tbath_offset)
      os.mkdir(parent_dir+dir_name)
      os.mkdir(parent_dir+dir_name+"/err")
      
      os.system("cp ./"+pdb_filename+" ./"+dir_name)
      os.system("cp ./"+mtz_filename+" ./"+dir_name)


      sub_filename = parent_dir+dir_name+"/ENSEMBLE_ref_tx_%s_ptls_%s_tbath_%s.sh" % (Tx, ptls, tbath_offset)
      sub_file = open(sub_filename, "w")

      sub_file_text = """#!/bin/bash
#SBATCH --account=def-rchica
#SBATCH --mail-user=niayeshzarifi@outlook.com
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name=PDBFile_Name
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=1000M
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
/home/niayesh/Apps/phenix/phenix-1.15.2-3472/build/bin/phenix.ensemble_refinement %s%s/%s %s%s/%s tx=%s ptls=%s wxray_coupled_tbath_offset=%s 

""" % (parent_dir, dir_name, pdb_filename, parent_dir, dir_name, mtz_filename, Tx, ptls, tbath_offset)

      sub_file.write(sub_file_text)

      qsub_line = "sbatch "+sub_filename
      qsub_lines.append(qsub_line)
      
      directories.append(parent_dir+dir_name)      

shell_script = open("submit_ensemble_ref_to_cluster.sh", "w")
shell_script.write("#!/bin/bash\n\n")
for i in range(len(qsub_lines)):
  shell_script.write("cd "+directories[i-1]+"\n")
  shell_script.write(qsub_lines[i-1]+"\n")
shell_script.write("cd ..")

os.system("chmod 777 submit_ensemble_ref_to_cluster.sh")
