#python2
#run to summarize ensemble refinement if started with "prepare_ensemble_ref-for_qb3cluster.py"
#parent_dir is hard-coded, copy from other script
#logfile_name is hard-coded, and will be the same for all of your individual runs

import string
#import re
import os

Tx_list = ["1"]
ptls_list = ["1.0", "0.9", "0.7", "0.6"]
tbath_offset_list = ["5"]


parent_dir = "./"



output_file = open(parent_dir+"review_ensemble_optimization.out", "w")

for Tx in Tx_list:
  for ptls in ptls_list:
    for tbath_offset in tbath_offset_list:

      parameters = "Tx = %s  ;  ptls = %s  ;  wxray_coupled_tbath_offset = %s \n" % (Tx, ptls, tbath_offset)

      output_file.write(parameters)

      dir_name = "ENSEMBLE_ref_tx_%s_ptls_%s_tbath_%s" % (Tx, ptls, tbath_offset)
      for filename in os.listdir(parent_dir+dir_name):
        if filename.endswith(".log"):
            logfile_name=filename
      logfile = open(parent_dir+dir_name+"/"+logfile_name, "r")
      for line in logfile.readlines():
        if line.startswith("FINAL Rwork = 0."):
          output_file.write(line)
        if line.startswith("Ensemble size :  "):
          output_file.write(line)

