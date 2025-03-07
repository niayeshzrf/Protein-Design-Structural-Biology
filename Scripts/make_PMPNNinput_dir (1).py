import os
import shutil
for direct in os.listdir("best10_RFdiff1st"):
    os.mkdir("P_Input/"+direct[:-4])
    shutil.copy("best10_RFdiff1st/"+direct, "P_Input/"+direct[:-4]+"/"+direct)
