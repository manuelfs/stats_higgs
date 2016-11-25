#!/usr/bin/env python

###### Script to send 58 jobs to the batch producing workspaces for each signal mass point
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

# Setting folders
model = "TChiHH"
ntu_date = "2016_08_10"
lumi = "36.2"

lumi_s = lumi.replace(".","p")
infolder  = "/net/cms2/cms2r0/babymaker/babies/"+ntu_date+"/"+model+"/merged_higmc_higtight/"
outfolder = "/net/cms2/cms2r0/babymaker/wspaces/"+ntu_date+"/"+model+"_lumi"+lumi_s+"/" 
runfolder = outfolder+"run/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

#input datasets
inputfiles = [i for i in os.listdir(infolder) if "SMS" in i]

os.system("JobSetup.csh")
njobs = 18
files_job = (len(inputfiles)+njobs-1)/njobs
ifile = 0
ijob = 0
for file in inputfiles:
  ifile += 1
  # Creating executable
  if ifile % files_job == 1 or files_job == 1:
    ijob += 1
    exename = runfolder+"/wspace_sig_"+str(ijob)+".sh"
    fexe = open(exename,"w")
    os.system("chmod u+x "+exename)
    fexe.write("#!/bin/bash\n\n")
  fexe.write("./run/wspace_sig.exe -f "+infolder+'/'+file+' -o '+outfolder+' -l '+lumi+'\n')
  if ifile % files_job == 0 or ifile == len(inputfiles): 
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh "+exename
    #print cmd
    os.system(cmd)

print "\nSubmitted "+str(ifile)+" files in "+str(ijob)+" jobs. Output goes to "+outfolder+"\n"
sys.exit(0)
