#!/usr/bin/env python

###### Script to send 58 jobs to the batch finding the limit for each signal mass point
import os, sys, subprocess
import pprint
import glob
import json
import string
import time

# Setting folders
model = "TChiHH"
ntu_date = "2017_06_26"
lumi = "35.9"

lumi_s = lumi.replace(".","p")
infolder  = "/net/cms27/cms27r0/babymaker/wspaces/"+ntu_date+"/"+model+"_lumi"+lumi_s+"/" 
# infolder  = "/net/cms2/cms2r0/babymaker/datacards/"+ntu_date+"/"+model+"/" 
runfolder = "batch_"+model+"/" 
if not os.path.exists(runfolder):
  os.system("mkdir -p "+runfolder)

do_cards = (not infolder.find("datacards")==-1)
if do_cards: 
  release = 'CMSSW_7_4_7'
else:
  release = 'CMSSW_7_4_14'
#input datasets
inputfiles = [i for i in os.listdir(infolder) if "xsecNom" in i]
if do_cards:
  inputfiles = [i for i in os.listdir(infolder) if "datacard" in i]

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
    exename = runfolder+"/find_limit_sig_"+str(ijob)+".sh"
    fexe = open(exename,"w")
    os.system("chmod u+x "+exename)
    fexe.write("#!/bin/bash\n\n")
    fexe.write(". /cvmfs/cms.cern.ch/cmsset_default.sh \n")
    if do_cards:
      fexe.write("cd ~/cmssw/CMSSW_7_4_7/src/ \n")
    else:
      fexe.write("cd ~/cmssw/stats/CMSSW_7_4_14/src/ \n")
    fexe.write("eval `scramv1 runtime -sh` \n")
    fexe.write("cd ~/code/stats_higgs ; \n\n")
  if do_cards:
    fexe.write("./run/scan_point.exe --cards -f "+infolder+'/'+file+' >> txt/limits_'+model+'_'+str(ijob)+'.txt\n')
  else: 
    # fexe.write("./python/run_combine.py "+infolder+'/'+file+' --full_fit --overwrite\n')
    fexe.write("./run/scan_point.exe -f "+infolder+'/'+file+' >> txt/limits_'+model+'_'+str(ijob)+'.txt\n')
  if ifile % files_job == 0 or ifile == len(inputfiles): 
    fexe.close()
    cmd = "JobSubmit.csh ./run/wrapper.sh "+release+" ./"+exename
    #print cmd
    os.system(cmd)

print "\nSubmitted "+str(ifile)+" files in "+str(ijob)+" jobs\n"
sys.exit(0)
