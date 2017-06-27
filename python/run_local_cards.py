#! /usr/bin/env python
import os, sys


masses = ['127']
masses.extend([str(i) for i in range(150,1001,25)])

os.system("./compile.sh")
for mass in masses:
    cardname = "datacard_SMS-TChiHH_mGluino-"+mass+"_mLSP-1_35p9ifb.txt"
    print "Processing", mass
    cmd = "./run/scan_point.exe -f cards/"+cardname + " >> limits_cards.txt"
    os.system(cmd)
