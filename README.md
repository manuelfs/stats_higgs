stats_higgs
=========
Statistics tools for the Higgsino->HH->bbbb analysis. Base code by Adam Dishaw.

Initial setup
-------------
stats_higgs relies on having CMSSW_7_4_14 or greater with the ROOT6 version of the Higgs combine tool
installed in ~/cmssw. Only CMSSW_7_4_14 is needed to produce workspaces and compute limits and
significances. If one also wishes to extract fitted values for the signal strength or other model parameters,
CMSSW_7_1_5 with the ROOT5 version of the Higgs combine tool is also needed. To perform the full setup with
both CMSSW versions, simply execute

    ./run/initial_setup.sh

