#!/bin/bash 

DIRECTORY=`pwd`
RELEASE=$1
SCRIPT=$2

cd /homes/aovcharova/cmssw/${RELEASE}/src/
. /net/cms2/cms2r0/babymaker/cmsset_default.sh
eval `scramv1 runtime -sh`
cd $DIRECTORY;

# Execute command
$SCRIPT