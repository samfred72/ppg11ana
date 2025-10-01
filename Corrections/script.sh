#!/bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

hostname

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo rsyncing from $this_dir
echo running: $this_script $*

PARTICLE=$1
TYPE=$2
SEG=$3

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/install/

#printenv 

echo "particle..."
echo ${PARTICLE}
echo "type..."
echo $TYPE
echo "segment..."
echo $SEG
root -b -l -q "MakeTruthHist.cpp(\"${PARTICLE}\",\"${TYPE}\",${SEG})"

echo all done
echo "script done"
