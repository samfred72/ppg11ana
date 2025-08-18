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

SECTION=$1
range=$2

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/install/

printenv 

echo "runnumber..."
echo ${SECTION}
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/pythia_MB_hadded/
hadd -f -k /sphenix/tg/tg01/jets/samfred/run25/pythia_MB_hadded/run28_${SECTION}.root /sphenix/tg/tg01/jets/samfred/run25/pythia_MB/outtree_queue_run28_v${range}.list

echo all done
echo "script done"
