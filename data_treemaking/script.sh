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

RUNNUMBER=$1

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/ppg11/data_treemaking/install/

printenv 

echo "runnumber..."
echo ${RUNNUMBER}
hadd -f -k /sphenix/tg/tg01/jets/samfred/run25/ppg11ana_hadded/run${RUNNUMBER}.root /sphenix/tg/tg01/jets/samfred/run25/ppg11ana/*${RUNNUMBER}*

echo all done
echo "script done"
