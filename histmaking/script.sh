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
TRIGGER=$2

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/install/

printenv 

echo "runnumber..."
echo ${SECTION}
echo "trigger..."
echo $TRIGGER
root "truthhistmaker.C(${SECTION},\"${TRIGGER}\")"

echo all done
echo "script done"
