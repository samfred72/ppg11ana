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

TRIGGER=$1
LIST=$2
ISTEST=0
[[ -n $3 ]] && ISTEST=$3


source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/ppg11/treemaking/install/

#printenv 

echo "input files..."
echo ${TRIGGER} ${LIST}
root -l -q -b "Fun4All_macro_mc.C(\"${LIST}\",\"${TRIGGER}\",$ISTEST)"

echo all done
echo "script done"
