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

SECTION="_${1}"
range=$2
TRIGGER="MB"
if [[ $SECTION == "_-1" ]]; then
  TRIGGER=$2
  SECTION=""
  range="*"
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/install/

printenv 

echo "inputs..."
echo ${TRIGGER}${SECTION}
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/pythia_${TRIGGER}_hadded/
hadd -f -k /sphenix/tg/tg01/jets/samfred/run25/pythia_${TRIGGER}_hadded/run28${SECTION}.root /sphenix/tg/tg01/jets/samfred/run25/pythia_$TRIGGER/outtree_queue_run28_v${range}.list

echo all done
echo "script done"
