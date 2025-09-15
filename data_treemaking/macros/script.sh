#!/bin/bash
export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

LIST=$1
TEST=0
[ -n "$2" ] && TEST=$2

#if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
#  cd $_CONDOR_SCRATCH_DIR
#else
#  echo condor scratch NOT set
#  exit -1
#fi
#FULLPATH=`psql FileCatalog -t -c "select full_file_path from files where lfn = '${1}';"`

hostname

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
#echo rsyncing from $this_dir
echo running: $this_script $*

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/samfred/run25/ppg11/data_treemaking/install/

#printenv 
LIST=$1
echo "input files..."
root -l -q -b "Fun4All_macro.C(\"${LIST}\")"

echo all done
echo "script done"
