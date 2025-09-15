#/bin/bash
runnumber=""
runlist=""
prod=""
cdb=""
dst=""

show_help() {

  echo
  echo OPTIONS:
  echo
  echo -r     runnumber \(must choose runnumber xor runlist\)
  echo -l     runlist \(must choose runnumber xor runlist\)
  echo -p     production tag, e.g. new \(required\)
  echo -c     cdb tag, e.g. newcdbtag_v005 \(required\)
  echo -d     DST type, e.g. DST_CALO \(required\)
  echo -D     Target directory, e.g. filelists/initialruns \(required\)
  echo
  exit
}

while getopts "hr:l:" flag; do
  case $flag in
    h) show_help ;;
    r) runnumber=${OPTARG} ;;
    l) runlist=${OPTARG};;
  esac
done

if [[ -z $runnumber && -z $runlist ]]; then
  echo
  echo "Must supply either run number or run list!"
  show_help
  exit
fi
if [[ -n $runnumber && -n $runlist ]]; then
  echo
  echo "Cannot supply both run number or run list!"
  show_help
  exit
fi

runs=$runnumber
splitnumber=1
if [[ -n $runlist ]]; then
  runs=`cat $runlist`
  splitnumber=5
fi

rm queue.list
touch queue.list

for run in $runs; do
  bash filelists/split.sh $run $splitnumber
  ls filelists/queue/*$run*  | xargs -n 1 basename >> queue.list
done
