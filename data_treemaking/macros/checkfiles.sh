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

while getopts "hr:l:p:c:d:D:" flag; do
  case $flag in
    h) show_help ;;
    r) runnumber=${OPTARG} ;;
    l) runlist=${OPTARG};;
    p) prod=${OPTARG} ;;
    c) cdb=${OPTARG};;
    d) dst=${OPTARG};;
    D) DIR=${OPTARG};;
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
dir=""
if [[ $dst == "DST_CALO" ]]; then
  dir="calocalib"
elif [[ $dst == "DST_CALOFITTING" ]]; then
  dir="calofitting"
else
  echo
  echo "Wrong DST type. Figure it out..."
  show_help
  exit
fi
echo "Files incoming..."
if [[ -n $runnumber ]]; then
  bash filelists/createdstlist.sh -r $runnumber -p $prod -c $cdb -d $dst
elif
  [[ -n $runlist ]]; then
  bash filelists/createdstlist.sh -l $runlist -p $prod -c $cdb -d $dst
fi
