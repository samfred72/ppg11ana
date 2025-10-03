#!/bin/bash
rm log/*MB*
bash submit.sh
for item in {0..19}; do
  while ! grep -q "script done" "log/${item}MB.out" 2>&1; do
    echo "Waiting on segment $item..." 
    tail -n 1 "log/${item}MB.out"
    sleep 5
  done
  echo "Segment $item complete"
done
bash dohadd.sh
#root alphafitmass.C
