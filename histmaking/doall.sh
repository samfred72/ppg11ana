#!/bin/bash
rm log/*MB*
bash submit.sh
for item in {0..9}; do
  while ! grep -q "script done" "log/${item}MB.out" 2>&1; do
    echo "Waiting on segment $item..." 
    sleep 5
  done
  echo "Segment $item complete"
done
bash dohadd.sh
root -b -l -q "onefitmass.C(\"pi0\",\"MC\",\"MB\",2)"
root -b -l -q "onefitmass.C(\"pi0\",\"MC\",\"MB\",3)"
root -b -l -q "onefitmass.C(\"pi0\",\"MC\",\"MB\",4)"
root -b -l -q "onefitmass.C(\"pi0\",\"MC\",\"MB\",5)"
root -b -l -q "onefitmass.C(\"pi0\",\"MC\",\"MB\",6)"
root -b -l -q "onefitmass.C(\"pi0\",\"MC\",\"MB\",7)"
cd ..
bash plotdata.sh pi0
bash plotcross.sh pi0 1
cd --
