#!/bin/bash
rm log/*
bash submit_data.sh
for item in `cat runlist.list`; do
  while ! grep -q "script done" "log/${item}.out" 2>&1; do
    output=`tail -n 1 "log/${item}.out"`
    echo "Waiting on segment ${item}... ${output}" 
    sleep 5
  done
  echo "Segment $item complete"
done
rm data_hists/mass_data.root
hadd data_hists/mass_data.root data_hists/*
