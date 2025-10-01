#!/bin/bash
run=$1  
split -l $2 filelists/dst_calofitting-000${run}.list temp_
n=0
for file in temp_*; do
  m=`printf "%04d\n" $n`
  mv "$file" "filelists/queue/queue_${run}_v${m}.list"
  ((n++))
done
