#!/bin/bash
dir=$1
mkdir -p $dir
split -l 10 combined_dstlist.list temp_
n=0
for file in temp_*; do
  m=`printf "%05d\n" $n`
  mv "$file" "$dir/queue_run28_v${m}.list"
  ((n++))
done
