#!/bin/bash
name=$1
particle=$2
type=$3
trigger=$4
pt=$5
val=$6
line=`grep -n "$name\s*\[$particle\]\[$type\]\[$trigger\]\[$pt\]" onefitmass.h | cut -d : -f 1`
echo "editing line `cat onefitmass.h | head -n $line | tail -n 1`"
sed -i -E   "s/$name\s*\[$particle\]\[$type\]\[$trigger\]\[$pt\] = -?0?.?(\w+);/$name \[$particle\]\[$type\]\[$trigger\]\[$pt\] = $val;/g" onefitmass.h
echo "line now     `cat onefitmass.h | head -n $line | tail -n 1`"
