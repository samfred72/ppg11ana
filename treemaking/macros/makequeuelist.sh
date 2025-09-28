#!/bin/bash
rm queue.list
touch queue.list
for file in `ls --color=no filelists/pythia_MB/`; do
  echo "MB $file" >> queue.list
done
for file in `ls --color=no filelists/pythia_Jet10/`; do
  echo "Jet10 $file" >> queue.list
done
for file in `ls --color=no filelists/pythia_Jet20/`; do
  echo "Jet20 $file" >> queue.list
done
for file in `ls --color=no filelists/pythia_Jet30/`; do
  echo "Jet30 $file" >> queue.list
done
