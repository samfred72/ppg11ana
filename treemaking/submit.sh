#!/bin/bash

mkdir -p /tmp/samfred
rm -r /sphenix/tg/tg01/jets/samfred/run25/pythia_MB_hadded
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/pythia_MB_hadded
rm -r /sphenix/tg/tg01/jets/samfred/run25/pythia_Jet10_hadded
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/pythia_Jet10_hadded
rm -r /sphenix/tg/tg01/jets/samfred/run25/pythia_Jet20_hadded
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/pythia_Jet20_hadded
rm -r /sphenix/tg/tg01/jets/samfred/run25/pythia_Jet30_hadded
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/pythia_Jet30_hadded
condor_submit condor.job
