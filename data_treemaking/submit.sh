#!/bin/bash

mkdir -p /tmp/samfred
mkdir -p /sphenix/tg/tg01/jets/samfred/run25/ppg11ana_hadded
rm /sphenix/tg/tg01/jets/samfred/run25/ppg11ana_hadded/*
condor_submit condor.job
