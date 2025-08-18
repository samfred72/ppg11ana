#!/bin/bash
rm queue.list
ls --color=no filelists/pythia_$1/ > queue.list$1
