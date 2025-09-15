#!/bin/bash
rm queue.list$1
ls --color=no filelists/pythia_$1/ > queue.list$1
