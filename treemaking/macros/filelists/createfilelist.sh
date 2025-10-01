#!/bin/bash
CreateFileList.pl -run 28 -type 26 -nop DST_MBD_EPD DST_CALO_CLUSTER DST_GLOBAL G4Hits
bash combine_dstlists.sh
bash split.sh "MB"
CreateFileList.pl -run 28 -type 12 -nop DST_MBD_EPD DST_CALO_CLUSTER DST_GLOBAL G4Hits
bash combine_dstlists.sh
bash split.sh "Jet10"
CreateFileList.pl -run 28 -type 21 -nop DST_MBD_EPD DST_CALO_CLUSTER DST_GLOBAL G4Hits
bash combine_dstlists.sh
bash split.sh "Jet20"
CreateFileList.pl -run 28 -type 11 -nop DST_MBD_EPD DST_CALO_CLUSTER DST_GLOBAL G4Hits
bash combine_dstlists.sh
bash split.sh "Jet30"
CreateFileList.pl -run 28 -type 36 -nop DST_MBD_EPD DST_CALO_CLUSTER DST_GLOBAL G4Hits
bash combine_dstlists.sh
bash split.sh "Jet5"
