#!/bin/bash
if [[ $2 == "" ]]; then
  echo "Needs a second arg for smear 1/0!"
  exit
fi

root -b -l -q "histmaking/draw_yield.C(\"${1}\",\"MC\",$2)"
root -b -l -q "histmaking/draw_yield.C(\"${1}\",\"data\")"
root -b -l -q "Corrections/drawEff.cpp(\"${1}\",$2)"
root -b -l -q "CrossSectionClosure/CompCrossData.C(\"${1}\")"
display CrossSectionClosure/plots/cross_data_${1}.pdf
