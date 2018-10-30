#!/bin/bash

if [ $# -ne 3 ]; then
     echo "" 
     echo "Script updates macroPeriodsLoop.sh with input paramters"
     echo ""
     exit 1
fi     

Mtype=$1
production=$2
whichCuts=$3

changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Cuts/Scripts/LoopScripts/macroPeriodsLoop.sh
sed -i.bak "20,30s/Mtype=.*/Mtype=\"${Mtype}\"/" ${changeFile}
sed -i.bak "20,30s/production=.*/production=\"${production}\"/" ${changeFile}
sed -i.bak "20,30s/whichCuts=.*/whichCuts=\"${whichCuts}\"/" ${changeFile}
