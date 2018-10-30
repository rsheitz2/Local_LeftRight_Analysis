#!/bin/bash

if [ $# -ne 4 ]; then
     echo "" 
     echo "Script updates input python script with input paramters"
     echo ""
     exit 1
fi     

changeFile=$1
mass_type=$2
production=$3
which_cuts=$4

if [ ! -f ${changeFile} ]; then #Check full path entered

    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Cuts/${changeFile}
    if [ ! -f ${changeFile} ]; then #Try looking for file in local path
	echo "${changeFile}  File does not exist"
	exit 1
    fi
fi

sed -i.bak "5,15s/mass_type=\".*\"/mass_type=\"${mass_type}\"/" ${changeFile}
sed -i.bak "5,15s/production=\".*\"/production=\"${production}\"/" ${changeFile}
sed -i.bak "5,15s/which_cuts=\".*\"/which_cuts=\"${which_cuts}\"/" ${changeFile}
