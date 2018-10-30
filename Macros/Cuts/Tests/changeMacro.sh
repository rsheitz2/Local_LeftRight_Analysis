#!/bin/bash

if [ $# -ne 5 ]; then
    echo "" 
    echo "Script updates input macro name with input paramters"
    echo "   This is the most generic change.sh file specific for macros only"
    echo ""
    echo "Enter full path or macro name in cuts path "
    echo "   of macro to change"
    echo ""
    exit 1
fi

changeFile=$1
period=$2
Mtype=$3
production=$4
whichCuts=$5


if [ ! -f ${changeFile} ]; then #Full path entered

    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Cuts/${changeFile}
    if [ ! -f ${changeFile} ]; then #Try looking for file in local path
	echo "${changeFile}  File does not exist"
	exit 1
    fi
fi

sed -i.bak "s/TString period =.*;/TString period =\"${period}\";/" ${changeFile}
sed -i.bak "s/TString Mtype =.*;/TString Mtype =\"${Mtype}\";/" ${changeFile}
sed -i.bak "s/TString production =.*;/TString production =\"${production}\";/" ${changeFile}
sed -i.bak "s/TString whichCuts =.*;/TString whichCuts =\"${whichCuts}\";/" ${changeFile}

sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
