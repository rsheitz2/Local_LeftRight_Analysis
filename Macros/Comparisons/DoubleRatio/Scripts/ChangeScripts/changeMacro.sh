#!/bin/bash

if [ $# -ne 8 ]; then
     echo "" 
     echo "Script updates input script with input parameters"
     echo ""
     
else
    changeFile=$1
    nBins=$2
    nHbins=$3
    period=$4
    Mtype=$5
    physBinned=$6
    production=$7
    whichTSA=$8

    if [ ! -f ${changeFile} ]; then
	echo "$changeFile  does not exist"
	exit 1
    fi
    
    sed -i.bak "s/const Int_t nBins =*.;/const Int_t nBins =${nBins};/" ${changeFile}
    sed -i.bak "s/const Int_t nHbins =*.;/const Int_t nHbins =${nHbins};/" ${changeFile}
    sed -i.bak "s/TString period =.*;/TString period =\"${period}\";/" ${changeFile}
    sed -i.bak "s/TString Mtype =.*;/TString Mtype =\"${Mtype}\";/" ${changeFile}
    sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${changeFile}
    sed -i.bak "s/TString production =.*;/TString production =\"${production}\";/" ${changeFile}
    sed -i.bak "s/TString whichTSA =.*;/TString whichTSA =\"${whichTSA}\";/" ${changeFile}
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
fi
