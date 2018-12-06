#!/bin/bash

if [ $# -ne 6 ]; then
     echo "" 
     echo "Script updates Macro.C with input paramters"
     echo ""
     
else
    changeFile=$1
    phiScut=$2
    A_siv=$3
    additionalCuts=$4
    N_gen=$5
    alphaScale=$6
    
    if [ ! -f $changeFile ]; then
	echo "$changeFile"
	echo " does not exist"
	exit 1
    fi

    sed -i.bak "s/Double_t phiScut =.*;/Double_t phiScut =${phiScut};/" ${changeFile}
    sed -i.bak "s/Double_t A_siv =.*;/Double_t A_siv =${A_siv};/" ${changeFile}
    sed -i.bak "s/TString additionalCuts =.*;/TString additionalCuts =\"${additionalCuts}\";/" ${changeFile}
    sed -i.bak "s/Int_t N_gen =.*;/Int_t N_gen =${N_gen};/" ${changeFile}    
    sed -i.bak "s/Bool_t alphaScale =.*;/Bool_t alphaScale =${alphaScale};/" ${changeFile}
	
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}

fi
