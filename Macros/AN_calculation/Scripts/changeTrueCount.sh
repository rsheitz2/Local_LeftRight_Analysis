#!/bin/bash

if [ $# -ne 7 ]; then
     echo "" 
     echo "Script updates trueCount.C with input paramters"
     echo ""
     
else
    nBins=$1
    period=$2
    fitMrangeType=$3
    fitMmin=$4
    fitMmax=$5
    physBinned=$6
    process=$7
    
    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/trueCount.C
    sed -i.bak "s/const Int_t nBins =*.;/const Int_t nBins =${nBins};/" ${changeFile}
    sed -i.bak "s/TString period_Mtype =.*;/TString period_Mtype =\"${period}_${fitMrangeType}\";/" ${changeFile}
    sed -i.bak "s/Double_t Mmin =.*;/Double_t Mmin =${fitMmin};/" ${changeFile}
    sed -i.bak "s/Double_t Mmax =.*;/Double_t Mmax =${fitMmax};/" ${changeFile}
    sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${changeFile}
    sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
fi
