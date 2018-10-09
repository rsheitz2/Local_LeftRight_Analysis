#!/bin/bash

if [ $# -ne 11 ]; then
    echo "" 
    echo "Script updates fullTargSysFunctMFit.C with input paramters"
    echo ""
    exit 1
fi

nBins=$1
period_Mtype=$2
hbins=$3
physBinned=$4
process=$5
LR_Mmin=$6
LR_Mmax=$7
Mmin=$8
Mmax=$9
whichFit=${10}
binRange=${11}

changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/fullTargSysFunctMFit.C
sed -i.bak "s/const Int_t nBins =.*;/const Int_t nBins =${nBins};/" ${changeFile}
sed -i.bak "s/TString period_Mtype =.*;/TString period_Mtype =\"${period_Mtype}\";/" ${changeFile}
sed -i.bak "s/Int_t hbins =.*;/Int_t hbins =${hbins};/" ${changeFile}
sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${changeFile}
sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
sed -i.bak "s/Double_t LR_Mmin =.*;/Double_t LR_Mmin =${LR_Mmin};/" ${changeFile}
sed -i.bak "s/Double_t LR_Mmax =.*;/Double_t LR_Mmax =${LR_Mmax};/" ${changeFile}
sed -i.bak "s/Double_t Mmin =.*;/Double_t Mmin =${Mmin};/" ${changeFile}
sed -i.bak "s/Double_t Mmax =.*;/Double_t Mmax =${Mmax};/" ${changeFile}
sed -i.bak "s/TString whichFit =.*;/TString whichFit =\"${whichFit}\";/" ${changeFile}
sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
sed -i.bak "s/TString binRange =.*;/TString binRange =\"${binRange}\";/" ${changeFile}
