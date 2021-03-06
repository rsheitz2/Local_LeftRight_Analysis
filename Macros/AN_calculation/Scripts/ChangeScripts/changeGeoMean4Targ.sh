#!/bin/bash

if [ $# -ne 14 ]; then
     echo "" 
     echo "Script updates GeoMean4Targ.C with input paramters"
     echo ""
     
else
    nBins=$1
    period=$2
    fitMrangeType=$3
    hbins=$4
    physBinned=$5
    process=$6
    LR_Mmin=$7
    LR_Mmax=$8
    fitMmin=$9
    fitMmax=${10}
    whichFit=${11}
    binRange=${12}
    production=${13}
    additionalCuts=${14}
    
    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/GeoMean4Targ.C
    sed -i.bak "s/const Int_t nBins =.*;/const Int_t nBins =${nBins};/" ${changeFile}
    sed -i.bak "s/TString period_Mtype =.*;/TString period_Mtype =\"${period}_${fitMrangeType}\";/" ${changeFile}
    sed -i.bak "s/Int_t hbins =.*;/Int_t hbins =${hbins};/" ${changeFile}
    sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${changeFile}
    sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
    sed -i.bak "s/TString lrMrange =.*;/TString lrMrange =\"${LR_Mmin}_${LR_Mmax}\";/" ${changeFile}
    sed -i.bak "s/TString fitMrange =.*;/TString fitMrange =\"${fitMmin}_${fitMmax}\";/" ${changeFile}
    sed -i.bak "s/TString whichFit =.*;/TString whichFit =\"${whichFit}\";/" ${changeFile}
    sed -i.bak "s/TString binRange =.*;/TString binRange =\"${binRange}\";/" ${changeFile}
    sed -i.bak "s/TString production =.*;/TString production =\"${production}\";/" ${changeFile}
    sed -i.bak "s/TString additionalCuts =.*;/TString additionalCuts =\"${additionalCuts}\";/" ${changeFile}
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
fi
