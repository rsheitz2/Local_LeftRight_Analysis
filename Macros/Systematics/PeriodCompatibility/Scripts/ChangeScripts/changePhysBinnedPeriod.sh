#!/bin/bash

if [ $# -ne 8 ]; then
     echo "" 
     echo "Script updates physBinnedPeriod.C with input paramters"
     echo ""
else
    nBins=$1
    fitMrangeType=$2
    hbins=$3
    physBinned=$4
    process=$5
    lrMrange=$6
    fitMrange=$7
    whichFit=$8

    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/physBinnedPeriod.C
    sed -i.bak "s/const Int_t nBins =.*;/const Int_t nBins =${nBins};/" ${changeFile}
    sed -i.bak "s/TString fitMrangeType =\".*\";/TString fitMrangeType =\"${fitMrangeType}\";/" ${changeFile}
    sed -i.bak "s/Int_t hbins =.*;/Int_t hbins =${hbins};/" ${changeFile}
    sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${changeFile}
    sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
    sed -i.bak "s/TString lrMrange =.*;/TString lrMrange =\"${lrMrange}\";/" ${changeFile}
    sed -i.bak "s/TString fitMrange =.*;/TString fitMrange =\"${fitMrange}\";/" ${changeFile}
    sed -i.bak "s/TString whichFit =.*;/TString whichFit =\"${whichFit}\";/" ${changeFile}
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}

fi
