#!/bin/bash

if [ $# -ne 14 ]; then
     echo "" 
     echo "Script updates pullDist.C with input paramters"
     echo ""
else
    nBins=$1
    fitMrangeType=$2
    hbins=$3
    process=$4
    lrMrange=$5
    fitMrange=$6
    
    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/pullDist.C
    sed -i.bak "s/const Int_t nBins =.*;/const Int_t nBins =${nBins};/" ${changeFile}
    sed -i.bak "s/TString fitMrangeType =\".*\";/TString fitMrangeType =\"${fitMrangeType}\";/" ${changeFile}
    sed -i.bak "s/Int_t hbins =.*;/Int_t hbins =${hbins};/" ${changeFile}
    sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
    sed -i.bak "s/TString lrMrange =.*;/TString lrMrange =\"${lrMrange}\";/" ${changeFile}
    sed -i.bak "s/TString fitMrange =.*;/TString fitMrange =\"${fitMrange}\";/" ${changeFile}

    #Assuming 4 physics binning variables at the moment
    sed -i.bak "s/TString physBinned\[nPhysBinned\] =.*;/TString physBinned\[nPhysBinned\] =\{\"${7}\", \"${8}\", \"${9}\", \"${10}\"\};/" ${changeFile}
    sed -i.bak "s/TString whichFit\[nPhysBinned\] =.*;/TString whichFit\[nPhysBinned\] =\{\"${11}\", \"${12}\", \"${13}\", \"${14}\"\};/" ${changeFile}
    
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
fi
