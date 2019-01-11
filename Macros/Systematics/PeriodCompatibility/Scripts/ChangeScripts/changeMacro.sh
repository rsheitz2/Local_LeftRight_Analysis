#!/bin/bash

if [ $# -ne 12 ]; then
     echo "" 
     echo "Script updates input Macro.C with input paramters"
     echo ""
else
    changeFile=$1
    nBins=$2
    fitMrangeType=$3
    hbins=$4
    physBinned=$5
    process=$6
    lrMrange=$7
    fitMrange=$8
    binRange=$9
    whichFit=${10}
    production=${11}
    additionalCuts=${12}

    if [ ! -f $changeFile ]; then
	echo "$changeFile does not exist"
	exit 1
    fi
    
    sed -i.bak "s/const Int_t nBins =.*;/const Int_t nBins =${nBins};/" ${changeFile}
    sed -i.bak "s/TString fitMrangeType =\".*\";/TString fitMrangeType =\"${fitMrangeType}\";/" ${changeFile}
    sed -i.bak "s/Int_t hbins =.*;/Int_t hbins =${hbins};/" ${changeFile}
    sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${changeFile}
    sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
    sed -i.bak "s/TString lrMrange =.*;/TString lrMrange =\"${lrMrange}\";/" ${changeFile}
    sed -i.bak "s/TString fitMrange =.*;/TString fitMrange =\"${fitMrange}\";/" ${changeFile}
    sed -i.bak "s/TString binRange =.*;/TString binRange =\"${binRange}\";/" ${changeFile}
    sed -i.bak "s/TString whichFit =.*;/TString whichFit =\"${whichFit}\";/" ${changeFile}
    sed -i.bak "s/TString production =.*;/TString production =\"${production}\";/" ${changeFile}
    sed -i.bak "s/TString additionalCuts =.*;/TString additionalCuts =\"${additionalCuts}\";/" ${changeFile}
    
    sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
fi
