#!/bin/bash

if [ $# -ne 10 ]; then
    echo "" 
    echo "Script updates allSysErrorFA.C with input paramters"
    echo ""
    exit 1
fi

nBins=$1
period_Mtype=$2
hbins=$3
process=$4
lrMrange=$5
fitMrange=$6

changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/allSysErrorFA.C
sed -i.bak "s/const Int_t nBins =.*;/const Int_t nBins =${nBins};/" ${changeFile}
sed -i.bak "s/TString period_Mtype =.*;/TString period_Mtype =\"${period_Mtype}\";/" ${changeFile}
sed -i.bak "s/Int_t hbins =.*;/Int_t hbins =${hbins};/" ${changeFile}
sed -i.bak "s/TString process =.*;/TString process =\"${process}\";/" ${changeFile}
sed -i.bak "s/TString lrMrange =.*;/TString lrMrange =\"${lrMrange}\";/" ${changeFile}
sed -i.bak "s/TString fitMrange =.*;/TString fitMrange =\"${fitMrange}\";/" ${changeFile}
sed -i.bak "s/TString whichFit\[\] =.*;/TString whichFit\[\] =\{\"${7}\", \"${8}\", \"${9}\", \"${10}\"\};/" ${changeFile}
sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${changeFile}
