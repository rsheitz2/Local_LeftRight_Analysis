#!/bin/bash

if [ $# -ne 11 ]; then
     echo "" 
     echo "Script updates FApipeline.C with input paramters"
     echo ""
     exit 1
fi     

period=$1
fitMrangeType=$2
nBins=$3
hbins=$4
fitMmin=$5
fitMmax=$6
physBinned=$7
process=$8
LR_Mmin=$9
LR_Mmax=${10}
whichFit=${11}

changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Scripts/FApipeline.sh
sed -i.bak "30,60s/period=\".*\"/period=\"${period}\"/" ${changeFile}
sed -i.bak "30,60s/fitMrangeType=\".*\"/fitMrangeType=\"${fitMrangeType}\"/" ${changeFile}
sed -i.bak "30,60s/nBins=.*/nBins=${nBins}/" ${changeFile}
sed -i.bak "30,60s/hbins=.*/hbins=${hbins}/" ${changeFile}
sed -i.bak "30,60s/fitMmin=.*/fitMmin=${fitMmin}/" ${changeFile}
sed -i.bak "30,60s/fitMmax=.*/fitMmax=${fitMmax}/" ${changeFile}
sed -i.bak "30,60s/physBinned=.*/physBinned=\"${physBinned}\"/" ${changeFile}
sed -i.bak "30,60s/process=.*/process=\"${process}\"/" ${changeFile}
sed -i.bak "30,60s/LR_Mmin=.*/LR_Mmin=${LR_Mmin}/" ${changeFile}
sed -i.bak "30,60s/LR_Mmax=.*/LR_Mmax=${LR_Mmax}/" ${changeFile}
sed -i.bak "30,60s/whichFit=.*/whichFit=\"${whichFit}\"/" ${changeFile}
