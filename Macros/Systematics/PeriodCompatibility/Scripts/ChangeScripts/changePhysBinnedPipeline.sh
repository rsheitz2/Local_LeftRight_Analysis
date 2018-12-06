#!/bin/bash

if [ $# -ne 10 ]; then
     echo "" 
     echo "Script updates pipeline.C with input paramters"
     echo ""
     
else
    fitMrangeType=$1
    nBins=$2
    hbins=$3
    physBinned=$4
    process=$5
    LR_Mmin=$6
    LR_Mmax=$7
    fitMmin=$8
    fitMmax=$9
    whichFit=${10}
    
    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Scripts/physBinnedPipeline.sh
    sed -i.bak "30,60s/fitMrangeType=\".*\"/fitMrangeType=\"${fitMrangeType}\"/" ${changeFile}
    sed -i.bak "30,60s/nBins=.*/nBins=${nBins}/" ${changeFile}
    sed -i.bak "30,60s/hbins=.*/hbins=${hbins}/" ${changeFile}
    sed -i.bak "30,60s/physBinned=.*/physBinned=\"${physBinned}\"/" ${changeFile}
    sed -i.bak "30,60s/process=.*/process=\"${process}\"/" ${changeFile}
    sed -i.bak "30,60s/LR_Mmin=.*/LR_Mmin=${LR_Mmin}/" ${changeFile}
    sed -i.bak "30,60s/LR_Mmax=.*/LR_Mmax=${LR_Mmax}/" ${changeFile}
    sed -i.bak "30,60s/fitMmin=.*/fitMmin=${fitMmin}/" ${changeFile}
    sed -i.bak "30,60s/fitMmax=.*/fitMmax=${fitMmax}/" ${changeFile}
    sed -i.bak "30,60s/whichFit=.*/whichFit=\"${whichFit}\"/" ${changeFile}

fi
