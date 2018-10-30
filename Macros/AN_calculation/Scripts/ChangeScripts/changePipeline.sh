#!/bin/bash

if [ $# -ne 13 ]; then
     echo "" 
     echo "Script updates pipeline.C with input paramters"
     echo ""
     
else
    period=$1
    fitMrangeType=$2
    nBins=$3
    hbins=$4
    physBinned=$5
    process=$6
    LR_Mmin=$7
    LR_Mmax=$8
    fitMmin=$9
    fitMmax=${10}
    whichFit=${11}
    binRange=${12}
    binFile=${13}

    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Scripts/pipeline.sh
    sed -i.bak "25,50s/period=\".*\"/period=\"${period}\"/" ${changeFile}
    sed -i.bak "25,50s/fitMrangeType=\".*\"/fitMrangeType=\"${fitMrangeType}\"/" ${changeFile}
    sed -i.bak "25,50s/nBins=.*/nBins=${nBins}/" ${changeFile}
    sed -i.bak "25,50s/hbins=.*/hbins=${hbins}/" ${changeFile}
    sed -i.bak "25,50s/physBinned=.*/physBinned=\"${physBinned}\"/" ${changeFile}
    sed -i.bak "25,50s/process=.*/process=\"${process}\"/" ${changeFile}
    sed -i.bak "25,50s/LR_Mmin=.*/LR_Mmin=${LR_Mmin}/" ${changeFile}
    sed -i.bak "25,50s/LR_Mmax=.*/LR_Mmax=${LR_Mmax}/" ${changeFile}
    sed -i.bak "25,50s/fitMmin=.*/fitMmin=${fitMmin}/" ${changeFile}
    sed -i.bak "25,50s/fitMmax=.*/fitMmax=${fitMmax}/" ${changeFile}
    sed -i.bak "25,50s/whichFit=.*/whichFit=\"${whichFit}\"/" ${changeFile}
    sed -i.bak "25,50s/binRange=.*/binRange=\"${binRange}\"/" ${changeFile}
    sed -i.bak "25,50s%binFile=.*%binFile=\"${binFile}\"%" ${changeFile}
fi
