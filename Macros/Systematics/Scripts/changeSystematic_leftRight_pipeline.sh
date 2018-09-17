#!/bin/bash

if [ $# -ne 6 ]; then
     echo "" 
     echo "Script updates systematic_leftRight_pipeline.C with input paramters"
     echo ""
     
else
    period=$1
    fitMrangeType=$2
    nBins=$3
    hbins=$4
    fitMmin=$5
    fitMmax=$6

    changeFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/Scripts/systematic_leftRight_pipeline.sh
    sed -i.bak "20,40s/period=\".*\"/period=\"${period}\"/" ${changeFile}
    sed -i.bak "20,40s/fitMrangeType=\".*\"/fitMrangeType=\"${fitMrangeType}\"/" ${changeFile}
    sed -i.bak "20,40s/nBins=.*/nBins=${nBins}/" ${changeFile}
    sed -i.bak "20,40s/hbins=.*/hbins=${hbins}/" ${changeFile}
    sed -i.bak "20,40s/fitMmin=.*/fitMmin=${fitMmin}/" ${changeFile}
    sed -i.bak "20,40s/fitMmax=.*/fitMmax=${fitMmax}/" ${changeFile}
fi
