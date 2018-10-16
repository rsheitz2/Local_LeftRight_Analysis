#!/bin/bash

if [ $# -ne 14 ]; then
     echo "" 
     echo "Script updates input file with input paramters"
     echo ""
     exit 1
fi     

inFile=$1
period=$2
fitMrangeType=$3
nBins=$4
hbins=$5
fitMmin=$6
fitMmax=$7
physBinned=$8
process=$9
LR_Mmin=${10}
LR_Mmax=${11}
whichFit=${12}
binRange=${13}
binFile=${14}


changeFile=${PWD}/${inFile}
if [ ! -f ${changeFile} ]; then
    changeFile=${inFile}
fi
if [ ! -f ${changeFile} ]; then #try again
    echo "$changeFile does not exist"
    exit 1
fi

sed -i.bak "25,60s/period=\".*\"/period=\"${period}\"/" ${changeFile}
sed -i.bak "25,60s/fitMrangeType=\".*\"/fitMrangeType=\"${fitMrangeType}\"/" ${changeFile}
sed -i.bak "25,60s/nBins=.*/nBins=${nBins}/" ${changeFile}
sed -i.bak "25,60s/hbins=.*/hbins=${hbins}/" ${changeFile}
sed -i.bak "25,60s/fitMmin=.*/fitMmin=${fitMmin}/" ${changeFile}
sed -i.bak "25,60s/fitMmax=.*/fitMmax=${fitMmax}/" ${changeFile}
sed -i.bak "25,60s/physBinned=.*/physBinned=\"${physBinned}\"/" ${changeFile}
sed -i.bak "25,60s/process=.*/process=\"${process}\"/" ${changeFile}
sed -i.bak "25,60s/LR_Mmin=.*/LR_Mmin=${LR_Mmin}/" ${changeFile}
sed -i.bak "25,60s/LR_Mmax=.*/LR_Mmax=${LR_Mmax}/" ${changeFile}
sed -i.bak "25,60s/whichFit=.*/whichFit=\"${whichFit}\"/" ${changeFile}
sed -i.bak "25,60s/binRange=.*/binRange=\"${binRange}\"/" ${changeFile}
sed -i.bak "25,60s%binFile=.*%binFile=\"${binFile}\"%" ${changeFile}
