#!/bin/bash

if [ $# -ne 17 ]; then
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
production=${15}
phiPhotonCut=${16}
additionalCuts=${17}


changeFile=${PWD}/${inFile}
if [ ! -f ${changeFile} ]; then
    changeFile=${inFile}
fi
if [ ! -f ${changeFile} ]; then #try again
    echo "$changeFile does not exist"
    exit 1
fi

sed -i.bak "25,70s/period=\".*\"/period=\"${period}\"/" ${changeFile}
sed -i.bak "25,70s/fitMrangeType=\".*\"/fitMrangeType=\"${fitMrangeType}\"/" ${changeFile}
sed -i.bak "25,70s/nBins=.*/nBins=${nBins}/" ${changeFile}
sed -i.bak "25,70s/hbins=.*/hbins=${hbins}/" ${changeFile}
sed -i.bak "25,70s/fitMmin=.*/fitMmin=${fitMmin}/" ${changeFile}
sed -i.bak "25,70s/fitMmax=.*/fitMmax=${fitMmax}/" ${changeFile}
sed -i.bak "25,70s/physBinned=.*/physBinned=\"${physBinned}\"/" ${changeFile}
sed -i.bak "25,70s/process=.*/process=\"${process}\"/" ${changeFile}
sed -i.bak "25,70s/LR_Mmin=.*/LR_Mmin=${LR_Mmin}/" ${changeFile}
sed -i.bak "25,70s/LR_Mmax=.*/LR_Mmax=${LR_Mmax}/" ${changeFile}
sed -i.bak "25,70s/whichFit=.*/whichFit=\"${whichFit}\"/" ${changeFile}
sed -i.bak "25,70s/binRange=.*/binRange=\"${binRange}\"/" ${changeFile}
sed -i.bak "25,70s%binFile=.*%binFile=\"${binFile}\"%" ${changeFile}
sed -i.bak "25,70s%production=.*%production=\"${production}\"%" ${changeFile}
sed -i.bak "25,70s%phiPhotonCut=.*%phiPhotonCut=\"${phiPhotonCut}\"%" ${changeFile}
sed -i.bak "25,70s%additionalCuts=.*%additionalCuts=\"${additionalCuts}\"%" ${changeFile}
