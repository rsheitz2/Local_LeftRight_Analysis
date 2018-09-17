#!/bin/bash

if [ $# -lt 1 ]; then
    echo "" 
    echo "This script gets all the data from systematic_leftRight"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" to see the current settings"
    echo "     Or enter a number great than 0 to run the script (i.e. 1)"
    echo ""
    exit
fi

#General variables
Steps=$1
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant



##Setup___ first setup search line
##########
##Step ONE settings
period="WAll"
fitMrangeType="HMDY"
nBins=3
hbins=150
fitMmin=4.30 #true fit mass range
fitMmax=8.50 #true fit mass range











##Setup ends, last setup search line

binFile=${analysisPath}/Presents/DATA/RealData/${fitMrangeType}/BinValues/WAll_${fitMrangeType}_${nBins}bins.txt
InputData=${analysisPath}/Presents/DATA/RealData/
if [ ${fitMrangeType} == "HMDY" ]; then #Speed optimization for lower data set HMDY
    InputData+=${fitMrangeType}/${period}_${fitMrangeType}.root
else
    InputData+=LowM_AMDY/${period}_LowM_AMDY.root
fi
echo ""
echo "______Step ONE settings____"
echo "Period:   ${period}"
echo "Fit mass range type:  ${fitMrangeType}"
echo "Number of kinematic bins:   ${nBins}"
echo "Number of histogram bins in M distribution:  ${hbins}"
echo "Binning file:"
echo "    ${binFile}"
echo "Input data:"
echo "    ${InputData}"
echo "Min left/right mass range:        ${fitMmin}"
echo "Max left/right mass range:        ${fitMmax}"
echo " "

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Run script
echo "_______Step ONE_____"
echo "      systematic_leftRight"

#Start checks
pathOne=${analysisPath}/Local_LeftRight_Analysis
if [ ! -x ${pathOne}/systematic_leftRight ]; then
    echo "systematic_leftRight does not exist"
    exit 1
fi
if [ ! -f ${binFile} ] || [ ! -f ${InputData} ]; then
    echo "Step ONE file does not exist:"
    echo "${binFile}"
    echo "${InputData}"
    exit 1
fi

#Execute Step ONE if files don't already exist
stepOne_Out=${pathOne}"/Data/systematic_leftRight_"${period}"_"${fitMrangeType}${fitMmin}_${fitMmax}"_"${nBins}"bins_"${hbins}"hbin.root"
if [ ! -f ${stepOne_OutData} ]; then
    echo "Making systematic_leftRight data:"
    echo ""
    echo ""
    echo ""
    ${pathOne}/systematic_leftRight -i${fitMmin} -a${fitMmax} -Q${stepOne_Out} -b${binFile} -Z${hbins} -f${InputData} 
else
    echo "systematic_leftRight file already exist"
fi
