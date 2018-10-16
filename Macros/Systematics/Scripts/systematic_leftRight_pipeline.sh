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
fitMrangeType="LowM_AMDY"
nBins=5
binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/WAll_JPsi25_43_${nBins}bins.txt
hbins=150
fitMmin=2.00
fitMmax=7.50
binRange="25_43"
whichFit="ten"








##Setup ends, last setup search line (40)
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
if [ $whichFit != "true" ]; then
    fitMmin=1.00
    fitMmax=8.50
    echo "Mass ranged updated to  $fitMmin - $fitMmax "
fi

stepOne_Out=${pathOne}"/Data/systematic_leftRight_"${period}"_"${fitMrangeType}${fitMmin}_${fitMmax}"_"${nBins}"bins"${binRange}"_"${hbins}"hbin.root"
if [ ! -f ${stepOne_Out} ]; then
    echo "Making systematic_leftRight one data:"
    echo ""
    echo ""
    echo ""
    ${pathOne}/systematic_leftRight -i${fitMmin} -a${fitMmax} -Q${stepOne_Out} -b${binFile} -Z${hbins} -f${InputData} 
else
    echo "systematic_leftRight one file already exist"
fi


stepTwo_Out=${pathOne}"/Data/systematic_leftRight_"${period}"_"${fitMrangeType}${fitMmin}_${fitMmax}"_"${nBins}"bins"${binRange}"_"${hbins}"hbin_optionER.root"
if [ ! -f ${stepTwo_Out} ]; then
    echo "Making systematic_leftRight two data:"
    echo ""
    echo ""
    echo ""
    ${pathOne}/systematic_leftRight -i${fitMmin} -a${fitMmax} -Q${stepTwo_Out} -b${binFile} -E -R -Z${hbins} -f${InputData} 
else
    echo "systematic_leftRight file two already exist"
fi
