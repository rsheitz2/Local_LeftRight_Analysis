#!/bin/bash

if [ $# -ne 2 ]; then
    echo "" 
    echo -n "This script loops over physics binning for input pipeline file"
    echo ""
    echo ""
    echo "To run this script provide as an argument:"
    echo "     pipeline file to loop over (i.e. Scripts/FA2targ_ratioCals_pipeline.sh)"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number of steps to take (i.e. 1-3)"
    echo ""
    exit 1
fi

#General variables
loopFile=$1
Steps=$2
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant




##Setup___  first line (30) to seach setup
##########
####Additional settings
#production="slot1"
#phiPhotonCut=0.53 #HMDY=0.1866, #LowM_AMDY=0.195
###Step ONE settings ###### DY
#period="WAll"
#fitMrangeType="HMDY"
#nBins=3
#binFile=${analysisPath}/Presents/DATA/RealData/HMDY/BinValues/WAll_HMDY_${nBins}bins.txt
#hbins=150
#fitMmin=4.30  #true fit mass range
#fitMmax=8.50  #true fit mass range
#binRange="43_85"
###Step TWO settings
#process="DY"
#LR_Mmin=4.30
#LR_Mmax=8.50
#physBinned=("xN" "xPi" "xF" "pT")
#whichFit=("true" "true" "true" "true")
###Step THREE settings

###Additional settings
production="slot1"
phiPhotonCut=0.53 #HMDY=0.187, #LowM_AMDY=0.195
##Step ONE settings  ########JPsi
period="WAll"
fitMrangeType="LowM_AMDY"
nBins=5
binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/WAll_JPsi25_43_${nBins}bins.txt
hbins=150
fitMmin=2.00  #true fit mass range
fitMmax=7.50  #true fit mass range
binRange="25_43"
##Step TWO settings
process="JPsi"
LR_Mmin=2.00
LR_Mmax=5.00
physBinned=("xN" "xPi" "xF" "pT")
whichFit=("eight" "eight" "eight" "eight")

additionalCuts=phiS$phiPhotonCut #add and new cuts here.  This should include all cuts used

##Setup___ last line (60) to search setup
lrMrange="${LR_Mmin}_${LR_Mmax}"
fitMrange="${fitMmin}_${fitMmax}"
period_Mtype="${period}_${fitMrangeType}"
echo ""
echo "______Step ONE settings____"
echo "Period:   ${period}"
echo "Fit mass range type:  ${fitMrangeType}"
echo "Number of kinematic bins:   ${nBins}"
echo "Number of histogram bins in M distribution:  ${hbins}"
echo "Min left/right mass range:        ${fitMmin}"
echo "Max left/right mass range:        ${fitMmax}"
echo " "
echo "______Step TWO settings____"
echo "Integrated physics process:   ${process}"
if [ ${whichFit[0]} == "true" ]; then
    echo "Min left/right mass range:        ${fitMmin}"
    echo "Max left/right mass range:        ${fitMmax}"
    LR_Mmin=${fitMmin}
    LR_Mmax=${fitMmax}
    lrMrange="${LR_Mmin}_${LR_Mmax}"
else
    echo "Min integration range:            ${LR_Mmin}"
    echo "Max integration range:            ${LR_Mmax}"
    echo "Minimum fit mass range:                          ${fitMmin}"
    echo "Maximum fit mass range:                          ${fitMmax}"	
fi
echo " "
echo "______Step THREE settings____"
echo " "
echo "Kinematic binning types:       ${physBinned[@]}"
echo "Fits considered:         ${whichFit[@]}"

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Setup
HOME=${analysisPath}/Local_LeftRight_Analysis/Macros
aNPath=${HOME}/AN_calculation

#Intial save files to be changed
cp ${PWD}/${loopFile} ${PWD}/${loopFile}_tmp

for i in `seq 0 3`; do
    echo ""
    echo "Physics Binned ${physBinned[$i]}"
    echo ""
    
    #pipeline changes
    ${aNPath}/Scripts/changeGenericPipeline.sh ${loopFile} ${period} $fitMrangeType $nBins $hbins $fitMmin $fitMmax ${physBinned[$i]} $process $LR_Mmin $LR_Mmax \
	      ${whichFit[$i]} ${binRange} ${binFile} ${production} ${phiPhotonCut} $additionalCuts
    #Execute
    ${PWD}/${loopFile} ${Steps} >> ${PWD}/${loopFile}_log.txt
    if [ $? != 0 ]; then
	echo "${loopFile}.sh did not execute well"
	mv ${PWD}/${loopFile} ${PWD}/${loopFile}.bak
	mv ${PWD}/${loopFile}_tmp ${PWD}/${loopFile}
	exit 1
    else
	rm ${PWD}/${loopFile}_log.txt
    fi
done

#Clean up
mv ${PWD}/${loopFile}_tmp ${PWD}/${loopFile}
rm ${PWD}/${loopFile}.bak
