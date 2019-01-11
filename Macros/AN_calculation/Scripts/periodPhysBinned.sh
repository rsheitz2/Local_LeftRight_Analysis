#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script loops over all physics periods for all physics binnings"
    echo ""
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number of steps to take (i.e. 1-3)"
    echo ""
    exit 1
fi

#General variables
Steps=$1
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant







##Setup___  first line (25) to seach setup
##########
###Additional settings
production="slot1"
phiPhotonCut=0.0 #0.0, 0.044, 0.088, 0.17, 0.36, 0.53, 0.71, 0.88, 1.07 #HMDY=0.1866, #LowM_AMDY=0.195
##Step ONE settings ###### DY
fitMrangeType="HMDY"
nBins=1
binFile=${analysisPath}/Presents/DATA/RealData/HMDY/BinValues/slot1WAll_HMDY_${nBins}bins.txt
#binFile=${analysisPath}/Presents/DATA/RealData/HMDY/BinValues/WAll_HMDY_${nBins}bins.txt
hbins=150
fitMmin=4.30  #true fit mass range
fitMmax=8.50  #true fit mass range
binRange="43_85"
##Step TWO settings
process="DY"
LR_Mmin=4.30
LR_Mmax=8.50
whichFit="true"


###Additional settings
#production="slot1"
#phiPhotonCut=0.53 #HMDY=0.187, #LowM_AMDY=0.195
###Step ONE settings  ########JPsi
#fitMrangeType="LowM_AMDY"
#nBins=5
#hbins=150
#fitMmin=2.00  #true fit mass range
#fitMmax=8.50  #true fit mass range
##fitMmax: xN=6.00, xPi=6.30, xF=pT=8.50
#binRange="25_43"
#binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/slot1WAll_JPsi${binRange}_${nBins}bins.txt
###Step TWO settings
#process="JPsi"
#LR_Mmin=2.00
#LR_Mmax=5.00
#whichFit="thirteen"


additionalCuts=phiS$phiPhotonCut #add and new cuts here.  This should include all cuts used




##Setup___ last line (70) to search setup
period=("W07" "W08" "W09" "W10" "W11" "W12" "W13" "W14" "W15" "WAll")
physBinned=("xN" "xPi" "xF" "pT" "M")
loopFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Scripts/Src/pipeline.sh

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
if [ ${whichFit} == "true" ]; then
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
echo "______Additional settings____"
echo "production used:   $production"
echo "phi photon target frame cut:    $phiPhotonCut"
echo "All cuts used:     $additionalCuts"
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
cp ${loopFile} ${loopFile}_tmp

for per in ${period[@]}; do
    echo ""
    echo "Period $per"
    echo ""

    for phy in ${physBinned[@]}; do
	echo ""
	echo "Physics Binned $phy"
	echo ""
    
	#pipeline changes
	${aNPath}/Scripts/ChangeScripts/changeGenericPipeline.sh ${loopFile} ${per} $fitMrangeType $nBins $hbins $fitMmin $fitMmax $phy $process $LR_Mmin $LR_Mmax \
		 ${whichFit} ${binRange} ${binFile} ${production} ${phiPhotonCut} ${additionalCuts}
	#Execute
	${loopFile} ${Steps} >> ${loopFile}_log.txt
	if [ $? != 0 ]; then
	    echo "${loopFile}.sh did not execute well"
	    mv ${loopFile} ${loopFile}.bak
	    mv ${loopFile}_tmp ${loopFile}
	    exit 1
	else
	    rm ${loopFile}_log.txt
	fi
    done #physics binning
done #period binning

#Clean up
mv ${loopFile}_tmp ${loopFile}
rm ${loopFile}.bak
