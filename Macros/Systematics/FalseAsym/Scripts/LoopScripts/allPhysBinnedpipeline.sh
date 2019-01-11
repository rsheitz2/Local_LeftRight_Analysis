#!/bin/bash

if [ $# -ne 2 ]; then
    echo "" 
    echo -n "This script loops over physics binning for input pipeline file"
    echo ""
    echo ""
    echo "To run this script provide as an argument:"
    echo "     pipeline file to loop over (i.e. Scripts/FA2targ_ratioCals_pipeline.sh)"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo "     Or enter a number greater than 10 to only run specific macro from input file (i.e. 11)"
    echo ""
    exit 1
fi

#General variables
loopFile=$1
Steps=$2
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant




##Setup___  first line (25) to seach setup
##########
##Additional settings
production="slot1"
phiPhotonCut="0.0"
##Step ONE settings ###### DY
period="W10"
fitMrangeType="HMDY"
nBins=3
binFile="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/HMDY/BinValues/slot1WAll_HMDY_3bins.txt"
hbins=150
fitMmin=4.30
fitMmax=8.50
binRange="43_85"
##Step TWO settings
process="DY"
LR_Mmin=4.30
LR_Mmax=8.50
whichFit="true"
##Step THREE settings

###Additional settings
#production="slot1"
#phiPhotonCut="0.0"
####Step ONE settings  ########JPsi
#period="W10"
#fitMrangeType="HMDY"
#nBins=3
#binFile="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/HMDY/BinValues/slot1WAll_HMDY_3bins.txt"
#hbins=150
#fitMmin=4.30
#fitMmax=8.50
#binRange="43_85"
###Step TWO settings
#process="DY"
#LR_Mmin=4.30
#LR_Mmax=8.50
#whichFit="true"





additionalCuts="phiS0.0"

##Setup___ last line (70) to search setup
lrMrange="${LR_Mmin}_${LR_Mmax}"
fitMrange="${fitMmin}_${fitMmax}"
period_Mtype="${period}_${fitMrangeType}"
physBinned=("xN" "xPi" "xF" "pT" "M")
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
echo "______Step THREE settings____"
echo " "
echo "Kinematic binning types:       ${physBinned[@]}"
echo "Fit considered:         ${whichFit}"

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Checks
if [ ! -f ${loopFile} ]; then #Full path not given for input file

    tmp_file=${loopFile}
    loopFile=${PWD}/${loopFile}
    if [ ! -f ${loopFile} ]; then
	echo "loopFile does not exist"
	echo ${tmp_file}
	exit 1
    fi
fi

#Basic Setup

HOME=${analysisPath}/Local_LeftRight_Analysis/Macros
aNPath=${HOME}/AN_calculation
sysPath=${HOME}/Systematics/FalseAsym
loopFileBase=$(basename ${loopFile})
loopFileBaseName=${loopFileBase%.sh}

#Intial save files to be changed
cp ${loopFile} ${sysPath}/Scripts/LoopScripts/${loopFileBase}_tmp

for phys in ${physBinned[@]}; do
    echo ""
    echo "Physics Binned ${phys}"
    echo ""
    
    #pipeline changes
    ${sysPath}/Scripts/ChangeScripts/changePipeline.sh ${loopFile} ${period} $fitMrangeType $nBins $hbins $fitMmin $fitMmax ${phys} $process $LR_Mmin $LR_Mmax \
	      ${whichFit} ${binRange} ${binFile} ${production} ${phiPhotonCut} ${additionalCuts}
    #Execute
    ${loopFile} ${Steps} >> ${sysPath}/Scripts/LoopScripts/${loopFileBaseName}_log.txt
    if [ $? != 0 ]; then
	echo "${loopFile}.sh did not execute well"
	mv ${loopFile} ${loopFile}.bak
	mv ${sysPath}/Scripts/LoopScripts/${loopFileBase}_tmp ${loopFile}
	exit 1
    else
	rm ${sysPath}/Scripts/LoopScripts/${loopFileBaseName}_log.txt
    fi
done


echo ""
echo "Integrated:  Physics Binned xN"
echo ""

#pipeline changes
${sysPath}/Scripts/ChangeScripts/changePipeline.sh ${loopFile} ${period} $fitMrangeType 1 $hbins $fitMmin $fitMmax xN $process $LR_Mmin $LR_Mmax \
	  ${whichFit} ${binRange} ${binFile} ${production} ${phiPhotonCut} ${additionalCuts}
#Execute
${loopFile} ${Steps} >> ${sysPath}/Scripts/LoopScripts/${loopFileBaseName}_log.txt
if [ $? != 0 ]; then
    echo "${loopFile}.sh did not execute well"
    mv ${loopFile} ${loopFile}.bak
    mv ${sysPath}/Scripts/LoopScripts/${loopFileBase}_tmp ${loopFile}
    exit 1
else
    rm ${sysPath}/Scripts/LoopScripts/${loopFileBaseName}_log.txt
fi


#Clean up
mv ${sysPath}/Scripts/LoopScripts/${loopFileBase}_tmp ${loopFile}
rm ${loopFile}.bak
