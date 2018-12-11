#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script loops over all physics periods for all physics binnings"
    echo ""
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number of steps to take (i.e. 1)"
    echo ""
    exit 1
fi

#General variables
Steps=$1



##Setup___  first line (20) to seach setup
##########
nBins=1
nHbins=16
Mtype="HMDY"
production="slot1"




##Setup___ last line (30) to search setup
period=("W07" "W08" "W09" "W10" "W11" "W12" "W13" "W14" "W15" "WAll")
physBinned=("xN" "xPi" "xF" "pT" "M")

echo ""
echo "______Step ONE settings____"
echo "Number of kinematic bins:   ${nBins}"
echo "Fit in histogram bins of:  ${nHbins}"
echo "Mass type considered:    $Mtype"
echo "production used:   $production"
echo " "
echo "Kinematic binning types:       ${physBinned[@]}"
echo "Periods considered:         ${period[@]}"

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Setup
HOME=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio

#Intial save files to be changed
loopFile=${HOME}/doubleRatio.C
cp $loopFile ${loopFile}_tmp

for per in ${period[@]}; do
    echo ""
    echo "Period $per"
    echo ""

    for phy in ${physBinned[@]}; do
	echo ""
	echo "Physics Binned $phy"
	echo ""

	#pipeline changes
	${HOME}/Scripts/ChangeScripts/changeMacro.sh ${loopFile} $nBins $nHbins $per $Mtype $phy ${production}
	#Execute
	root -l -b -q "${loopFile}(1)" >> ${loopFile}_log.txt
	if [ $? != 0 ]; then
	    echo "${loopFile} did not execute well"
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
