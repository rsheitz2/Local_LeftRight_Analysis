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
nBins=3
Mtype="HMDY"
production="slot1"
minimizer="MINOS" #"MIGRAD", "MINOS"




##Setup___ last line (30) to search setup
period=("W07" "W08" "W09" "W10" "W11" "W12" "W13" "W14" "W15" "WAll")
physBinned=("xN" "xPi" "xF" "pT" "M")

echo ""
echo "______Step ONE settings____"
echo "Number of kinematic bins:   ${nBins}"
echo "Mass type considered:    $Mtype"
echo "production used:   $production"
echo " "
echo "Kinematic binning types:       ${physBinned[@]}"
echo "Periods considered:         ${period[@]}"

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Setup
HOME=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting

#Intial save files to be changed
cp ${HOME}/weightingCounts.C ${HOME}/weightingCounts_tmp.C
cp ${HOME}/fitEventWeighting.C ${HOME}/fitEventWeighting_tmp.C
cp ${HOME}/wAvg.C ${HOME}/wAvg_tmp.C

for phy in ${physBinned[@]}; do
    echo ""
    echo "Physics Binned $phy"
    echo ""
    
    for per in ${period[@]}; do
	echo ""
	echo "Period $per"
	echo ""

	#pipeline changes
	${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/weightingCounts.C $nBins $per $Mtype $phy ${production} "blank"
	#Execute
	root -l -b -q ${HOME}/weightingCounts.C >> ${HOME}/weightingCounts_log.txt
	if [ $? != 0 ]; then
	    echo "${HOME}/weightingCounts.C did not execute well"
	    mv ${HOME}/weightingCounts.C ${HOME}/weightingCounts.C.bak
	    mv ${HOME}/weightingCounts_tmp.C ${HOME}/weightingCounts.C
	    exit 1
	else
	    rm ${HOME}/weightingCounts_log.txt
	fi

	#pipeline changes
	${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/fitEventWeighting.C $nBins $per $Mtype $phy ${production} ${minimizer}
	#Execute
	root -l -b -q ${HOME}/fitEventWeighting.C >> ${HOME}/fitEventWeighting_log.txt
	if [ $? != 0 ]; then
	    echo "${HOME}/fitEventWeighting.C did not execute well"
	    mv ${HOME}/fitEventWeighting.C ${HOME}/fitEventWeighting.C.bak
	    mv ${HOME}/fitEventWeighting_tmp.C ${HOME}/fitEventWeighting.C
	    exit 1
	else
	    rm ${HOME}/fitEventWeighting_log.txt
	fi
    done #period binning

    #pipeline changes
    ${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/wAvg.C $nBins $nHbins "blank" $Mtype $phy ${production} ${minimizer}
    #Execute
    root -l -b -q ${HOME}/wAvg.C >> ${HOME}/wAvg_log.txt
    if [ $? != 0 ]; then
    	echo "wAvg.C did not execute well"
    	mv ${HOME}/wAvg.C ${HOME}/wAvg.C.bak
    	mv ${HOME}/wAvg_tmp.C ${HOME}/wAvg.C
    	exit 1
    else
    	rm ${HOME}/wAvg_log.txt
    fi
done #physics binning

#Clean up
mv ${HOME}/weightingCounts_tmp.C ${HOME}/weightingCounts.C
mv ${HOME}/fitEventWeighting_tmp.C ${HOME}/fitEventWeighting.C
mv ${HOME}/wAvg_tmp.C ${HOME}/wAvg.C
rm ${HOME}/*.bak
