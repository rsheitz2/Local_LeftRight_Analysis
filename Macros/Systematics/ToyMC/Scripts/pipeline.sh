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



##Setup___  first line (20) to seach setup
##########
phiScut=(0.0 0.044 0.088 0.17 0.26 0.36 0.53 0.71 0.88 1.07 1.24 1.44)
A_siv=(0.0 0.005 0.01 0.02 0.05 0.1)
N_gen=1000
alphaScale=true




##Setup___ last line (30) to search setup

echo ""
echo "______Setup____"
echo "phiScuts:   ${phiScut[@]}"
echo "Sivers amplitude:   ${A_siv[@]}"
echo "Number of generated events per target:  ${N_gen}"
echo "To generate different numbers of events per target?  ${alphaScale}"
echo ""

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Setup
HOME=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/ToyMC

#Intial save files to be changed
cp ${HOME}/generator.C ${HOME}/generator_tmp.C
cp ${HOME}/GeoMean4Targ.C ${HOME}/GeoMean4Targ_tmp.C

for phiS in ${phiScut[@]}; do
    echo ""
    echo "phiS   $phiS"
    echo ""

    for A in ${A_siv[@]}; do
	echo ""
	echo "Sivers amplitude  $A"
	echo ""

	additionalCuts=phiS${phiS}

	#Generator changes
	${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/generator.C $phiS $A $additionalCuts $N_gen $alphaScale
	
	#Execute
	root -l -q -b "${HOME}/generator.C" >> ${HOME}/generator_log.txt
	if [ $? != 0 ]; then
	    echo "generator.C did not execute well"
	    mv ${HOME}/generator.C ${HOME}/generator.C.bak
	    mv ${HOME}/generator_tmp.C ${HOME}/generator.C
	    exit 1
	fi

	
	#GeoMean changes
	${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/GeoMean4Targ.C $phiS $A $additionalCuts $N_gen $alphaScale
	
	#Execute
	root -l -q -b "${HOME}/GeoMean4Targ.C(1)" >> ${HOME}/GeoMean4Targ_log.txt
	if [ $? != 0 ]; then
	    echo "GeoMean4Targ.C did not execute well"
	    mv ${HOME}/GeoMean4Targ.C ${HOME}/GeoMean4Targ.C.bak
	    mv ${HOME}/GeoMean4Targ_tmp.C ${HOME}/GeoMean4Targ.C
	    exit 1
	fi
	
    done #Sivers binning
done #phiS binning

#Clean up
mv ${HOME}/generator_tmp.C ${HOME}/generator.C
mv ${HOME}/GeoMean4Targ_tmp.C ${HOME}/GeoMean4Targ.C
rm ${HOME}/*.bak
