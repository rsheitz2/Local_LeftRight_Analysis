#!/bin/bash

if [ $# -lt 2 ]; then
    echo "" 
    echo "This script runs input macros over all periods"
    echo "All period considered:"
    echo "\"W07\" \"W08\" \"W09\" \"W10\" \"W11\" \"W12\" \"W13\" \"W14\" \"W15\" \"WAll\""
    echo ""
    echo "To run this script provide as an argument:"
    echo "     Macro to loop over"
    echo "     Enter a number g.t. 0 to proceed (i.e. 1, or \"h\" to see options)"
    echo "     Enter a number g.t. 10 to skip settings output (i.e. 11)"
    echo ""
    exit 1
fi
Macro=$1
Steps=$2



##Setup___ (20) first setup search line 
##########
period=("W07" "W08" "W09" "W10" "W11" "W12" "W13" "W14" "W15" "WAll")
Mtype="HMDY" #"HMDY", "LowM_AMDY"
production="t3" #"t3"=t3, "slot1"=t5
whichCuts="FinalCuts" #"FinalCuts", "PhastCuts"




##Setup ends, (30) last setup search line
if [ ${Steps} == "h" ] || [ ${Steps} -lt 10 ]; then
    echo ""
    echo "______Settings____"
    echo "Mass range type:         ${Mype}"
    echo "Production considered:   ${production}  (\"\"=t3, slot1=t5)"
    echo "Which cuts considered:   $whichCuts     (\"FinalCuts\", \"FinalCuts\")"
    echo " "
    echo " "
fi
if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic checks
if [ $production != "t3" ] && [ $production !="slot1" ]; then
    echo "Invalid production option:   $production"
    exit 1
fi
if [ $whichCuts != "FinalCuts" ] && [ $whichCuts != "PhastCuts" ]; then
    echo "Invalid cuts option:   $whichCuts"
    exit 1
fi

#General variables
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis
cutsPath=${analysisPath}/Macros/Cuts
macroBase=$(basename ${Macro})
macroBaseName=${macroBase%.C}

#Intial save files to be changed
cp ${cutsPath}/${macroBase} ${cutsPath}/${macroBaseName}_tmp.C

for per in ${period[@]}; do
    echo ""
    echo "Period ${per}  for ${macroBase}"
    echo ""

    #macro changes
    ${cutsPath}/Scripts/ChangeScripts/changeMacro.sh ${cutsPath}/${macroBase} ${per} $Mtype $production $whichCuts
    
    #Execute
    root -l -b -q "${cutsPath}/${macroBase}(1)" >> ${cutsPath}/${macroBaseName}_log.txt
    if [ $? != 0 ]; then
	echo "${cutsPath}/${macroBase} did not execute well"
	mv ${cutsPath}/${macroBase} ${cutsPath}/${macroBase}.bak
	mv ${cutsPath}/${macroBaseName}_tmp.C ${cutsPath}/${macroBase}
	exit 1
    else
	rm ${cutsPath}/${macroBaseName}_log.txt
    fi
done

#Clean up
mv ${cutsPath}/${macroBaseName}_tmp.C ${cutsPath}/${macroBase}
rm ${cutsPath}/${macroBase}.bak
