#!/bin/bash

if [ $# -lt 1 ]; then
    echo "" 
    echo "This script gets all data ready to make a table of cuts"
    echo ""
    echo "Steps in this script"
    echo "Step ONE:    tableOfCuts.C for multiple periods"
    echo "Step TWO:    combine_tables.py for all periods from previous step"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" to see the current settings"
    echo "     Or enter the number of steps to proceed through (i.e. 1, 2)"
    echo ""
    exit 1
fi
Steps=$1


##Setup___ (20) first setup search line 
##########
Mtype="LowM_AMDY"
production="t3" #"t3"=t3, "slot1"=t5
whichCuts="FinalCuts"




##Setup ends, (30) last setup search line
echo ""
echo "______Settings____"
echo "Mass range type:         ${Mype}"
echo "Production considered:   ${production}  (\"\"=t3, slot1=t5)"
echo "Which cuts considered:   $whichCuts     (\"FinalCuts\", \"FinalCuts\")"
echo " "
echo " "

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
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant
cutsPath=${analysisPath}/Macros/Cuts


#Step ONE
echo "_______Step ONE_____"
echo "      tableOfCuts.C"

#Intial save files to be changed
cp ${cutsPath}/Scripts/macroPeriodsLoop.sh ${cutsPath}/Scripts/macroPeriodsLoop_tmp.sh

#macro changes
${cutsPath}/Scripts/changePeriodsLoop.sh $Mtype $production $whichCuts
    
#Execute
${cutsPath}/Scripts/macroPeriodsLoop.sh ${cutsPath}/tableOfCuts.C 1
if [ $? != 0 ]; then
    echo "${cutsPath}/Scripts/macroPeriodsLoop.sh did not execute well"
    mv ${cutsPath}/Scripts/macroPeriodsLoop.sh ${cutsPath}/Scripts/macroPeriodsLoop.sh.bak
    mv ${cutsPath}/Scripts/macroPeriodsLoop_tmp.sh ${cutsPath}/Scripts/macroPeriodsLoop.sh
    exit 1
fi

#Clean up
mv ${cutsPath}/Scripts/macroPeriodsLoop_tmp.sh ${cutsPath}/Scripts/macroPeriodsLoop.sh
rm ${cutsPath}/Scripts/macroPeriodsLoop.sh.bak 


#Step TWO
echo "_______Step TWO_____"
echo "      combine_tables.py"

#Intial save files to be changed
cp ${cutsPath}/combine_tables.py ${cutsPath}/combine_tables_tmp.py

#script changes
${cutsPath}/Scripts/changePython.sh ${cutsPath}/combine_tables.py $Mtype $production $whichCuts
    
#Execute
python ${cutsPath}/combine_tables.py 1
if [ $? != 0 ]; then
    echo "${cutsPath}/combine_tables.py did not execute well"
    mv ${cutsPath}/combine_tables.py ${cutsPath}/combine_tables.py.bak
    mv ${cutsPath}/combine_tables_tmp.py ${cutsPath}/combine_tables.py
    exit 1
fi

#Clean up
mv ${cutsPath}/combine_tables_tmp.py ${cutsPath}/combine_tables.py
rm ${cutsPath}/combine_tables.py.bak 
