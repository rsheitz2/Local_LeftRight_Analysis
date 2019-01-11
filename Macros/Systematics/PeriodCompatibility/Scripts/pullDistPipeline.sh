#!/bin/bash

if [ $# -ne 1 ]; then
    echo ""
    echo "Pipeline for by period AN"
    echo "This script works for:  physBinnedPeriod.C, pullDist.C, allPhysBinned.C"
    echo ""
    echo "To do:  falseApullDist.C, acceptPullDist.C"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo ""
    exit 1
fi    

Steps=$1



##Setup___
##########
###HMDY
fitMrangeType="HMDY"
nBins=3
hbins=150
process="DY"
lrMrange="4.30_8.50" #0.01 precision  #does nothing for true fit
fitMrange="4.30_8.50" #0.01 precision
binRange="43_85"
whichFit="true"
production="slot1" #"t3", "slot1"
additionalCuts="phiS0.0" 

##Setup___
##########

#General variables/Basic setup
HOME=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros
thisSysPath=${HOME}/Systematics/PeriodCompatibility
physBinned=("xN" "xPi" "xF" "pT" "M")

echo ""
echo "______Step ONE settings____"
echo "Fit mass range type:  ${fitMrangeType}"
echo "Number of kinematic bins:   ${nBins}"
echo "Number of histogram bins in M distribution:  ${hbins}"
echo "Min/Max left/right mass range:        ${fitMrange}"
echo "Integrated physics process:   ${process}"
echo "Min/Max integration range:        ${lrMrange}"
echo "Binned in which mass range:       $binRange"
echo "WhichFit per physics binning       ${whichFit}"
echo "Which production considered:       ${production}"
echo "Which other cuts considered:       ${additionalCuts}"
echo " "
echo "Physics binnings considered    ${physBinned[*]}"
echo " "

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic checks
if [ ${whichFit} == "true" ]; then
    if [ ${lrMrange} != ${fitMrange} ]; then
	echo "set lrMrange equal to fitMrange ust in case..."
	exit 1
    fi
fi

echo " "
echo "physBinnedPeriod.C"
echo " "
cp ${thisSysPath}/physBinnedPeriod.C ${thisSysPath}/tmp_physBinnedPeriod.C
if [ -f ${thisSysPath}/Scripts/Logs/physBinnedPeriod_log.txt ]; then
    rm ${thisSysPath}/Scripts/Logs/physBinnedPeriod_log.txt
fi

for phys in ${physBinned[@]}; do
    #Changes
    ${thisSysPath}/Scripts/ChangeScripts/changeMacro.sh ${thisSysPath}/physBinnedPeriod.C $nBins $fitMrangeType $hbins $phys $process $lrMrange $fitMrange $binRange \
		  $whichFit $production $additionalCuts
    #Execute 
    root -l -b -q "${thisSysPath}/physBinnedPeriod.C(1)" >> ${thisSysPath}/Scripts/Logs/physBinnedPeriod_log.txt
    if [ $? != 0 ]; then
	echo "physBinnedPeriod.C did not execute well"
	mv ${thisSysPath}/tmp_physBinnedPeriod.C ${thisSysPath}/physBinnedPeriod.C
	exit 1
    fi
done
mv ${thisSysPath}/tmp_physBinnedPeriod.C ${thisSysPath}/physBinnedPeriod.C

echo " "
echo "pullDist.C"
echo " "
cp ${thisSysPath}/pullDist.C ${thisSysPath}/tmp_pullDist.C
if [ -f ${thisSysPath}/Scripts/Logs/pullDist_log.txt ]; then
    rm ${thisSysPath}/Scripts/Logs/pullDist_log.txt
fi

#Changes
${thisSysPath}/Scripts/ChangeScripts/changeMacro.sh ${thisSysPath}/pullDist.C $nBins $fitMrangeType $hbins $phys $process $lrMrange $fitMrange $binRange \
		  $whichFit $production $additionalCuts
    
#Execute 
root -l -b -q "${thisSysPath}/pullDist.C(1)" >> ${thisSysPath}/Scripts/Logs/pullDist_log.txt
if [ $? != 0 ]; then
    echo "pullDist.C did not execute well"
    mv ${thisSysPath}/tmp_pullDist.C ${thisSysPath}/pullDist.C
    exit 1
fi
    
#Final Cleanup
mv ${thisSysPath}/tmp_pullDist.C ${thisSysPath}/pullDist.C


echo " "
echo "allPhysBinned.C"
echo " "
cp ${thisSysPath}/allPhysBinned.C ${thisSysPath}/tmp_allPhysBinned.C
if [ -f ${thisSysPath}/Scripts/Logs/allPhysBinned_log.txt ]; then
    rm ${thisSysPath}/Scripts/Logs/allPhysBinned_log.txt
fi

#Changes
${thisSysPath}/Scripts/ChangeScripts/changeMacro.sh ${thisSysPath}/allPhysBinned.C $nBins $fitMrangeType $hbins $phys $process $lrMrange $fitMrange $binRange \
		  $whichFit $production $additionalCuts
    
#Execute 
root -l -b -q "${thisSysPath}/allPhysBinned.C(1)" >> ${thisSysPath}/Scripts/Logs/allPhysBinned_log.txt
if [ $? != 0 ]; then
    echo "allPhysBinned.C did not execute well"
    mv ${thisSysPath}/tmp_allPhysBinned.C ${thisSysPath}/allPhysBinned.C
    exit 1
fi
    
#Final Cleanup
mv ${thisSysPath}/tmp_allPhysBinned.C ${thisSysPath}/allPhysBinned.C
rm ${thisSysPath}/{allPhysBinned.C.bak,pullDist.C.bak,physBinnedPeriod.C.bak}	


