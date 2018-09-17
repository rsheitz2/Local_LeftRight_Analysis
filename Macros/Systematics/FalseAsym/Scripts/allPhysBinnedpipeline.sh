#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by physics binning AN from"
    echo " a FALSE 4 target geometric mean calculation"
    echo "Scrip works by calling FalseAsym/Scripts/FApipeline.sh"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:    leftRight_byTarget"
    echo "Step TWO:    functMFit.C or trueCount.C"
    echo "     Determine AN by target"
    echo "Step THREE:   falseGeoMean4Targ_targFlip.C && falseGeoMean_splitTarg.C"
    echo "     Make false asymmetries"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo ""
    exit 1
fi

#General variables
Steps=$1
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant




##Setup___  first line to seach setup
##########
##Step ONE settings
period="WAll"
fitMrangeType="LowM_AMDY"
nBins=5
hbins=150
fitMmin=1.0  #true fit mass range
fitMmax=8.50  #true fit mass range
##Step TWO settings
process="JPsi"
LR_Mmin=2.80
LR_Mmax=3.50
physBinned=("xN" "xPi" "xF" "pT")
whichFit=("six" "six" "seven" "seven")
##Step THREE settings














##Setup___ last line to search setup
lrMrange="${LR_Mmin}_${LR_Mmax}"
fitMrange="${fitMmin}_${fitMmax}"

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
sysPath=${HOME}/Systematics/FalseAsym

#Intial save files to be changed
cp ${sysPath}/Scripts/FApipeline.sh ${sysPath}/Scripts/tmp_FApipeline.sh

for i in `seq 0 3`; do
    #Step THREE falseAsymmetry setup/checks    
    #FApipeline changes
    ${sysPath}/Scripts/changeFApipeline.sh ${period} $fitMrangeType $nBins $hbins $fitMmin $fitMmax ${physBinned[$i]} $process $LR_Mmin $LR_Mmax ${whichFit[$i]}
    #Execute
    ${sysPath}/Scripts/FApipeline.sh 1 >> ${sysPath}/Scripts/log_FApipeline.txt
    if [ $? != 0 ]; then
	echo "FApipeline did not execute well"
	mv ${sysPath}/Scripts/FApipeline.sh ${sysPath}/Scripts/FApipeline.sh.bak
	mv ${sysPath}/Scripts/tmp_FApipeline.sh ${sysPath}/Scripts/FApipeline.sh
	exit 1
    else
	rm ${sysPath}/Scripts/log_FApipeline.txt
    fi
done

#Clean up
mv ${sysPath}/Scripts/tmp_FApipeline.sh ${sysPath}/Scripts/FApipeline.sh
rm ${sysPath}/Scripts/FApipeline.sh.bak

#Step FIVE allSysErrorFA.C setup/checks
echo " "
echo "allSysErrorFA.C"
echo " "
stepFivePath=${sysPath}/Data/allSysErrorFA
stepFiveFile=${stepFivePath}/allSysErrorFA_${whichFit}
if [ ${whichFit} == "true" ]; then
    stepFiveFile+=_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned[$i]}${nBins}.root
else
    stepFiveFile+=${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned[$i]}${nBins}.root
fi

if [ ! -f ${stepFiveFile} ]; then
    #allSysErrorFA.C changes
    cp ${sysPath}/allSysErrorFA.C ${sysPath}/tmp_allSysErrorFA.C
    ${sysPath}/Scripts/changeAllSysErrorFA.sh $nBins ${period}_${fitMrangeType} $hbins $process $lrMrange $fitMrange ${whichFit[@]}

    #Execute allSysErrorFA.C
    root -l -b -q "${sysPath}/allSysErrorFA.C(1)" >> ${sysPath}/log_allSysErrorFA.txt
    if [ $? != 0 ]; then
	echo "allSysErrorFA.C did not execute well"
	mv ${sysPath}/allSysErrorFA.C ${sysPath}/allSysErrorFA.C.bak
	mv ${sysPath}/tmp_allSysErrorFA.C ${sysPath}/allSysErrorFA.C
	exit 1
    else
	#cleanup allSysErrorFA.C
	mv ${sysPath}/tmp_allSysErrorFA.C ${sysPath}/allSysErrorFA.C
	rm ${sysPath}/allSysErrorFA.C.bak
	rm ${sysPath}/log_allSysErrorFA.txt
    fi
else
    echo "File  $stepFiveFile  file already exist"
fi
