#!/bin/bash

if [ $# -lt 1 ]; then
    echo "" 
    echo "This script gets all the data for and runs the macros in AN_calculation folder"
    echo ""
    echo "Steps in this script"
    echo "Step ONE:    leftRight_byTarget"
    echo "Step TWO:    functMFit.C or trueCounts.C"
    echo "     Determine AN by target by mass fitting or true counts"
    echo "Step THREE:  GeoMean4Targ.C"
    echo "     Calculate acceptance free AN"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" to see the current settings"
    echo "     Or enter the number of steps to proceed through (i.e. 1, 2, 3)"
    echo ""
    exit 1
fi
Steps=$1

#General variables
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant

##Setup___ (25) first setup search line 
##########
##Step ONE settings
period="W07"
fitMrangeType="LowM_AMDY"
nBins=3
binFile="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/JPsi/BinValues/Wall_JPsi25_43_3bins.txt"
hbins=150
fitMmin=2.87
fitMmax=3.38
binRange="25_43"
##Step TWO settings
physBinned="xN"
process="JPsi"
LR_Mmin=2.87
LR_Mmax=3.38
whichFit="true"
##Additional settings
production="slot1"
phiPhotonCut="0.0"

additionalCuts="phiS0.0"



##Setup ends, (50) last setup search line

#Default setup
if [ ${whichFit} == "true" ]; then
    hbins=150
fi

InputData=${analysisPath}/Presents/DATA/RealData/
if [ ${fitMrangeType} == "HMDY" ]; then #Speed optimization for lower data set HMDY
    InputData+=${fitMrangeType}/
    if [ ${production} == "slot1" ]; then
	InputData+=${production}${period}_${fitMrangeType}.root
    else
	InputData+=${period}_${fitMrangeType}.root
    fi
else
    InputData+=LowM_AMDY/
    if [ $production == "slot1" ]; then
	InputData+=${production}${period}_LowM_AMDY.root
    else
	InputData+=${period}_LowM_AMDY.root
    fi
fi

Mmin=1.00
Mmax=8.50

#Echo out setup
echo ""
echo "______Step ONE settings____"
echo "Period:   ${period}"
echo "Fit mass range type:  ${fitMrangeType}"
echo "Number of kinematic bins:   ${nBins}"
echo "Number of histogram bins in M distribution:  ${hbins}"
echo "Binning file:"
echo "    ${binFile}"
echo "Input data:"
echo "    ${InputData}"
echo "Min left/right mass range:        ${fitMmin}"
echo "Max left/right mass range:        ${fitMmax}"
echo " "
echo "______Step TWO settings____"
echo "Integrated physics process:   ${process}"
if [ ${whichFit} == "true" ]; then
    echo "Min left/right mass range:        ${fitMmin}"
    echo "Max left/right mass range:        ${fitMmax}"
    Mmin=${fitMmin}
    Mmax=${fitMmax}
else
    echo "Min integration range:        ${LR_Mmin}"
    echo "Max integration range:        ${LR_Mmax}"
    echo "Minimum fit mass range:                          ${fitMmin}"
    echo "Maximum fit mass range:                          ${fitMmax}"	
fi
echo "Kinematic binning type:       ${physBinned}"
echo "Fit considered:         ${whichFit}"
echo " "
echo "______Additional settings____"
echo "production used:   $production"
echo "phi photon target frame cut:    $phiPhotonCut"
echo "All cuts used:     $additionalCuts"
echo " "
echo " "

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi
#Step ONE
echo "_______Step ONE_____"
echo "      leftRight_byTarget"

#Step ONE checks
pathOne=${analysisPath}/Local_LeftRight_Analysis
if [ ! -x ${pathOne}/leftRight_byTarget ]; then
    echo "leftRight_byTarget does not exist"
    exit 1
fi
if [ ! -f ${binFile} ] || [ ! -f ${InputData} ]; then
    echo "Step ONE file does not exist:"
    echo "${binFile}"
    echo "${InputData}"
    exit 1
fi

#Execute Step ONE if files don't already exist
stepOne_OutData=${pathOne}"/Data/leftRight_byTarget_" 
stepOne_OutData+=${period}"_"${fitMrangeType}${Mmin}_${Mmax}"_"${nBins}"bins"${binRange}"_"${hbins}"hbin"_${production}_${additionalCuts}".root"
if [ ! -f ${stepOne_OutData} ]; then
    echo "Making stepOneData polarization correct data:"
    echo ""
    echo ""
    echo ""
    ${pathOne}/leftRight_byTarget -i${Mmin} -a${Mmax} -Q${stepOne_OutData} -b${binFile} -Z${hbins} -E${phiPhotonCut} -f${InputData} \
	      >> ${pathOne}/log_leftRight_byTarget.txt
    if [ $? != 0 ]; then
	echo "leftRight_byTarget did not execute well"
	exit 1
    else
	#clean up leftRight_byTarget
	rm ${pathOne}/log_leftRight_byTarget.txt
    fi
else
    echo "StepOneData polarization correct already exist"
fi
echo " "

if [ ${Steps} -lt 2 ]; then 
    exit 0
fi
#Step TWO
echo "_______Step TWO_____"
pathTwo=${pathOne}/Macros/AN_calculation
stepTwo_OutData=${pathTwo}/Data
if [ ${whichFit} == "true" ]; then
    echo "trueCount.C"
    stepTwo_OutData+=/trueCount/trueCount_${period}_${fitMrangeType}_${process}${fitMmin}_${fitMmax}_${binRange}${physBinned}${nBins}
elif [ ${whichFit} == "MC" ]; then
    echo "mcMFit.C"
    echo "mcMFit not setup correctly at the moment"
    exit 1
    stepTwo_OutData+=/mcMFit/mcMFit_${whichFit}${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}_${hbins}hbin
else
    echo "functMFit.C"
    stepTwo_OutData+=/functMFit/functMFit_${whichFit}${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${binRange}${physBinned}${nBins}_${hbins}hbin
fi
stepTwo_OutData+=_${production}_${additionalCuts}"_corr.root"

echo "Making stepTwoData"
echo ""
echo ""
if [ ${whichFit} == "true" ]; then
    #Prepare trueCount.C macro settings
    cp ${pathTwo}/trueCount.C ${pathTwo}/tmpTrueTwo.C
    ${pathTwo}/Scripts/ChangeScripts/changeTrueCount.sh $nBins $period $fitMrangeType $fitMmin $fitMmax $physBinned $process $binRange $production $additionalCuts
    
    #Execute pol corrected and pol unCorr
    root -l -b -q "${pathTwo}/trueCount.C(1)" >> ${pathTwo}/log_trueCount.txt
    if [ $? != 0 ]; then
	echo "trueCount.C did not execute well"
	mv ${pathTwo}/trueCount.C ${pathTwo}/trueCount.C.bak	    
	mv ${pathTwo}/tmpTrueTwo.C ${pathTwo}/trueCount.C
	exit 1
    else
	#clean up trueCount.C
	mv ${pathTwo}/tmpTrueTwo.C ${pathTwo}/trueCount.C
	rm ${pathTwo}/trueCount.C.bak
	rm ${pathTwo}/log_trueCount.txt
    fi
elif [ ${whichFit} == "MC" ];then
    echo "Currently MC fit is not supported"
    exit 1
    #Prepare mcMFit.C macro settings
    cp ${pathTwo}/mcMFit.C ${pathTwo}/tmpMcTwo.C
    ${pathTwo}/Scripts/ChangeScripts/changeMcMFit.sh $nBins $period $fitMrangeType $hbins $physBinned $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax $whichFit $binRange
    
    #Execute pol corrected and pol unCorr
    root -l -b -q "${pathTwo}/mcMFit.C(1)" >> ${pathTwo}/log_mcMFit.txt
    if [ $? != 0 ]; then
	echo "mcMFit.C did not execute well"
	mv ${pathTwo}/mcMFit.C ${pathTwo}/mcMFit.C.bak	    
	mv ${pathTwo}/tmpMcTwo.C ${pathTwo}/mcMFit.C
	exit 1
    else
	#clean up mcMFit.C
	mv ${pathTwo}/tmpMcTwo.C ${pathTwo}/mcMFit.C
	rm ${pathTwo}/mcMFit.C.bak
	rm ${pathTwo}/log_mcMFit.txt
    fi
else
    #Prepare functMFit.C macro settings
    cp ${pathTwo}/functMFit.C ${pathTwo}/tmpFunctTwo.C
    ${pathTwo}/Scripts/ChangeScripts/changeFunctMFit.sh $nBins $period $fitMrangeType $hbins $physBinned $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax $whichFit $binRange \
	      $production $additionalCuts

    #Execute pol corrected and pol unCorr
    root -l -b -q "${pathTwo}/functMFit.C(1)" >> ${pathTwo}/log_functMFit.txt
    if [ $? != 0 ]; then
	echo "functMFit.C did not execute well"
	mv ${pathTwo}/functMFit.C ${pathTwo}/functMFit.C.bak	    
	mv ${pathTwo}/tmpFunctTwo.C ${pathTwo}/functMFit.C
	exit 1
    else
	#clean up functMFit.C
	mv ${pathTwo}/tmpFunctTwo.C ${pathTwo}/functMFit.C
	rm ${pathTwo}/functMFit.C.bak
	rm ${pathTwo}/log_functMFit.txt
    fi
fi
echo " "


if [ ${Steps} -lt 3 ]; then 
    exit 0
fi
#Step THREE
echo "_______Step THREE_____"
echo "GeoMean4Targ.C"

#Step THREE checks
if [ ! -f ${stepTwo_OutData} ]; then
    echo "stepTwoData was not created well!"
    exit 1
fi

stepThree_OutData=${pathTwo}/Data/GeoMean4Targ/GeoMean4Targ_${whichFit}${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${binRange}${physBinned}${nBins}
stepThree_OutData+="_${production}_${additionalCuts}.root"
echo "Making stepThreeData"
echo ""
echo ""

#Prepare GeoMean4Targ.C macro settings
cp ${pathTwo}/GeoMean4Targ.C ${pathTwo}/tmpThree.C
${pathTwo}/Scripts/ChangeScripts/changeGeoMean4Targ.sh $nBins $period $fitMrangeType $hbins $physBinned $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax $whichFit $binRange \
	  $production $additionalCuts

#Execute 
root -l -b -q "${pathTwo}/GeoMean4Targ.C(1)" >> ${pathTwo}/log_GeoMean4Targ.txt
if [ $? != 0 ]; then
    echo "GeoMean4Targ.C did not execute well"
    mv ${pathTwo}/GeoMean4Targ.C ${pathTwo}/GeoMean4Targ.C.bak	    
    mv ${pathTwo}/tmpThree.C ${pathTwo}/GeoMean4Targ.C
    exit 1
else
    #clean up GeoMean4Targ.C
    mv ${pathTwo}/tmpThree.C ${pathTwo}/GeoMean4Targ.C
    rm ${pathTwo}/GeoMean4Targ.C.bak
    rm ${pathTwo}/log_GeoMean4Targ.txt
fi
