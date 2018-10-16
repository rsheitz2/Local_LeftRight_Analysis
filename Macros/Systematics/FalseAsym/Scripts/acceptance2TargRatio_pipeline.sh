#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by period for acceptance2TargRatio.C"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:    systematic_leftRight_pipeline.sh"
    echo "Step TWO:    FA2targ_ratioCals.C (fullTargSysFunctMFit.C) "
    echo "Step THREE:  acceptance2TargRatio.C"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo -n "     Or enter a number greater than 0 and less than 11 (i.e. 1)"
    echo " to get all data needed for this macro"
    echo "      a number greater than 10 to only run macro (i.e. 11)"
    echo ""
    exit 1
fi

#General variables
Steps=$1
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant




##Setup___  first line (25) to seach setup
##########
##Step ONE settings
period="WAll"
fitMrangeType="LowM_AMDY"
nBins=5
binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/WAll_JPsi25_43_${nBins}bins.txt
hbins=150
fitMmin=2.00
fitMmax=7.50
binRange="25_43"
##Step TWO settings
physBinned="xN"
process="JPsi"
LR_Mmin=2.90
LR_Mmax=3.30
whichFit="ten"


















##Setup___ last line (60) to search setup
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
echo "Kinematic binning type:       ${physBinned}"
echo "Fit considered:         ${whichFit}"
echo " "

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Setup
HOME=${analysisPath}/Local_LeftRight_Analysis/Macros
sysPath=${HOME}/Systematics
sysFApath=${sysPath}/FalseAsym

if [ ${Steps} -lt 11 ]; then #get dependency files for macro
    #Step One FA2targ_ratioCals_pipeline.sh
    ########
    echo " "
    echo "FA2targ_ratioCals_pipeline.sh"
    echo " "

    #FA2targ_ratioCals Pipeline changes
    cp ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh ${sysFApath}/Scripts/tmp_FA2targ_ratioCals_pipeline.sh
    ${sysFApath}/Scripts/changeGenericPipeline.sh ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh ${period} $fitMrangeType $nBins $hbins $fitMmin $fitMmax \
		${physBinned[$i]} $process $LR_Mmin $LR_Mmax ${whichFit[$i]} ${binRange} ${binFile}
    
    #Execute
    ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh 1 >> ${sysFApath}/Scripts/log_FA2targ_ratioCals_pipeline.txt
    if [ $? != 0 ]; then
	echo "FA2targ_ratioCals_pipeline.sh did not execute well"
	mv ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh.bak
	mv ${sysFApath}/Scripts/tmp_FA2targ_ratioCals_pipeline.sh ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh
	exit 1
    else
	#clean up FA2targ_ratioCals_pipeline
	mv ${sysFApath}/Scripts/tmp_FA2targ_ratioCals_pipeline.sh ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh
	rm ${sysFApath}/Scripts/FA2targ_ratioCals_pipeline.sh.bak
	rm ${sysFApath}/Scripts/log_FA2targ_ratioCals_pipeline.txt
    fi
fi

#Step Two acceptance2TargRatio.C
#############
#acceptance2TargRatio changes
echo " "
echo "acceptance2TargRatio.C"
echo " "
cp ${sysFApath}/acceptance2TargRatio.C ${sysFApath}/tmp_acceptance2TargRatio.C
${sysFApath}/Scripts/changeGeneric.sh acceptance2TargRatio ${nBins} $period_Mtype $hbins ${physBinned} $process $lrMrange $fitMrange $binRange ${whichFit}

#Execute
root -l -b -q "${sysFApath}/acceptance2TargRatio.C(1)" >> ${sysFApath}/log_acceptance2TargRatio.txt
if [ $? != 0 ]; then
    echo "acceptance2TargRatio.C did not execute well"
    mv ${sysFApath}/acceptance2TargRatio.C ${sysFApath}/acceptance2TargRatio.C.bak
    mv ${sysFApath}/tmp_acceptance2TargRatio.C ${sysFApath}/acceptance2TargRatio.C
    exit 1
else
    #clean up acceptance2TargRatio
    mv ${sysFApath}/tmp_acceptance2TargRatio.C ${sysFApath}/acceptance2TargRatio.C
    rm ${sysFApath}/acceptance2TargRatio.C.bak
    rm ${sysFApath}/log_acceptance2TargRatio.txt
fi
