#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by period for physics binning acceptanceFourTargRatio.C"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:    falseGeoMean4Targ_targFlips.C"
    echo "Step TWO:    acceptanceFourTargRatio.C"
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
fitMrangeType="HMDY"
nBins=3
binFile=${analysisPath}/Presents/DATA/RealData/HMDY/BinValues/WAll_HMDY_${nBins}bins.txt
hbins=150
fitMmin=4.30
fitMmax=8.50
binRange="43_85"
##Step TWO settings
physBinned="xPi"
process="DY"
LR_Mmin=4.30
LR_Mmax=8.50
whichFit="true"
















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
    #Step One systematic_leftRight setup/checks
    ########
    echo " "
    echo "FApipeline"
    echo " "
    
    #systematic_leftRight Pipeline changes
    cp ${sysFApath}/Scripts/FApipeline.sh ${sysFApath}/Scripts/tmp_FApipeline.sh
    ${sysFApath}/Scripts/changeFApipeline.sh ${period} $fitMrangeType $nBins $hbins $fitMmin $fitMmax $physBinned $process $LR_Mmin $LR_Mmax $whichFit $binRange $binFile

    #Execute
    ${sysFApath}/Scripts/FApipeline.sh 1 >> ${sysFApath}/Scripts/log_FApipeline.txt
    if [ $? != 0 ]; then
	echo "FApipeline.sh did not execute well"
	mv ${sysFApath}/Scripts/FApipeline.sh ${sysFApath}/Scripts/FApipeline.sh.bak
	mv ${sysFApath}/Scripts/tmp_FApipeline.sh ${sysFApath}/Scripts/FApipeline.sh
	exit 1
    else
	#clean up FApipeline
	mv ${sysFApath}/Scripts/tmp_FApipeline.sh ${sysFApath}/Scripts/FApipeline.sh
	rm ${sysFApath}/Scripts/FApipeline.sh.bak
	rm ${sysFApath}/Scripts/log_FApipeline.txt
    fi
fi

#Step Two acceptanceFourTargRatio.C (fullTargSysFunctMFit.C)
#############
#acceptanceFourTargRatio changes
echo " "
echo "acceptanceFourTargRatio.C"
echo " "
cp ${sysFApath}/acceptanceFourTargRatio.C ${sysFApath}/tmp_acceptanceFourTargRatio.C
${sysFApath}/Scripts/changeGeneric.sh acceptanceFourTargRatio ${nBins} $period_Mtype $hbins ${physBinned} $process $lrMrange $fitMrange $binRange ${whichFit}

#Execute
root -l -b -q "${sysFApath}/acceptanceFourTargRatio.C(1)" >> ${sysFApath}/log_acceptanceFourTargRatio.txt
if [ $? != 0 ]; then
    echo "acceptanceFourTargRatio.C did not execute well"
    mv ${sysFApath}/acceptanceFourTargRatio.C ${sysFApath}/acceptanceFourTargRatio.C.bak
    mv ${sysFApath}/tmp_acceptanceFourTargRatio.C ${sysFApath}/acceptanceFourTargRatio.C
    exit 1
else
    #clean up acceptanceFourTargRatio
    mv ${sysFApath}/tmp_acceptanceFourTargRatio.C ${sysFApath}/acceptanceFourTargRatio.C
    rm ${sysFApath}/acceptanceFourTargRatio.C.bak
    rm ${sysFApath}/log_acceptanceFourTargRatio.txt
fi
