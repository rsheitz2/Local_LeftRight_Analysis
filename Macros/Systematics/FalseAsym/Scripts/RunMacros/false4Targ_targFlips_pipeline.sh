#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by period for physics binning AN from"
    echo " a FALSE 4 target geometric mean calculation"
    echo "Scrip works by calling AN_calculation/Scripts/pipeline.sh"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:     pipeline.sh"
    echo "Step TWO:     falseGeoMean4Targ_targFlip.C"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo "     Or enter a number greater than 10 to skip pipeline.sh step (i.e. 11)"
    echo ""
    exit 1
fi

#General variables
Steps=$1
analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant

##Setup___  first line (25) to seach setup
##########
###Additional settings
production="slot1"
phiPhotonCut=0.53 #HMDY=0.1866, #LowM_AMDY=0.195
##Step ONE settings
period="WAll"  #HMDY
fitMrangeType="HMDY"
nBins=3
binFile=${analysisPath}/Presents/DATA/RealData/HMDY/BinValues/WAll_HMDY_${nBins}bins.txt
hbins=150
fitMmin=4.30  #true fit mass range
fitMmax=8.50  #true fit mass range
binRange="43_85"
##Step TWO settings
physBinned="xN"
process="DY" 
LR_Mmin=4.30
LR_Mmax=8.50
whichFit="true"
##Step THREE settings


###Additional settings
#production="slot1"
#phiPhotonCut=0.53 #HMDY=0.1866, #LowM_AMDY=0.195
###Step ONE settings  #LowM_AMDY
#period="WAll"
#fitMrangeType="LowM_AMDY"
#nBins=5
#binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/WAll_JPsi25_43_${nBins}bins.txt
#hbins=150
#fitMmin=2.90  #true fit mass range
#fitMmax=3.30  #true fit mass range
#binRange="25_43"
###Step TWO settings
#physBinned="xN"
#process="JPsi" 
#LR_Mmin=2.90
#LR_Mmax=3.30
#whichFit="eight"
###Step THREE settings

additionalCuts=phiS$phiPhotonCut #add and new cuts here.  This should include all cuts used

##Setup___ last line (70) to search setup
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
echo "______Step THREE settings____"
echo " "

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Basic Setup
HOME=${analysisPath}/Local_LeftRight_Analysis/Macros
aNPath=${HOME}/AN_calculation
sysPath=${HOME}/Systematics
sysFApath=${sysPath}/FalseAsym

if [ ${Steps} -lt 10 ]; then #Do all previous steps
    #Step One systematic_leftRight setup/checks
    ########
    #Pipeline changes
    echo " "
    echo "Pipeline"
    echo " "
    cp ${aNPath}/Scripts/pipeline.sh ${aNPath}/Scripts/tmp_pipeline.sh
    ${aNPath}/Scripts/ChangeScripts/changePipeline.sh ${period} $fitMrangeType $nBins $hbins ${physBinned} $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax ${whichFit} \
	     ${binRange} ${binFile} ${production} ${phiPhotonCut} ${additionalCuts}

    #Execute
    ${aNPath}/Scripts/pipeline.sh 2 >> ${aNPath}/Scripts/log_pipeline.txt
    if [ $? != 0 ]; then
	echo "pipeline.sh did not execute well"
	mv ${aNPath}/Scripts/pipeline.sh ${aNPath}/Scripts/pipeline.sh.bak
	mv ${aNPath}/Scripts/tmp_pipeline.sh ${aNPath}/Scripts/pipeline.sh
	exit 1
    else
	#clean up pipeline
	mv ${aNPath}/Scripts/tmp_pipeline.sh ${aNPath}/Scripts/pipeline.sh
	rm ${aNPath}/Scripts/pipeline.sh.bak
	rm ${aNPath}/Scripts/log_pipeline.txt
    fi
fi


#Step TWO
#############
#falseAsymmetry_targFlips setup/checks    
echo " "
echo "falseGeoMean4Targ_targFlips.C"
echo " "
stepThreePath=${sysFApath}/Data/TargFlip
stepThreeFile=${stepThreePath}/falseGeoMean4Targ_${whichFit}
if [ ${whichFit} == "true" ]; then
    stepThreeFile+=_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}.root
else
    stepThreeFile+=${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}.root
fi

#falseGeoMean4Targ_targFlips.C changes
cp ${sysFApath}/falseGeoMean4Targ_targFlips.C ${sysFApath}/tmp_falseGeoMean4Targ_targFlips.C
${sysFApath}/Scripts/ChangeScripts/changeMacro.sh falseGeoMean4Targ_targFlips $nBins ${period}_${fitMrangeType} $hbins $physBinned $process $lrMrange $fitMrange \
	    $binRange ${whichFit} ${production} ${additionalCuts}

#Execute falseGeoMean4Targ_targFlips.C
root -l -b -q "${sysFApath}/falseGeoMean4Targ_targFlips.C(1)" >> ${sysFApath}/log_falseGeoMean4Targ_targFlips.txt
if [ $? != 0 ]; then
    echo "falseGeoMean4Targ_targFlips.C did not execute well"
    mv ${sysFApath}/falseGeoMean4Targ_targFlips.C ${sysFApath}/falseGeoMean4Targ_targFlips.C.bak
    mv ${sysFApath}/tmp_falseGeoMean4Targ_targFlips.C ${sysFApath}/falseGeoMean4Targ_targFlips.C
    exit 1
else
    #cleanup falseGeoMean4Targ_targFlips.C
    mv ${sysFApath}/tmp_falseGeoMean4Targ_targFlips.C ${sysFApath}/falseGeoMean4Targ_targFlips.C
    rm ${sysFApath}/falseGeoMean4Targ_targFlips.C.bak
    rm ${sysFApath}/log_falseGeoMean4Targ_targFlips.txt
fi
