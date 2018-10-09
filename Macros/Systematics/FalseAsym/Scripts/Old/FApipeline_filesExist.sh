#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by period for physics binning AN from"
    echo " a FALSE 4 target geometric mean calculation"
    echo "Scrip works by calling AN_calculation/Scripts/pipeline.sh"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:    leftRight_byTarget, systematic_leftRight"
    echo "Step TWO:    functMFit.C or trueCount.C, "
    echo "     Determine AN by target"
    echo "Step THREE:   falseGeoMean4Targ_targFlip.C, falseGeoMean4Targ_splitTarg.C"
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
binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/Wall_JPsi25_43_5bins.txt
hbins=150
fitMmin=2.00  #true fit mass range
fitMmax=7.50  #true fit mass range
binRange="25_43"
##Step TWO settings
physBinned="pT"
process="JPsi"
LR_Mmin=2.90
LR_Mmax=3.30
whichFit="ten"
##Step THREE settings












##Setup___ last line (60) to search setup
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

#Step One systematic_leftRight setup/checks
########
echo " "
echo "systematic_leftRight"
echo " "
sysLeftRight_Out=${analysisPath}"/Local_LeftRight_Analysis/Data/systematic_leftRight_"${period}"_"${fitMrangeType}${fitMmin}_${fitMmax}"_"
sysLeftRight_Out+=${nBins}"bins"${binRange}"_"${hbins}"hbin.root"

#systematic_leftRight Pipeline changes
cp ${sysPath}/Scripts/systematic_leftRight_pipeline.sh ${sysPath}/Scripts/tmp_systematic_leftRight_pipeline.sh
${sysPath}/Scripts/changeSystematic_LeftRight_Pipeline.sh ${period} $fitMrangeType $nBins $hbins $fitMmin $fitMmax $binFile $binRange

#Execute
${sysPath}/Scripts/systematic_leftRight_pipeline.sh 1 >> ${sysPath}/Scripts/log_systematic_leftRight_pipeline.txt
if [ $? != 0 ]; then
    echo "systematic_leftRight_pipeline.sh did not execute well"
    mv ${sysPath}/Scripts/systematic_leftRight_pipeline.sh ${sysPath}/Scripts/systematic_leftRight_pipeline.sh.bak
    mv ${sysPath}/Scripts/tmp_systematic_leftRight_pipeline.sh ${sysPath}/Scripts/systematic_leftRight_pipeline.sh
    exit 1
else
    #clean up systematic_leftRight_pipeline
    mv ${sysPath}/Scripts/tmp_systematic_leftRight_pipeline.sh ${sysPath}/Scripts/systematic_leftRight_pipeline.sh
    rm ${sysPath}/Scripts/systematic_leftRight_pipeline.sh.bak
    rm ${sysPath}/Scripts/log_systematic_leftRight_pipeline.txt
fi

#Step One/Two leftRight_byTarget/funcMFit.C or trueCount.C
#############
#Pipeline changes
echo " "
echo "Pipeline"
echo " "
cp ${aNPath}/Scripts/pipeline.sh ${aNPath}/Scripts/tmp_pipeline.sh
${aNPath}/Scripts/changePipeline.sh ${period} $fitMrangeType $nBins $hbins ${physBinned} $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax ${whichFit} ${binRange} ${binFile}

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
exit 1 #cleanup
#Step THREE
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

if [ ! -f ${stepThreeFile} ]; then
    #falseGeoMean4Targ_targFlips.C changes
    cp ${sysFApath}/falseGeoMean4Targ_targFlips.C ${sysFApath}/tmp_falseGeoMean4Targ_targFlips.C
    ${sysFApath}/Scripts/changeFalseTargFlips.sh $nBins ${period}_${fitMrangeType} $hbins $physBinned $process $lrMrange $fitMrange ${whichFit}

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
else
    echo "File  $stepThreeFile  file already exist"
fi
echo "Script stops early!!!"
exit 1 #cleanup
#Step THREE falseAsymmetry_splitTarg setup/checks
#To do sysFunctMFit
echo " "
echo "falseGeoMean4Targ_splitTarg.C"
echo " "
FA_splitTarg_out=${sysFApath}/Data/SplitTarg_
if [ ${whichFit} == "true" ]; then
    FA_splitTarg_out+=true_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}.root
else
    FA_splitTarg_out+=${whichFit}${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}.root
fi

if [ ! -f ${FA_splitTarg_out} ]; then
    #falseGeoMean4Targ_splitTarg.C changes
    cp ${sysFApath}/falseGeoMean4Targ_splitTarg.C ${sysFApath}/tmp_falseGeoMean4Targ_splitTarg.C
    ${sysFApath}/Scripts/changeFalseSplitTarg.sh $nBins ${period}_${fitMrangeType} $hbins $physBinned $process $lrMrange $fitMrange ${whichFit} $binRange

    #Execute falseGeoMean4Targ_splitTarg.C
    root -l -b -q "${sysFApath}/falseGeoMean4Targ_splitTarg.C(1)"  >> ${sysFApath}/log_falseGeoMean4Targ_splitTarg.txt
    if [ $? != 0 ]; then
	echo "falseGeoMean4Targ_splitTarg.C did not execute well"
	mv ${sysFApath}/falseGeoMean4Targ_splitTarg.C ${sysFApath}/falseGeoMean4Targ_splitTarg.C.bak
	mv ${sysFApath}/tmp_falseGeoMean4Targ_splitTarg.C ${sysFApath}/falseGeoMean4Targ_splitTarg.C
	exit 1
    else
	#clean up falseGeoMean4Targ_splitTarg
	mv ${sysFApath}/tmp_falseGeoMean4Targ_splitTarg.C ${sysFApath}/falseGeoMean4Targ_splitTarg.C
	rm ${sysFApath}/falseGeoMean4Targ_splitTarg.C.bak
	rm ${sysFApath}/log_falseGeoMean4Targ_splitTarg.txt
    fi
else
    echo "File  $FA_splitTarg_out  file already exist"
fi

##Step FOUR
#############
#sysErrorFA.C setup/checks
echo " "
echo "sysErrorFA.C"
echo " "
stepFourPath=${sysFApath}/Data/sysError
stepFourFile=${stepFourPath}/sysErrorFA_${whichFit}
if [ ${whichFit} == "true" ]; then
    stepFourFile+=_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}.root
else
    stepFourFile+=${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}.root
fi

if [ ! -f ${stepFourFile} ]; then
    #sysErrorFA.C changes
    cp ${sysFApath}/sysErrorFA.C ${sysFApath}/tmp_sysErrorFA.C
    ${sysFApath}/Scripts/changeSysError.sh $nBins ${period}_${fitMrangeType} $hbins $physBinned $process $lrMrange $fitMrange ${whichFit}

    #Execute sysErrorFA.C
    root -l -b -q "${sysFApath}/sysErrorFA.C(1)" >> ${sysFApath}/log_sysErrorFA.txt
    if [ $? != 0 ]; then
	echo "sysErrorFA.C did not execute well"
	mv ${sysFApath}/sysErrorFA.C ${sysFApath}/sysErrorFA.C.bak
	mv ${sysFApath}/tmp_sysErrorFA.C ${sysFApath}/sysErrorFA.C
	exit 1
    else
	#clean up sysErrorFA
	mv ${sysFApath}/tmp_sysErrorFA.C ${sysFApath}/sysErrorFA.C
	rm ${sysFApath}/sysErrorFA.C.bak
	rm ${sysFApath}/log_sysErrorFA.txt
    fi
else
    echo "File  $stepFourFile  file already exist"
fi
