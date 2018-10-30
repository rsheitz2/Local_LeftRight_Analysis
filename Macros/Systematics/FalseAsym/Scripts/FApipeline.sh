#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo "This is the final script to run all macros in this folder "
    echo ""
    echo "Script gets all the data for all physics binning setup to get the systematics from acceptance changes"
    echo "    this is done for input period only"
    echo ""
    echo "Final output is in Data/acceptanceFourTargRatio/SystematicError"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo "     Or enter a number greater than 10 to skip steps outside of this directory (i.e. 11)"
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
process="DY"
LR_Mmin=4.30
LR_Mmax=8.50
whichFit="true"
##Step THREE settings


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
#process="JPsi" 
#LR_Mmin=2.90
#LR_Mmax=3.30
#whichFit="eight"
###Step THREE settings











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


#############
#falseAsymmetry_targFlips setup/checks    
echo " "
echo "falseGeoMean4Targ_targFlips_pipeline.sh"
echo " "

#allPhysBinnedpipeline changes
cp ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh
${sysFApath}/Scripts/ChangeScripts/changePipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${period} $fitMrangeType $nBins $hbins \
	    $fitMmin $fitMmax null_physBinned $process $LR_Mmin $LR_Mmax $whichFit $binRange $binFile

#Execute 
${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/RunMacros/false4Targ_targFlips_pipeline.sh $Steps \
	    >> ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline_log.txt
if [ $? != 0 ]; then
    echo "false4Targ_targFlips_pipeline.sh did not execute well"
    mv ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh.bak
    mv ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh
    exit 1
else
    rm ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline_log.txt
fi


#############
#acceptanceFourTargRatio setup/checks    
echo " "
echo "acceptanceFourTargRatio_pipeline.sh"
echo " "

#allPhysBinnedpipeline changes
cp ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh
${sysFApath}/Scripts/ChangeScripts/changePipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${period} $fitMrangeType $nBins $hbins \
	    $fitMmin $fitMmax null_physBinned $process $LR_Mmin $LR_Mmax $whichFit $binRange $binFile

#Execute 
${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/RunMacros/acceptanceFourTargRatio_pipeline.sh 11 \
	    >> ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline_log.txt
if [ $? != 0 ]; then
    echo "acceptanceFourTargRatio_pipeline.sh did not execute well"
    mv ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh.bak
    mv ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh
    exit 1
else
    rm ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline_log.txt
fi


#Final clean up pipeline
mv ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh
rm ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh.bak
