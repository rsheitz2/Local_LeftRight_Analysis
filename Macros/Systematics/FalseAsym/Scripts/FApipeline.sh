#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo "This is the final script to run all macros in this folder "
    echo "This Script loops over all periods and all physics binnings"
    echo "Recommended to skip steps outside of this directory for the time being..."
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
###Additional settings
production="slot1"
phiPhotonCut=0.53 #HMDY=0.1866, #LowM_AMDY=0.195
##Step ONE settings
fitMrangeType="HMDY"
nBins=3
binFile=${analysisPath}/Presents/DATA/RealData/HMDY/BinValues/slot1WAll_HMDY_${nBins}bins.txt
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

###Additional settings
#production="slot1"
#phiPhotonCut=0.53 #HMDY=0.187, #LowM_AMDY=0.195
###Step ONE settings  #LowM_AMDY
#fitMrangeType="LowM_AMDY"
#nBins=5
#binFile=${analysisPath}/Presents/DATA/RealData/JPsi/BinValues/slot1WAll_JPsi${binRange}_${nBins}bins.txt
#hbins=150
#fitMmin=2.00  #true fit mass range
#fitMmax=8.50  #true fit mass range
#binRange="25_43"
###Step TWO settings
#process="JPsi" 
#LR_Mmin=2.00
#LR_Mmax=5.00
#whichFit="thirteen"
###Step THREE settings


additionalCuts=phiS$phiPhotonCut #add and new cuts here.  This should include all cuts used



##Setup___ last line (70) to search setup
lrMrange="${LR_Mmin}_${LR_Mmax}"
fitMrange="${fitMmin}_${fitMmax}"
periods=("W07" "W08" "W09" "W10" "W11" "W12" "W13" "W14" "W15" "WAll")

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
echo "Fit considered:         ${whichFit}"
echo " "
echo "______Step THREE settings____"
echo " "

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi
if [ ${Steps} -lt 10 ]; then 
    echo "Recommended to skip steps outside of this directory for the time being..."
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

for p in ${periods[@]}; do
    echo ""
    echo "Period $p"
    echo ""
    
    ${sysFApath}/Scripts/ChangeScripts/changePipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${p} $fitMrangeType $nBins $hbins \
		$fitMmin $fitMmax null_physBinned $process $LR_Mmin $LR_Mmax $whichFit $binRange $binFile ${production} ${phiPhotonCut} ${additionalCuts}

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
done

#Final clean up pipeline
mv ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh
rm ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh.bak


#############
#acceptanceFourTargRatio setup/checks    
echo " "
echo "acceptanceFourTargRatio_pipeline.sh"
echo " "

#allPhysBinnedpipeline changes
cp ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh

for p in ${periods[@]}; do
    echo ""
    echo "Period $p"
    echo ""
    
    ${sysFApath}/Scripts/ChangeScripts/changePipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh ${p} $fitMrangeType $nBins $hbins \
		$fitMmin $fitMmax null_physBinned $process $LR_Mmin $LR_Mmax $whichFit $binRange $binFile ${production} ${phiPhotonCut} ${additionalCuts}

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
done

#Final clean up pipeline
mv ${sysFApath}/Scripts/LoopScripts/tmp_allPhysBinnedpipeline.sh ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh
rm ${sysFApath}/Scripts/LoopScripts/allPhysBinnedpipeline.sh.bak
