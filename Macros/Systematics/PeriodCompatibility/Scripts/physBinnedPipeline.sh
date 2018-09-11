#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by period for physic binning AN from"
    echo " a 4 target geometric mean calculation"
    echo "Scrip works by calling AN_calculation/Scripts/pipeline.sh"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:    leftRight_byTarget"
    echo "Step TWO:    functMFit.C or trueCount.C"
    echo "     Determine AN by target"
    echo "Step THREE:  GeoMean4Targ.C"
    echo "     Calculate acceptance free AN"
    echo "Step FOUR:   physBinnedPeriod.C"
    echo "   plot kinematices by period"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo ""
    
else
    Steps=$1

    #General variables
    analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant


    ##Setup___  first line to seach setup
    ##########
    ##Step ONE settings
    fitMrangeType="HMDY"
    nBins=3
    hbins=150
    fitMmin=4.30
    fitMmax=8.50
    ##Step TWO settings
    physBinned="pT"
    process="DY"
    LR_Mmin=2.80
    LR_Mmax=3.50
    whichFit="true"
    ##Step THREE settings
    
    ##Step FOUR settings
    lrMrange="${LR_Mmin}_${LR_Mmax}"
    fitMrange="${fitMmin}_${fitMmax}"











    ##Setup___ last line to search setup

    echo ""
    echo "______Step ONE settings____"
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
    echo "Kinematic binning type:       ${physBinned}"
    echo " "
    echo "______Step THREE settings____"
    echo " "
    echo "______Step FOUR settings____"
    echo "Min/Max integration range:         ${lrMrange}"
    echo "Min/Max fit range:         ${fitMrange}"
    echo " "

    if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
	exit 1
    fi

    #Basic Setup
    HOME=${analysisPath}/Local_LeftRight_Analysis/Macros
    aNPath=${HOME}/AN_calculation
    sysPath=${HOME}/Systematics/PeriodCompatibility
    
    #Period loop
    periods=(07 08 09 10 11 12 13 14 15)
    cp ${aNPath}/Scripts/pipeline.sh ${aNPath}/Scripts/tmp_pipeline.sh
    for p in ${periods[*]} ; do
    	echo "Performing period  W$p"
    	echo " "
    
    	#Pipeline changes
    	${aNPath}/Scripts/changePipeline.sh W$p $fitMrangeType $nBins $hbins $physBinned $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax $whichFit
	
    	#Execute
    	${aNPath}/Scripts/pipeline.sh 3 | tee ${sysPath}/Scripts/OutputBackup/pipeline_W${p}_output.txt
    	wait
    	
    done
    
    #clean up pipeline
    mv ${aNPath}/Scripts/tmp_pipeline.sh ${aNPath}/Scripts/pipeline.sh

    echo " "
    echo "physBinnedPeriod.C"
    echo " "
    #Step FOUR checks
    stepThreePath=${aNPath}/Data/GeoMean4Targ
    stepThreeFileBase=${stepThreePath}/GeoMean4Targ_${whichFit}
    stepFourFileBase=${sysPath}/physBinnedPeriod_${whichFit}
    if [ ${whichFit} != "true" ]; then
	stepThreeFileBase += ${fitMmin}_${fitMmax}
	stepFourFileBase += ${fitMrange}
    fi
    
    for p in ${periods[*]} ; do
	if [ ! -f ${stepThreeFileBase}_W${p}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}.root ]; then
	    echo ${stepThreeFileBase}_W${p}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}.root
	    echo "GeoMean4Targ data not created for period:    W$p"
	    exit 1
	fi	
    done

    if [ ! -f ${stepFourFileBase}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}_${hbins}hbin.root ]; then
	#Step FOUR changes
	cp ${sysPath}/physBinnedPeriod.C ${sysPath}/tmp_phyBinnedPeriod.C
	${sysPath}/Scripts/changePhysBinnedPeriod.sh $nBins $fitMrangeType $hbins $physBinned $process $lrMrange $fitMrange $whichFit
	
	#Execute 
	root -l -b -q "${sysPath}/physBinnedPeriod.C(1)"

	#Cleanup folder
	mv ${sysPath}/tmp_phyBinnedPeriod.C ${sysPath}/physBinnedPeriod.C
	rm ${sysPath}/physBinnedPeriod.C.bak	
    else
	echo -n "physBinnedPeriod_${whichFit}${fitMrange}_${fitMrangeType}_${process}${lrMrange}_${physBinned}${nBins}_${hbins}hbin.root "
	echo "already exist"
    fi

    
    
fi
