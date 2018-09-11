#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo -n "This script gets all the data setup by period and by physics binning for AN from"
    echo " a 4 target geometric mean calculation"
    echo "Script works by calling physBinnedPipeline.sh"
    echo ""
    echo "Steps in this script (from AN_calculation folder)"
    echo "Step ONE:    leftRight_byTarget"
    echo "Step TWO:    functMFit.C or trueCount.C"
    echo "     Determine AN by target"
    echo "Step THREE:  GeoMean4Targ.C"
    echo "     Calculate acceptance free AN"
    echo "Step FOUR:   pullDist.C"
    echo " "
    echo "To Do:"
    echo "Make mass binning work"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number greater than 0 (i.e. 1)"
    echo ""
    
else
    Steps=$1

    #General variables
    analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant

    ##Setup___
    ##########
    ##Step ONE settings
    fitMrangeType="HMDY"
    nBins=3
    hbins=150
    fitMmin=4.30 #0.01 precision  
    fitMmax=8.50 #0.01 precision
    ##Step TWO settings
    process="DY"
    LR_Mmin=4.30 #0.01 precision  #does nothing for true fit
    LR_Mmax=8.50 #0.01 precision  #does noting for true fit
    ##Step THREE settings
    
    ##Step FOUR settings
    physBinned=("xN" "xPi" "xF" "pT")
    fits=("true" "true" "true" "true")

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
    echo "Min integration range:        ${LR_Mmin}"
    echo "Max integration range:        ${LR_Mmax}"
    echo " "
    echo "______Step THREE settings____"
    echo " "
    echo "______Step FOUR settings____"
    echo "Physics binnings considered    ${physBinned[*]}"
    echo "Fits per physics binning       ${fits[*]}"
    echo " "

    if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
	exit 1
    fi

    #Basic Setup
    HOME=${analysisPath}/Local_LeftRight_Analysis/Macros
    sysPath=${HOME}/Systematics/PeriodCompatibility
    if [ ${fits[0]} == "true" ]; then
	if [ ${LR_Mmin} != ${fitMmin} ] || [ ${LR_Mmax} != ${LR_Mmax} ]; then
	    echo "set LR_Mmin equal to fitMmin  &&   LR_Mmin equal to fitMmin"
	    echo "just in case"
	    exit 1
	fi
    fi
    lrMrange=${LR_Mmin}_${LR_Mmax}
    fitMrange=${fitMmin}_${fitMmax}
    
    #Physics Binning loop
    #cp ${sysPath}/Scripts/physBinnedPipeline.sh ${sysPath}/Scripts/tmp_physBinnedPipeline.sh 
    #for i in `seq 0 3`; do
    #	echo "Performing physic binning    ${physBinned[$i]}"
    #	echo "Fit type for  ${physBinned[$i]}    is     ${fits[$i]}"
    #	echo " "
    #
    #	##physBinnedPipeline changes
    #	${sysPath}/Scripts/changePhysBinnedPipeline.sh $fitMrangeType $nBins $hbins ${physBinned[$i]} $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax ${fits[$i]}
    #	
    #	#Execute
    #	${sysPath}/Scripts/physBinnedPipeline.sh 1 
    #	wait
    #	
    #done
    #
    ##clean up pipeline
    #mv ${sysPath}/Scripts/tmp_physBinnedPipeline.sh ${sysPath}/Scripts/physBinnedPipeline.sh
    #rm ${sysPath}/Scripts/physBinnedPipeline.sh.bak

    echo " "
    echo "pullDist.C"
    echo " "
    #Step FOUR checks
    physBinnedPeriod_fileBase=${sysPath}/Data/physBinned/physBinnedPeriod_
    for i in `seq 0 3`; do
	physBinnedPeriodData=${physBinnedPeriod_fileBase}
	if [ ${fits[$i]} == "true" ]; then
	    physBinnedPeriodData+=true_${fitMrangeType}_${process}${fitMmin}_${fitMmax}_${physBinned[$i]}${nBins}.root
	else
	    physBinnedPeriodData+=${fits[$i]}${fitMmin}_${fitMmax}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned[$i]}${nBins}_${hbins}hbin.root
	fi
	
	if [ ! -f ${physBinnedPeriodData} ]; then
	    echo ${physBinnedPeriodData}
	    echo "physBinnedPeriod.C data does not exist for    ${physBinned[$i]}"
	    exit 1
	fi
    done

    pullDistData=${sysPath}/Data/pullDist/pullDist_
    if [ ${fits[0]} == "true" ]; then
	pullDistData+=${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${nBins}bins.root
    else
	pullDistData+=${fitMmin}_${fitMmax}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${nBins}bins_${hbins}hbin.root
    fi
    if [ ! -f ${pullDistData} ]; then
    	#Step FOUR changes
    	cp ${sysPath}/pullDist.C ${sysPath}/tmp_pullDist.C
    	${sysPath}/Scripts/changePullDist.sh $nBins $fitMrangeType $hbins $process $lrMrange $fitMrange ${physBinned[@]} ${fits[@]}
    	
    	#Execute 
    	root -l -b -q "${sysPath}/pullDist.C(1)"
    
    	#Cleanup folder
    	mv ${sysPath}/tmp_pullDist.C ${sysPath}/pullDist.C
    	rm ${sysPath}/pullDist.C.bak	
    else
    	echo "${pullDistData} already exist"
    fi
    
fi
