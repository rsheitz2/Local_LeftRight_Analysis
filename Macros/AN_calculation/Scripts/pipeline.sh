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
    
else
    Steps=$1
    
    #General variables
    analysisPath=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant
    
    ##Setup___ first setup search line
    ##########
    ##Step ONE settings
    period="WAll"
    fitMrangeType="HMDY"
    nBins=3
    hbins=150
    fitMmin=4.30
    fitMmax=8.50
    ##Step TWO settings
    physBinned="pT"
    process="DY"
    LR_Mmin=4.30 #does nothing with whichFit==true
    LR_Mmax=8.50 #does nothing with whichFit==true
    whichFit="true"
    ##Step THREE settings









    ##Setup ends, last setup search line
    
    if [ ${whichFit} == "true" ]; then
	hbins=150 #this should always be the case
    fi
    
    binFile=${analysisPath}/Presents/DATA/RealData/${fitMrangeType}/BinValues/WAll_${fitMrangeType}_${nBins}bins.txt
    InputData=${analysisPath}/Presents/DATA/RealData/
    if [ ${fitMrangeType} == "HMDY" ]; then
	InputData+=${fitMrangeType}/${period}_${fitMrangeType}.root
    else
	InputData+=LowM_AMDY/${period}_LowM_AMDY.root
    fi
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
    else
	echo "Min integration range:        ${LR_Mmin}"
	echo "Max integration range:        ${LR_Mmax}"
	echo "Minimum fit mass range:                          ${fitMmin}"
	echo "Maximum fit mass range:                          ${fitMmax}"	
    fi
    echo "Kinematic binning type:       ${physBinned}"
    echo "Fit considered:         ${whichFit}"
    echo " "
    echo "______Step THREE settings____"
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
    stepOne_OutnoCorr=${stepOne_OutData}${period}"_"${fitMrangeType}${fitMmin}_${fitMmax}"_"${nBins}"bins_noCorr.root"
    stepOne_OutData+=${period}"_"${fitMrangeType}${fitMmin}_${fitMmax}"_"${nBins}"bins_"${hbins}"hbin.root"
    if [ ! -f ${stepOne_OutData} ]; then
	echo "Making stepOneData polarization correct data:"
	echo ""
	echo ""
	echo ""
	${pathOne}/leftRight_byTarget -i${fitMmin} -a${fitMmax} -Q${stepOne_OutData} -b${binFile} -N${nBins} -Z${hbins} -f${InputData} 
    else
	echo "StepOneData polarization correct already exist"
    fi
    if [ ! -f ${stepOne_OutnoCorr} ]; then
	echo "Making stepOneData polarization unCorr data:"
	echo ""
	echo ""
	echo ""
	${pathOne}/leftRight_byTarget -i${fitMmin} -a${fitMmax} -Q${stepOne_OutnoCorr} -P -b${binFile} -N${nBins} -Z${hbins} -f${InputData}
    else
	echo "StepOneData pol unCorr already exist"
    fi
    echo " "
    

    if [ ${Steps} -lt 2 ]; then 
	exit 1
    fi
    #Step TWO
    echo "_______Step TWO_____"
    pathTwo=${pathOne}/Macros/AN_calculation
    stepTwo_OutData=${pathTwo}/Data
    if [ ${whichFit} == "true" ]; then
	echo "trueCount.C"
	stepTwo_OutData+=/trueCount/trueCount_${period}_${fitMrangeType}_${process}${fitMmin}_${fitMmax}_${physBinned}${nBins}
    else
	echo "functMFit.C"
	stepTwo_OutData+=/functMFit/functMFit_${whichFit}${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}_${hbins}hbin
    fi
    stepTwo_noCorrData=${stepTwo_OutData}"_noCorr.root"
    stepTwo_OutData+="_corr.root"
    
    #Step TWO checks
    if [ ! -f ${stepTwo_OutData} ]; then
	echo "Making stepTwoData"
	echo ""
	echo ""

	if [ ${whichFit} == "true" ]; then
	    #Prepare functMFit.C macro settings
	    cp ${pathTwo}/trueCount.C ${pathTwo}/tmpTrueTwo.C
	    ${pathTwo}/Scripts/changeTrueCount.sh $nBins $period $fitMrangeType $fitMmin $fitMmax $physBinned $process 
	    
	    #Execute pol corrected and pol unCorr
	    root -l -b -q "${pathTwo}/trueCount.C(true, 1)"
	    root -l -b -q "${pathTwo}/trueCount.C(false, 1)"

	    #Cleanup folder
	    rm ${pathTwo}/trueCount.C.bak
	    mv ${pathTwo}/tmpTrueTwo.C ${pathTwo}/trueCount.C
	else
	    #Prepare functMFit.C macro settings
	    cp ${pathTwo}/functMFit.C ${pathTwo}/tmpFunctTwo.C
	    ${pathTwo}/Scripts/changeFunctMFit.sh $nBins $period $fitMrangeType $hbins $physBinned $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax $whichFit

	    #Execute pol corrected and pol unCorr
	    root -l -b -q "${pathTwo}/functMFit.C(true, 1)"
	    root -l -b -q "${pathTwo}/functMFit.C(false, 1)"

	    #Cleanup folder
	    rm ${pathTwo}/functMFit.C.bak
	    mv ${pathTwo}/tmpFunctTwo.C ${pathTwo}/functMFit.C
	fi
    else
	echo "StepTwoData already exist"
    fi
    echo " "

    
    if [ ${Steps} -lt 3 ]; then 
	exit 1
    fi
    #Step THREE
    echo "_______Step THREE_____"
    echo "GeoMean4Targ.C"
    
    #Step THREE checks
    if [ ! -f ${stepTwo_OutData} ]; then
	echo "stepTwoData was not created well!"
	exit 1
    fi

    stepThree_OutData=${pathTwo}/Data/GeoMean4Targ/GeoMean4Targ_${whichFit}${fitMmin}_${fitMmax}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}
    stepThree_OutData+=".root"
    if [ ! -f ${stepThree_OutData} ]; then
	echo "Making stepThreeData"
	echo ""
	echo ""

	#Prepare functMFit.C macro settings
	cp ${pathTwo}/GeoMean4Targ.C ${pathTwo}/tmpThree.C
	${pathTwo}/Scripts/changeGeoMean4Targ.sh $nBins $period $fitMrangeType $hbins $physBinned $process $LR_Mmin $LR_Mmax $fitMmin $fitMmax $whichFit
	
	#Execute 
	root -l -b -q "${pathTwo}/GeoMean4Targ.C(1)"

	#Cleanup folder
	rm ${pathTwo}/GeoMean4Targ.C.bak
	mv ${pathTwo}/tmpThree.C ${pathTwo}/GeoMean4Targ.C
    else
	echo "StepThreeData already exist"
    fi
    
fi