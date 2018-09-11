#!/bin/bash

if [ $# -lt 1 ]; then
    echo "" 
    echo "This script gets all the data for and runs the macros in AN_calculation folder"
    echo ""
    echo "Steps in this script"
    echo "Step ONE:    leftRight_byTarget"
    echo "Step TWO:    functMFit.C"
    echo "     Determine AN by target"
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
    
    ##Setup___
    ##########
    ##Step ONE settings
    period="W10"
    fitMrangeType="LowM_AMDY"
    nBins=5
    hbins=150
    fitMmin=1.0
    fitMmax=8.5
    binFile=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/AMDY/BinValues/WAll_AMDY_5bins.txt
    InputData=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/W10_LowM_AMDY.root
    ##Step TWO settings
    process="JPsi"
    LR_Mmin=2.8
    LR_Mmax=3.5
    physBinned="xF"
    ##Step THREE settings
    whichFit="25_85_six"
    
    echo ""
    echo "______Step ONE settings____"
    echo "Period:   ${period}"
    echo "Fit mass range type:  ${fitMrangeType}"
    echo "Number of kinematic bins:   ${nBins}"
    echo "Number of histogram bins in M distribution:  ${hbins}"
    echo "Minimum mass range:                          ${fitMmin}"
    echo "Maximum mass range:                          ${fitMmax}"
    echo "Binning file:"
    echo "    ${binFile}"
    echo "Input data:"
    echo "    ${InputData}"
    echo " "
    echo "______Step TWO settings____"
    echo "Integrated physics process:   ${process}"
    echo "Min integration range:        ${LR_Mmin}"
    echo "Max integration range:        ${LR_Mmax}"
    echo "Kinematic binning type:       ${physBinned}"
    echo "    Remember to check what fitting is being used"
    echo "    Fit settings are not input in this script"
    echo " "
    echo "______Step THREE settings____"
    echo "Fit considered:         ${whichFit}"
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
    stepOne_OutnoCorr=${stepOne_OutData}${period}"_"${fitMrangeType}"_"${nBins}"bins_noCorr.root"
    stepOne_OutData+=${period}"_"${fitMrangeType}"_"${nBins}"bins_"${hbins}"hbin.root"
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
    echo "functMFit.C"

    #Step TWO checks
    pathTwo=${pathOne}/Macros/AN_calculation
    stepTwo_OutData=${pathTwo}/Data/functMFit/functMFit_${whichFit}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}_${hbins}hbin
    stepTwo_noCorrData=${stepTwo_OutData}"_noCorr.root"
    stepTwo_OutData+="_corr.root"
    if [ ! -f ${stepTwo_OutData} ]; then
	echo "Making stepTwoData"
	echo ""
	echo ""

	#Prepare functMFit.C macro settings
	cp ${pathTwo}/functMFit.C ${pathTwo}/tmpTwo.C
	sed -i.bak "s/nBins =[1-255]/nBins =${nBins}/" ${pathTwo}/functMFit.C
	sed -i.bak "s/period_Mtype =.*;/period_Mtype =\"${period}_${fitMrangeType}\";/" ${pathTwo}/functMFit.C
	sed -i.bak "s/hbins =.*;/hbins =${hbins};/" ${pathTwo}/functMFit.C
	sed -i.bak "s/TString physBinned =.*;/TString physBinned =\"${physBinned}\";/" ${pathTwo}/functMFit.C
	sed -i.bak "s/TString process = .*;/TString process = \"${process}\";/" ${pathTwo}/functMFit.C
	sed -i.bak "s/Double_t LR_Mmin =.*, LR_Mmax =.*;/Double_t LR_Mmin =${LR_Mmin}, LR_Mmax =${LR_Mmax};/" ${pathTwo}/functMFit.C
	sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${pathTwo}/functMFit.C

	#Execute pol corrected and pol unCorr
	root -l -b -q "${pathTwo}/functMFit.C(true, 1)"
	root -l -b -q "${pathTwo}/functMFit.C(false, 1)"

	#Cleanup folder
	rm ${pathTwo}/functMFit.C.bak
	mv ${pathTwo}/tmpTwo.C ${pathTwo}/functMFit.C
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

    stepThree_OutData=${pathTwo}/Data/GeoMean4Targ/GeoMean4Targ_${whichFit}_${period}_${fitMrangeType}_${process}${LR_Mmin}_${LR_Mmax}_${physBinned}${nBins}
    stepThree_OutData+=".root"
    if [ ! -f ${stepThree_OutData} ]; then
	echo "Making stepThreeData"
	echo ""
	echo ""

	#Prepare functMFit.C macro settings
	cp ${pathTwo}/GeoMean4Targ.C ${pathTwo}/tmpThree.C
	sed -i.bak "s/nBins =[1-255]/nBins =${nBins}/" ${pathTwo}/GeoMean4Targ.C
	sed -i.bak "s/period_Mtype =.*;/period_Mtype =\"${period}_${fitMrangeType}\";/" ${pathTwo}/GeoMean4Targ.C
	sed -i.bak "s/hbins =.*;/hbins =${hbins};/" ${pathTwo}/GeoMean4Targ.C
	sed -i.bak "s/TString physBinned = .*;/TString physBinned = \"${physBinned}\";/" ${pathTwo}/GeoMean4Targ.C
	sed -i.bak "s/TString process = .*;/TString process = \"${process}\";/" ${pathTwo}/GeoMean4Targ.C
	sed -i.bak "s/TString lrMrange = .*;/TString lrMrange = \"${LR_Mmin}_${LR_Mmax}\";/" ${pathTwo}/GeoMean4Targ.C
	sed -i.bak "s/Bool_t toWrite =.*;/Bool_t toWrite =true;/" ${pathTwo}/GeoMean4Targ.C

	#Execute 
	root -l -b -q "${pathTwo}/GeoMean4Targ.C(1)"

	#Cleanup folder
	rm ${pathTwo}/GeoMean4Targ.C.bak
	mv ${pathTwo}/tmpThree.C ${pathTwo}/GeoMean4Targ.C
    else
	echo "StepThreeData already exist"
    fi
    
fi
