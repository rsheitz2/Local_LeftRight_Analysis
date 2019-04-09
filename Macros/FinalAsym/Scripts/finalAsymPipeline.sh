#!/bin/bash

if [ $# -ne 1 ]; then
    echo "" 
    echo "Pipeline script for \"integrated\" and \"all physics\" binning in:"
    echo "   wAvg.C,  physBinnedData.C"
    echo "   Integrated is already included"
    echo ""
    echo "To do:"
    echo "Include other macros sysWphysBinnedData.C and sysOverStat.C"
    echo ""
    echo "To run this script provide as an argument:"
    echo "     \"h\" or \"0\" to see the current settings"
    echo "     Or enter a number of steps to take (i.e. 1)"
    echo ""
    exit 1
fi

#General variables
Steps=$1



##Setup___  first line (20) to seach setup
##########
####Additional settings
#production="slot1"
#phiPhotonCut=0.0 
#Mtype="HMDY"
#nBins=3 #integrated taken care of
#hbins=150
#lrMrange="4.30_8.50"  #true fit mass range
#binRange="43_85"
#process="DY"
#fitMrange="4.30_8.50"
#whichFit="true"


###JPsi
production="slot1"
phiPhotonCut=0.0 #HMDY=0.187, #LowM_AMDY=0.195
Mtype="LowM_AMDY"
nBins=3 #integrated taken care of
hbins=150
lrMrange="2.87_3.38" #"3.08_3.17" #"2.87_3.38" 
binRange=#"29_34" #"31_32" #"29_34" #25_43
process="JPsi"
fitMrange=$lrMrange
whichFit="true"


additionalCuts=phiS$phiPhotonCut #add and new cuts here.  This should include all cuts used

##Setup___ last line (50) to search setup

#Basic Setup
physBinned=("xN" "xPi" "xF" "pT" "M")
HOME=/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym

echo ""
echo "______settings____"
echo "Fit mass range type:  ${Mtype}"
echo "Number of kinematic bins:   ${nBins}"
echo "Number of histogram bins in M distribution:  ${hbins}"
echo "Left/right mass range:        ${lrMrange}"
echo "Integrated physics process:   ${process}"
echo "Fit considered:   ${whichFit}" 
echo "Mass integration range:            ${fitMrange}"
echo ""
echo "Kinematic binning types:       ${physBinned[@]}"
echo ""

if [ ${Steps} == "h" ] || [ ${Steps} -lt 1 ]; then #Help option, output settings
    exit 1
fi

#Intial save files to be changed
cp ${HOME}/wAvg.C ${HOME}/wAvg_tmp.C

#Integrated
${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/wAvg.C 1 $Mtype $hbins "xN" $process $lrMrange $fitMrange $binRange $whichFit $production $additionalCuts

#Execute
root -q -b -l ${HOME}/wAvg.C >> ${HOME}/wAvg_log.txt
if [ $? != 0 ]; then
    echo "wAvg.C did not execute well"
    mv ${HOME}/wAvg.C ${HOME}/wAvg.C.bak
    mv ${HOME}/wAvg_tmp.C ${HOME}/wAvg.C
    exit 1
else
    rm ${HOME}/wAvg_log.txt
fi

#Loop over phys binning
for phys in ${physBinned[@]}; do
    echo ""
    echo "Physics Binned $phys"
    echo ""
    
    #Macro changes
    ${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/wAvg.C $nBins $Mtype $hbins $phys $process $lrMrange $fitMrange $binRange $whichFit $production $additionalCuts
    #Execute
    root -q -b -l ${HOME}/wAvg.C >> ${HOME}/wAvg_log.txt
    if [ $? != 0 ]; then
	echo "wAvg.C did not execute well"
	mv ${HOME}/wAvg.C ${HOME}/wAvg.C.bak
	mv ${HOME}/wAvg_tmp.C ${HOME}/wAvg.C
	exit 1
    else
	rm ${HOME}/wAvg_log.txt
    fi
    
done

#Clean up
mv ${HOME}/wAvg_tmp.C ${HOME}/wAvg.C
rm ${HOME}/wAvg.C.bak


#Final plot
cp ${HOME}/physBinnedData.C ${HOME}/physBinnedData_tmp.C
${HOME}/Scripts/ChangeScripts/changeMacro.sh ${HOME}/physBinnedData.C $nBins $Mtype $hbins $phys $process $lrMrange $fitMrange $binRange $whichFit $production \
       $additionalCuts
#Execute
root -b -q -l "${HOME}/physBinnedData.C(1)" >> ${HOME}/physBinnedData_log.txt
if [ $? != 0 ]; then
    echo "physBinnedData.C did not execute well"
    mv ${HOME}/physBinnedData.C ${HOME}/physBinnedData.C.bak
    mv ${HOME}/physBinnedData_tmp.C ${HOME}/physBinnedData.C
    exit 1
else
    rm ${HOME}/physBinnedData_log.txt
fi
    
mv ${HOME}/physBinnedData_tmp.C ${HOME}/physBinnedData.C
rm ${HOME}/physBinnedData.C.bak
