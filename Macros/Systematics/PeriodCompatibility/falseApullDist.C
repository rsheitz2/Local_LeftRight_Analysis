#include "include/helperFunctions.h"

void falseApullDist(TString start=""){
  //Setup_______________
  /*const Int_t nBins =3;//HMDY
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  const Int_t nBins =4;//JPsi
  TString fitMrangeType ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi
  TString lrMrange ="2.87_3.38";
  TString fitMrange ="2.87_3.38";
  TString binRange ="29_34"; //"25_43";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  //To write change value in /include/oneGraph.C and edit for writting
  //Setup_______________

  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/TargFlip";
  TString pathPull = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility";

  if (start==""){
    cout <<"\nThis Script computes the pull distributions for an L/R asymmetry";
    cout <<"\nThe asymmetry is determined ";
    cout << "(# periods)*(# physics kinematics)*(# physics kinematic bins)\n";
    cout << "Input data should be a TGraphErrors for each period and";
    cout << "physical kinematic";
    cout << "\n\nCurrent Setup:" << endl;
    cout << "Data comes from:                 " << pathFA << endl;
    cout << "Number of physics kinematic bins " << nBins << endl;
    cout << "Mass type considered:            " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Binned in Mass range:       " << binRange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Fit considered:     " << whichFit << endl;
    cout << "Production considered:   " << production << endl;
    cout << "additional cuts considered:   " << additionalCuts << endl;
    cout << "\nUsage:" << endl;
    cout <<"root \'PullDist(1)\'\n" << endl;
    exit(EXIT_FAILURE);
  }
  
  if (whichFit == "true" ){
    if (lrMrange != fitMrange){
      cout << "Error: lrMrange != fitMrange with whichFit==true" << endl;
      exit(EXIT_FAILURE);
    }
  }
      
  TString fStart = Form("%s/falseGeoMean4Targ_true_", pathFA.Data());
  TString fMiddle =
    Form("_%s_%s%s_%s", fitMrangeType.Data(), process.Data(), fitMrange.Data(),
	 binRange.Data());
  TString fEnd = Form("%i_%s_%s.root", nBins, production.Data(),
			   additionalCuts.Data());
  cout << "\n" << endl;
  
  /*//2Targ pol false asymmetry
  TString faPol = "falseAN_pol";
  cout << "2 Targ pol false asymmetry" << endl;
  cout << "Jura/Saleve acceptance" << endl;
  gROOT->ProcessLine(Form(".x %s/include/makePull.C(\"%s\", \"%s\", \"%s\", \"%s\", %i)",
			  pathPull.Data(), fStart.Data(), fMiddle.Data(),
			  fEnd.Data(), faPol.Data(), nBins));//*/

  /*//2Targ subper false asymmetry
  cout << "2 Targ subper false asymmetry" << endl;
    cout << "upstream Saleve/Jura, downstream Jura/Saleve" << endl;
  TString faSubper = "falseAN_subper";
  gROOT->ProcessLine(Form(".x %s/include/makePull.C(\"%s\", \"%s\", \"%s\", \"%s\", %i)",
			  pathPull.Data(), fStart.Data(), fMiddle.Data(),
			  fEnd.Data(), faSubper.Data(), nBins));//*/

  //upstream false asymmetry
  TString faUpS = "falseAN_2Targ_upS";
  cout << "Normal geomean upstream false asymmetry" << endl;
  gROOT->ProcessLine(Form(".x %s/include/makePull.C(\"%s\", \"%s\", \"%s\", \"%s\", %i)",
			  pathPull.Data(), fStart.Data(), fMiddle.Data(),
			  fEnd.Data(), faUpS.Data(), nBins));//*/

  /*//downstream false asymmetry
  cout << "Normal geomean downstream false asymmetry" << endl;
  TString faDownS = "falseAN_2Targ_downS";
  gROOT->ProcessLine(Form(".x %s/include/makePull.C(\"%s\", \"%s\", \"%s\", \"%s\", %i)",
			  pathPull.Data(), fStart.Data(), fMiddle.Data(),
			  fEnd.Data(), faDownS.Data(), nBins));//*/

  cout << "\n" << endl;
  
}
