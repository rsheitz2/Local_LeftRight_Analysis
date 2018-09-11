#include "include/helperFunctions.h"


void CalAbsSys(TGraphErrors *gFA_1, TGraphErrors *gFA_2, Double_t *systematics){
  Double_t *yFA_1 = gFA_1->GetY();
  Double_t *yFA_2 = gFA_2->GetY();

  Double_t *e_yFA_1 = gFA_1->GetEY();
  Double_t *e_yFA_2 = gFA_2->GetEY();

  if (gFA_1->GetN() != gFA_2->GetN() ) {
    cout << "Binning number inconsistency" << endl;
    exit(EXIT_FAILURE);
  }

  for (Int_t i=0; i<gFA_1->GetN(); i++) {
    Double_t sys_1 = yFA_1[i]/e_yFA_1[i];
    if (sys_1 < 0) sys_1 *= -1.0;

    Double_t sys_2 = yFA_2[i]/e_yFA_2[i];
    if (sys_2 < 0) sys_2 *= -1.0;
    
    systematics[i] = (sys_1 + sys_2)/2.0;
  }
}


void sysErrorFA(TString start =""){
  //Setup_______________
  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/TargFlip/";

  const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString physBinned ="xF";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit ="true";

  Bool_t toWrite =false;
  //Setup_______________
  
  if (start==""){//Basic info
    cout << "\nScript calculates systematic error due false AN asymmetries\n";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C or trueCount.C ";
    cout <<"-> fasleGeoMean4Targ_targFlips.C  -> sysErrorFA.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'sysErrorFA(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "False asymmetry data coming from:            " << pathFA << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Period and Mass type considered:   " << period_Mtype << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    if (whichFit == "true"){
      if (lrMrange != fitMrange){
	cout << "L/R mass range does not equal fit mass range for true fit\n";
	exit(EXIT_FAILURE);
      }
      cout << "L/R mass range:     " << lrMrange << endl;
    }
    else{
      cout << "LR integral mass range:     " << lrMrange << endl;
      cout << "Fit mass range:     " << fitMrange << endl;
    }
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //File name setup && get file
  TString inFileFA;
  if (whichFit=="true"){
    inFileFA =
      Form("falseGeoMean4Targ_true_%s_%s%s_%s%i.root", period_Mtype.Data(),
	   process.Data(), lrMrange.Data(), physBinned.Data(), nBins);
  }
  else{
    inFileFA =
      Form("falseGeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbin.root", whichFit.Data(),
	   fitMrange.Data(),period_Mtype.Data(), process.Data(),lrMrange.Data(),
	   physBinned.Data(), nBins, hbins);
  }
  
  TFile *fFA = TFile::Open(pathFA+inFileFA);
  if (!fFA){
    cout << "False asym file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data from files
  TGraphErrors *gFA_pol = (TGraphErrors*)fFA->Get("falseAN_pol");
  TGraphErrors *gFA_subper = (TGraphErrors*)fFA->Get("falseAN_subper");

  //Make FA/statisics average
  Double_t FA_sys[nBins];
  Double_t *xvals = gFA_pol->GetX();
  CalAbsSys(gFA_pol, gFA_subper, FA_sys);
    
  //Draw systematic error
  TGraph *gSys = new TGraphErrors(nBins, xvals, FA_sys);
  SetUp(gSys);

  TCanvas* c1 = new TCanvas();
  gSys->Draw("AP");
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/sysError";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =
      Form("%s/sysErrorFA_true_%s_%s%s_%s%i.root", thisDirPath.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   physBinned.Data(), nBins);
  }
  else{
    fOutput =
      Form("%s/sysErrorFA_%s%s_%s_%s%s_%s%i.root", thisDirPath.Data(),
	   whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
	   process.Data(), lrMrange.Data(), physBinned.Data(), nBins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    gSys->Write("gSys");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "False asymmetry data from:            " << pathFA << endl;
  cout << "False asymmetry file used:                   " << inFileFA << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Period and Mass type considered:   " << period_Mtype << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
  cout << "AN physical process:        " << process << endl;
  if (whichFit == "true"){
    cout << "L/R mass range:     " << lrMrange << endl;
  }
  else{
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
  }
  cout << "Which fit considered:       " << whichFit << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
