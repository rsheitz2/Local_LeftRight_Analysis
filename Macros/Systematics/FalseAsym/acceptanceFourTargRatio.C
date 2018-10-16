#include "include/helperFunctions.h"


void Cal_AccRatioAndError(TGraphErrors *g, Double_t *alpha, Double_t *error){
  //Calculate acceptance ratio without polarization values in asymmetry
  Double_t *yval = g->GetY();
  Double_t *eYval = g->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    alpha[i] = (1+yval[i])/(1-yval[i]);

    error[i] = 2.0*eYval[i]/( (1-yval[i])*(1-yval[i]) );
  }
}


void Cal_AccRatioAndError(TGraphErrors *g, TGraph *g_Pol,
			  Double_t *alpha, Double_t *error){
  Double_t *yval = g->GetY();
  Double_t *eYval = g->GetEY();
  Double_t *yPol = g_Pol->GetY();
  for (Int_t i=0; i<g->GetN(); i++) {
    alpha[i] = (1+yval[i]*yPol[i])/(1-yval[i]*yPol[i]);

    error[i] = 2.0*eYval[i]*yPol[i]/( (1-yval[i]*yPol[i])*(1-yval[i]*yPol[i]) );
  }
}


void acceptanceFourTargRatio(TString start =""){
  //Setup_______________
  const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit ="true";//*/
  
  /*const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.90_3.30";
  TString fitMrange ="2.00_7.50";
  //TString fitMrange ="2.90_3.30";
  TString whichFit ="ten";//*/

  Bool_t toWrite =false;
  //Setup_______________
  TString binRange ="25_34";//not used at the moment
  
  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
TargFlip/";

  if (start==""){//Basic info
    cout << "\nScript calculates acceptance ratios for 4targ geomean" << endl;
    cout << "Calculation input comes from 2targ geomean false asymmetries\n";
    cout << "\nInput needed = TGraphErrors of alpha/beta acceptance ratios\n";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "systematic_leftRight -> falseGeoMean4Targ_targFlips.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'acceptanceRatio(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << pathFA << endl;
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
  
  //File name setups && get files
  TString inputFiles;
  if (whichFit=="true"){
    inputFiles = Form("falseGeoMean4Targ_true_%s_%s%s_%s%i.root",
		      period_Mtype.Data(), process.Data(), fitMrange.Data(),
		      physBinned.Data(), nBins);
  }
  else{
    inputFiles = Form("falseGeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbins.root",
		      whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
		      process.Data(), lrMrange.Data(), physBinned.Data(),
		      nBins, hbins);
  }
  
  TFile *fFAs = TFile::Open(pathFA+inputFiles);
  if (!fFAs ){
    cout << "Ratios file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data from files/Draw FAs
  TGraphErrors *g_subper =   (TGraphErrors*)fFAs->Get("falseAN_subper");
  TGraph *g_Pol = (TGraph*)fFAs->Get("Polarization");

  if (!g_subper || !g_Pol){
    cout << "Graphs in FA file do not all exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  TCanvas* c1 = new TCanvas("False Asym", "False Asym");
  g_subper->Draw("AP"); DrawLine(g_subper, 0.0);
  g_subper->SetTitle("False_Asym");
  g_subper->GetYaxis()->SetRangeUser(-0.15, 0.15);
  
  //Calculation alpha/beta ratios and Draw them
  Double_t a_subper[nBins], eA_subper[nBins];
  Cal_AccRatioAndError(g_subper, g_Pol, a_subper, eA_subper);
  
  Double_t *xvals = g_subper->GetX();
  Double_t ex[nBins] = {0.0};
  TGraphErrors *g_alpha_subper = new TGraphErrors(nBins, xvals, a_subper, ex,
						 eA_subper);
  SetUp(g_alpha_subper); 

  TCanvas* c2 = new TCanvas("Acceptance Ratio", "Acceptance Ratio");
  g_alpha_subper->Draw("AP"); g_alpha_subper->SetTitle("Acc_ratio");
  g_alpha_subper->SetMarkerColor(kBlue); 
  g_alpha_subper->GetYaxis()->SetRangeUser(0.6, 1.4);
  DrawLine(g_alpha_subper, 1.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
acceptanceFourTargRatio";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/acceptanceFourTargRatio_true_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), period_Mtype.Data(), process.Data(),
		  lrMrange.Data(), physBinned.Data(), nBins);
  }
  else{
    fOutput =Form("%s/acceptanceFourTargRatio_%s%s_%s_%s%s_%s%i_%ihbins.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_subper->Write("falseAN_subper");
    g_Pol->Write("Polarization");
    g_alpha_subper->Write("alpha_subper");
    
    fResults->Close();
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathFA << endl;
  cout << "Input data:        " << inputFiles << endl;
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
