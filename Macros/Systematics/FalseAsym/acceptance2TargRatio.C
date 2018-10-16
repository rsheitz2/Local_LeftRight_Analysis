#include "include/helperFunctions.h"


void Cal_AccRatioAndError(TGraphErrors *g, Double_t *alpha, Double_t *error){
  Double_t *yval = g->GetY();
  Double_t *eYval = g->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    alpha[i] = (1+yval[i])/(1-yval[i]);

    error[i] = 2.0*eYval[i]/( (1-yval[i])*(1-yval[i]) );
  }
}


void acceptance2TargRatio(TString start =""){
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
    TString whichFit ="ten";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
FA2targ_ratioCals/";
  
  if (start==""){//Basic info
    cout << "\nScript calculates acceptance ratios for 4targ geomean" << endl;
    cout << "Calculation input comes from 2targ geomean false asymmetries\n";
    cout << "\nInput needed = TGraphErrors of alpha/beta acceptance ratios\n";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "systematic_leftRight -> FA2targ_ratioCals.C " << endl;
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
    inputFiles = Form("FA2targ_ratioCals_true_%s_%s%s_%s%i.root",
		      period_Mtype.Data(), process.Data(), fitMrange.Data(),
		      physBinned.Data(), nBins);
  }
  else{
    inputFiles = Form("FA2targ_ratioCals_%s%s_%s_%s%s_%s%i_%ihbins.root",
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
  TGraphErrors *g_upSup =   (TGraphErrors*)fFAs->Get("falseAN_upSup");
  TGraphErrors *g_upSdown = (TGraphErrors*)fFAs->Get("falseAN_upSdown");

  TCanvas* c1 = new TCanvas();
  g_upSup->Draw("AP"); g_upSdown->Draw("Psame");
  DrawLine(g_upSup, 0.0);
  
  //Calculation alpha/beta ratios and Draw them
  Double_t a_upSup[nBins], a_upSdown[nBins];
  Double_t eA_upSup[nBins], eA_upSdown[nBins];
  Cal_AccRatioAndError(g_upSup, a_upSup, eA_upSup);
  Cal_AccRatioAndError(g_upSdown, a_upSdown, eA_upSdown);

  Double_t *xvals = g_upSup->GetX();
  Double_t ex[nBins] = {0.0};
  TGraphErrors *g_alpha_upSup = new TGraphErrors(nBins, xvals, a_upSup, ex,
						 eA_upSup);
  TGraphErrors *g_alpha_upSdown = new TGraphErrors(nBins, xvals, a_upSdown, ex,
						   eA_upSdown);
  SetUp(g_alpha_upSup); SetUp(g_alpha_upSdown);

  TCanvas* c2 = new TCanvas();
  g_alpha_upSup->Draw("AP"); g_alpha_upSdown->Draw("Psame");
  g_alpha_upSup->SetMarkerColor(kBlue); g_alpha_upSdown->SetMarkerColor(kRed);
  g_alpha_upSup->GetYaxis()->SetRangeUser(0.9, 1.1);
  DrawLine(g_alpha_upSup, 1.0);
  
  //Calculate and draw final ratio kappa
  Double_t final_a[nBins], eFinal_a[nBins];
  for (Int_t i=0; i<nBins; i++) {
    final_a[i] = TMath::Sqrt( a_upSup[i]*a_upSdown[i] );

    Double_t error = eA_upSup[i]*eA_upSup[i]/(a_upSup[i]*a_upSup[i]);
    error += eA_upSdown[i]*eA_upSdown[i]/(a_upSdown[i]*a_upSdown[i]);
    eFinal_a[i] = 0.5*TMath::Sqrt(final_a[i])*TMath::Sqrt(error);
  }
  
  TGraphErrors* g_final_a = new TGraphErrors(nBins, xvals, final_a, ex, eFinal_a);
  SetUp(g_final_a);

  TCanvas* c3 = new TCanvas();
  g_final_a->Draw("AP"); g_final_a->GetYaxis()->SetRangeUser(0.8, 1.2);
  DrawLine(g_final_a, 1.0);

    //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
acceptanceRatio";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/acceptanceRatio_true_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), period_Mtype.Data(), process.Data(),
		  lrMrange.Data(), physBinned.Data(), nBins);
  }
  else{
    fOutput =Form("%s/acceptanceRatio_%s%s_%s_%s%s_%s%i_%ihbins.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_upSup->Write("falseAN_upSup");
    g_upSdown->Write("falseAN_upSdown");
    g_final_a->Write("final_alpha");
    g_alpha_upSup->Write("alpha_upSup");
    g_alpha_upSdown->Write("beta_upSdown");
    
    fResults->Close();
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathFA << endl;
  cout << "Input P corrected data:        " << inputFiles[0] << endl;
  cout << "Input unCorr P data:           " << inputFiles[1] << endl;
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
