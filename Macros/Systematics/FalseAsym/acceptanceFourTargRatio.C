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
  /*const Int_t nBins =3;//HMDY
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString physBinned ="xPi";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/
  
  const Int_t nBins =5;//JPsi
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xPi";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.80_3.44";
  TString fitMrange ="2.80_3.44";
  TString binRange ="25_43";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";

  Bool_t toWrite =false;
  //Setup_______________
  
  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
TargFlip/";

  if (start==""){//Basic info
    cout << "\nScript calculates acceptance ratios for 4targ geomean" << endl;
    cout << "Calculation input comes from 4targ geomean false asymmetries\n";
    cout << "\nInput needed = TGraphErrors of alpha acceptance ratios\n";
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
    inputFiles = Form("falseGeoMean4Targ_%s_%s_%s%s_%s%s%i_%s_%s.root",
		      whichFit.Data(), period_Mtype.Data(), process.Data(),
		      lrMrange.Data(), binRange.Data(),physBinned.Data(), nBins,
		      production.Data(), additionalCuts.Data());
  }
  else{
    inputFiles =
      Form("falseGeoMean4Targ_%s%s_%s_%s%s_%s%s%i_%ihbins_%s_%s.root",
	   whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
	   process.Data(), lrMrange.Data(), binRange.Data(), physBinned.Data(),
	   nBins, hbins, production.Data(), additionalCuts.Data());
  }
  
  TFile *fFAs = TFile::Open(pathFA+inputFiles);
  if (!fFAs ){
    cout << "Ratios file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data from files/Draw FAs
  TGraphErrors *g_FAsubper =   (TGraphErrors*)fFAs->Get("falseAN_subper");
  TGraph *g_Pol = (TGraph*)fFAs->Get("Polarization");

  //Jura/Saleve acceptance
  TGraphErrors *g_FApol =   (TGraphErrors*)fFAs->Get("falseAN_pol");
  TGraphErrors *g_FAupS =   (TGraphErrors*)fFAs->Get("falseAN_2Targ_upS");
  TGraphErrors *g_FAdownS =   (TGraphErrors*)fFAs->Get("falseAN_2Targ_downS");

  if (!g_FAsubper || !g_FApol){
    cout << "Graphs in FA file do not all exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  TCanvas* c1 = new TCanvas("False Asym", "False Asym");
  g_FAsubper->Draw("AP"); DrawLine(g_FAsubper, 0.0);
  g_FAsubper->SetTitle("False_Asym");
  
  //Calculation alpha ratio and Draw them
  Double_t a_subper[nBins], eA_subper[nBins];
  Cal_AccRatioAndError(g_FAsubper, g_Pol, a_subper, eA_subper);

  //Jura/Saleve acceptance
  Double_t a_pol[nBins], eA_pol[nBins];
  Cal_AccRatioAndError(g_FApol, g_Pol, a_pol, eA_pol);
  Double_t a_upS[nBins], eA_upS[nBins];
  Cal_AccRatioAndError(g_FAupS, g_Pol, a_upS, eA_upS);
  Double_t a_downS[nBins], eA_downS[nBins];
  Cal_AccRatioAndError(g_FAdownS, g_Pol, a_downS, eA_downS);
  
  Double_t *xvals = g_FAsubper->GetX();
  Double_t ex[nBins] = {0.0};
  TGraphErrors *g_alpha_subper = new TGraphErrors(nBins, xvals, a_subper, ex,
						 eA_subper);
  TGraphErrors *g_alpha_pol = new TGraphErrors(nBins, xvals, a_pol, ex,
						 eA_pol);
  TGraphErrors *g_alpha_upS = new TGraphErrors(nBins, xvals, a_upS, ex,
						 eA_upS);
  TGraphErrors *g_alpha_downS = new TGraphErrors(nBins, xvals, a_downS, ex,
						 eA_downS);
  SetUp(g_alpha_subper); SetUp(g_alpha_pol);
  SetUp(g_alpha_upS); SetUp(g_alpha_downS); 

  TCanvas* c2 = new TCanvas("Acceptance Ratio", "Acceptance Ratio");
  c2->Divide(1,2);
  c2->cd(1);
  g_alpha_subper->Draw("AP"); g_alpha_subper->SetTitle("Acc_ratio");
  g_alpha_subper->SetMarkerColor(kBlue);
  if (process=="DY")
    g_alpha_subper->GetYaxis()->SetRangeUser(0.8, 1.2);
  else
    g_alpha_subper->GetYaxis()->SetRangeUser(0.93, 1.07);
  DrawLine(g_alpha_subper, 1.0);
  
  c2->cd(2);
  g_alpha_upS->Draw("AP"); g_alpha_upS->SetTitle("Acc_ratio_upS");
  g_alpha_upS->SetMarkerColor(kRed);
  g_alpha_downS->Draw("Psame"); g_alpha_downS->SetTitle("Acc_ratio_downS");
  g_alpha_downS->SetMarkerColor(kBlue);
  g_alpha_pol->Draw("Psame");
  g_alpha_upS->GetYaxis()->SetRangeUser(0.93, 1.07);
  DrawLine(g_alpha_upS, 1.0);

  //Calculate final systematic error
  Double_t sysErr[nBins], sysErr_statErr[nBins];
  Double_t *ey_subper = g_FAsubper->GetEY();
  Double_t nSigma =1.0;
  for (Int_t i=0; i<nBins; i++) {
    sysErr[i] = 1.0-a_subper[i];
    if (sysErr[i] < 0) sysErr[i] *= -1.0;
    sysErr[i] += nSigma*eA_subper[i];
    sysErr[i] /= 2.0;

    sysErr_statErr[i] = sysErr[i]/ey_subper[i];
  }
  TGraph *g_sys = new TGraph(nBins, xvals, sysErr);
  SetUp(g_sys); g_sys->SetTitle("Systematic_Error");
  TGraph *g_sys_stat = new TGraph(nBins, xvals, sysErr_statErr);
  SetUp(g_sys_stat); g_sys_stat->SetTitle("SysError_StatError");
  
  TCanvas* c3 = new TCanvas(); c3->Divide(2);
  c3->cd(1); g_sys->Draw("AP");
  c3->cd(2); g_sys_stat->Draw("AP");
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
acceptanceFourTargRatio";
  TString fOutput, fSystematics;
  if (whichFit=="true"){
    fOutput =Form("%s/acceptanceFourTargRatio_true_%s_%s%s_%s%s%i_%s_%s.root",
		  thisDirPath.Data(), period_Mtype.Data(), process.Data(),
		  lrMrange.Data(), binRange.Data(), physBinned.Data(), nBins, 
		  production.Data(), additionalCuts.Data());
    fSystematics =
      Form("%s/SystematicError/accSys4TargRatio_true_%s_%s%s_%s%s%i_%s_%s.root",
	   thisDirPath.Data(), period_Mtype.Data(), process.Data(),
	   lrMrange.Data(), binRange.Data(), physBinned.Data(), nBins,
	   production.Data(), additionalCuts.Data());
  }
  else{
    fOutput =
      Form("%s/acceptanceFourTargRatio_%s%s_%s_%s%s_%s%s%i_%ihbins_%s_%s.root",
	   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   binRange.Data(), physBinned.Data(), nBins, hbins, production.Data(),
	   additionalCuts.Data());
    fSystematics =
      Form("%s/SystematicError/\
accSys4TargRatio_%s%s_%s_%s%s_%s%s%i_%ihbins_%s_%s.root",
	   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   binRange.Data(), physBinned.Data(), nBins, hbins, production.Data(),
	   additionalCuts.Data());
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_FAsubper->Write("falseAN_subper");
    g_Pol->Write("Polarization");
    g_alpha_subper->Write("alpha_subper");
    g_alpha_pol->Write("alpha_pol");
    g_alpha_upS->Write("alpha_upS");
    g_alpha_downS->Write("alpha_downS");
    fResults->Close();

    TFile *fSys = new TFile(fSystematics, "RECREATE");
    g_sys->Write("gSys");
    g_sys_stat->Write("gSys_Stat");
    fSys->Close();
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
