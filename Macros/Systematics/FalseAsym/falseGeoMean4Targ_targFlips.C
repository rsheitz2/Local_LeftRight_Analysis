#include "include/helperFunctions.h"

Double_t FalseAmp_pol(Double_t NL[][4], Double_t NR[][4],
		      Double_t Pol, Int_t bi){
  //False asymmetry flipping polarized down L/R counts
  //NL[][4] => 4 = nTargPol
  
  Double_t Lup, Rup;
  Lup = NL[bi][upS_up]*NR[bi][upS_down];
  Rup = NR[bi][upS_up]*NL[bi][upS_down];

  Double_t Ldown, Rdown;
  Ldown = NL[bi][downS_up]*NR[bi][downS_down];
  Rdown = NR[bi][downS_up]*NL[bi][downS_down];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.5);
  R = TMath::Power(R, 0.5);

  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;

  return A;
}


Double_t FalseAmp_subper(Double_t NL[][4], Double_t NR[][4],
			 Double_t Pol, Int_t bi){
  //False asymmetry flipping sub period 1 L/R counts
  //NL[][4] => 4 = nTargPol
  
  Double_t Lup, Rup;
  Lup = NR[bi][upS_up]*NL[bi][upS_down];
  Rup = NL[bi][upS_up]*NR[bi][upS_down];

  Double_t Ldown, Rdown;
  Ldown = NL[bi][downS_up]*NR[bi][downS_down];
  Rdown = NR[bi][downS_up]*NL[bi][downS_down];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.5);
  R = TMath::Power(R, 0.5);

  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;

  return A;
}


Double_t e_Amp_pol(Double_t NL[][4], Double_t NR[][4],
		   Double_t e_NL[][4], Double_t e_NR[][4],
		   Double_t Pol, Int_t bi){
  //NL[][4] => 4 = nTargPol
  Int_t nTargPol=4;
  
  Double_t LinvSum2=0.0, RinvSum2=0.0;
  for (Int_t tr=0; tr<nTargPol; tr++) {
    LinvSum2 += e_NL[bi][tr]*e_NL[bi][tr]/(NL[bi][tr]*NL[bi][tr] );
    RinvSum2 += e_NR[bi][tr]*e_NR[bi][tr]/(NR[bi][tr]*NR[bi][tr] );
  }
  
  Double_t Lup, Rup;
  Lup = NL[bi][upS_up]*NR[bi][upS_down];
  Rup = NR[bi][upS_up]*NL[bi][upS_down];

  Double_t Ldown, Rdown;
  Ldown = NL[bi][downS_up]*NR[bi][downS_down];
  Rdown = NR[bi][downS_up]*NL[bi][downS_down];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.5);
  R = TMath::Power(R, 0.5);
    
  Double_t e = L*R/( (L + R)*(L + R) );
  Double_t error = e*TMath::Sqrt(LinvSum2 + RinvSum2)/Pol;

  return error;
}


Double_t e_Amp_subper(Double_t NL[][4], Double_t NR[][4],
		   Double_t e_NL[][4], Double_t e_NR[][4],
		   Double_t Pol, Int_t bi){
  //NL[][4] => 4 = nTargPol
  Int_t nTargPol=4;
  
  Double_t LinvSum2=0.0, RinvSum2=0.0;
  for (Int_t tr=0; tr<nTargPol; tr++) {
    LinvSum2 += e_NL[bi][tr]*e_NL[bi][tr]/(NL[bi][tr]*NL[bi][tr] );
    RinvSum2 += e_NR[bi][tr]*e_NR[bi][tr]/(NR[bi][tr]*NR[bi][tr] );
  }

  Double_t Lup, Rup;
  Lup = NR[bi][upS_up]*NL[bi][upS_down];
  Rup = NL[bi][upS_up]*NR[bi][upS_down];

  Double_t Ldown, Rdown;
  Ldown = NL[bi][downS_up]*NR[bi][downS_down];
  Rdown = NR[bi][downS_up]*NL[bi][downS_down];
  
  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.5);
  R = TMath::Power(R, 0.5);

  Double_t e = L*R/( (L + R)*(L + R) );
  Double_t error = e*TMath::Sqrt(LinvSum2 + RinvSum2)/Pol;

  return error;
}


void CalAmp_AmpErr(Double_t *Fasym, Double_t *e_Fasym,
		   Double_t NL[][4], Double_t NR[][4],
		   Double_t e_NL[][4], Double_t e_NR[][4],
		   Double_t *Pol, Int_t nBins, TString whichFasym){

  for (Int_t bi=0; bi<nBins; bi++) {
    if (whichFasym=="pol"){
      Fasym[bi] = FalseAmp_pol(NL, NR, Pol[bi], bi);
      e_Fasym[bi] = e_Amp_pol(NL, NR, e_NL, e_NR, Pol[bi], bi);
    }
    else if (whichFasym=="subper"){
      Fasym[bi] = FalseAmp_subper(NL, NR, Pol[bi], bi);
      e_Fasym[bi] = e_Amp_subper(NL, NR, e_NL, e_NR, Pol[bi], bi);
    }
  }
}



void falseGeoMean4Targ_targFlips(TString start =""){
  //Setup_______________
  const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xF";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit ="true";

  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/";
  
  if (start==""){//Basic info
    cout << "\nScript calculates false AN asymmetries using the 4 target";
    cout << " geometric mean method" << endl;
    cout << "False asymmetries are determined by flipping L/R counts of ";
    cout << "polarization down targ" << endl;
    cout << "    and by flipping L/R counts of sub period 1 ";
    cout << "(defined as upS_up && downS_down)" << endl;
    cout << "\nInput needed is a TGraphErrors of AN per target/pol";
    cout<<" with polarization corrections and without polarization corrections";
    cout << "\nand a TGraphErrors of Left and Right counts per target/pol";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C or trueCount.C ";
    cout <<"-> falseGeoMean4Targ_targFlips.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'falseGeoMean4Targ_polRevs(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << pathAN << endl;
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
  TString inputFiles[2];//{corrPath, noCorrPath}
  if (whichFit=="true"){
    inputFiles[0] = Form("trueCount_%s_%s%s_%s%i_corr.root",
			 period_Mtype.Data(), process.Data(), lrMrange.Data(),
			 physBinned.Data(), nBins);
    inputFiles[1] = Form("trueCount_%s_%s%s_%s%i_noCorr.root",
			 period_Mtype.Data(), process.Data(), lrMrange.Data(),
			 physBinned.Data(), nBins);
    pathAN += "trueCount/";
  }
  else{
    inputFiles[0] = Form("functMFit_%s%s_%s_%s%s_%s%i_%ihbin_corr.root",
			 whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
			 process.Data(), lrMrange.Data(), physBinned.Data(),
			 nBins, hbins);
    inputFiles[1] = Form("functMFit_%s%s_%s_%s%s_%s%i_%ihbin_noCorr.root",
			 whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
			 process.Data(), lrMrange.Data(), physBinned.Data(),
			 nBins, hbins);
    pathAN += "functMFit/";
  }
  
  TFile *fAN = TFile::Open(pathAN+inputFiles[0]);
  TFile *fAN_noCorr = TFile::Open(pathAN+inputFiles[1]);

  if (!fAN || !fAN_noCorr ){
    cout << "RD or RD_noCorr file does not exist " << endl;
    exit(EXIT_FAILURE);
  }
  
  //File data names setup
  TString inputAN[4];//{upS_up, upS_down, downS_up, downS_down}
  TString inputLeft[4];//{upS_up, upS_down, downS_up, downS_down}
  TString inputRight[4];//{upS_up, upS_down, downS_up, downS_down}
  
  TString targName[4] = {"upS_up", "upS_down", "downS_up", "downS_down"};  
  for (Int_t i=0; i<4; i++){
    inputAN[i] = Form("AN_%s", targName[i].Data() );
    inputLeft[i] = Form("Left_%s", targName[i].Data() );
    inputRight[i] = Form("Right_%s", targName[i].Data() );
  }

  //Get Data from files
  const Int_t nTargPol =4;
  TGraphErrors *g_corr[nTargPol], *g_noCorr[nTargPol];
  TGraphErrors *g_Left[nTargPol], *g_Right[nTargPol];
  
  for (Int_t tr=0; tr<nTargPol; tr++) {
    g_corr[tr] = (TGraphErrors*)fAN->Get(inputAN[tr]);
    g_noCorr[tr] = (TGraphErrors*)fAN_noCorr->Get(inputAN[tr]);
    
    g_Left[tr] = (TGraphErrors*)fAN->Get(inputLeft[tr]);
    g_Right[tr] = (TGraphErrors*)fAN->Get(inputRight[tr]);
  }
  
  //Get L/R counts 
  Double_t LeftCounts[nBins][nTargPol], RightCounts[nBins][nTargPol];
  Double_t e_LeftCounts[nBins][nTargPol], e_RightCounts[nBins][nTargPol];
  
  for (Int_t tr=0; tr<nTargPol; tr++) {
    Double_t *y_Left = g_Left[tr]->GetY();
    Double_t *y_Right = g_Right[tr]->GetY();
    Double_t *ey_Left = g_Left[tr]->GetEY();
    Double_t *ey_Right = g_Right[tr]->GetEY();
    for (Int_t bi=0; bi<nBins; bi++) {
      LeftCounts[bi][tr] = y_Left[bi];
      RightCounts[bi][tr] = y_Right[bi];
      e_LeftCounts[bi][tr] = ey_Left[bi];
      e_RightCounts[bi][tr] = ey_Right[bi];
    }
  }

  //Get polarization values
  Double_t Pol[nBins];
  Double_t *y_corr = g_corr[0]->GetY();
  Double_t *y_noCorr = g_noCorr[0]->GetY();
  GetPolarization(y_noCorr, y_corr, Pol, nBins);

  //Make 4Targ False Asymmetries
  Double_t fAN_pol[nBins], e_fAN_pol[nBins];
  Double_t fAN_subper[nBins], e_fAN_subper[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_corr[0]->GetX();

  CalAmp_AmpErr(fAN_pol, e_fAN_pol, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, Pol, nBins,
		"pol");
  CalAmp_AmpErr(fAN_subper, e_fAN_subper, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, Pol, nBins,
		"subper");
    
  //Draw false AN
  TGraphErrors *g_fAN_pol =
    new TGraphErrors(nBins, xvals, fAN_pol, ex, e_fAN_pol);
  SetUp(g_fAN_pol);
  TGraphErrors *g_fAN_subper
    = new TGraphErrors(nBins, xvals, fAN_subper, ex, e_fAN_subper);
  SetUp(g_fAN_subper);

  TCanvas* c1 = new TCanvas();
  Double_t yMax = 0.4;
  g_fAN_pol->Draw("AP"); g_fAN_pol->SetMarkerColor(kBlue);
  g_fAN_subper->Draw("Psame"); g_fAN_subper->SetMarkerColor(kRed);
  g_fAN_pol->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_fAN_pol, 0.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/Data/TargFlip/falseGeoMean4Targ_%s_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), whichFit.Data(), period_Mtype.Data(),
		  process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins);
  }
  else{
    fOutput =Form("%s/Data/TargFlip/falseGeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbins.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_fAN_pol->Write("falseAN_pol");
    g_fAN_subper->Write("falseAN_subper");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathAN << endl;
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
