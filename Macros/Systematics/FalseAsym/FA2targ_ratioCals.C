#include "include/helperFunctions.h"

//2Targ Geomean
Double_t AsymCal(Double_t L, Double_t R){
  L = TMath::Sqrt(L);
  R = TMath::Sqrt(R);
  
  return (L - R)/(L + R);
}


Double_t AsymErrorCal(Double_t L, Double_t R,
		      Double_t LinvSum2, Double_t RinvSum2){
  Double_t epsilon = AsymCal(L, R);
  
  LinvSum2 = TMath::Sqrt(LinvSum2);
  RinvSum2 = TMath::Sqrt(RinvSum2);
  
  L = TMath::Sqrt(L);
  R = TMath::Sqrt(R);
  Double_t dL = 0.5*L*LinvSum2;
  Double_t dR = 0.5*R*RinvSum2;
  
  Double_t error =
    (1 - epsilon)*(1-epsilon)*dL*dL + (1 + epsilon)*(1 + epsilon)*dR*dR;
  error = TMath::Sqrt(error);
  error *= 1.0/(L + R);
  
  return error;
}


Double_t FalseAmp_upSup(Double_t NL[][4], Double_t NR[][4], Int_t bi){
  Double_t L, R;
  L = NL[bi][upS_up]*NR[bi][downS_down]; //FA 
  R = NR[bi][upS_up]*NL[bi][downS_down];
  //sqrt(a^Pup_uL*a^Pdown_dL / a^Pup_uR*a^Pdown_dR )
  
  Double_t A = AsymCal(L, R);

  return A;
}


Double_t e_Amp_upSup(Double_t NL[][4], Double_t NR[][4],
		     Double_t e_NL[][4], Double_t e_NR[][4], Int_t bi){
  Double_t LinvSum2=0.0;
  LinvSum2 += e_NL[bi][upS_up]*e_NL[bi][upS_up]/(NL[bi][upS_up]*NL[bi][upS_up]);
  LinvSum2 +=
    e_NR[bi][downS_down]*e_NR[bi][upS_down]/
    (NR[bi][downS_down]*NR[bi][downS_down]);
  Double_t RinvSum2=0.0;
  RinvSum2 += e_NL[bi][upS_up]*e_NL[bi][upS_up]/(NL[bi][upS_up]*NL[bi][upS_up]);
  RinvSum2 +=
    e_NR[bi][downS_down]*e_NR[bi][downS_down]/
    (NR[bi][downS_down]*NR[bi][downS_down]);

  Double_t L = NL[bi][upS_up]*NR[bi][downS_down];
  Double_t R = NR[bi][upS_up]*NL[bi][downS_down];
  
  Double_t error = AsymErrorCal(L, R, LinvSum2, RinvSum2);

  return error;
}


Double_t FalseAmp_upSdown(Double_t NL[][4], Double_t NR[][4], Int_t bi){
  Double_t L, R;
  L = NR[bi][upS_down]*NL[bi][downS_up]; //FA
  R = NL[bi][upS_up]*NR[bi][downS_up];
  //sqrt(a^Pdown_uL*a^Pdup_dL / a^Pdown_uR*a^Pup_dR )
  
  Double_t A = AsymCal(L, R);

  return A;
}

Double_t e_Amp_upSdown(Double_t NL[][4], Double_t NR[][4],
		       Double_t e_NL[][4], Double_t e_NR[][4], Int_t bi){
  Double_t LinvSum2=0.0;
  LinvSum2 +=
    e_NR[bi][upS_down]*e_NR[bi][upS_down]/(NR[bi][upS_down]*NR[bi][upS_down]);
  LinvSum2 +=
    e_NL[bi][downS_up]*e_NL[bi][downS_up]/
    (NL[bi][downS_up]*NL[bi][downS_up]);
  Double_t RinvSum2=0.0;
  RinvSum2 +=
    e_NL[bi][upS_down]*e_NL[bi][upS_down]/(NL[bi][upS_down]*NL[bi][upS_down]);
  RinvSum2 +=
    e_NR[bi][downS_up]*e_NR[bi][downS_up]/
    (NR[bi][downS_up]*NR[bi][downS_up]);
  
  Double_t L = NR[bi][upS_down]*NL[bi][downS_up];
  Double_t R = NL[bi][upS_up]*NR[bi][downS_up];

  Double_t error = AsymErrorCal(L, R, LinvSum2, RinvSum2);

  return error;
}


void CalAmp_AmpErr(Double_t *Fasym, Double_t *e_Fasym,
		   Double_t NL[][4], Double_t NR[][4],
		   Double_t e_NL[][4], Double_t e_NR[][4],
		   Int_t nBins, TString whichFasym){
  for (Int_t bi=0; bi<nBins; bi++) {
    if (whichFasym=="Pup"){
      Fasym[bi] = FalseAmp_upSup(NL, NR, bi);
      e_Fasym[bi] = e_Amp_upSup(NL, NR, e_NL, e_NR, bi);
    }
    else if (whichFasym=="Pdown"){
      Fasym[bi] = FalseAmp_upSdown(NL, NR, bi);
      e_Fasym[bi] = e_Amp_upSdown(NL, NR, e_NL, e_NR, bi);
    }
  }
}


void FA2targ_ratioCals(TString start =""){
  //Setup_______________
  /*const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";//*/
  
  const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.90_3.30";
  TString fitMrange ="2.00_7.50";
  TString binRange ="25_43";
  TString whichFit ="ten";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathLR = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/";
  
  if (start==""){//Basic info
    cout << "\nScript calculates false AN asymmetries using the a two target";
    cout << " geometric mean method" << endl;
    cout << "False asymmetries are determined by looking a Jura-Saleve ";
    cout << "for upstream polarized up or upstream polarized down" << endl;
    cout << "\nInput needed = TGraphErrors of Left and Right counts ";
    cout << "per target/polarization";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "(true fit) systematic_leftRight -> FA2targ_ratioCals.C " << endl;
    cout << "(other fit) systematic_leftRight -> fullTargSysFunctMFit.C -> ";
    cout << "FA2targ_ratioCals.C " << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'FA2targ_ratioCals(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << pathLR << endl;
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
    inputFiles = Form("systematic_leftRight_%s%s_%ibins%s_%ihbin.root",
		      period_Mtype.Data(), fitMrange.Data(), nBins,
		      binRange.Data(), hbins);
    pathLR += "Data/";
  }
  else{
    inputFiles = Form("fullTargSysFunctMFit_%s%s_%s_%s%s_%s%i_%ihbin.root",
		      whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
		      process.Data(), lrMrange.Data(), physBinned.Data(),
		      nBins, hbins);
    pathLR += "Macros/Systematics/FalseAsym/Data/fullTargSysFunctMFit/";
  }

  TFile *fAN = TFile::Open(pathLR+inputFiles);
  if (!fAN ){
    cout << "RD file does not exist " << endl;
    exit(EXIT_FAILURE);
  }
  
  //File data names setup
  const Int_t nTargPol =4;
  TString inputLeft[nTargPol];//{upS_up, upS_down, downS_up, downS_down}
  TString inputRight[nTargPol];//{upS_up, upS_down, downS_up, downS_down}
  
  TString targName[nTargPol] = {"upstream_up", "upstream_down",
				"downstream_up", "downstream_down"};  
  for (Int_t i=0; i<nTargPol; i++){
    inputLeft[i] = Form("%s_left_%s", physBinned.Data(), targName[i].Data() );
    inputRight[i] = Form("%s_right_%s", physBinned.Data(), targName[i].Data() );
  }
  
  //Get Data from files
  TGraphErrors *g_Left[nTargPol], *g_Right[nTargPol];
  for (Int_t tr=0; tr<nTargPol; tr++) {
    g_Left[tr] = (TGraphErrors*)fAN->Get(inputLeft[tr]);
    g_Right[tr] = (TGraphErrors*)fAN->Get(inputRight[tr]);
  }
  
  //Get L/R counts 
  Double_t LeftCounts[nBins][nTargPol], RightCounts[nBins][nTargPol];
  Double_t e_LeftCounts[nBins][nTargPol], e_RightCounts[nBins][nTargPol];

  for (Int_t tr=0; tr<nTargPol; tr++) {
    Double_t *y_Left = g_Left[tr]->GetY();
    Double_t *y_Right = g_Right[tr]->GetY();
    Double_t *ey_Left = (whichFit=="true") ? NULL : g_Left[tr]->GetEY();
    Double_t *ey_Right = (whichFit=="true") ? NULL : g_Right[tr]->GetEY();
      
    for (Int_t bi=0; bi<nBins; bi++) {
      LeftCounts[bi][tr] = y_Left[bi];
      RightCounts[bi][tr] = y_Right[bi];
      if (whichFit == "true"){
	e_LeftCounts[bi][tr] = TMath::Sqrt(y_Left[bi] );
	e_RightCounts[bi][tr] = TMath::Sqrt(y_Right[bi] );
      }
      else{
	e_LeftCounts[bi][tr] = ey_Left[bi];
	e_RightCounts[bi][tr] = ey_Right[bi];
      }
    }
  }
  
  //Make 4Targ False Asymmetries
  Double_t fAN_upSup[nBins], e_fAN_upSup[nBins];
  Double_t fAN_upSdown[nBins], e_fAN_upSdown[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_Left[0]->GetX();
  
  CalAmp_AmpErr(fAN_upSup, e_fAN_upSup, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, nBins, "Pup");
  CalAmp_AmpErr(fAN_upSdown, e_fAN_upSdown, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, nBins, "Pdown");
  
  //Draw false AN
  TGraphErrors *g_fAN_upSup =
    new TGraphErrors(nBins, xvals, fAN_upSup, ex, e_fAN_upSup);
  SetUp(g_fAN_upSup);
  TGraphErrors *g_fAN_upSdown
    = new TGraphErrors(nBins, xvals, fAN_upSdown, ex, e_fAN_upSdown);
  SetUp(g_fAN_upSdown);

  TCanvas* c1 = new TCanvas();
  Double_t yMax = 0.08;
  g_fAN_upSup->Draw("AP"); g_fAN_upSup->SetMarkerColor(kBlue);
  g_fAN_upSdown->Draw("Psame"); g_fAN_upSdown->SetMarkerColor(kRed);
  g_fAN_upSup->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_fAN_upSup, 0.0);

  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
FA2targ_ratioCals";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/FA2targ_ratioCals_true_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), period_Mtype.Data(), process.Data(),
		  lrMrange.Data(), physBinned.Data(), nBins);
  }
  else{
    fOutput =Form("%s/FA2targ_ratioCals_%s%s_%s_%s%s_%s%i_%ihbins.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_fAN_upSup->Write("falseAN_upSup");
    g_fAN_upSdown->Write("falseAN_upSdown");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathLR << endl;
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

