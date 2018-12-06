#include "include/helperFunctions.h"

Double_t AsymCal(Double_t L, Double_t R){
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);
  
  return (L - R)/(L + R);
}


Double_t AsymErrorCal(Double_t L, Double_t R,
		      Double_t LinvSum2, Double_t RinvSum2){
  Double_t epsilon = AsymCal(L, R);
  
  LinvSum2 = TMath::Sqrt(LinvSum2);
  RinvSum2 = TMath::Sqrt(RinvSum2);
  
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);
  Double_t dL = 0.25*L*LinvSum2;
  Double_t dR = 0.25*R*RinvSum2;
  
  Double_t error =
    (1 - epsilon)*(1-epsilon)*dL*dL + (1 + epsilon)*(1 + epsilon)*dR*dR;
  error = TMath::Sqrt(error);
  error *= 1.0/(L + R);
  
  return error;
}


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
  Double_t A = AsymCal(L, R);

  return A/Pol;
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
  Double_t A = AsymCal(L, R);

  return A/Pol;
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
  Double_t error = AsymErrorCal(L, R, LinvSum2, RinvSum2);

  return error/Pol;
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
  Double_t error = AsymErrorCal(L, R, LinvSum2, RinvSum2);
    
  return error/Pol;
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
  const Int_t nBins =3;//HMDY
  TString period_Mtype ="W12_HMDY";
  Int_t hbins =150;
  TString physBinned ="xPi";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.53";
  
  //const Int_t nBins =5;//JPsi
  //TString period_Mtype ="W12_LowM_AMDY";
  //Int_t hbins =150;
  //TString physBinned ="xPi";//xN, xPi, xF, pT, M
  //TString process ="JPsi";//JPsi, psi, DY
  //TString lrMrange ="2.00_5.00";
  //TString fitMrange ="2.00_7.50";
  //TString binRange ="25_43";
  //TString whichFit ="thirteen";
  //TString production ="slot1";//"t3", "slot1"
  //TString additionalCuts ="phiS0.195";
  
  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/";
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  
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
  TString inputFiles, RDfile;
  if (whichFit=="true"){
    inputFiles =
      inputFiles =
      Form("trueCount_%s_%s%s_%s%s%i_%s_%s_corr.root", period_Mtype.Data(),
	   process.Data(), fitMrange.Data(), binRange.Data(), physBinned.Data(),
	   nBins, production.Data(), additionalCuts.Data());
    pathAN += "trueCount/";

    RDfile =Form("leftRight_byTarget_%s%s_%ibins%s_%ihbin_%s_%s.root",
		 period_Mtype.Data(), lrMrange.Data(), nBins, binRange.Data(),
		 hbins, production.Data(), additionalCuts.Data());
  }
  else{
    inputFiles = Form("functMFit_%s%s_%s_%s%s_%s%s%i_%ihbin_%s_%s_corr.root",
		      whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
		      process.Data(), lrMrange.Data(), binRange.Data(),
		      physBinned.Data(), nBins, hbins, production.Data(),
		      additionalCuts.Data());
    pathAN += "functMFit/";

    RDfile =
      Form("leftRight_byTarget_%s1.00_8.50_%ibins%s_%ihbin_%s_%s.root",
	   period_Mtype.Data(), nBins, binRange.Data(), hbins,
	   production.Data(), additionalCuts.Data() );
  }
  
  TFile *fAN = TFile::Open(pathAN+inputFiles);
  TFile *fRD  = TFile::Open(pathRD + RDfile);

  if (!fAN || !fRD ){
    cout << "fAN or fRD file does not exist " << endl;
    cout << inputFiles << "\n";
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
  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  TGraph* g_Dil =(TGraph*)fRD->Get(Form("%s_Dil", physBinned.Data()));
  
  const Int_t nTargPol =4;
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
    Double_t *ey_Left = g_Left[tr]->GetEY();
    Double_t *ey_Right = g_Right[tr]->GetEY();
    for (Int_t bi=0; bi<nBins; bi++) {
      LeftCounts[bi][tr] = y_Left[bi];
      RightCounts[bi][tr] = y_Right[bi];
      e_LeftCounts[bi][tr] = ey_Left[bi];
      e_RightCounts[bi][tr] = ey_Right[bi];
    }
  }
  
  //Make 4Targ False Asymmetries
  Double_t fAN_pol[nBins], e_fAN_pol[nBins];
  Double_t fAN_subper[nBins], e_fAN_subper[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_Pol->GetX();

  //Get polarization values
  Double_t Pol[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);
  TGraph *g_PolDil = new TGraph(nBins, xvals, Pol); SetUp(g_PolDil);

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
  Double_t yMax = (nBins==3) ? 0.75 : 0.25;
  g_fAN_pol->Draw("AP"); g_fAN_pol->SetMarkerColor(kRed);
  g_fAN_subper->Draw("Psame"); g_fAN_subper->SetMarkerColor(kBlue);
  g_fAN_pol->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_fAN_pol, 0.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/TargFlip";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/falseGeoMean4Targ_%s_%s_%s%s_%s%s%i_%s_%s.root",
		  thisDirPath.Data(), whichFit.Data(), period_Mtype.Data(),
		  process.Data(), lrMrange.Data(), binRange.Data(),
		  physBinned.Data(), nBins, production.Data(),
		  additionalCuts.Data());
  }
  else{
    fOutput =Form("%s/falseGeoMean4Targ_%s%s_%s_%s%s_%s%s%i_%ihbins_%s_%s.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  binRange.Data(), physBinned.Data(), nBins, hbins,
		  production.Data(), additionalCuts.Data());
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_fAN_pol->Write("falseAN_pol");
    g_fAN_subper->Write("falseAN_subper");
    g_PolDil->Write("Polarization");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathAN << endl;
  cout << "Input P corrected data:        " << inputFiles << endl;
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
