#include "include/helperFunctions.h"

Double_t Amp(Double_t NL[][4], Double_t NR[][4],
	     Double_t Pol, Int_t bi){
  //NL[][4] => 4 = nTargPol
  
  Double_t Lup, Rup;
  Lup = NL[bi][upS_up]*NL[bi][upS_down];
  Rup = NR[bi][upS_up]*NR[bi][upS_down];

  Double_t Ldown, Rdown;
  Ldown = NL[bi][downS_up]*NL[bi][downS_down];
  Rdown = NR[bi][downS_up]*NR[bi][downS_down];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);
  
  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;
  
  return A;
}


Double_t e_Amp(Double_t NL[][4], Double_t NR[][4],
	       Double_t e_NL[][4], Double_t e_NR[][4],
	       Double_t Pol, Int_t bi){
  //NL[][4] => 4 = nTargPol
  Int_t nTargPol=4;
  
  Double_t L=1.0, R=1.0;
  Double_t LinvSum2=0.0, RinvSum2=0.0;
  for (Int_t tr=0; tr<nTargPol; tr++) {
    L *= NL[bi][tr]; R *= NR[bi][tr];
    LinvSum2 += e_NL[bi][tr]*e_NL[bi][tr]/(NL[bi][tr]*NL[bi][tr] );
    RinvSum2 += e_NR[bi][tr]*e_NR[bi][tr]/(NR[bi][tr]*NR[bi][tr] );
  }
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  Double_t e = L*R/( 2*(L + R)*(L + R) );
  Double_t error = e*TMath::Sqrt(LinvSum2 + RinvSum2)/Pol;

  return error;
}


void GeoMean4Targ(TString start =""){
  //Setup_______________
  Double_t phiScut =1.44;
  TString additionalCuts ="phiS0.0";
  Double_t A_siv =0.1;
  Int_t N_gen =1000;
  Bool_t alphaScale =true;

  TString physBinned ="xN";//"xF", "pT"
  const Int_t nBins =1;//# of physBinned bins
      
  Bool_t toWrite =true;
  //Setup_______________
  
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Macros/Systematics/ToyMC/Data/";
  if (start==""){//Basic info
    cout << "\nScript calculates AN using the 4 target geometric mean method";
    cout << "\nInput needed is a TGraphErrors of AN per target/pol";
    cout<<" with polarization corrections and without polarization corrections";
    cout << "\nand a TGraphErrors of Left and Right counts per target/pol";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C or trueCount.C ";
    cout <<"-> GeoMean4Targ.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'AN_4TargGeoMean_fromMFit(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << pathRD << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //File name setups && get files
  TString RDfile = Form("Gen/generator_%0.3f_%s_gen%i_alpha%i.root",
			A_siv, additionalCuts.Data(), N_gen, alphaScale);
  TFile *fRD  = OpenFile(pathRD + RDfile);
  
  //Get data from files
  TString inputLeft[4];//{upS_up, upS_down, downS_up, downS_down}
  TString inputRight[4];//{upS_up, upS_down, downS_up, downS_down}
  
  TString targName[4] = {"upS_up", "upS_down", "downS_up", "downS_down"};  
  for (Int_t i=0; i<4; i++){
    inputLeft[i] = Form("Left_%s", targName[i].Data() );
    inputRight[i] = Form("Right_%s", targName[i].Data() );
  }

  const Int_t nTargPol =4;
  TGraphErrors *g_Left[nTargPol], *g_Right[nTargPol];
  for (Int_t tr=0; tr<nTargPol; tr++) {
    g_Left[tr] = (TGraphErrors*)fRD->Get(inputLeft[tr]);
    g_Right[tr] = (TGraphErrors*)fRD->Get(inputRight[tr]);
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

  //Get Polarization
  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  Double_t *Pol = g_Pol->GetY();
  
  //Make 4Targ Asym
  Double_t AN_4targ[nBins], e_AN_4targ[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_Pol->GetX();
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_4targ[bi] = Amp(LeftCounts, RightCounts, Pol[bi], bi);
    e_AN_4targ[bi] = e_Amp(LeftCounts, RightCounts,
			   e_LeftCounts, e_RightCounts, Pol[bi], bi);
    if (e_AN_4targ[bi] < 10e-9) {
      cout << "Error: AN error way too small" << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Draw AN
  TGraphErrors *g_AN = new TGraphErrors(nBins, xvals, AN_4targ, ex, e_AN_4targ);
  SetUp(g_AN);

  Double_t yMax =0.15;
  TCanvas* c1 = new TCanvas();
  g_AN->Draw("AP"); g_AN->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN, 0.0);
  
  //Write output/final settings
  TString thisDirPath = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/ToyMC/Data";
  TString fOutput = Form("%s/GeoMean/GeoMean_%0.3f_%s_gen%i_alpha%i.root",
			 thisDirPath.Data(), A_siv, additionalCuts.Data(),
			 N_gen, alphaScale);
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_AN->Write("AN");
  }

  cout << "\nSettings________" << endl;
  cout << "Data coming from:            " << pathRD << endl;
  cout << "Input P corrected data:        " << RDfile << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
  cout << " " << endl;
}
