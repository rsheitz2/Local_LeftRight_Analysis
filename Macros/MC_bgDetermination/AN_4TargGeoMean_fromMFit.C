#include "helperFunctions.h"

//Setup_______________
const Int_t nBins =5;

TString physAmp = "psi";
TString physBinned = "pT";//"xN" "xPi" "xF" "pT"
Bool_t toWrite =true;
//Setup_______________  


Double_t Amp(Double_t NL[][nBins], Double_t NR[][nBins],
	     Double_t Pol, Int_t bi){

  Double_t Lup, Rup;
  Lup = NL[0][bi]*NL[1][bi]; Rup = NR[0][bi]*NR[1][bi];

  Double_t Ldown, Rdown;
  Ldown = NL[2][bi]*NL[3][bi]; Rdown = NR[2][bi]*NR[3][bi];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;

  return A;
}


Double_t e_Amp(Double_t NL[][nBins], Double_t NR[][nBins],
	       Double_t e_NL[][nBins], Double_t e_NR[][nBins],
	       Double_t Pol, Int_t bi){
  Int_t nTarg=4;

  Double_t L=1.0, R=1.0;
  Double_t LinvSum2=0.0, RinvSum2=0.0;
  for (Int_t tr=0; tr<nTarg; tr++) {
    L *= NL[tr][bi]; R *= NR[tr][bi];
    LinvSum2 += e_NL[tr][bi]*e_NL[tr][bi]/(NL[tr][bi]*NL[tr][bi] );
    RinvSum2 += e_NR[tr][bi]*e_NR[tr][bi]/(NR[tr][bi]*NR[tr][bi] );
  }
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  Double_t e = L*R/( 2*(L + R)*(L + R) );
  Double_t error = e*( TMath::Sqrt(LinvSum2) + TMath::Sqrt(RinvSum2) )/Pol;
  
  return error;
}


void AN_4TargGeoMean_fromMFit(){
  TString pathFitAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/MC_bgDetermination/Data";
  TFile *fAN = TFile::Open(Form("%s/FitMass_%s_2.50_8.50_corr_%s.root",
				pathFitAN.Data(), physBinned.Data(),
				physAmp.Data() )  );
  TFile *fAN_noCorr = TFile::Open(Form("%s/FitMass_%s_2.50_8.50_noCorr_%s.root",
				       pathFitAN.Data(), physBinned.Data(),
				       physAmp.Data() ) ); 
  if (!fAN || !fAN_noCorr ){
    cout << "RD or RD_noCorr file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nTargPol =4;
  TString targName[nTargPol] = {"upS_up", "upS_down", "downS_up", "downS_down"};
  TGraphErrors *g_corr[nTargPol], *g_noCorr[nTargPol];
  TGraphErrors *g_Left[nTargPol], *g_Right[nTargPol];
  
  for (Int_t tr=0; tr<nTargPol; tr++) {
    g_corr[tr] = (TGraphErrors*)fAN->Get(Form("%s_%s", physAmp.Data(),
					      targName[tr].Data()) );
    g_noCorr[tr] = (TGraphErrors*)fAN_noCorr->Get(Form("%s_%s", physAmp.Data(),
						       targName[tr].Data()) );

    g_Left[tr] = (TGraphErrors*)fAN->Get(Form("%s_Left_%s", physAmp.Data(),
					      targName[tr].Data()) );
    g_Right[tr] = (TGraphErrors*)fAN->Get(Form("%s_Right_%s", physAmp.Data(),
					       targName[tr].Data()) );
  }

  if (g_corr[0]->GetN() != nBins){
    cout << "nBins not defined well!!!" << endl;
    exit(EXIT_FAILURE);
  }

  Double_t LeftCounts[nTargPol][nBins], RightCounts[nTargPol][nBins];
  Double_t e_LeftCounts[nTargPol][nBins], e_RightCounts[nTargPol][nBins];
  Double_t Pol[nBins];
  for (Int_t tr=0; tr<nTargPol; tr++) {
    Double_t *y_corr = g_corr[tr]->GetY();
    Double_t *y_noCorr = g_noCorr[tr]->GetY();
    
    Double_t *y_Left = g_Left[tr]->GetY();
    Double_t *y_Right = g_Right[tr]->GetY();
    Double_t *ey_Left = g_Left[tr]->GetEY();
    Double_t *ey_Right = g_Right[tr]->GetEY();
    for (Int_t bi=0; bi<nBins; bi++) {
      if (tr==0) Pol[bi] = y_noCorr[bi]/y_corr[bi];

      LeftCounts[tr][bi] = y_Left[bi];
      RightCounts[tr][bi] = y_Right[bi];
      e_LeftCounts[tr][bi] = ey_Left[bi];
      e_RightCounts[tr][bi] = ey_Right[bi];
    }
  }


  Double_t AN_4targ[nBins], e_AN_4targ[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_corr[0]->GetX();
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_4targ[bi] = Amp(LeftCounts, RightCounts, Pol[bi], bi);
    e_AN_4targ[bi] = e_Amp(LeftCounts, RightCounts,
			   e_LeftCounts, e_RightCounts, Pol[bi], bi);
  }

  TGraphErrors *g_AN = new TGraphErrors(nBins, xvals, AN_4targ, ex, e_AN_4targ);
  SetUpTGraph(g_AN);

  TCanvas* c1 = new TCanvas();
  g_AN->Draw("AP");
  DrawLine(g_AN, 0.0);
  

  TString fOutput = Form("AN_%s_%s.root", physBinned.Data(), physAmp.Data() );
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_AN->Write("AN");
  }

  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "AN for physics type:  " << physAmp << endl;
  cout << "nBins considered:  " << nBins << endl;
  cout << "Binned in:  " << physBinned << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
