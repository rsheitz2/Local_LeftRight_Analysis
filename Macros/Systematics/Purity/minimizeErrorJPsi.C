#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

void minimizeErrorJPsi(){

  //Basic Setup
  TString pathAn="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ/";
  TString pathSys="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/Purity/Data/Purity/";

  const Int_t nMassRanges =4;
  TString massRanges[] = {"2.52_3.72", "2.72_3.52", "2.92_3.32",
			  "3.02_3.22"};

  //Error determination 
  Double_t delta[nMassRanges];
  Double_t xSigmas[] = {3.0, 2.0, 1.0, 0.5};
  Double_t deltaAn[nMassRanges], deltaPurity[nMassRanges];
  
  for (Int_t m=0; m<nMassRanges; m++) {
    //Get AN errors
    TString anName =
      Form("GeoMean4Targ_true_WAll_LowM_AMDY_JPsi%s_25_43xN1_slot1_phiS0.0.root",
	   massRanges[m].Data() );

    TFile *fAN = OpenFile(pathAn + anName);
    TGraphErrors *gAN =(TGraphErrors*)fAN->Get("AN");

    //Get purity
    TString purityName =
      Form("mcMFit_MC2.00_8.50_WAll_LowM_AMDY_JPsi%s_25_43pT5_150hbin.root",
	   massRanges[m].Data() );

    TFile *fPurity = OpenFile(pathSys + purityName);
    TH1D *hPurity =(TH1D*)fPurity->Get("h_purity");
    
    Double_t p = hPurity->GetMean();

    //Plot errors
    Double_t *eAN = gAN->GetEY();

    delta[m] = eAN[0]*eAN[0];
    delta[m] *= (p*p + (1-p)*(1-p))/(p*p);
    delta[m] = TMath::Sqrt(delta[m]);

    deltaAn[m] = eAN[0]*eAN[0];
    deltaPurity[m] = (p*p + (1-p)*(1-p))/(p*p);

    cout << "An factor: " << deltaAn[m]
	 << "   purity factor: " << deltaPurity[m]
	 << "   purity: " << p
	 << " +/- " << hPurity->GetRMS() << endl;
  }

  //Graphs
  TGraph *gMinimize = new TGraph(nMassRanges, xSigmas, delta);
  SetUp(gMinimize);

  TGraph *g_dAn = new TGraph(nMassRanges, xSigmas, deltaAn);
  TGraph *g_dPurity = new TGraph(nMassRanges, xSigmas, deltaPurity);
  SetUp(g_dAn); SetUp(g_dPurity);

  TCanvas* c1 = new TCanvas(); c1->Divide(2);
  c1->cd(1); gMinimize->Draw("AP");

  c1->cd(2);
  g_dPurity->Draw("AP");
  g_dAn->Draw("Psame");
  g_dPurity->SetMarkerColor(kBlue);
  g_dPurity->SetMarkerStyle(25);
}
