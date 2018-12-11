#include "include/helperFunctions.h"

void drawOneGeoMean(){
  //Setup_______________
  const Int_t nPhiScut =12;
  TString phiScut[nPhiScut] ={"0.0", "0.044", "0.088", "0.17", "0.26", "0.36",
			      "0.53", "0.71", "0.88", "1.07", "1.24", "1.44"};

  Double_t A_siv =0.0;

  Int_t N_gen =1000;
  Bool_t alphaScale =false;
  //Setup_______________

  TString thisDirPath = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/ToyMC/Data";

  //Aesthetics
  TCanvas* c1 = new TCanvas();
  TCanvas* cTarg = new TCanvas(); cTarg->Divide(2, nPhiScut/2);
  Double_t scaleYaxis =3.0;

  //Loop
  for (Int_t phiS=0; phiS<nPhiScut; phiS++) {
    TString input =
      Form("%s/GeoMean/GeoMean_%0.3f_phiS%s_gen%i_alpha%i.root",
	   thisDirPath.Data(), A_siv, phiScut[phiS].Data(),
	   N_gen, alphaScale);
    TFile *f = OpenFile(input);

    TGraphErrors *g = (TGraphErrors*)f->Get("AN");
    g->SetMarkerColor(phiS+1);

    c1->cd();
    Offset(g, 4.0*phiScut[phiS].Atof()-1 );
    if (phiS==0) {
      Double_t yMin = g->GetY()[0] - g->GetEY()[0]*scaleYaxis;
      Double_t yMax = g->GetY()[0] + g->GetEY()[0]*scaleYaxis;
      g->GetYaxis()->SetRangeUser(yMin, yMax);
      g->SetTitle(Form("Asiv %0.3f", A_siv));
      g->Draw("AP");
      g->GetXaxis()->SetLimits(-0.2, 2*TMath::Pi() );
      DrawLine(g, 0.0);
      DrawLine(g, A_siv, 2);
    }
    else {
      g->Draw("Psame");
    }

    //Individual targets
    TGraphErrors* g_AN_upS_up = (TGraphErrors*)f->Get("AN_upS_up");
    TGraphErrors* g_AN_upS_down = (TGraphErrors*)f->Get("AN_upS_down");
    TGraphErrors* g_AN_downS_up = (TGraphErrors*)f->Get("AN_downS_up");
    TGraphErrors* g_AN_downS_down = (TGraphErrors*)f->Get("AN_downS_down");

    cTarg->cd(phiS+1);
    g_AN_upS_up->GetYaxis()->SetRangeUser(-0.3, 0.3);
    g_AN_upS_up->Draw("AP"); DrawLine(g_AN_upS_up, A_siv);
    g_AN_upS_down->Draw("Psame");
    g_AN_downS_up->Draw("Psame");
    g_AN_downS_down->Draw("Psame");
  }//nPhiScut loop


  //Settings
  cout << "\nSettings Used_______________" << endl;
  cout << "\nPhiS cuts:" << endl;
  for (Int_t i=0; i<nPhiScut; i++) {
    cout << phiScut[i] << "   ";
  }
  cout << "\n\nSivers amplitudes simulated:" << endl;
  cout << A_siv << "   ";
  cout << "\n\nNumber of events Generated:  " << N_gen << endl;
  cout << "Number of events in each target are different:  " <<alphaScale<<endl;
  cout << " " << endl;
}
