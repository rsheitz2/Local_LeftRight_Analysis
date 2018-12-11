#include "include/helperFunctions.h"

void wAvgAll(){
  //Setup_______________
  const Int_t nPhiScut =12;
  TString phiScut[nPhiScut] ={"0.0", "0.044", "0.088", "0.17", "0.26", "0.36",
			      "0.53", "0.71", "0.88", "1.07", "1.24", "1.44"};

  const Int_t nAsiv =6;
  Double_t A_siv[nAsiv] ={0.0, 0.005, 0.01, 0.02, 0.05, 0.1};

  Int_t N_gen =10000;
  Bool_t alphaScale =false;
  //Setup_______________

  TString thisDirPath = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/ToyMC/Data";

  //Aesthetics
  TCanvas* c1 = new TCanvas(); c1->Divide(2, nAsiv/2);
  Double_t scaleYaxis =3.0;

  //Loops
  for (Int_t siv=0; siv<nAsiv; siv++) {
    for (Int_t phiS=0; phiS<nPhiScut; phiS++) {
      TString input =
	Form("%s/GeoMean/GeoMean_%0.3f_phiS%s_gen%i_alpha%i.root",
	     thisDirPath.Data(), A_siv[siv], phiScut[phiS].Data(),
	     N_gen, alphaScale);
      TFile *f = OpenFile(input);

      TGraphErrors *g = (TGraphErrors*)f->Get("Wavg");
      g->SetMarkerColor(phiS+1);

      c1->cd(siv+1);
      Offset(g, 4.0*phiScut[phiS].Atof()-1 );
      if (phiS==0) {
	Double_t yMin = g->GetY()[0] - g->GetEY()[0]*scaleYaxis;
	Double_t yMax = g->GetY()[0] + g->GetEY()[0]*scaleYaxis;
	g->GetYaxis()->SetRangeUser(yMin, yMax);
	g->SetTitle(Form("Asiv %0.3f", A_siv[siv]));
	g->Draw("AP");
	g->GetXaxis()->SetLimits(-0.2, 2*TMath::Pi() );
	DrawLine(g, 0.0);
	DrawLine(g, 2.0*A_siv[siv]/TMath::Pi(), 2);
      }
      else {
	g->Draw("Psame");
      }
    }//nPhiScut loop
  }//siv loop

  //Settings
  cout << "\nSettings Used_______________" << endl;
  cout << "\nPhiS cuts:" << endl;
  for (Int_t i=0; i<nPhiScut; i++) {
    cout << phiScut[i] << "   ";
  }
  cout << "\n\nSivers amplitudes simulated:" << endl;
  for (Int_t i=0; i<nAsiv; i++) {
    cout << A_siv[i] << "   ";
  }
  cout << "\n\nNumber of events Generated:  " << N_gen << endl;
  cout << "Number of events in each target are different:  " <<alphaScale<<endl;
  cout << " " << endl;
}
