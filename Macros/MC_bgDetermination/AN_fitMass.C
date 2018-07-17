#include "helperFunctions.h"

//Setup_______________
const Int_t nBins =5;

TString physBinned = "pT";//"xN" "xPi" "xF" "pT"
Bool_t toWrite=false; Bool_t PolCorr =true;
const Double_t Mmin=2.5, Mmax=8.5;

const Int_t nFits = 4;
TH1D *hFits[nFits];
TString fitTypes[nFits] = {"JPsi", "psi", "OC", "AMDY"};
//Setup_______________  


Double_t FitMCs(Double_t *x, Double_t *par){
  Double_t nJPsi = hFits[0]->GetBinContent(hFits[0]->FindBin(x[0] ) );
  Double_t npsi = hFits[1]->GetBinContent(hFits[1]->FindBin(x[0] ) );
  Double_t nOC = hFits[2]->GetBinContent(hFits[2]->FindBin(x[0] ) );
  Double_t nAMDY = hFits[3]->GetBinContent(hFits[3]->FindBin(x[0] ) );
    
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}


Double_t FitGetPars(TH1D* h, Int_t bin,
		    Double_t *JPsi, Double_t *psi, Double_t *OC, Double_t *AMDY,
		    Double_t *e_JPsi, Double_t *e_psi, Double_t *e_OC,
		    Double_t *e_AMDY){

  TF1 *fitFunc = new TF1("fitFunc", FitMCs, Mmin, Mmax, nFits);
  fitFunc->SetParLimits(0, 0, 2e6); fitFunc->SetParLimits(1, 0, 2e6);
  fitFunc->SetParLimits(2, 0, 2e6); fitFunc->SetParLimits(3, 0, 2e6);
  
  TFitResultPtr status = h->Fit("fitFunc", "RWLS", "", Mmin, Mmax);
  if (status->Status() ){
    cout << "Fit failed!!" << endl;
    exit(EXIT_FAILURE);
  }

  JPsi[bin] = fitFunc->GetParameter(0);
  psi[bin] = fitFunc->GetParameter(1);
  OC[bin] = fitFunc->GetParameter(2);
  AMDY[bin] = fitFunc->GetParameter(3);

  e_JPsi[bin] = fitFunc->GetParError(0);
  e_psi[bin] = fitFunc->GetParError(1);
  e_OC[bin] = fitFunc->GetParError(2);
  e_AMDY[bin] = fitFunc->GetParError(3);

  return fitFunc->GetChisquare();
}


Double_t MakeAsym(Double_t L, Double_t R, Double_t P){
  Double_t A = L - R;
  A /= ( L + R );
  A /= P;

  return A;
}


//Only errors from fit
/*Double_t MakeAsymError(Double_t L, Double_t R, Double_t e_L, Double_t e_R){
  Double_t e = e_L*e_L/( L*L ) + e_R*e_R/( R*R );
  e = TMath::Sqrt( e );
  e *= 2.0*L*R/( (L+R)*(L+R) );

  return e;
}//*/


//Only statistical errors
Double_t MakeAsymError(Double_t L, Double_t R, Double_t e_L, Double_t e_R,
		       Double_t P){
  Double_t e = 1.0/L  + 1.0/R;
  e = TMath::Sqrt( e );
  e *= 2.0*L*R/( (L+R)*(L+R) );
  e /= P;
  
  return e;
}


void AN_fitMass(){
    
  TFile *fFits[nFits]; 
  for (Int_t f=0; f<nFits; f++) {
    fFits[f] = TFile::Open(Form("FitDist/MC_%s_%.2f_%.2f.root",
				fitTypes[f].Data(), Mmin,Mmax));
    if (!fFits[f]) {
      cout << "File:   " << f << "   does not exist " << endl;
      exit(EXIT_FAILURE);
    }

    hFits[f] = (TH1D*)fFits[f]->Get(Form("h_%s", fitTypes[f].Data() ));
  }

  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";

  TFile *fRD = TFile::Open(pathRD+"leftRight_byTarget_WAll_AMDY_5bins.root");
  TFile *fRD_noCorr =
    TFile::Open(pathRD+"leftRight_byTarget_WAll_AMDY_5bins_noCorr.root");
  if (!fRD || !fRD_noCorr ){
    cout << "RD or RD_noCorr file does not exist " << endl;
    exit(EXIT_FAILURE);
  }
  
  const Int_t nTargPol =4;
  TCanvas* c1[nTargPol];
  for (Int_t c=0; c<nTargPol; c++) {
    c1[c] = new TCanvas(); c1[c]->Divide(2, nBins); }

  TH1D *hRD_upS_up[2*nBins], *hRD_upS_down[2*nBins];
  TH1D *hRD_downS_up[2*nBins], *hRD_downS_down[2*nBins];

  //{ {Left}, {Right} }
  Double_t JPsi_upS_up[2][nBins], e_JPsi_upS_up[2][nBins];
  Double_t psi_upS_up[2][nBins], e_psi_upS_up[2][nBins];
  Double_t OC_upS_up[2][nBins], e_OC_upS_up[2][nBins];
  Double_t AMDY_upS_up[2][nBins], e_AMDY_upS_up[2][nBins];

  Double_t JPsi_upS_down[2][nBins], e_JPsi_upS_down[2][nBins];
  Double_t psi_upS_down[2][nBins], e_psi_upS_down[2][nBins];
  Double_t OC_upS_down[2][nBins], e_OC_upS_down[2][nBins];
  Double_t AMDY_upS_down[2][nBins], e_AMDY_upS_down[2][nBins];

  Double_t JPsi_downS_up[2][nBins], e_JPsi_downS_up[2][nBins];
  Double_t psi_downS_up[2][nBins], e_psi_downS_up[2][nBins];
  Double_t OC_downS_up[2][nBins], e_OC_downS_up[2][nBins];
  Double_t AMDY_downS_up[2][nBins], e_AMDY_downS_up[2][nBins];

  Double_t JPsi_downS_down[2][nBins], e_JPsi_downS_down[2][nBins];
  Double_t psi_downS_down[2][nBins], e_psi_downS_down[2][nBins];
  Double_t OC_downS_down[2][nBins], e_OC_downS_down[2][nBins];
  Double_t AMDY_downS_down[2][nBins], e_AMDY_downS_down[2][nBins];
  
  for (Int_t bi=0, lr=0; bi<nBins; bi++, lr+=2) {

    hRD_upS_up[lr] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s%i", 
					  physBinned.Data(), bi) );
    hRD_upS_up[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_up_%s%i",
					    physBinned.Data(), bi) );

    hRD_upS_down[lr] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_down_%s%i",
					    physBinned.Data(), bi) );
    hRD_upS_down[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_down_%s%i",
					      physBinned.Data(), bi) );

    hRD_downS_up[lr] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_up_%s%i",
					    physBinned.Data(), bi) );
    hRD_downS_up[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_downstream_up_%s%i",
					      physBinned.Data(), bi) );

    hRD_downS_down[lr] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_down_%s%i",
					      physBinned.Data(), bi) );
    hRD_downS_down[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_downstream_down_%s%i",
						physBinned.Data(), bi) );
    
    c1[0]->cd(2*bi+1);
    FitGetPars(hRD_upS_up[lr], bi,
	       JPsi_upS_up[0], psi_upS_up[0], OC_upS_up[0], AMDY_upS_up[0],
	       e_JPsi_upS_up[0], e_psi_upS_up[0], e_OC_upS_up[0], e_AMDY_upS_up[0]);
    c1[0]->cd(2*bi+2);
    FitGetPars(hRD_upS_up[lr+1], bi,
	       JPsi_upS_up[1], psi_upS_up[1], OC_upS_up[1], AMDY_upS_up[1],
	       e_JPsi_upS_up[1], e_psi_upS_up[1], e_OC_upS_up[1], e_AMDY_upS_up[1]);

    c1[1]->cd(2*bi+1);
    FitGetPars(hRD_upS_down[lr], bi,
	       JPsi_upS_down[0], psi_upS_down[0], OC_upS_down[0], AMDY_upS_down[0],
	       e_JPsi_upS_down[0], e_psi_upS_down[0], e_OC_upS_down[0], e_AMDY_upS_down[0]);
    c1[1]->cd(2*bi+2);
    FitGetPars(hRD_upS_down[lr+1], bi,
	       JPsi_upS_down[1], psi_upS_down[1], OC_upS_down[1], AMDY_upS_down[1],
	       e_JPsi_upS_down[1], e_psi_upS_down[1], e_OC_upS_down[1], e_AMDY_upS_down[1]);

    c1[2]->cd(2*bi+1);
    FitGetPars(hRD_downS_up[lr], bi,
	       JPsi_downS_up[0], psi_downS_up[0], OC_downS_up[0], AMDY_downS_up[0],
	       e_JPsi_downS_up[0], e_psi_downS_up[0], e_OC_downS_up[0], e_AMDY_downS_up[0]);
    c1[2]->cd(2*bi+2);
    FitGetPars(hRD_downS_up[lr+1], bi,
	       JPsi_downS_up[1], psi_downS_up[1], OC_downS_up[1], AMDY_downS_up[1],
	       e_JPsi_downS_up[1], e_psi_downS_up[1], e_OC_downS_up[1], e_AMDY_downS_up[1]);
    
    c1[3]->cd(2*bi+1);
    FitGetPars(hRD_downS_down[lr], bi,
	       JPsi_downS_down[0], psi_downS_down[0], OC_downS_down[0], AMDY_downS_down[0],
	       e_JPsi_downS_down[0], e_psi_downS_down[0], e_OC_downS_down[0], e_AMDY_downS_down[0]);
    c1[3]->cd(2*bi+2);
    FitGetPars(hRD_downS_down[lr+1], bi,
	       JPsi_downS_down[1], psi_downS_down[1], OC_downS_down[1], AMDY_downS_down[1],
	       e_JPsi_downS_down[1], e_psi_downS_down[1], e_OC_downS_down[1], e_AMDY_downS_down[1]);
  }


  Double_t ex[nBins] = {0.0};
  TGraphErrors* g_asym
    =(TGraphErrors*)fRD->Get(Form("%s_asym",physBinned.Data()));
  TGraphErrors* g_asym_noCorr
    =(TGraphErrors*)fRD_noCorr->Get(Form("%s_asym",physBinned.Data()));
  if (g_asym->GetN() != nBins){
    cout << "nBins not defined well!!!" << endl;
    exit(EXIT_FAILURE);
  }
  Double_t *xvals = g_asym->GetX();
  Double_t *yvals =g_asym->GetY();
  Double_t *yvals_noCorr =g_asym_noCorr->GetY();
  Double_t Pol[nBins];
  if (PolCorr){
    for (Int_t bi=0; bi<nBins; bi++) Pol[bi] = yvals_noCorr[bi]/yvals[bi];
  }
  else {
    for (Int_t bi=0; bi<nBins; bi++) Pol[bi] = 1.0;
  }
  
  Double_t AN_JPsi_upS_up[nBins], e_AN_JPsi_upS_up[nBins];
  Double_t AN_JPsi_upS_down[nBins], e_AN_JPsi_upS_down[nBins];
  Double_t AN_JPsi_downS_up[nBins], e_AN_JPsi_downS_up[nBins];
  Double_t AN_JPsi_downS_down[nBins], e_AN_JPsi_downS_down[nBins];
  
  Double_t AN_psi_upS_up[nBins], e_AN_psi_upS_up[nBins];
  //Double_t AN_psi_upS_down[nBins], e_AN_psi_upS_down[nBins];
  
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_JPsi_upS_up[bi] = MakeAsym(JPsi_upS_up[0][bi], JPsi_upS_up[1][bi],
				  Pol[bi]);
    e_AN_JPsi_upS_up[bi] =
      MakeAsymError(JPsi_upS_up[0][bi], JPsi_upS_up[1][bi],
		    e_JPsi_upS_up[0][bi], e_JPsi_upS_up[1][bi], Pol[bi]);
    AN_JPsi_upS_down[bi] = MakeAsym(JPsi_upS_down[0][bi], JPsi_upS_down[1][bi],
				    Pol[bi]);
    e_AN_JPsi_upS_down[bi] =
      MakeAsymError(JPsi_upS_down[0][bi], JPsi_upS_down[1][bi],
		    e_JPsi_upS_down[0][bi], e_JPsi_upS_down[1][bi], Pol[bi]);
    AN_JPsi_downS_up[bi] = MakeAsym(JPsi_downS_up[0][bi], JPsi_downS_up[1][bi],
				    Pol[bi]);
    e_AN_JPsi_downS_up[bi] =
      MakeAsymError(JPsi_downS_up[0][bi], JPsi_downS_up[1][bi],
		    e_JPsi_downS_up[0][bi], e_JPsi_downS_up[1][bi], Pol[bi]);
    AN_JPsi_downS_down[bi] = MakeAsym(JPsi_downS_down[0][bi],
				      JPsi_downS_down[1][bi], Pol[bi]);
    e_AN_JPsi_downS_down[bi] =
      MakeAsymError(JPsi_downS_down[0][bi], JPsi_downS_down[1][bi],
		    e_JPsi_downS_down[0][bi], e_JPsi_downS_down[1][bi],Pol[bi]);
    

    AN_psi_upS_up[bi] = MakeAsym(psi_upS_up[0][bi], psi_upS_up[1][bi], Pol[bi]);
    e_AN_psi_upS_up[bi] =
      MakeAsymError(psi_upS_up[0][bi], psi_upS_up[1][bi],
		    e_psi_upS_up[0][bi], e_psi_upS_up[1][bi], Pol[bi]);
  }

  TGraphErrors *g_AN_JPsi_upS_up = new TGraphErrors(nBins, xvals,AN_JPsi_upS_up,
						    ex, e_AN_JPsi_upS_up);
  TGraphErrors *g_AN_JPsi_upS_down =
    new TGraphErrors(nBins, xvals,AN_JPsi_upS_down, ex, e_AN_JPsi_upS_down);
  TGraphErrors *g_AN_JPsi_downS_up =
    new TGraphErrors(nBins, xvals,AN_JPsi_downS_up, ex, e_AN_JPsi_downS_up);
  TGraphErrors *g_AN_JPsi_downS_down =
    new TGraphErrors(nBins, xvals,AN_JPsi_downS_down, ex, e_AN_JPsi_downS_down);
  
  TGraphErrors *g_AN_psi_upS_up = new TGraphErrors(nBins, xvals,AN_psi_upS_up,
						    ex, e_AN_psi_upS_up);
  SetUpTGraph(g_AN_JPsi_upS_up); SetUpTGraph(g_AN_JPsi_upS_down);
  SetUpTGraph(g_AN_JPsi_downS_up); SetUpTGraph(g_AN_JPsi_downS_down); 
  
  SetUpTGraph(g_AN_psi_upS_up, 3, 0.002, nBins);

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
  TCanvas* cAN = new TCanvas(); cAN->Divide(2, 2);
  Double_t yMax = 0.5;
  cAN->cd(1);
  g_AN_JPsi_upS_up->Draw("AP");
  g_AN_JPsi_upS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
  g_AN_psi_upS_up->Draw("Psame");
  DrawLine(g_AN_JPsi_upS_up, 0.0);

  cAN->cd(2);
  g_AN_JPsi_upS_down->Draw("AP");
  g_AN_JPsi_upS_down->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_JPsi_upS_down, 0.0);

  cAN->cd(3);
  g_AN_JPsi_downS_up->Draw("AP");
  g_AN_JPsi_downS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_JPsi_downS_up, 0.0);

  cAN->cd(4);
  g_AN_JPsi_downS_down->Draw("AP");
  g_AN_JPsi_downS_down->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_JPsi_downS_down, 0.0);


  TString fOutput = Form("FitMass_%s_%.2f_%.2f_",
			   physBinned.Data(), Mmin, Mmax);
  if (PolCorr) fOutput += "corr.root";
  else fOutput += "noCorr.root";
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_AN_JPsi_upS_up->Write("JPsi_upS_up");
    g_AN_JPsi_upS_down->Write("JPsi_upS_down");
    g_AN_JPsi_downS_up->Write("JPsi_downS_up");
    g_AN_JPsi_downS_down->Write("JPsi_downS_down");    
  }

  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Mass Range is: " << Mmin << "  -  " << Mmax << endl;
  cout << "Polarization was performed:  " << PolCorr << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
