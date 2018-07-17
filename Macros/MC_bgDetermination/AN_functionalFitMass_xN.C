#include "helperFunctions.h"

//Setup_______________
const Int_t nBins =5;

TString physBinned = "xPi";//"xN" "xPi" "xF" "pT"
Bool_t toWrite =false;
Bool_t PolCorr =true;
const Double_t Mmin=2.5, Mmax=8.5;

const Int_t nPar = 10;
Int_t Failure=0;
//Setup_______________


//Fit results/Monitoring_______________
Int_t iter=0;

Double_t r_Chi2[nBins*8], r_NDF[nBins*8], r_RedChi2[nBins*8];
Double_t r_Xpoints[nBins*8];

Double_t r_JPsi_M[nBins*8], r_JPsi_W[nBins*8];
Double_t r_psi_M[nBins*8], r_psi_W[nBins*8];
Double_t r_DY_S[nBins*8], r_OC_S[nBins*8], r_CombBg_S[nBins*8];
Double_t r_psi_JPsi_M[nBins*8], r_psiJPsi_M[nBins*8];
Double_t r_psi_JPsi_W[nBins*8], r_psiJPsi_W[nBins*8];

Double_t er_JPsi_M[nBins*8], er_JPsi_W[nBins*8];
Double_t er_psi_M[nBins*8], er_psi_W[nBins*8];
Double_t er_DY_S[nBins*8], er_OC_S[nBins*8], er_CombBg_S[nBins*8];
Double_t er_psi_JPsi_M[nBins*8], er_psiJPsi_M[nBins*8];
Double_t er_psi_JPsi_W[nBins*8], er_psiJPsi_W[nBins*8];

TH1D *r_hRD_upS_up[2*nBins], *r_hRD_upS_down[2*nBins];
TH1D *r_hRD_downS_up[2*nBins], *r_hRD_downS_down[2*nBins];
//Fit results/Monitoring_______________


/*Double_t FitMCs(Double_t *x, Double_t *par){

  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[5] )/par[6];
  Double_t JPsi = par[4]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[8] )/par[9];
  Double_t psi = par[7]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[10]*TMath::Exp( par[11]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
  }*/


Double_t FitMCs(Double_t *x, Double_t *par){

  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
    
  Double_t arg_JPsi = ( x[0] - par[3] )/par[4];
  Double_t JPsi = par[2]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[6] )/par[7];
  Double_t psi = par[5]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );

  return CombBg + JPsi + psi + DY;
}


Double_t GaussInt(Double_t A, Double_t sigma){
  
  return TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
}


Double_t GaussIntError(Double_t A, Double_t sigma,
		       Double_t eA, Double_t eSigma){
  Double_t e = eA*eA/( A*A ) + eSigma*eSigma/( sigma*sigma );
  e = TMath::Sqrt( e );
  e *= TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
  
  return e;
}


Double_t ExpoInt(Double_t A, Double_t b){
  if (b > 0){
    cout << "Cannot do exponential integral" << endl;
    exit(EXIT_FAILURE);
  }
  
  return -1.0*A/b;
}


Double_t FitGetPars(TH1D* h, Int_t bin,
		    Double_t *JPsi, Double_t *psi, Double_t *OC, Double_t *AMDY,
		    Double_t *e_JPsi, Double_t *e_psi, Double_t *e_OC,
		    Double_t *e_AMDY){

  TF1 *fitFunc = new TF1("fitFunc", FitMCs, Mmin, Mmax, nPar);
  fitFunc->SetParameter(0, 8.5e2); fitFunc->SetParameter(1, -3.5);//CombBg
  fitFunc->SetParameter(2, 6.2e3); fitFunc->SetParameter(3, 3.1);//JPsi
  fitFunc->SetParameter(4, 0.15);
  fitFunc->SetParameter(5, 1e2); fitFunc->SetParameter(6, 3.6);//psi
  fitFunc->SetParameter(7, 0.15);
  fitFunc->SetParameter(8, 8.5e2); fitFunc->SetParameter(9, -0.9);//DY
  
  fitFunc->SetParLimits(0, 100, 850); fitFunc->SetParLimits(1, -8.0, -1);//Bg
  
  fitFunc->SetParLimits(2, 1e3, 2e5); fitFunc->SetParLimits(3, 2.5, 3.6);//JPsi
  fitFunc->SetParLimits(4, 0, 1.0);
  
  fitFunc->SetParLimits(5, 10, 2e3); //psi
  fitFunc->SetParLimits(6, 3.5, 4.1); fitFunc->SetParLimits(7, 0.1, 0.5);

  fitFunc->SetParLimits(8, 400, 8e4); fitFunc->SetParLimits(9, -5.0, 0);//Bg
  

  h->Sumw2();
  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ){
    cout << "Fit failed!!" << endl;
    cout << h->GetTitle() << "   bin: " << bin << endl;
    Failure++;
    //exit(EXIT_FAILURE);
  }

  Double_t *pars = fitFunc->GetParameters();
  JPsi[bin] = GaussInt(pars[2], pars[4]);
  psi[bin] = GaussInt(pars[5], pars[7]);
  OC[bin] = 0.0;
  AMDY[bin] = 0.0;

  const Double_t *ePars = fitFunc->GetParErrors();
  e_JPsi[bin] = GaussIntError(pars[2], pars[4], ePars[2], ePars[4] );
  e_psi[bin] = GaussIntError(pars[5], pars[7], ePars[5], ePars[7] );
  e_OC[bin] = 0.0;
  e_AMDY[bin] = 0.0;


  TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
  f_CombBg->SetParameters(pars[0], pars[1] );
  f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

  TF1 *f_JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2])  )",
	    0, Mmax);
  f_JPsi->SetParameters(pars[2], pars[3], pars[4]);
  f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2])  )",
	    0, Mmax);
  f_psi->SetParameters(pars[5], pars[6], pars[7] );
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");

  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
  f_DY->SetParameters(pars[8], pars[9] );
  f_DY->SetLineColor(kBlue); f_DY->Draw("same");

  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();

  //Fit results
  r_Chi2[iter] =Chi2; r_NDF[iter] =ndf; r_RedChi2[iter] =Chi2/ndf;
  if (!iter) r_Xpoints[iter] = 1;
  else if ( !(iter%8) ) r_Xpoints[iter] = r_Xpoints[iter-1]+4;
  else r_Xpoints[iter] = r_Xpoints[iter-1]+1;

  r_JPsi_M[iter] = pars[3]; r_JPsi_W[iter] = pars[4];
  r_psi_M[iter] = pars[6]; r_psi_W[iter] = pars[7];
  er_JPsi_M[iter] = ePars[3]; er_JPsi_W[iter] = ePars[4];
  er_psi_M[iter] = ePars[3]; er_psi_W[iter] = ePars[7];

  r_psi_JPsi_M[iter] =r_psi_M[iter] - r_JPsi_M[iter];
  r_psiJPsi_M[iter] =r_psi_M[iter]/r_JPsi_M[iter];
  r_psi_JPsi_W[iter] =r_psi_W[iter] - r_JPsi_W[iter];
  r_psiJPsi_W[iter] =r_psi_W[iter]/r_JPsi_W[iter];
  er_psi_JPsi_M[iter] = DiffError(er_psi_M[iter], er_JPsi_M[iter]);
  er_psiJPsi_M[iter] = RatioError(r_psi_M[iter], r_JPsi_M[iter],
				  er_psi_M[iter], er_JPsi_M[iter] );
  er_psi_JPsi_W[iter] = DiffError(er_psi_W[iter], er_JPsi_W[iter]);
  er_psiJPsi_W[iter] = RatioError(r_psi_W[iter], r_JPsi_W[iter],
				  er_psi_W[iter], er_JPsi_W[iter] );

  if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_left_upstream_up_%s%i",physBinned.Data(),bin),  21) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_upS_up[2*bin]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_right_upstream_up_%s%i",physBinned.Data(),bin),  22) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_upS_up[2*bin+1]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_left_upstream_down_%s%i",physBinned.Data(),bin),  21) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_upS_down[2*bin]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_right_upstream_down_%s%i",physBinned.Data(),bin),  22) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_upS_down[2*bin+1]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_left_downstream_up_%s%i",physBinned.Data(),bin),  24) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_downS_up[2*bin]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_right_downstream_up_%s%i",physBinned.Data(),bin),  24) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_downS_up[2*bin+1]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_left_downstream_down_%s%i",physBinned.Data(),bin),  26) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_downS_down[2*bin]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
  else if (strncmp(Form("%s", h->GetTitle() ), Form("MuMu_right_downstream_down_%s%i",physBinned.Data(),bin),  26) == 0 ){
    for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
      Double_t x[] = { h->GetBinCenter(bi) };
      r_hRD_downS_down[2*bin+1]->SetBinContent(bi, FitMCs(x, pars));
    }
  }
    
  iter++;
  //Fit results
  
  return Chi2/ndf;
}


Double_t MakeAsym(Double_t L, Double_t R, Double_t P){
  Double_t A = L - R;
  A /= ( L + R );
  A /= P;

  return A;
}


Double_t MakeAsymError(Double_t L, Double_t R, Double_t e_L, Double_t e_R,
		       Double_t P){
  Double_t dL2 = L + e_L*e_L;
  Double_t dR2 = R + e_R*e_R;


  if ( (L< e_L*e_L) || (R< e_R*e_R) ){
    cout << "Fit errors close to Stat errors   ";
    cout << "L " << L << "  e_L " << e_L << "   R " << R <<"   e_R "<<e_R<<endl;
  }

  Double_t e = dL2/( L*L )  + dR2/( R*R );
  e = TMath::Sqrt( e );
  e *= 2.0*L*R/( (L+R)*(L+R) );
  e /= P;
  
  return e;
}


void AN_functionalFitMass_xN(){
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data";

  TFile *fRD = TFile::Open(Form("%s/leftRight_byTarget_WAll_AMDY_%ibins.root",
				pathRD.Data(), nBins) );
  TFile *fRD_noCorr =
    TFile::Open(Form("%s/leftRight_byTarget_WAll_AMDY_%ibins_noCorr.root",
		     pathRD.Data(), nBins) );
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

  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 30, 0, 6);
  TH1D* hChi2_upS = new TH1D("hChi2_upS", "hChi2_upS", 30, 0, 6);
  TH1D* hChi2_downS = new TH1D("hChi2_downS", "hChi2_downS", 30, 0, 6);
  SetUpHist(hChi2); SetUpHist(hChi2_upS); SetUpHist(hChi2_downS);
  
  for (Int_t bi=0, lr=0; bi<nBins; bi++, lr+=2) {
    hRD_upS_up[lr] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s%i", 
					  physBinned.Data(), bi) );
    r_hRD_upS_up[lr] = new TH1D(Form("r_hRD_upS_up_%i", lr),
				Form("r_hRD_upS_up_%i", lr),
				hRD_upS_up[lr]->GetNbinsX(), Mmin, Mmax);
    hRD_upS_up[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_up_%s%i",
					    physBinned.Data(), bi) );
    r_hRD_upS_up[lr+1] = new TH1D(Form("r_hRD_upS_up_%i", lr+1),
				  Form("r_hRD_upS_up_%i", lr+1),
				  hRD_upS_up[lr+1]->GetNbinsX(), Mmin, Mmax);

    hRD_upS_down[lr] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_down_%s%i",
					    physBinned.Data(), bi) );
    r_hRD_upS_down[lr] = new TH1D(Form("r_hRD_upS_down_%i", lr),
				  Form("r_hRD_upS_down_%i", lr),
				  hRD_upS_down[lr]->GetNbinsX(), Mmin, Mmax);
    hRD_upS_down[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_down_%s%i",
					      physBinned.Data(), bi) );
    r_hRD_upS_down[lr+1] = new TH1D(Form("r_hRD_upS_down_%i", lr+1),
				  Form("r_hRD_upS_down_%i", lr+1),
				  hRD_upS_down[lr+1]->GetNbinsX(), Mmin, Mmax);

    hRD_downS_up[lr] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_up_%s%i",
					    physBinned.Data(), bi) );
    r_hRD_downS_up[lr] = new TH1D(Form("r_hRD_downS_up_%i", lr),
				  Form("r_hRD_downS_up_%i", lr),
				  hRD_downS_up[lr]->GetNbinsX(), Mmin, Mmax);
    hRD_downS_up[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_downstream_up_%s%i",
					      physBinned.Data(), bi) );
    r_hRD_downS_up[lr+1] = new TH1D(Form("r_hRD_downS_up_%i", lr+1),
				    Form("r_hRD_downS_up_%i", lr+1),
				    hRD_downS_up[lr+1]->GetNbinsX(), Mmin, Mmax);
    

    hRD_downS_down[lr] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_down_%s%i",
					      physBinned.Data(), bi) );
    r_hRD_downS_down[lr] = new TH1D(Form("r_hRD_downS_down_%i", lr),
				  Form("r_hRD_downS_down_%i", lr),
				  hRD_downS_down[lr]->GetNbinsX(), Mmin, Mmax);
    hRD_downS_down[lr+1] = (TH1D*)fRD->Get(Form("MuMu_right_downstream_down_%s%i",
						physBinned.Data(), bi) );
    r_hRD_downS_down[lr+1] = new TH1D(Form("r_hRD_downS_down_%i", lr+1),
				  Form("r_hRD_downS_down_%i", lr+1),
				  hRD_downS_down[lr+1]->GetNbinsX(), Mmin, Mmax);

    SetUpHist(hRD_upS_up[lr]); SetUpHist(hRD_upS_up[lr+1]);
    SetUpHist(hRD_upS_down[lr]); SetUpHist(hRD_upS_down[lr+1]);
    SetUpHist(hRD_downS_up[lr]); SetUpHist(hRD_downS_up[lr+1]);
    SetUpHist(hRD_downS_down[lr]); SetUpHist(hRD_downS_down[lr+1]); 
    
    c1[0]->cd(2*bi+1); gPad->SetLogy();
    Double_t Chi =FitGetPars(hRD_upS_up[lr], bi,
			     JPsi_upS_up[0], psi_upS_up[0], OC_upS_up[0],
			     AMDY_upS_up[0],
			     e_JPsi_upS_up[0], e_psi_upS_up[0], e_OC_upS_up[0],
			     e_AMDY_upS_up[0]);
    hChi2->Fill(Chi); hChi2_upS->Fill(Chi);
        
    c1[0]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_upS_up[lr+1], bi,
		    JPsi_upS_up[1], psi_upS_up[1], OC_upS_up[1],
		    AMDY_upS_up[1],
		    e_JPsi_upS_up[1], e_psi_upS_up[1], e_OC_upS_up[1],
		    e_AMDY_upS_up[1]);
    hChi2->Fill(Chi); hChi2_upS->Fill(Chi);

    c1[1]->cd(2*bi+1); gPad->SetLogy();
    Chi =FitGetPars(hRD_upS_down[lr], bi,
		    JPsi_upS_down[0], psi_upS_down[0], OC_upS_down[0],
		    AMDY_upS_down[0],
		    e_JPsi_upS_down[0], e_psi_upS_down[0], e_OC_upS_down[0],
		    e_AMDY_upS_down[0]);
    hChi2->Fill(Chi); hChi2_upS->Fill(Chi);
    
    c1[1]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_upS_down[lr+1], bi,
		    JPsi_upS_down[1], psi_upS_down[1], OC_upS_down[1],
		    AMDY_upS_down[1],
		    e_JPsi_upS_down[1], e_psi_upS_down[1], e_OC_upS_down[1],
		    e_AMDY_upS_down[1]);
    hChi2->Fill(Chi); hChi2_upS->Fill(Chi);

    c1[2]->cd(2*bi+1); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_up[lr], bi,
		    JPsi_downS_up[0], psi_downS_up[0], OC_downS_up[0],
		    AMDY_downS_up[0],
		    e_JPsi_downS_up[0], e_psi_downS_up[0], e_OC_downS_up[0],
		    e_AMDY_downS_up[0]);
    hChi2->Fill(Chi); hChi2_downS->Fill(Chi);
    
    c1[2]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_up[lr+1], bi,
		    JPsi_downS_up[1], psi_downS_up[1], OC_downS_up[1],
		    AMDY_downS_up[1],
		    e_JPsi_downS_up[1], e_psi_downS_up[1], e_OC_downS_up[1],
		    e_AMDY_downS_up[1]);
    hChi2->Fill(Chi); hChi2_downS->Fill(Chi);
    
    c1[3]->cd(2*bi+1); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_down[lr], bi,
		    JPsi_downS_down[0], psi_downS_down[0], OC_downS_down[0],
		    AMDY_downS_down[0],
		    e_JPsi_downS_down[0], e_psi_downS_down[0], e_OC_downS_down[0],
		    e_AMDY_downS_down[0]);
    hChi2->Fill(Chi); hChi2_downS->Fill(Chi);
    
    c1[3]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_down[lr+1], bi,
		    JPsi_downS_down[1], psi_downS_down[1], OC_downS_down[1],
		    AMDY_downS_down[1],
		    e_JPsi_downS_down[1], e_psi_downS_down[1], e_OC_downS_down[1],
		    e_AMDY_downS_down[1]);
    hChi2->Fill(Chi); hChi2_downS->Fill(Chi);
  }

  TCanvas* c2[nTargPol];
  for (Int_t c=0; c<nTargPol; c++) {
    c2[c] = new TCanvas(); c2[c]->Divide(2, nBins); }

  for (Int_t bi=0, lr=0; bi<nBins; bi++, lr+=2) {
    c2[0]->cd(2*bi+1); r_hRD_upS_up[2*bi]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_upS_up[2*bi]->Divide(hRD_upS_up[2*bi]);
    SetUpHist(r_hRD_upS_up[2*bi]); r_hRD_upS_up[2*bi]->Draw();
    DrawLine(r_hRD_upS_up[2*bi], 1.0);

    c2[0]->cd(2*bi+2); r_hRD_upS_up[2*bi+1]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_upS_up[2*bi+1]->Divide(hRD_upS_up[2*bi+1]);
    SetUpHist(r_hRD_upS_up[2*bi+1]); r_hRD_upS_up[2*bi+1]->Draw();
    DrawLine(r_hRD_upS_up[2*bi+1], 1.0);

    c2[1]->cd(2*bi+1); r_hRD_upS_down[2*bi]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_upS_down[2*bi]->Divide(hRD_upS_down[2*bi]);
    SetUpHist(r_hRD_upS_down[2*bi]); r_hRD_upS_down[2*bi]->Draw();
    DrawLine(r_hRD_upS_down[2*bi], 1.0);

    c2[1]->cd(2*bi+2); r_hRD_upS_down[2*bi+1]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_upS_down[2*bi+1]->Divide(hRD_upS_down[2*bi+1]);
    SetUpHist(r_hRD_upS_down[2*bi+1]); r_hRD_upS_down[2*bi+1]->Draw();
    DrawLine(r_hRD_upS_down[2*bi+1], 1.0);

    c2[2]->cd(2*bi+1); r_hRD_downS_up[2*bi]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_downS_up[2*bi]->Divide(hRD_downS_up[2*bi]);
    SetUpHist(r_hRD_downS_up[2*bi]); r_hRD_downS_up[2*bi]->Draw();
    DrawLine(r_hRD_downS_up[2*bi], 1.0);

    c2[2]->cd(2*bi+2); r_hRD_downS_up[2*bi+1]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_downS_up[2*bi+1]->Divide(hRD_downS_up[2*bi+1]);
    SetUpHist(r_hRD_downS_up[2*bi+1]); r_hRD_downS_up[2*bi+1]->Draw();
    DrawLine(r_hRD_downS_up[2*bi+1], 1.0);

    c2[3]->cd(2*bi+1); r_hRD_downS_down[2*bi]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_downS_down[2*bi]->Divide(hRD_downS_down[2*bi]);
    SetUpHist(r_hRD_downS_down[2*bi]); r_hRD_downS_down[2*bi]->Draw();
    DrawLine(r_hRD_downS_down[2*bi], 1.0);

    c2[3]->cd(2*bi+2); r_hRD_downS_down[2*bi+1]->GetYaxis()->SetRangeUser(0.3, 1.7);
    r_hRD_downS_down[2*bi+1]->Divide(hRD_downS_down[2*bi+1]);
    SetUpHist(r_hRD_downS_down[2*bi+1]); r_hRD_downS_down[2*bi+1]->Draw();
    DrawLine(r_hRD_downS_down[2*bi+1], 1.0);
  }

  
  TGraph* gChi2 = new TGraph(nBins*8, r_Xpoints, r_Chi2);
  SetUpTGraph(gChi2); gChi2->GetYaxis()->SetTitle("Chi2");
  TGraph* gNDF = new TGraph(nBins*8, r_Xpoints, r_NDF);
  SetUpTGraph(gNDF); gNDF->GetYaxis()->SetTitle("NDF");
  TGraph* gRedChi2 = new TGraph(nBins*8, r_Xpoints, r_RedChi2);
  SetUpTGraph(gRedChi2); gRedChi2->GetYaxis()->SetTitle("RedChi2");

  TCanvas* cChi2 = new TCanvas(); cChi2->Divide(2, 2);
  cChi2->cd(1); hChi2->Draw();
  hChi2->SetLineColor(3); hChi2_upS->Draw("sames");
  hChi2->SetLineColor(7); hChi2_downS->Draw("sames");
  cChi2->cd(2); gChi2->Draw("AP");
  cChi2->cd(3); gNDF->Draw("AP");
  cChi2->cd(4); gRedChi2->Draw("AP");


  Double_t ex8nBins[nBins*8] = {0.0};
  TGraphErrors* gJPsi_M
    = new TGraphErrors(nBins*8, r_Xpoints, r_JPsi_M, ex8nBins, er_JPsi_M);
  TGraphErrors* gJPsi_W
    = new TGraphErrors(nBins*8, r_Xpoints, r_JPsi_W, ex8nBins, er_JPsi_W);
  TGraphErrors* gpsi_M
    = new TGraphErrors(nBins*8, r_Xpoints, r_psi_M, ex8nBins, er_psi_M);
  TGraphErrors* gpsi_W
    = new TGraphErrors(nBins*8, r_Xpoints, r_psi_W, ex8nBins, er_psi_W);
  SetUpTGraph(gJPsi_M); SetUpTGraph(gJPsi_W);
  SetUpTGraph(gpsi_M); SetUpTGraph(gpsi_W);

  TCanvas* cPsi = new TCanvas(); cPsi->Divide(2, 2);
  cPsi->cd(1); gJPsi_M->Draw("AP"); gJPsi_M->SetTitle("JPsi Mass");
  cPsi->cd(2); gJPsi_W->Draw("AP"); gJPsi_M->SetTitle("JPsi W");
  cPsi->cd(3); gpsi_M->Draw("APsame"); gpsi_M->SetMarkerColor(kBlue);
  gpsi_M->SetTitle("psi Mass");
  cPsi->cd(4); gpsi_W->Draw("APsame"); gpsi_W->SetMarkerColor(kBlue);
  gpsi_W->SetTitle("psi W");

  //tmp
  /*for (Int_t i=0, e=0; i<nBins*8; i++) {
    if ( e > 0 ){
      cout << i << endl;
      r_psi_JPsi_M[i] = r_psi_JPsi_M[i-4];
      r_psiJPsi_M[i] = r_psiJPsi_M[i-4];
      
      er_psi_JPsi_M[i] = er_psi_JPsi_M[i-4];
      er_psiJPsi_M[i] = er_psiJPsi_M[i-4];

      if (e == 4) e = 0;
      else e++;
    }
    else if ( !( (i+1)%4 ) ) e = 1;
      
  }//*/
  //tmp
  
  TGraphErrors* gpsi_JPsi_M
    = new TGraphErrors(nBins*8, r_Xpoints, r_psi_JPsi_M,ex8nBins,er_psi_JPsi_M);
  TGraphErrors* gpsiJPsi_M
    = new TGraphErrors(nBins*8, r_Xpoints, r_psiJPsi_M, ex8nBins, er_psiJPsi_M);
  TGraphErrors* gpsi_JPsi_W
    = new TGraphErrors(nBins*8, r_Xpoints, r_psi_JPsi_W,ex8nBins,er_psi_JPsi_W);
  TGraphErrors* gpsiJPsi_W
    = new TGraphErrors(nBins*8, r_Xpoints, r_psiJPsi_W, ex8nBins, er_psiJPsi_W);
  SetUpTGraph(gpsi_JPsi_M); SetUpTGraph(gpsiJPsi_M);
  SetUpTGraph(gpsi_JPsi_W); SetUpTGraph(gpsiJPsi_W);
  
  TCanvas *cDiffPsi = new TCanvas(); cDiffPsi->Divide(2, 2);
  cDiffPsi->cd(1); gpsi_JPsi_M->Draw("AP");gpsi_JPsi_M->SetTitle("psi M shift");
  cDiffPsi->cd(2); gpsiJPsi_M->Draw("AP"); gpsiJPsi_M->SetTitle("psi M ratio");
  cDiffPsi->cd(3); gpsi_JPsi_W->Draw("AP");gpsi_JPsi_W->SetTitle("psi W shift");
  cDiffPsi->cd(4); gpsiJPsi_W->Draw("AP"); gpsiJPsi_W->SetTitle("psi W ratio");
  

  TGraphErrors* gCombBg_S
    = new TGraphErrors(nBins*8, r_Xpoints, r_CombBg_S, ex8nBins, er_CombBg_S);
  SetUpTGraph(gCombBg_S); 
  TCanvas* cSlope = new TCanvas();
  gCombBg_S->Draw("AP"); gCombBg_S->SetTitle("CombBg_slope");

  
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
  
  //Change to psi/JPsi below
  Double_t AN_JPsi_upS_up[nBins], e_AN_JPsi_upS_up[nBins];
  Double_t AN_JPsi_upS_down[nBins], e_AN_JPsi_upS_down[nBins];
  Double_t AN_JPsi_downS_up[nBins], e_AN_JPsi_downS_up[nBins];
  Double_t AN_JPsi_downS_down[nBins], e_AN_JPsi_downS_down[nBins];
  
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
  }


  TGraphErrors* g_JPsi_Left_upS_up =
    new TGraphErrors(nBins, xvals, JPsi_upS_up[0], ex, e_JPsi_upS_up[0]);
  TGraphErrors* g_JPsi_Right_upS_up =
    new TGraphErrors(nBins, xvals, JPsi_upS_up[1], ex, e_JPsi_upS_up[1]);
  TGraphErrors* g_JPsi_Left_upS_down =
    new TGraphErrors(nBins, xvals, JPsi_upS_down[0], ex, e_JPsi_upS_down[0]);
  TGraphErrors* g_JPsi_Right_upS_down =
    new TGraphErrors(nBins, xvals, JPsi_upS_down[1], ex, e_JPsi_upS_down[1]);
  TGraphErrors* g_JPsi_Left_downS_up =
    new TGraphErrors(nBins, xvals, JPsi_downS_up[0], ex, e_JPsi_downS_up[0]);
  TGraphErrors* g_JPsi_Right_downS_up =
    new TGraphErrors(nBins, xvals, JPsi_downS_up[1], ex, e_JPsi_downS_up[1]);
  TGraphErrors* g_JPsi_Left_downS_down =
    new TGraphErrors(nBins, xvals, JPsi_downS_down[0], ex,
		     e_JPsi_downS_down[0]);
  TGraphErrors* g_JPsi_Right_downS_down =
    new TGraphErrors(nBins, xvals, JPsi_downS_down[1], ex,
		     e_JPsi_downS_down[1]);
  
  SetUpTGraph(g_JPsi_Left_upS_up); SetUpTGraph(g_JPsi_Right_upS_up);
  SetUpTGraph(g_JPsi_Left_upS_down); SetUpTGraph(g_JPsi_Right_upS_down);
  SetUpTGraph(g_JPsi_Left_downS_up); SetUpTGraph(g_JPsi_Right_downS_up);
  SetUpTGraph(g_JPsi_Left_downS_down); SetUpTGraph(g_JPsi_Right_downS_down);

  TCanvas* cLR = new TCanvas(); cLR->Divide(2, 2);
  cLR->cd(1); g_JPsi_Left_upS_up->Draw("AP");//Upstream
  g_JPsi_Right_upS_up->Draw("Psame"); g_JPsi_Right_upS_up->SetMarkerColor(kRed);
  cLR->cd(2); g_JPsi_Left_upS_down->Draw("AP");
  g_JPsi_Right_upS_down->Draw("Psame");
  g_JPsi_Right_upS_down->SetMarkerColor(kRed);
  cLR->cd(3); g_JPsi_Left_downS_up->Draw("AP");//Downstream
  g_JPsi_Right_downS_up->Draw("Psame");
  g_JPsi_Right_downS_up->SetMarkerColor(kRed);
  cLR->cd(4); g_JPsi_Left_downS_down->Draw("AP");
  g_JPsi_Right_downS_down->Draw("Psame");
  g_JPsi_Right_downS_down->SetMarkerColor(kRed);
  

  TGraphErrors *g_AN_JPsi_upS_up = new TGraphErrors(nBins, xvals,AN_JPsi_upS_up,
						    ex, e_AN_JPsi_upS_up);
  TGraphErrors *g_AN_JPsi_upS_down =
    new TGraphErrors(nBins, xvals,AN_JPsi_upS_down, ex, e_AN_JPsi_upS_down);
  TGraphErrors *g_AN_JPsi_downS_up =
    new TGraphErrors(nBins, xvals,AN_JPsi_downS_up, ex, e_AN_JPsi_downS_up);
  TGraphErrors *g_AN_JPsi_downS_down =
    new TGraphErrors(nBins, xvals,AN_JPsi_downS_down, ex, e_AN_JPsi_downS_down);
  
  SetUpTGraph(g_AN_JPsi_upS_up); SetUpTGraph(g_AN_JPsi_upS_down);
  SetUpTGraph(g_AN_JPsi_downS_up); SetUpTGraph(g_AN_JPsi_downS_down);
  

  TCanvas* cAN = new TCanvas(); cAN->Divide(2, 2);
  Double_t yMax = 0.5;
  cAN->cd(1);
  g_AN_JPsi_upS_up->Draw("AP");
  g_AN_JPsi_upS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
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

    g_JPsi_Left_upS_up->Write("JPsi_Left_upS_up");
    g_JPsi_Left_upS_down->Write("JPsi_Left_upS_down");
    g_JPsi_Left_downS_up->Write("JPsi_Left_downS_up");
    g_JPsi_Left_downS_down->Write("JPsi_Left_downS_down");

    g_JPsi_Right_upS_up->Write("JPsi_Right_upS_up");
    g_JPsi_Right_upS_down->Write("JPsi_Right_upS_down");
    g_JPsi_Right_downS_up->Write("JPsi_Right_downS_up");
    g_JPsi_Right_downS_down->Write("JPsi_Right_downS_down");
  }

  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Mass Range is: " << Mmin << "  -  " << Mmax << endl;
  cout << "Polarization was performed:  " << PolCorr << endl;
  cout << "Physics binning is:  " << physBinned << endl;
  cout << "Fit failed:  " << Failure << "  times" << endl;
  cout << "Number of histogram bins used:  "<< hRD_upS_up[0]->GetNbinsX()<<endl;
  cout << "Outputting asymmetry for    JPsi" << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
  
}
