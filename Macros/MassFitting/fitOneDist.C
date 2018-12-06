#include "include/helperFunctions.h"

Double_t Get_twelve_Ratio(TString targ){
  if (targ =="UpS") return 1.1485;
  else if (targ =="DownS") return 1.1418;
  else{
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t twelve_CrystalBall(Double_t *x, Double_t *par){
  //f(x; alpha, n, xBar, sigma, Amp)
  //   Amp = par[0], xBar = par[1], sigma = par[2]
  //   alpha = par[3], n = par[4]
  Double_t A = TMath::Power(par[4]/par[3], par[4])*
    TMath::Exp(-par[3]*par[3]/2.0);
  
  Double_t B = par[4]/par[3] - par[3];

  Double_t C
    =(par[4]/par[3])*(1.0/(par[4]-1.0))*TMath::Exp(-par[3]*par[3]/2.0);
  
  Double_t D
    =TMath::Sqrt( TMath::Pi()/2.0 )*(1+TMath::Erf( par[3]/TMath::Sqrt(2) ) );

  Double_t Norm = 1.0/(par[2]*(C+D) );

  Double_t arg = (x[0] - par[1])/par[2];
  if (arg > -par[3] ) return par[0]*Norm*TMath::Exp(-0.5*arg*arg);
  else return par[0]*Norm*A*TMath::Power((B - arg), -par[4]);
}


Double_t Fit_twelve_UpS(Double_t *x, Double_t *par){
  //Double_t xShift = x[0]-par[11];
  //Double_t xShift = x[0]-2.00;
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_twelve_Ratio("UpS");

  Double_t par_JPsi[] = {par[0], par[1], par[2], par[3], par[4]};
  Double_t JPsi = twelve_CrystalBall(x, par_JPsi);
  
  Double_t par_psi[] =
  {par[5], par[1]*ratioPsi, par[2]*ratioPsi, par[3]*ratioPsi, par[4]};
  //{par[5], par[1]*ratioPsi, par[2]*ratioPsi, par[3], par[4]};
  Double_t psi = twelve_CrystalBall(x, par_psi);

  //Double_t CombBg = par[7]*TMath::Exp( par[8]*xShift );
  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  //Double_t DY = par[9]*TMath::Exp( par[10]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_twelve(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax, Bool_t hIsUpS){
  //Normal parameter setup
  //     Fixed parameter setup
  Double_t massJPsi =3.131, widthJPsi=0.17;
  Double_t alphaJPsi =1.5, nJPsi =5.0;
  Double_t Bg_slope =-1.85;
  Double_t DY_slope = -0.90, DY_Mmin = 4.5;
    
  //     variable parameter setup
  Double_t ratioPsi = (hIsUpS) ?
    Get_twelve_Ratio("UpS") : Get_twelve_Ratio("DownS");
  
  Double_t A_JPsi = h->GetBinContent(h->FindBin(massJPsi) );
  Double_t A_psi = h->GetBinContent(h->FindBin(massJPsi*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(DY_Mmin) );

  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));
  A_Bg = A_Bg - A_DY;

  Double_t C
    =(nJPsi/alphaJPsi)*(1.0/(nJPsi-1.0))*TMath::Exp(-alphaJPsi*alphaJPsi/2.0);
  Double_t D
    =TMath::Sqrt( TMath::Pi()/2.0 )*(1+TMath::Erf( alphaJPsi/TMath::Sqrt(2) ) );
  Double_t Norm = 1.0/(widthJPsi*(C+D) );

  A_JPsi /= Norm;
  A_psi /= Norm;

  if (A_JPsi == 0) A_JPsi = 200.;
  if (A_psi == 0) A_psi = 200.;
  if (A_Bg <= 0) A_Bg = 100.;
  if (A_DY == 0) A_DY = 200.;

  Int_t ipar =0;
  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, massJPsi); ipar++;
  fitFunc->SetParameter(ipar, widthJPsi); ipar++;
  fitFunc->SetParameter(ipar, alphaJPsi); ipar++;
  fitFunc->SetParameter(ipar, nJPsi); ipar++;
  
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  //fitFunc->SetParameter(ipar, alphaJPsi); ipar++;
  
  fitFunc->SetParameter(ipar, A_Bg); ipar++;//CombBg
  fitFunc->SetParameter(ipar, Bg_slope); ipar++;
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, DY_slope); ipar++;
  fitFunc->SetParameter(ipar, Mmin); ipar++;//Mmin
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }

  /*ipar=0;
  Double_t factor =5.0;
  //fitFunc->SetParLimits(0, A_JPsi/(factor), A_JPsi*factor); ipar++;//JPsi
  //fitFunc->SetParLimits(1, 3.0, 3.6); ipar++;
  //fitFunc->SetParLimits(2, 0.1, 0.22); ipar++;

  fitFunc->SetParLimits(3, 0.5, 5.0); ipar++;
  fitFunc->SetParLimits(4, 1.0, 10.0); ipar++;
  
  fitFunc->SetParLimits(5, A_psi/(factor), A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(6, 0.5, 5.0); ipar++;
  //
  //fitFunc->SetParLimits(7, 0.0, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(8, -7.0, -1.2); ipar++;
  //fitFunc->SetParLimits(9, 0.0, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(10, -2.0, 0.0); ipar++;//*/
  
  fitFunc->SetParLimits(10, Mmin, Mmin); ipar++;//Mmin
}


TF1* ComponentFuncts_twelve(Double_t *pars, Double_t Mmin, Double_t Mmax,
			   TString process, Bool_t hIsUpS){
  Double_t psiMW = (hIsUpS) ? Get_twelve_Ratio("UpS") : Get_twelve_Ratio("DownS");
  
  TF1 *JPsi = new TF1("f_JPsi", twelve_CrystalBall, Mmin, Mmax, 5);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
  
  TF1 *psi = new TF1("f_psi", twelve_CrystalBall, Mmin, Mmax, 5);
  //psi->SetParameters(pars[5], pars[1]*psiMW, pars[2]*psiMW, pars[6], pars[4]);
  psi->SetParameters(pars[5], pars[1]*psiMW, pars[2]*psiMW, pars[3]*psiMW, pars[4]);

  TF1 *Bg = new TF1("f_Bg", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  //Bg->SetParameters(pars[7], pars[8], Mmin);
  Bg->SetParameters(pars[6], pars[7], Mmin);
  
  TF1 *DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  //DY->SetParameters(pars[9], pars[10], Mmin);
  DY->SetParameters(pars[8], pars[9], Mmin);

  JPsi->SetLineColor(kGreen); 
  psi->SetLineColor(kGreen);
  Bg->SetLineColor(kGreen); 
  DY->SetLineColor(kGreen); 

  JPsi->Draw("same");
  psi->Draw("same");
  Bg->Draw("same");
  DY->Draw("same");
  
  if (process == "JPsi") {
    JPsi->SetLineColor(kBlack);
    return JPsi;
  }
  else if (process == "psi") {
    psi->SetLineColor(kBlack);
    return psi;
  }
  else if (process == "DY") {
    DY->SetLineColor(kBlack);
    return DY;
  }
  else {
    cout << "Wrong input process Fit_twelve" << endl;
    exit(EXIT_FAILURE);
  }
}


void fitOneDist(){
  //Setup______
  TString period_Mtype ="W07_LowM_AMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;
  TString binRange ="25_43";
  TString physBinned ="pT";//"xN", "xPi", "xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.00;
  Double_t LR_Mmax =5.00;//L/R counts mass range
  Double_t Mmin =2.00;//Fit Mass minimum
  Double_t Mmax =7.50;//Fit Mass maximum
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.195"; //HMDY=0.1866, #LowM_AMDY=0.195
  TString whichFit ="twelve";
  //Setup______

  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =
    Form("leftRight_byTarget_%s1.00_8.50_%ibins%s_%ihbin_%s_%s.root",
	 period_Mtype.Data(), nBins, binRange.Data(), hbins, production.Data(),
	 additionalCuts.Data());
  TFile *fRD  = OpenFile(pathRD + RDfile);

  TCanvas* c1 = new TCanvas(); gPad->SetLogy();
  TH1D *hRD =
    (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s1", physBinned.Data()) );
  hRD->Sumw2();

  Int_t nPar =11;
  TF1 *fitFunc = new TF1("fitFunc", Fit_twelve_UpS, Mmin, Mmax, nPar);
  Paras_twelve(hRD, fitFunc, nPar, Mmin, Mmax, true);
  hRD->Fit(fitFunc, "RL", "", Mmin, Mmax);

  Double_t *pars = fitFunc->GetParameters();
  TF1 *f_LR =ComponentFuncts_twelve(pars, Mmin, Mmax, "JPsi", true);
}
