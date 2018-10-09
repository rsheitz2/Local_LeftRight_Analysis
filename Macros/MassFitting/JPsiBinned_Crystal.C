#include "include/helperFunctions.h"

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_eight.h"

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//1 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_nine.h"


void DoFit(TH1D *h, TMatrixDSym &cov, Double_t Mmin, Double_t Mmax){
  h->Sumw2(); 
  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ){
    cout << "Fit failed!!" << endl;
    cout << h->GetTitle() << endl;
    //exit(EXIT_FAILURE);
  }

  cov = status->GetCovarianceMatrix();
}


Double_t FitGetPars(TH1D **h, TH1D **hRatio, Int_t bin, Double_t *LR,
		    Double_t *e_LR, Double_t LR_Mmin, Double_t LR_Mmax,
		    TString whichFit, Double_t Mmin, Double_t Mmax,
		    Double_t MWpars[][2], Double_t e_MWpars[][2]){
  
  Double_t processPars[8] = {0.0};//{aJPsi,mJPsi,wJPsi,Apsi,Mpsi,Wpsi,aDY,cDY}
  Double_t e_processPars[8] = {0.0};//{aJPsi,mJPsi,wJPsi,Apsi,Mpsi,Wpsi,aDY,cDY}
  Double_t LR_cov[25];
  
  Bool_t hIsUpS =false;
  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_upstream", 18) == 0)
    hIsUpS =true;
  else if (strncmp(Form("%s", h[bin]->GetTitle() ),"MuMu_right_upstream",19)==0)
    hIsUpS =true;

  //Fit Setup_____
  TF1 *fitFunc = NULL;
  Int_t nPar;
  if (whichFit =="eight"){
    fitFunc = SetupFunc_eight(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_eight(fitFunc, processPars, LR_cov, cov, "JPsi", nPar, hIsUpS,
		    e_processPars);
  }
  else if (whichFit =="nine"){
    fitFunc = SetupFunc_nine(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_nine(fitFunc, processPars, LR_cov, cov, "JPsi", nPar, hIsUpS,
		     e_processPars);
  }
  else{
    cout << "Int valid fit type:   " << whichFit << endl;
    exit(EXIT_FAILURE);
  }

  //Draw physics processes
  Double_t *pars = fitFunc->GetParameters();
  TF1 *f_JPsi =
    new TF1("f_JPsi", f_CrystalBall, Mmin, Mmax, 5);
  f_JPsi->SetParameters(processPars[0], processPars[1], processPars[2],
			pars[3], pars[4]);
  f_JPsi->SetLineColor(kGreen); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", f_CrystalBall, Mmin, Mmax, 5);
  f_psi->SetParameters(processPars[3], processPars[4], processPars[5],
			pars[3], pars[4]);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");

  Double_t DY_pars[] = {processPars[6], processPars[7], Mmin};
  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  f_DY->SetParameters(DY_pars[0], DY_pars[1], DY_pars[2]);
  f_DY->SetLineColor(kGreen); f_DY->Draw("same");

  //TF1 *f_Bg = new TF1("f_Bg", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  //f_Bg->SetParameters(pars[4], pars[5], Mmin);
  //f_Bg->SetLineColor(kGreen); f_Bg->Draw("same");

  //Integrate for L/R counts
  f_JPsi->SetLineColor(kBlack);

  LR_Mmin = processPars[1]-0.15;
  LR_Mmax = processPars[1]+0.15;
  cout << LR_Mmin << " " << LR_Mmax << endl;
  LR[bin] = f_JPsi->Integral(LR_Mmin, LR_Mmax);
  e_LR[bin] = f_JPsi->IntegralError(LR_Mmin, LR_Mmax, processPars, LR_cov);
  MWpars[bin][0] = processPars[1];
  MWpars[bin][1] = processPars[2];
  e_MWpars[bin][0] = e_processPars[1];
  e_MWpars[bin][1] = e_processPars[2];
  /*if (process =="JPsi"){
    f_JPsi->SetLineColor(kBlack);
    
    LR[bin] = f_JPsi->Integral(LR_Mmin, LR_Mmax);
    e_LR[bin] = f_JPsi->IntegralError(LR_Mmin, LR_Mmax, processPars, LR_cov);
  }
  else if (process =="psi"){
    f_psi->SetLineColor(kBlack); 
  
    LR[bin] = f_psi->Integral(LR_Mmin, LR_Mmax);
    e_LR[bin] = f_psi->IntegralError(LR_Mmin, LR_Mmax, &processPars[3], LR_cov);
    }*/

  TGraphErrors *gError = new TGraphErrors(hRatio[bin]->GetNbinsX() );
  for (Int_t i=1; i<hRatio[bin]->GetNbinsX()+1; i++) 
    gError->SetPoint(i, h[bin]->GetBinCenter(i), 0);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gError);
  
  for (Int_t bi=1; bi<h[bin]->GetNbinsX()+1; bi++) {
    Double_t numerator = hRatio[bin]->GetBinContent(bi);
    Double_t num_error = hRatio[bin]->GetBinError(bi);
    
    Double_t fitValue = gError->Eval(h[bin]->GetBinCenter(bi) );
    Double_t dem_error = gError->GetErrorY(bi);

    Double_t error = RatioError(numerator, fitValue, num_error, dem_error);
    Double_t ratio = numerator/fitValue;
    
    hRatio[bin]->SetBinContent(bi, ratio);
    hRatio[bin]->SetBinError(bi, error);
  }
  
  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();
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


void JPsiBinned_Crystal(TString start=""){
  //Setup___
  const Int_t nBins =5;
  
  //Double_t fitMmin[nBins] ={2.5, 2.5, 2.5, 2.5, 2.5};
  Double_t comFitMin = 2.2;
  Double_t fitMmin[nBins] ={comFitMin, comFitMin, comFitMin, comFitMin, comFitMin};

  Double_t comFitMax = 6.0;
  Double_t fitMmax[nBins] ={comFitMax, comFitMax, comFitMax, comFitMax, comFitMax};
  //Double_t fitMmax[nBins] ={4.50, 5.0, 6.0, 6.50, 8.50};
  //Double_t fitMmax[nBins] ={8.50, 8.50, 8.50, 8.50, 8.50};
  
  Double_t LR_Mmin =2.9, LR_Mmax =3.3;
  TString physBinned ="pT";//xN, xPi, xF, pT
  TString fname ="leftRight_byTarget_WAll_LowM_AMDY1.00_8.50_5bins25_43_150hbin.root";
  //TString whichFit[nBins] ={"seven", "seven", "seven", "six", "six"};
  TString whichFit[nBins] ={"eight", "eight", "eight", "eight", "eight"};
  //TString whichFit[nBins] ={"nine", "nine", "nine", "nine", "nine"};
  
  Bool_t toWrite =false;
  //Setup___
    
  TString pathLR = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";

  if (start==""){
    cout << "Script is for testing functional mass fitting binned in different";
    cout << " physics kinematics";
    cout << "\nUsage:" << endl;
    cout << "root \'JPsiBinned(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Data coming from:            " << pathLR << endl;
    cout << "Data file is:  " << fname << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "physBinned in dy variable:   " << physBinned << endl;
    cout << "Left/Right integration range:  " << fitMmin << " - " << fitMmax;
    cout << "\nWhich fits, fitMmin,  fitMmax:" << endl;
    for (Int_t i=0; i<nBins; i++) {
      cout << whichFit[i] << " " << fitMmin[i] << " " << fitMmax[i] << endl;;
    }
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  //Get and check for data file
  TFile *f_LR =TFile::Open(pathLR+fname);
  if (!f_LR ){
    cout << "LR file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get/fit histogram, get fit parameters, get chi2 values
  gStyle->SetOptStat(1110);  gStyle->SetOptFit(1111);
  TCanvas* cMuMu = new TCanvas(); cMuMu->Divide(1, nBins);
  TCanvas* cRatio = new TCanvas(); cRatio->Divide(1, nBins);
  TH1D *h_M[nBins], *hRatio_M[nBins];
  Double_t c_LR[nBins], e_LR[nBins];
  Double_t pars[nBins][2], e_pars[nBins][2];
  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 20, 0, 10);
  TGraphErrors *g_asym =
    (TGraphErrors*)f_LR->Get(Form("%s_asym", physBinned.Data()));
  Double_t *xvals = g_asym->GetX();
  for (Int_t bi=0; bi<nBins; bi++) {
    h_M[bi] =(TH1D*)f_LR->Get(Form("MuMu_left_upstream_up_%s%i",
				   physBinned.Data(), bi));
    SetUpHist(h_M[bi]);
    hRatio_M[bi] = (TH1D*)h_M[bi]->Clone();
    
    cMuMu->cd(bi+1); gPad->SetLogy();
    h_M[bi]->Draw();
  
    Double_t redChi =FitGetPars(h_M, hRatio_M, bi, c_LR, e_LR, LR_Mmin, LR_Mmax,
				whichFit[bi], fitMmin[bi], fitMmax[bi],
				pars, e_pars);
    hChi2->Fill(redChi);
    cout << "redChi:  " << redChi << "  " << physBinned<<"  "<<xvals[bi]<<"\n";

    cRatio->cd(bi+1);
    hRatio_M[bi]->Draw("e1");
    hRatio_M[bi]->GetYaxis()->SetRangeUser(0.8, 1.2);
    hRatio_M[bi]->GetXaxis()->SetRangeUser(fitMmin[bi], fitMmax[bi]);
    DrawLine(hRatio_M[bi], 1.0);
  }
  
  TCanvas* cChi2 = new TCanvas();
  hChi2->Draw();

  Double_t ex[nBins] = {0.0};
  
  TGraphErrors* g_LR = new TGraphErrors(nBins, xvals, c_LR, ex, e_LR);
  SetUp(g_LR);
  TCanvas* cLR = new TCanvas();
  g_LR->Draw("AP"); g_LR->SetTitle("Left/Right counts");
  g_LR->Fit("pol0");

  Double_t JPsi_M[nBins], e_JPsi_M[nBins];
  Double_t JPsi_W[nBins], e_JPsi_W[nBins];
  for (Int_t bi=0; bi<nBins; bi++) {
    JPsi_M[bi] = pars[bi][0];
    e_JPsi_M[bi] = e_pars[bi][0];

    JPsi_W[bi] = pars[bi][1];
    e_JPsi_W[bi] = e_pars[bi][1];
  }

  TGraphErrors* g_M = new TGraphErrors(nBins, xvals, JPsi_M, ex, e_JPsi_M);
  SetUp(g_M);  g_M->SetTitle("JPsi Mass");
  TGraphErrors* g_W = new TGraphErrors(nBins, xvals, JPsi_W, ex, e_JPsi_W);
  SetUp(g_W); g_W->SetTitle("JPsi Width");
  TCanvas* cPars = new TCanvas(); cPars->Divide(2);
  cPars->cd(1); g_M->Draw("AP");
  cPars->cd(2); g_W->Draw("AP");

  TString fOutput = Form("FitMass_%s_%.2f_%.2f_",
			 physBinned.Data(), fitMmin[0], fitMmax[0]);
  /*fOutput += "corr.root";
  else fOutput += "noCorr.root";
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_AN_JPsi_upS_up->Write("JPsi_upS_up");

    g_JPsi_Left_upS_up->Write("JPsi_Left_upS_up");

    g_JPsi_Right_upS_up->Write("JPsi_Right_upS_up");
    }*/

  cout << "\nSettings !!!!" << endl;
  cout << "Data coming from:            " << pathLR << endl;
  cout << "Data file is:  " << fname << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "physBinned in dy variable:   " << physBinned << endl;
  cout << "\nWhich fits, fitMmin,  fitMmax:" << endl;
  for (Int_t i=0; i<nBins; i++) {
    cout << whichFit[i] << " " << fitMmin[i] << " " << fitMmax[i] << endl;;
  }
  /*if (toWrite){
    cout << "\nFile:  " << fOutput << "   was written" << endl;
  }
  else cout << "\nFile: " << fOutput << " was NOT written" << endl;*/
}
