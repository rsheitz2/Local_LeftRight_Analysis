#include "include/helperFunctions.h"

Double_t Mmin, Mmax;

#include "include/Fits_1_85/Fit_one.h"
#include "include/Fits_1_85/Fit_two.h"
#include "include/Fits_1_85/Fit_three.h"
#include "include/Fits_1_85/Fit_five.h"
#include "include/Fits_1_85/Fit_six.h"

#include "include/Fits_25_85/Fit_25_85_six.h"
#include "include/Fits_25_85/Fit_25_85_seven.h"

void FitOutput(TH1D *h){
  h->Sumw2();
  //TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  TFitResultPtr status = h->Fit("fitFunc", "RLS", "", Mmin, Mmax);
  if (status->Status() ){
    cout << "Fit failed!!" << endl;
    cout << h->GetTitle() << endl;
      
    //exit(EXIT_FAILURE);
  }
}


void DoFit(TH1D *h, Int_t iFit, TH1D *h_r,
	   Double_t *nicePars, Double_t *eNicePars){
  Mmin = h->GetXaxis()->GetXmin();
  Mmax = h->GetXaxis()->GetXmax();
  if (iFit == 256 || 257) Mmin = 2.5;
  Int_t nhbins = h->GetNbinsX();
  Int_t nBinsPhys =1;

  Bool_t hIsUpS =false;
  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0)
    hIsUpS=true;

  
  TF1* fitFunc;
  TGraphErrors *gError = new TGraphErrors(nhbins );
  if (iFit==1){
    Int_t nPar =10;
    fitFunc = new TF1("fitFunc", Fit_one, Mmin, Mmax, nPar);
    Paras_one(fitFunc, nPar);
    FitOutput(h);
    
    Double_t *pars = fitFunc->GetParameters();
    Double_t ratioPsi = Get_1_85_one_Ratio();
    Draw_2Gauss_3Expo(pars[0], pars[1], pars[2],
		      pars[3], pars[1]*ratioPsi, pars[2]*ratioPsi,
		      pars[4], pars[5],
		      pars[6], pars[7],
		      pars[8], pars[9],
		      Mmin, Mmax);
    NicePars_one(h, fitFunc, nicePars, eNicePars);
  }
  else if (iFit==2){
    Int_t nPar =10;
    if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_two_upS, Mmin, Mmax, nPar);
    else fitFunc = new TF1("fitFunc", Fit_two_downS, Mmin, Mmax, nPar);
    Paras_two(fitFunc, nPar);
    FitOutput(h);

    Double_t *pars = fitFunc->GetParameters();
    if (hIsUpS){
      Double_t ratioPsi_upS = Get_1_85_two_Ratio("UpS");
      Draw_2Gauss_3Expo(pars[0], pars[1], pars[2],
			pars[3], pars[1]*ratioPsi_upS, pars[2]*ratioPsi_upS,
			pars[4], pars[5],
			pars[6], pars[7],
			pars[8], pars[9],
			Mmin, Mmax);
    }
    else{
      Double_t ratioPsi_downS = Get_1_85_two_Ratio("DownS");
      Draw_2Gauss_3Expo(pars[0], pars[1], pars[2],
			pars[3], pars[1]*ratioPsi_downS, pars[2]*ratioPsi_downS,
			pars[4], pars[5],
			pars[6], pars[7],
			pars[8], pars[9],
			Mmin, Mmax);
    }
    NicePars_two(h, fitFunc, nicePars, eNicePars);
  }
  else if (iFit==3){
    Int_t nPar =12;
    fitFunc = new TF1("fitFunc", Fit_three, Mmin, Mmax, nPar);
    Paras_three(fitFunc, nPar);
    FitOutput(h);

    Double_t *pars = fitFunc->GetParameters();
    Draw_2Gauss_3Expo(pars, Mmin, Mmax);
    NicePars_three(h, fitFunc, nicePars, eNicePars);
  }
  else if (iFit==5){
    Int_t nPar =8;
    if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_five_upS, Mmin, Mmax, nPar);
    else fitFunc = new TF1("fitFunc", Fit_five_downS, Mmin, Mmax, nPar);
    Paras_five(fitFunc, nPar);
    FitOutput(h);

    Double_t *pars = fitFunc->GetParameters();
    if (hIsUpS){
      Double_t ratioPsi_upS_2exp = Get_1_85_five_Ratio("UpS");
      Draw_2Gauss_2Expo(pars[0], pars[1], pars[2],
			pars[3], pars[1]*ratioPsi_upS_2exp,
			pars[2]*ratioPsi_upS_2exp,
			pars[4], pars[5],
			pars[6], pars[7],
			Mmin, Mmax);
    }
    else{
      Double_t ratioPsi_downS_2exp = Get_1_85_five_Ratio("DownS");
      Draw_2Gauss_2Expo(pars[0], pars[1], pars[2],
			pars[3], pars[1]*ratioPsi_downS_2exp,
			pars[2]*ratioPsi_downS_2exp,
			pars[4], pars[5],
			pars[6], pars[7],
			Mmin, Mmax);
    }
    NicePars_five(h, fitFunc, nicePars, eNicePars);
  }
  else if (iFit==6){
    Int_t nPar =10;
    fitFunc = new TF1("fitFunc", Fit_six, Mmin, Mmax, nPar);
    Paras_six(fitFunc, nPar);
    FitOutput(h);

    Double_t *pars = fitFunc->GetParameters();
    Draw_2Gauss_2Expo(pars, Mmin, Mmax);
    NicePars_six(h, fitFunc, nicePars, eNicePars);
  }
  else if (iFit==256){
    Int_t nPar =8;
    if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_25_85_six_UpS, Mmin,Mmax,nPar);
    else fitFunc = new TF1("fitFunc", Fit_25_85_six_DownS, Mmin, Mmax, nPar);
    
    Paras_25_85_six(fitFunc, nPar, nBinsPhys);
    FitOutput(h);

    Double_t *pars = fitFunc->GetParameters();
    if (hIsUpS){
      Double_t ratioPsi_upS_2exp = Get_25_85_six_Ratio("UpS");
      Draw_2Gauss_2Expo(pars[0], pars[1], pars[2],
			pars[3], pars[1]*ratioPsi_upS_2exp,
			pars[2]*ratioPsi_upS_2exp,
			pars[4], pars[5],
			pars[6], pars[7],
			Mmin, Mmax);
    }
    else{
      Double_t ratioPsi_downS_2exp = Get_25_85_six_Ratio("DownS");
      Draw_2Gauss_2Expo(pars[0], pars[1], pars[2],
			pars[3], pars[1]*ratioPsi_downS_2exp,
			pars[2]*ratioPsi_downS_2exp,
			pars[4], pars[5],
			pars[6], pars[7],
			Mmin, Mmax);
    }
    NicePars_25_85_six(h, fitFunc, nicePars, eNicePars);
  }
  else if (iFit==257){
    Int_t nPar =10;
    fitFunc = new TF1("fitFunc", Fit_25_85_seven, Mmin,Mmax,nPar);
        
    Paras_25_85_seven(fitFunc, nPar, nBinsPhys);
    FitOutput(h);

    Double_t *pars = fitFunc->GetParameters();
    Double_t ratioPsi = Get_25_85_seven_Ratio();
    Draw_2Crystal_2Expo(pars[0], pars[1], pars[2],//JPsi
			pars[3], pars[4],
			pars[5], pars[1]*ratioPsi, pars[2]*ratioPsi,//psi
			pars[3], pars[4],
			pars[6], pars[7],//Bg
			pars[8], pars[9],//DY
			Mmin, Mmax);
    
    NicePars_25_85_seven(h, fitFunc, nicePars, eNicePars);
  }
  else{
    cout << "Incorrect iFit value:   " << iFit << " given to DoFit 1\n";
    exit(EXIT_FAILURE);
  }

  
  //FillResult
  for (Int_t i=1; i<nhbins+1; i++) 
    gError->SetPoint(i, h->GetBinCenter(i), 0);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gError);
  
  for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
    
    Double_t numerator = h_r->GetBinContent(bi);
    Double_t num_error = h_r->GetBinError(bi);
    
    Double_t fitValue = gError->Eval(h->GetBinCenter(bi) );
    //Double_t fitValue = fitFunc->Eval(h->GetBinContent(bi) );
    Double_t dem_error = gError->GetErrorY(bi);
    
    h_r->SetBinContent(bi, numerator/fitValue);
    h_r->SetBinError(bi,RatioError(numerator,fitValue,num_error,dem_error));
  }

  //Chi2
  Int_t ndf = fitFunc->GetNDF();
  Double_t chi2 = fitFunc->GetChisquare();
  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0)
    cout << "Reduced chi2 upstream:  " << chi2/(1.0*ndf) << endl;
  else cout << "Reduced chi2 downstream:  " << chi2/(1.0*ndf) << endl;
}


void fitMuMu(Int_t iFit =1, Bool_t toWrite =false, TString start=""){
  //Setup_____
  TString inputPath = "InputData/";
  TString inputData = "MuMu_1_85_notScaled.root";
  //Setup_____
  
  if (start==""){
    cout << "Script gives the Chi2, fit to data ratios,";
    cout << " and import physical fit parameters" << endl;
    cout << "This script is specifically for mass range of 1.0-8.5 GeV" << endl;
    cout <<"\nUsage:" << endl;
    cout << "root \'Mass1_85_fitMuMu(Int_t iFit=1,Bool_t toWrite=false, 1)\'\n";
    cout << "\nCurrent settings:" << endl;
    cout << "Input data path:   " << inputPath << endl;
    exit(EXIT_FAILURE);
  }
  
  TFile *f = TFile::Open(inputPath+inputData);
  if (!f){
    cout << "File   " << inputData << "  does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  TH1D *h_mumu_up = (TH1D*)f->Get("h_mumu_up");
  TH1D *h_mumu_down = (TH1D*)f->Get("h_mumu_down");

  TH1D *r_mumu_up = (TH1D*)h_mumu_up->Clone();
  TH1D *r_mumu_down = (TH1D*)h_mumu_down->Clone();

  Double_t nicePars[8] = {0.0};
  Double_t eNicePars[8] = {0.0};
  
  TCanvas* cFit = new TCanvas(); cFit->Divide(2);
  cFit->cd(1);
  gPad->SetLogy();
  DoFit(h_mumu_up, iFit, r_mumu_up, nicePars, eNicePars);

  cFit->cd(2);
  gPad->SetLogy();
  DoFit(h_mumu_down, iFit, r_mumu_down, nicePars, eNicePars);

  
  TCanvas* cRatio = new TCanvas(); 
  r_mumu_up->Draw();
  r_mumu_down->Draw("same");
  DrawLine(r_mumu_up, 1.0);

  //Useful Parameters
  Double_t M_JPsi_up =nicePars[0], M_JPsi_down =nicePars[1];
  Double_t W_JPsi_up =nicePars[2], W_JPsi_down =nicePars[3];
  Double_t M_psi_up =nicePars[4], M_psi_down =nicePars[5];
  Double_t W_psi_up =nicePars[6], W_psi_down =nicePars[7];

  Double_t eM_JPsi_up =eNicePars[0], eM_JPsi_down =eNicePars[1];
  Double_t eW_JPsi_up =eNicePars[2], eW_JPsi_down =eNicePars[3];
  Double_t eM_psi_up =eNicePars[4], eM_psi_down =eNicePars[5];
  Double_t eW_psi_up =eNicePars[6], eW_psi_down =eNicePars[7];

  Double_t min_y = 0.0;
  Double_t max_y = 2.5;
  DrawLineY(r_mumu_up, M_JPsi_up, min_y, max_y, 4);
  DrawLineY(r_mumu_up, M_psi_up, min_y, max_y, 4);
  DrawLineY(r_mumu_up, M_JPsi_down, min_y, max_y, 2);
  DrawLineY(r_mumu_up, M_psi_down, min_y, max_y, 2);


  Double_t xval_up[] = {0.9 + iFit/10.0}, ex[] = {0.0};
  Double_t xval_down[] = {0.9 + iFit/10.0+0.02};
  TGraphErrors* g_MJPsi_UpS
    = new TGraphErrors(1, xval_up, &M_JPsi_up, ex, &eM_JPsi_up);
  TGraphErrors* g_MJPsi_DownS
    = new TGraphErrors(1, xval_down, &M_JPsi_down, ex, &eM_JPsi_down);
  TGraphErrors* g_WJPsi_UpS
    = new TGraphErrors(1, xval_up, &W_JPsi_up, ex, &eW_JPsi_up);
  TGraphErrors* g_WJPsi_DownS
    = new TGraphErrors(1, xval_down, &W_JPsi_down, ex, &eW_JPsi_down);

  TGraphErrors* g_Mpsi_UpS
    = new TGraphErrors(1, xval_up, &M_psi_up, ex, &eM_psi_up);
  TGraphErrors* g_Mpsi_DownS
    = new TGraphErrors(1, xval_down, &M_psi_down, ex, &eM_psi_down);
  TGraphErrors* g_Wpsi_UpS
    = new TGraphErrors(1, xval_up, &W_psi_up, ex, &eW_psi_up);
  TGraphErrors* g_Wpsi_DownS
    = new TGraphErrors(1, xval_down, &W_psi_down, ex, &eW_psi_down);
  
  SetUp(g_MJPsi_UpS, 20, 1); SetUp(g_MJPsi_DownS, 21, 2);
  SetUp(g_WJPsi_UpS, 20, 1); SetUp(g_WJPsi_DownS, 21, 2);

  SetUp(g_Mpsi_UpS, 20, 1); SetUp(g_Mpsi_DownS, 21, 2);
  SetUp(g_Wpsi_UpS, 20, 1); SetUp(g_Wpsi_DownS, 21, 2);
  
  TCanvas* cPars = new TCanvas(); cPars->Divide(2,2);
  cPars->cd(1);
  g_MJPsi_UpS->Draw("AP"); g_MJPsi_DownS->Draw("Psame");
  g_MJPsi_UpS->SetTitle("JPsi Mass (GeV)");

  cPars->cd(2);
  g_WJPsi_UpS->Draw("AP"); g_WJPsi_DownS->Draw("Psame");
  g_WJPsi_UpS->SetTitle("JPsi Width (GeV)");

  cPars->cd(3);
  g_Mpsi_UpS->Draw("AP"); g_Mpsi_DownS->Draw("Psame");
  g_Mpsi_UpS->SetTitle("psi' Mass (GeV)");

  cPars->cd(4);
  g_Wpsi_UpS->Draw("AP"); g_Wpsi_DownS->Draw("Psame");
  g_Wpsi_UpS->SetTitle("psi' Width (GeV)");

  //End of Macro Output
  cout << " " << endl;
  cout << "Input data is:  " << inputData << endl;
  cout << "iFit:  " << iFit << endl;
  switch(iFit){
  case 1 : cout << "2 Gaussians Constrained (psi'=a*JPsi) && 3 Exponentials\n";
    break;
  case 2 : cout << "2 Gaussians Constrained by target (psi'=a*JPsi) ";
    cout << "&& 3 Exponentials\n";
    break;
  case 3 : cout << "2 Gaussians && 3 Exponentials\n";
    break;
  case 4 : cout << "2 Gaussians Constrained (psi'=a*JPsi) && 2 Exponentials\n";
    break;
  case 5 : cout << "2 Gaussians Constrained by target (psi'=a*JPsi) ";
    cout << "&& 2 Exponentials\n";
    break;
  case 6 : cout << "2 Gaussians && 2 Exponentials\n";
    break;
  case 17 : cout<<"2 Gaussians psi' M shifted psi' W*aJPsi && 3 Exponentials\n";
    break;
  case 18 : cout<<"2 Gaussians psi' M shifted psi' W*aJPsi by target";
    cout << " && 3 Exponentials\n";
    break;
  case 256 : cout << "2 Gaussians Constrained by target (psi'=a*JPsi) ";
    cout << "&& 2 Exponentials, Mass Range 2.5-8.5\n";
    break;
  case 257 : cout << "2 CrystalBalls Constrained (psi'=a*JPsi) ";
    cout << "&& 2 Exponentials, Mass Range 2.5-8.5\n";
    break;
  }
  cout << " " << endl;
 
  //Write to file
  TString outName = Form("OutputData/Fit_%i_%s", iFit, inputData.Data() );
  cout << "File   " << outName;
  if (toWrite){
    TFile *fOut = TFile::Open(outName, "RECREATE");
    
    h_mumu_up->Write();
    h_mumu_down->Write();
    r_mumu_up->Write("r_mumu_up");
    r_mumu_down->Write("r_mumu_down");

    g_MJPsi_UpS->Write("MJPsi_UpS");
    g_WJPsi_UpS->Write("WJPsi_UpS");
    g_MJPsi_DownS->Write("MJPsi_DownS");
    g_WJPsi_DownS->Write("WJPsi_DownS");

    g_Mpsi_UpS->Write("Mpsi_UpS");
    g_Wpsi_UpS->Write("Wpsi_UpS");
    g_Mpsi_DownS->Write("Mpsi_DownS");
    g_Wpsi_DownS->Write("Wpsi_DownS");

    cout << "   was written" << endl;
    fOut->Close();
  }
  else cout << "   was  NOT  written" << endl;
}
