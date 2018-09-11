#include "include/helperFunctions.h"

//Setup_______________
const Int_t nBins =5;//5, 3
const Int_t hbins =150;//100, 120, 150, 200, 300

TString physBinned = "xF";//"xN" "xPi" "xF" "pT"
Bool_t toWrite =false;
Bool_t PolCorr =true;
const Double_t Mmin=2.5, Mmax=8.5;

const Int_t whichFit =9;
Int_t Failure=0;
Int_t nPar;

//For whichFit =7 hbins =150
//Double_t ratioPsi = 1.1587;//1.62 //red Chi2 from hbins==150
//Double_t ratioPsi = 1.16;//1.63
//Double_t ratioPsi = 1.17;//1.73
//Double_t ratioPsi = 1.18;//1.89
//Double_t ratioPsi = 1.19;//2.10
//Double_t ratioPsi = 1.15;//1.58
//Double_t ratioPsi = 1.14;//1.575
//Double_t ratioPsi = 1.13;//1.62
//Double_t ratioPsi = 1.135;//1.59
Double_t ratioPsi = 1.145;//1.5699
//Double_t ratioPsi = 1.143;//1.5704
//Double_t ratioPsi = 1.146;//fail
//Double_t ratioPsi = 1.144;//fail

//For whichFit =7 hbins =150 2expos fits
//Double_t ratioPsi = 1.145;//1.70 
//Double_t ratioPsi = 1.1587;//1.70
//Double_t ratioPsi = 1.16;//
//Double_t ratioPsi = 1.17;//1.79
//Double_t ratioPsi = 1.18;
//Double_t ratioPsi = 1.19;//2.13
//Double_t ratioPsi = 1.15;
//Double_t ratioPsi = 1.14;
//Double_t ratioPsi = 1.13;//1.84
//Double_t ratioPsi = 1.135;

//Double_t ratioPsi_UpS =1.145;//1.793
Double_t ratioPsi_UpS =1.15;//1.785
//Double_t ratioPsi_UpS =1.16;//1.811
//Double_t ratioPsi_UpS =1.19;//2.237
//Double_t ratioPsi_UpS =1.155;//1.79
//Double_t ratioPsi_UpS =1.145;//1.793

//Double_t ratioPsi_DownS =1.145;//1.347
//Double_t ratioPsi_DownS =1.15;//1.37
Double_t ratioPsi_DownS =1.14;//1.335
//Double_t ratioPsi_DownS =1.13;//1.346
//Double_t ratioPsi_DownS =1.135;//1.335
//Setup_______________


//Fit results/Monitoring_______________
Int_t iter=0;

Double_t r_Chi2[nBins*8], r_NDF[nBins*8], r_RedChi2[nBins*8];
Double_t r_Xpoints[nBins*8];

TH1D *r_hRD_upS_up_left[nBins], *r_hRD_upS_up_right[nBins];
TH1D *r_hRD_upS_down_left[nBins], *r_hRD_upS_down_right[nBins];
TH1D *r_hRD_downS_up_left[nBins], *r_hRD_downS_up_right[nBins];
TH1D *r_hRD_downS_down_left[nBins], *r_hRD_downS_down_right[nBins];
//Fit results/Monitoring_______________

Double_t FitMCs(Double_t *x, Double_t *par){//whichFit =1
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[5] )/par[6];
  Double_t JPsi = par[4]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[8] )/par[9];
  Double_t psi = par[7]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[10]*TMath::Exp( par[11]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


Double_t Fit_1Gauss_3Expo(Double_t *x, Double_t *par){//whichFit =2
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[5] )/par[6];
  Double_t JPsi = par[4]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );

  Double_t DY = par[7]*TMath::Exp( par[8]*xShift );
    
  return CombBg + OC + JPsi + DY;
}


Double_t Fit_1Gauss_2Expo(Double_t *x, Double_t *par){//whichFit =3
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[3] )/par[4];
  Double_t JPsi = par[2]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );

  Double_t DY = par[5]*TMath::Exp( par[6]*xShift );
    
  return CombBg + JPsi + DY;
}


Double_t Fit_2Gauss_2Expo(Double_t *x, Double_t *par){//whichFit =4
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[3] )/par[4];
  Double_t JPsi = par[2]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[6] )/par[7];
  Double_t psi = par[5]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + DY;
}


Double_t Fit_2Wigner_3Expo(Double_t *x, Double_t *par){//whichFit =5
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );

  Double_t JPsi = par[4]*TMath::BreitWigner(x[0], par[5], par[6]);
  Double_t psi = par[7]*TMath::BreitWigner(x[0], par[8], par[9]);

  Double_t DY = par[10]*TMath::Exp( par[11]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


Double_t Fit_1Wigner_3Expo(Double_t *x, Double_t *par){//whichFit =6
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );

  Double_t JPsi = par[4]*TMath::BreitWigner(x[0], par[5], par[6]);
  
  Double_t DY = par[7]*TMath::Exp( par[8]*xShift );
    
  return CombBg + OC + JPsi + DY;
}


Double_t Fit_2GaussConstrain_3Expo(Double_t *x, Double_t *par){//whichFit =7
  
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[5] )/par[6];
  Double_t JPsi = par[4]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[5]*ratioPsi )/(par[6]*ratioPsi);
  Double_t psi = par[7]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


Double_t Fit_2GaussConstrain_2Expo(Double_t *x, Double_t *par){//whichFit =8
  
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
    
  Double_t arg_JPsi = ( x[0] - par[3] )/par[4];
  Double_t JPsi = par[2]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[3]*ratioPsi )/(par[4]*ratioPsi);
  Double_t psi = par[5]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[6]*TMath::Exp( par[7]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_2GaussConstrainUpS_3Expo(Double_t *x, Double_t *par){//whichFit =9
  
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[5] )/par[6];
  Double_t JPsi = par[4]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[5]*ratioPsi_UpS )/(par[6]*ratioPsi_UpS);
  Double_t psi = par[7]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


Double_t Fit_2GaussConstrainDownS_3Expo(Double_t *x, Double_t *par){//whichFit =9
  
  Double_t xShift = x[0]-Mmin;

  Double_t CombBg = par[0]*TMath::Exp( par[1]*xShift );
  Double_t OC = par[2]*TMath::Exp( par[3]*xShift );
  
  Double_t arg_JPsi = ( x[0] - par[5] )/par[6];
  Double_t JPsi = par[4]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[5]*ratioPsi_DownS )/(par[6]*ratioPsi_DownS);
  Double_t psi = par[7]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


void FillResult(TH1D* h, TH1D* hResult, Double_t *pars, Int_t iFit){
  for (Int_t bi=0; bi<h->GetNbinsX(); bi++) {
    Double_t numerator =hResult->GetBinContent(bi);

    Double_t fitValue;
    Double_t x[] = { h->GetBinCenter(bi) };
    if (iFit==1) fitValue =FitMCs(x, pars);
    else if (iFit==2) fitValue =Fit_1Gauss_3Expo(x, pars);
    else if (iFit==3) fitValue =Fit_1Gauss_2Expo(x, pars);
    else if (iFit==4) fitValue =Fit_2Gauss_2Expo(x, pars);
    else if (iFit==5) fitValue =Fit_2Wigner_3Expo(x, pars);
    else if (iFit==6) fitValue =Fit_1Wigner_3Expo(x, pars);
    else if (iFit==7) fitValue =Fit_2GaussConstrain_3Expo(x, pars);
    else if (iFit==8) fitValue =Fit_2GaussConstrain_2Expo(x, pars);
    else if (iFit==9) {
      if (strncmp(Form("%s", h->GetTitle() ), "MuMu_left_upstream", 18) == 0 ){
	fitValue =Fit_2GaussConstrainUpS_3Expo(x, pars);
      }
      else if (strncmp(Form("%s", h->GetTitle() ), "MuMu_right_upstream", 19)==0){
	fitValue =Fit_2GaussConstrainUpS_3Expo(x, pars);
      }
      else if (strncmp(Form("%s", h->GetTitle() ), "MuMu_left_downstream",20)==0){
	fitValue =Fit_2GaussConstrainDownS_3Expo(x, pars);
      }
      else if (strncmp(Form("%s", h->GetTitle() ),"MuMu_right_downstream",21)==0){
	fitValue =Fit_2GaussConstrainDownS_3Expo(x, pars);
      }
      else {
	cout << "Wrong hist name " << endl;
	exit(EXIT_FAILURE);
      }
    }
    else {
      cout << "Incorrect iFit value:   " << iFit << " given to FillResult\n";
      exit(EXIT_FAILURE);
    }

    hResult->SetBinContent(bi, numerator/fitValue);
    hResult->SetBinError(bi, RatioError(numerator, fitValue,
					TMath::Sqrt(numerator),
					TMath::Sqrt(fitValue) ) );
  }
}


Double_t FitGetPars(TH1D* h, Int_t bin, vector<Double_t> &counts,
		    vector<Double_t> &e_counts, Int_t iFit){
  if (whichFit==1) nPar =12;
  else if (whichFit==2) nPar =9;
  else if (whichFit==3) nPar =7;
  else if (whichFit==4) nPar =10;
  else if (whichFit==5) nPar =12;
  else if (whichFit==6) nPar =9;
  else if (whichFit==7) nPar =10;
  else if (whichFit==8) nPar =8;
  else if (whichFit==9) nPar =10;
  else cout << "Incorrect iFit value:  " << iFit << "   to  0  FitGetPars\n";
  
  TF1 *fitFunc;

  if (iFit==1){
    fitFunc = new TF1("fitFunc", FitMCs, Mmin, Mmax, nPar);
    fitFunc->SetParameter(0, 8.5e2); fitFunc->SetParameter(1, -3.5);//CombBg
    fitFunc->SetParameter(2, 8.5e2); fitFunc->SetParameter(3, -2.5);//OC
    fitFunc->SetParameter(4, 6.2e3); fitFunc->SetParameter(5, 3.1);//JPsi
    fitFunc->SetParameter(6, 0.15);
    fitFunc->SetParameter(7, 1e2); fitFunc->SetParameter(8, 3.6);//psi
    fitFunc->SetParameter(9, 0.15);
    fitFunc->SetParameter(10, 2e2); fitFunc->SetParameter(11,-0.9);//DY

    fitFunc->SetParLimits(0, 100, 850); fitFunc->SetParLimits(1, -4.0, -3.0);//Bg
    fitFunc->SetParLimits(2, 100, 850); fitFunc->SetParLimits(3, -3.0, -2.0);//OC
    fitFunc->SetParLimits(4, 1e3, 2e5); fitFunc->SetParLimits(5, 2.5, 3.6);//JPsi
    fitFunc->SetParLimits(6, 0, 1.0);
    fitFunc->SetParLimits(7, 10, 2e5); //psi
    fitFunc->SetParLimits(8, 3.5, 4.1); fitFunc->SetParLimits(9, 0, 0.25);
    fitFunc->SetParLimits(10, 100, 3e2); fitFunc->SetParLimits(11, -1.0, 0);//DY
  }
  else if (iFit==2){
    fitFunc = new TF1("fitFunc", Fit_1Gauss_3Expo, Mmin, Mmax, nPar);
    fitFunc->SetParameter(0, 8.5e2); fitFunc->SetParameter(1, -3.5);//CombBg
    fitFunc->SetParameter(2, 8.5e2); fitFunc->SetParameter(3, -2.5);//OC
    fitFunc->SetParameter(4, 6.2e3); fitFunc->SetParameter(5, 3.1);//JPsi
    fitFunc->SetParameter(6, 0.15);
    fitFunc->SetParameter(7, 2e2); fitFunc->SetParameter(8,-0.9);//DY

    fitFunc->SetParLimits(0, 100, 850); fitFunc->SetParLimits(1, -4.0, -3.0);//Bg
    fitFunc->SetParLimits(2, 100, 850); fitFunc->SetParLimits(3, -3.0, -2.0);//OC
    fitFunc->SetParLimits(4, 1e3, 2e5); fitFunc->SetParLimits(5, 2.5, 3.6);//JPsi
    fitFunc->SetParLimits(6, 0, 1.0);
    fitFunc->SetParLimits(7, 100, 3e2); fitFunc->SetParLimits(8, -1.0, 0);//DY
  }
  else if (iFit==3){
    fitFunc = new TF1("fitFunc", Fit_1Gauss_2Expo, Mmin, Mmax, nPar);
    fitFunc->SetParameter(0, 8.5e2); fitFunc->SetParameter(1, -3.5);//CombBg
    fitFunc->SetParameter(2, 6.2e3); fitFunc->SetParameter(3, 3.1);//JPsi
    fitFunc->SetParameter(4, 0.15);
    fitFunc->SetParameter(5, 2e2); fitFunc->SetParameter(6,-0.9);//DY

    fitFunc->SetParLimits(0, 100, 850); fitFunc->SetParLimits(1, -4.0, -3.0);//Bg
    fitFunc->SetParLimits(2, 1e3, 2e5); fitFunc->SetParLimits(3, 2.5, 3.6);//JPsi
    fitFunc->SetParLimits(4, 0, 1.0);
    fitFunc->SetParLimits(5, 100, 3e2); fitFunc->SetParLimits(6, -1.0, 0);//DY
  }
  else if (iFit==4){
    fitFunc = new TF1("fitFunc", Fit_2Gauss_2Expo, Mmin, Mmax, nPar);
    fitFunc->SetParameter(0, 8.5e2); fitFunc->SetParameter(1, -3.5);//CombBg
    fitFunc->SetParameter(2, 6.2e3); fitFunc->SetParameter(3, 3.1);//JPsi
    fitFunc->SetParameter(4, 0.15);
    fitFunc->SetParameter(5, 1e2); fitFunc->SetParameter(6, 3.6);//psi
    fitFunc->SetParameter(7, 0.15);
    fitFunc->SetParameter(8, 2e2); fitFunc->SetParameter(9,-0.9);//DY

    fitFunc->SetParLimits(0, 100, 850); fitFunc->SetParLimits(1, -4.0, -3.0);//Bg
    fitFunc->SetParLimits(2, 1e3, 2e5); fitFunc->SetParLimits(3, 2.5, 3.6);//JPsi
    fitFunc->SetParLimits(4, 0, 1.0);
    fitFunc->SetParLimits(5, 10, 2e5); //psi
    fitFunc->SetParLimits(6, 3.5, 4.1); fitFunc->SetParLimits(7, 0, 0.25);
    fitFunc->SetParLimits(8, 100, 3e2); fitFunc->SetParLimits(9, -1.0, 0);//DY
  }
  else if (iFit==5){
    fitFunc = new TF1("fitFunc", Fit_2Wigner_3Expo, Mmin, Mmax, nPar);
    fitFunc->SetParameter(0, 8.5e2); fitFunc->SetParameter(1, -3.5);//CombBg
    fitFunc->SetParameter(2, 8.5e2); fitFunc->SetParameter(3, -2.5);//OC
    fitFunc->SetParameter(4, 6.2e3); fitFunc->SetParameter(5, 3.1);//JPsi
    fitFunc->SetParameter(6, 0.15);
    fitFunc->SetParameter(7, 1e2); fitFunc->SetParameter(8, 3.6);//psi
    fitFunc->SetParameter(9, 0.15);
    fitFunc->SetParameter(10, 2e2); fitFunc->SetParameter(11,-0.9);//DY

    fitFunc->SetParLimits(0, 10, 850); fitFunc->SetParLimits(1, -4.0, -3.0);//Bg
    fitFunc->SetParLimits(2, 10, 850); fitFunc->SetParLimits(3, -3.0, -2.0);//OC
    fitFunc->SetParLimits(4, 10, 2e5); fitFunc->SetParLimits(5, 2.5, 3.6);//JPsi
    fitFunc->SetParLimits(6, 0, 0.2);
    fitFunc->SetParLimits(7, 10, 2e5); //psi
    fitFunc->SetParLimits(8, 3.5, 4.1); fitFunc->SetParLimits(9, 0.1, 1.0);
    fitFunc->SetParLimits(10, 100, 3e2); fitFunc->SetParLimits(11, -1.0, 0);//DY
  }
  else if (iFit==6){
    fitFunc = new TF1("fitFunc", Fit_1Wigner_3Expo, Mmin, Mmax, nPar);
    
    Int_t ipar =0;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//CombBg
    fitFunc->SetParameter(ipar, -3.5); ipar++;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//OC
    fitFunc->SetParameter(ipar, -2.5); ipar++;
    fitFunc->SetParameter(ipar, 6.2e3); ipar++;//JPsi
    fitFunc->SetParameter(ipar, 3.1); ipar++;
    fitFunc->SetParameter(ipar, 0.15); ipar++;
    fitFunc->SetParameter(ipar, 2e2); ipar++;//DY
    fitFunc->SetParameter(ipar, -0.9); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }

    ipar=0;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//Bg
    fitFunc->SetParLimits(ipar, -4.0, -3.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//OC
    fitFunc->SetParLimits(ipar, -3.0, -2.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//JPsi
    fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
    fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
    fitFunc->SetParLimits(ipar, 100, 3e2); ipar++;//DY
    fitFunc->SetParLimits(ipar, -1.0, 0); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else if (iFit==7){
    fitFunc = new TF1("fitFunc", Fit_2GaussConstrain_3Expo, Mmin, Mmax, nPar);
    
    Int_t ipar =0;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//CombBg
    fitFunc->SetParameter(ipar, -3.5); ipar++;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//OC
    fitFunc->SetParameter(ipar, -2.5); ipar++;
    fitFunc->SetParameter(ipar, 6.2e3); ipar++;//JPsi
    fitFunc->SetParameter(ipar, 3.1); ipar++;
    fitFunc->SetParameter(ipar, 0.15); ipar++;
    fitFunc->SetParameter(ipar, 1e2); ipar++;//psi
    fitFunc->SetParameter(ipar, 2e2); ipar++;//DY
    fitFunc->SetParameter(ipar, -0.9); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }

    ipar=0;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//Bg
    fitFunc->SetParLimits(ipar, -4.0, -3.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//OC
    fitFunc->SetParLimits(ipar, -3.0, -2.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//JPsi
    fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
    fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//psi
    fitFunc->SetParLimits(ipar, 100, 3e2); ipar++;//DY
    fitFunc->SetParLimits(ipar, -1.0, 0); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else if (iFit==8){
    fitFunc = new TF1("fitFunc", Fit_2GaussConstrain_2Expo, Mmin, Mmax, nPar);
    
    Int_t ipar =0;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//CombBg
    fitFunc->SetParameter(ipar, -3.5); ipar++;
    fitFunc->SetParameter(ipar, 6.2e3); ipar++;//JPsi
    fitFunc->SetParameter(ipar, 3.1); ipar++;
    fitFunc->SetParameter(ipar, 0.15); ipar++;
    fitFunc->SetParameter(ipar, 1e2); ipar++;//psi
    fitFunc->SetParameter(ipar, 2e2); ipar++;//DY
    fitFunc->SetParameter(ipar, -0.9); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }

    ipar=0;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//Bg
    fitFunc->SetParLimits(ipar, -4.0, -3.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//JPsi
    fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
    fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//psi
    fitFunc->SetParLimits(ipar, 100, 3e2); ipar++;//DY
    fitFunc->SetParLimits(ipar, -1.0, 0); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else if (iFit==9){
    if (strncmp(Form("%s", h->GetTitle() ), "MuMu_left_upstream", 18) == 0 ){
      fitFunc = new TF1("fitFunc", Fit_2GaussConstrainUpS_3Expo,Mmin,Mmax,nPar);
    }
    else if (strncmp(Form("%s", h->GetTitle() ), "MuMu_right_upstream", 19)==0){
      fitFunc = new TF1("fitFunc", Fit_2GaussConstrainUpS_3Expo,Mmin,Mmax,nPar);
    }
    else if (strncmp(Form("%s", h->GetTitle() ), "MuMu_left_downstream",20)==0){
      fitFunc =new TF1("fitFunc",Fit_2GaussConstrainDownS_3Expo,Mmin,Mmax,nPar);
    }
    else if (strncmp(Form("%s", h->GetTitle() ),"MuMu_right_downstream",21)==0){
      fitFunc =new TF1("fitFunc",Fit_2GaussConstrainDownS_3Expo,Mmin,Mmax,nPar);
    }
    else {
      cout << "Wrong hist name " << endl;
      exit(EXIT_FAILURE);
    }
    
    Int_t ipar =0;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//CombBg
    fitFunc->SetParameter(ipar, -3.5); ipar++;
    fitFunc->SetParameter(ipar, 8.5e2); ipar++;//OC
    fitFunc->SetParameter(ipar, -2.5); ipar++;
    fitFunc->SetParameter(ipar, 6.2e3); ipar++;//JPsi
    fitFunc->SetParameter(ipar, 3.1); ipar++;
    fitFunc->SetParameter(ipar, 0.15); ipar++;
    fitFunc->SetParameter(ipar, 1e2); ipar++;//psi
    fitFunc->SetParameter(ipar, 2e2); ipar++;//DY
    fitFunc->SetParameter(ipar, -0.9); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }

    ipar=0;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//Bg
    fitFunc->SetParLimits(ipar, -4.0, -3.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 850); ipar++;//OC
    fitFunc->SetParLimits(ipar, -3.0, -2.0); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//JPsi
    fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
    fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
    fitFunc->SetParLimits(ipar, 10, 2e5); ipar++;//psi
    fitFunc->SetParLimits(ipar, 100, 3e2); ipar++;//DY
    fitFunc->SetParLimits(ipar, -1.0, 0); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else cout << "Incorrect iFit value:  " << iFit << "   to  1  FitGetPars\n";

  
  h->Sumw2();
  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ){
    cout << "Fit failed!!" << endl;
    cout << h->GetTitle() << "   bin: " << bin << endl;
    Failure++;
    //exit(EXIT_FAILURE);
  }
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *ePars = fitFunc->GetParErrors();
  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();

  
  if (iFit==1){
    Double_t JPsi = GaussInt(pars[4], pars[6]);
    Double_t psi = GaussInt(pars[7], pars[9]);
    Double_t OC = ExpoInt(pars[2], pars[3]);
    Double_t AMDY = ExpoInt(pars[10], pars[11]);
    counts.push_back(JPsi);
    counts.push_back(psi);
    counts.push_back(OC);
    counts.push_back(AMDY);

    Double_t e_JPsi = GaussIntError(pars[4], pars[6], ePars[4], ePars[6] );
    Double_t e_psi = GaussIntError(pars[7], pars[9], ePars[7], ePars[9] );
    Double_t e_OC = RatioError(pars[2], pars[3], ePars[2], ePars[3] );
    Double_t e_AMDY = RatioError(pars[10], pars[11], ePars[10], ePars[11] );
    e_counts.push_back(e_JPsi);
    e_counts.push_back(e_psi);
    e_counts.push_back(e_OC);
    e_counts.push_back(e_AMDY);

    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] );
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_OC->SetParameters(pars[2], pars[3] );
    f_OC->SetLineColor(6); f_OC->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[4], pars[5], pars[6]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_psi =
      new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_psi->SetParameters(pars[7], pars[8], pars[9] );
    f_psi->SetLineColor(kGreen); f_psi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[10], pars[11]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==2){
    Double_t JPsi = GaussInt(pars[4], pars[6]);
    Double_t OC = ExpoInt(pars[2], pars[3]);
    Double_t AMDY = ExpoInt(pars[7], pars[8]);
    counts.push_back(JPsi);
    counts.push_back(OC);
    counts.push_back(AMDY);

    Double_t e_JPsi = GaussIntError(pars[4], pars[6], ePars[4], ePars[6] );
    Double_t e_OC = RatioError(pars[2], pars[3], ePars[2], ePars[3] );
    Double_t e_AMDY = RatioError(pars[7], pars[8], ePars[7], ePars[8] );
    e_counts.push_back(e_JPsi);
    e_counts.push_back(e_OC);
    e_counts.push_back(e_AMDY);

    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] );
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_OC->SetParameters(pars[2], pars[3] );
    f_OC->SetLineColor(6); f_OC->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[4], pars[5], pars[6]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[10], pars[11]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==3){
    Double_t JPsi = GaussInt(pars[2], pars[4]);
    Double_t AMDY = ExpoInt(pars[5], pars[6]);
    counts.push_back(JPsi);
    counts.push_back(AMDY);

    Double_t e_JPsi = GaussIntError(pars[2], pars[4], ePars[2], ePars[4] );
    Double_t e_AMDY = RatioError(pars[5], pars[6], ePars[5], ePars[6] );
    e_counts.push_back(e_JPsi);
    e_counts.push_back(e_AMDY);

    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] );
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[2], pars[3], pars[4]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[5], pars[6]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==4){
    Double_t JPsi = GaussInt(pars[2], pars[4]);
    Double_t psi = GaussInt(pars[5], pars[7]);
    Double_t AMDY = ExpoInt(pars[8], pars[9]);
    counts.push_back(JPsi);
    counts.push_back(psi);
    counts.push_back(AMDY);

    Double_t e_JPsi = GaussIntError(pars[2], pars[4], ePars[2], ePars[4] );
    Double_t e_psi = GaussIntError(pars[5], pars[7], ePars[5], ePars[7] );
    Double_t e_AMDY = RatioError(pars[8], pars[9], ePars[8], ePars[9] );
    e_counts.push_back(e_JPsi);
    e_counts.push_back(e_psi);
    e_counts.push_back(e_AMDY);

    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] );
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[2], pars[3], pars[4]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_psi =
      new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_psi->SetParameters(pars[5], pars[6], pars[7]);
    f_psi->SetLineColor(kGreen); f_psi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[8], pars[9]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==5){
    /*Double_t JPsi = GaussInt(pars[4], pars[6]);
    Double_t psi = GaussInt(pars[7], pars[9]);
    Double_t OC = ExpoInt(pars[2], pars[3]);
    Double_t AMDY = ExpoInt(pars[10], pars[11]);
    counts.push_back(JPsi);
    counts.push_back(psi);
    counts.push_back(OC);
    counts.push_back(AMDY);

    Double_t e_JPsi = GaussIntError(pars[4], pars[6], ePars[4], ePars[6] );
    Double_t e_psi = GaussIntError(pars[7], pars[9], ePars[7], ePars[9] );
    Double_t e_OC = RatioError(pars[2], pars[3], ePars[2], ePars[3] );
    Double_t e_AMDY = RatioError(pars[10], pars[11], ePars[10], ePars[11] );
    e_counts.push_back(e_JPsi);
    e_counts.push_back(e_psi);
    e_counts.push_back(e_OC);
    e_counts.push_back(e_AMDY);*/

    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] );
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_OC->SetParameters(pars[2], pars[3] );
    f_OC->SetLineColor(6); f_OC->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::BreitWigner(x, [1], [2])",
	      0, Mmax);
    f_JPsi->SetParameters(pars[4], pars[5], pars[6]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::BreitWigner(x, [1], [2])",
	      0, Mmax);
    f_psi->SetParameters(pars[7], pars[8], pars[9]);
    f_psi->SetLineColor(kGreen); f_psi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[10], pars[11]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==6){
    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] ); 
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_OC->SetParameters(pars[2], pars[3] );
    f_OC->SetLineColor(6); f_OC->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::BreitWigner(x, [1], [2])",
	      0, Mmax);
    f_JPsi->SetParameters(pars[4], pars[5], pars[6]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[7], pars[8]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==7){
    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] ); 
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_OC->SetParameters(pars[2], pars[3] );
    f_OC->SetLineColor(6); f_OC->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[4], pars[5], pars[6]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_psi =
      new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_psi->SetParameters(pars[7], pars[5]*ratioPsi, pars[6]*ratioPsi);
    f_psi->SetLineColor(kGreen); f_psi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[8], pars[9]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==8){
    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] ); 
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[2], pars[3], pars[4]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_psi =
      new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_psi->SetParameters(pars[5], pars[3]*ratioPsi, pars[4]*ratioPsi);
    f_psi->SetLineColor(kGreen); f_psi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[6], pars[7]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else if (iFit==9){
    TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-2.5) )", 0,Mmax);
    f_CombBg->SetParameters(pars[0], pars[1] ); 
    f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

    TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_OC->SetParameters(pars[2], pars[3] );
    f_OC->SetLineColor(6); f_OC->Draw("same");

    TF1 *f_JPsi =
      new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_JPsi->SetParameters(pars[4], pars[5], pars[6]);
    f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

    TF1 *f_psi =
      new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	      0, Mmax);
    f_psi->SetParameters(pars[7], pars[5]*ratioPsi, pars[6]*ratioPsi);
    f_psi->SetLineColor(kGreen); f_psi->Draw("same");

    TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-2.5) )", 0, Mmax);
    f_DY->SetParameters(pars[8], pars[9]);
    f_DY->SetLineColor(kBlue); f_DY->Draw("same");
  }
  else cout << "Incorrect iFit value:  " << iFit << "   to  2  FitGetPars\n";
  
  
  //Fit results
  r_Chi2[iter] =Chi2; r_NDF[iter] =ndf; r_RedChi2[iter] =Chi2/ndf;
  if (!iter) r_Xpoints[iter] = 1;
  else if ( !(iter%8) ) r_Xpoints[iter] = r_Xpoints[iter-1]+4;
  else r_Xpoints[iter] = r_Xpoints[iter-1]+1;

  if (strncmp(Form("%s", h->GetTitle() ),
	      Form("MuMu_left_upstream_up_%s%i",physBinned.Data(),bin),  21)
      == 0 ){
    FillResult(h, r_hRD_upS_up_left[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_right_upstream_up_%s%i",physBinned.Data(),bin),  22)
	   == 0 ){
    FillResult(h, r_hRD_upS_up_right[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_left_upstream_down_%s%i",physBinned.Data(),bin),  23)
	   == 0 ){
    FillResult(h, r_hRD_upS_down_left[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_right_upstream_down_%s%i",physBinned.Data(),bin),  24)
	   == 0 ){
    FillResult(h, r_hRD_upS_down_right[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_left_downstream_up_%s%i",physBinned.Data(),bin),  23)
	   == 0 ){
    FillResult(h, r_hRD_downS_up_left[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_right_downstream_up_%s%i",physBinned.Data(),bin),
		   24)
	   == 0 ){
    FillResult(h, r_hRD_downS_up_right[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_left_downstream_down_%s%i",physBinned.Data(),bin),
		   25)
	   == 0 ){
    FillResult(h, r_hRD_downS_down_left[bin], pars, iFit);
  }
  else if (strncmp(Form("%s", h->GetTitle() ),
		   Form("MuMu_right_downstream_down_%s%i",physBinned.Data(),bin),
		   26)
	   == 0 ){
    FillResult(h, r_hRD_downS_down_right[bin], pars, iFit);
  }
  else {
    cout << "Wrong hist name " << endl;
    exit(EXIT_FAILURE);
  }
    
  iter++;
  //Fit results*/
  
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


void BinnedFunctMFits(TString start=""){
  if (start==""){
    cout << "Script is used to try out many different Functional Mass Fitting";
    cout << " functions" << endl;
    cout << "Mass fitting is done in physics bins" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'BinnedFunctMFits()\'" << endl;
    exit(EXIT_FAILURE);
  }
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data";

  TFile *fRD =
    TFile::Open(Form("%s/leftRight_byTarget_WAll_AMDY_%ibins_%ihbin.root",
		     pathRD.Data(), nBins, hbins) );
  cout << Form("%s/leftRight_byTarget_WAll_AMDY_%ibins_%ihbin.root",
	       pathRD.Data(), nBins, hbins) << endl;
  if (!fRD ){
    cout << "RD or RD_noCorr file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nTargPol =4;
  TCanvas* c1[nTargPol];
  for (Int_t c=0; c<nTargPol; c++) {
    c1[c] = new TCanvas(); c1[c]->Divide(2, nBins); }

  TH1D *hRD_upS_up_left[nBins], *hRD_upS_up_right[nBins];
  TH1D *hRD_upS_down_left[nBins], *hRD_upS_down_right[nBins];
  TH1D *hRD_downS_up_left[nBins], *hRD_downS_up_right[nBins];
  TH1D *hRD_downS_down_left[nBins], *hRD_downS_down_right[nBins];
  for (Int_t bi=0; bi<nBins; bi++) {
    //upS
    hRD_upS_up_left[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s%i", 
					       physBinned.Data(), bi) );
    r_hRD_upS_up_left[bi] = (TH1D*) hRD_upS_up_left[bi]->Clone();

    hRD_upS_up_right[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_up_%s%i", 
						physBinned.Data(), bi) );
    r_hRD_upS_up_right[bi] = (TH1D*) hRD_upS_up_right[bi]->Clone();

    hRD_upS_down_left[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_down_%s%i",
						 physBinned.Data(), bi) );
    r_hRD_upS_down_left[bi] = (TH1D*) hRD_upS_down_left[bi]->Clone();

    hRD_upS_down_right[bi]=(TH1D*)fRD->Get(Form("MuMu_right_upstream_down_%s%i",
						physBinned.Data(), bi) );
    r_hRD_upS_down_right[bi] = (TH1D*) hRD_upS_down_right[bi]->Clone();

    //downS
    hRD_downS_up_left[bi] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_up_%s%i",
						 physBinned.Data(), bi) );
    r_hRD_downS_up_left[bi] = (TH1D*) hRD_downS_up_left[bi]->Clone();

    hRD_downS_up_right[bi]=(TH1D*)fRD->Get(Form("MuMu_right_downstream_up_%s%i",
						physBinned.Data(), bi) );
    r_hRD_downS_up_right[bi] = (TH1D*) hRD_downS_up_right[bi]->Clone();

    hRD_downS_down_left[bi]=
      (TH1D*)fRD->Get(Form("MuMu_left_downstream_down_%s%i", 
			   physBinned.Data(), bi) );
    r_hRD_downS_down_left[bi] = (TH1D*) hRD_downS_down_left[bi]->Clone();

    hRD_downS_down_right[bi]
      = (TH1D*)fRD->Get(Form("MuMu_right_downstream_down_%s%i",
			     physBinned.Data(), bi) );
    r_hRD_downS_down_right[bi] = (TH1D*) hRD_downS_down_right[bi]->Clone();

    
    SetUpHist(hRD_upS_up_left[bi]); SetUpHist(hRD_upS_up_right[bi]);
    SetUpHist(hRD_upS_down_left[bi]); SetUpHist(hRD_upS_down_right[bi]);
    SetUpHist(hRD_downS_up_left[bi]); SetUpHist(hRD_downS_up_right[bi]);
    SetUpHist(hRD_downS_down_left[bi]); SetUpHist(hRD_downS_down_right[bi]);
  }


  vector<Double_t> c_upS_up_left, e_upS_up_left;
  vector<Double_t> c_upS_down_left, e_upS_down_left;
  vector<Double_t> c_downS_up_left, e_downS_up_left;
  vector<Double_t> c_downS_down_left, e_downS_down_left;
  vector<Double_t> c_upS_up_right, e_upS_up_right;
  vector<Double_t> c_upS_down_right, e_upS_down_right;
  vector<Double_t> c_downS_up_right, e_downS_up_right;
  vector<Double_t> c_downS_down_right, e_downS_down_right;

  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 20, 0, 4);
  TH1D* hChi2_upS = new TH1D("hChi2_upS", "hChi2_upS", 20, 0, 4);
  TH1D* hChi2_downS = new TH1D("hChi2_downS", "hChi2_downS", 20, 0, 4);
  Double_t avgChi = 0.0;
  
  for (Int_t bi=0; bi<nBins; bi++) {
    c1[0]->cd(2*bi+1); gPad->SetLogy();//UpS
    Double_t Chi =FitGetPars(hRD_upS_up_left[bi], bi, c_upS_up_left,
			     e_upS_up_left, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_upS->Fill(Chi); 
    
    c1[0]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_upS_up_right[bi], bi, c_upS_up_right,
		    e_upS_up_right, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_upS->Fill(Chi); 

    c1[1]->cd(2*bi+1); gPad->SetLogy();
    Chi =FitGetPars(hRD_upS_down_left[bi], bi, c_upS_down_left,
		    e_upS_down_left, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_upS->Fill(Chi); 
    
    c1[1]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_upS_down_right[bi], bi, c_upS_down_right,
		    e_upS_down_right, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_upS->Fill(Chi); 

    
    c1[2]->cd(2*bi+1); gPad->SetLogy();//DownS
    Chi =FitGetPars(hRD_downS_up_left[bi], bi, c_downS_up_left,
		    e_downS_up_left, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_downS->Fill(Chi); 
    
    c1[2]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_up_right[bi], bi, c_downS_up_right,
		    e_downS_up_right, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_downS->Fill(Chi); 

    c1[3]->cd(2*bi+1); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_down_left[bi], bi, c_downS_down_left,
		    e_downS_down_left, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_downS->Fill(Chi); 
    
    c1[3]->cd(2*bi+2); gPad->SetLogy();
    Chi =FitGetPars(hRD_downS_down_right[bi], bi, c_downS_down_right,
		    e_downS_down_right, whichFit);
    hChi2->Fill(Chi); avgChi+= Chi;
    hChi2_downS->Fill(Chi); 
  }

  avgChi /= (nTargPol*nBins*2);
  cout << "\nAverage reduced chi2 is:   " << avgChi << "\n\n";
  
  TGraph* gChi2 = new TGraph(nBins*8, r_Xpoints, r_Chi2);
  SetUpTGraph(gChi2); gChi2->SetTitle("Chi2");
  TGraph* gNDF = new TGraph(nBins*8, r_Xpoints, r_NDF);
  SetUpTGraph(gNDF); gNDF->SetTitle("NDF");
  TGraph* gRedChi2 = new TGraph(nBins*8, r_Xpoints, r_RedChi2);
  SetUpTGraph(gRedChi2); gRedChi2->SetTitle("RedChi2");

  gStyle->SetOptStat(111111);
  TCanvas* cChi2 = new TCanvas(); cChi2->Divide(2, 2);
  cChi2->cd(1); SetUpHist(hChi2);
  
  SetUpHist(hChi2_upS); SetUpHist(hChi2_downS);
  hChi2->SetLineColor(kBlack);
  hChi2_upS->SetLineColor(kGreen); hChi2_downS->SetLineColor(kBlue);
  TF1 *fChi2 = new TF1("fChi2", RedChiSquareDistr, 0.00, 4, 2); 
  Double_t hChi2_Int = hChi2_downS->Integral()/( hChi2_downS->GetNbinsX()/( hChi2_downS->GetXaxis()->GetXmax() - hChi2_downS->GetXaxis()->GetXmin() ) );
  Double_t ndf = hbins-nPar;
  fChi2->SetParameter(0, ndf); fChi2->SetParameter(1, ndf*hChi2_Int);
  fChi2->Draw(); SetUpTF(fChi2);
  hChi2->Draw("sames"); hChi2_upS->Draw("sames"); hChi2_downS->Draw("sames");
  
  cChi2->cd(2); gChi2->Draw("AP");
  cChi2->cd(3); gNDF->Draw("AP");
  cChi2->cd(4); gRedChi2->Draw("AP");

  
  TCanvas* c2[nTargPol];
  for (Int_t c=0; c<nTargPol; c++) {
    c2[c] = new TCanvas(); c2[c]->Divide(2, nBins); }
  Double_t r_minY =0.3, r_maxY =1.7;
  
  for (Int_t bi=0; bi<nBins; bi++) {
    c2[0]->cd(2*bi+1);//upS
    r_hRD_upS_up_left[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_upS_up_left[bi]);
    r_hRD_upS_up_left[bi]->Draw("e"); DrawLine(r_hRD_upS_up_left[bi], 1.0);

    c2[0]->cd(2*bi+2);
    r_hRD_upS_up_right[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_upS_up_right[bi]);
    r_hRD_upS_up_right[bi]->Draw("e"); DrawLine(r_hRD_upS_up_right[bi], 1.0);

    c2[1]->cd(2*bi+1);
    r_hRD_upS_down_left[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_upS_down_left[bi]);
    r_hRD_upS_down_left[bi]->Draw("e"); DrawLine(r_hRD_upS_down_left[bi], 1.0);

    c2[1]->cd(2*bi+2);
    r_hRD_upS_down_right[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_upS_down_right[bi]);
    r_hRD_upS_down_right[bi]->Draw("e"); DrawLine(r_hRD_upS_down_right[bi],1.0);


    c2[2]->cd(2*bi+1);//downS
    r_hRD_downS_up_left[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_downS_up_left[bi]);
    r_hRD_downS_up_left[bi]->Draw("e"); DrawLine(r_hRD_downS_up_left[bi], 1.0);

    c2[2]->cd(2*bi+2);
    r_hRD_downS_up_right[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_downS_up_right[bi]);
    r_hRD_downS_up_right[bi]->Draw("e"); DrawLine(r_hRD_downS_up_right[bi],1.0);

    c2[3]->cd(2*bi+1);
    r_hRD_downS_down_left[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_downS_down_left[bi]);
    r_hRD_downS_down_left[bi]->Draw("e");DrawLine(r_hRD_downS_down_left[bi],1.0);

    c2[3]->cd(2*bi+2);
    r_hRD_downS_down_right[bi]->GetYaxis()->SetRangeUser(r_minY, r_maxY);
    SetUpHist(r_hRD_downS_down_right[bi]);
    r_hRD_downS_down_right[bi]->Draw("e");DrawLine(r_hRD_downS_down_right[bi],1.0);
  }


  TH1D *r_hRD_upS_up_left_cpy[nBins], *r_hRD_upS_up_right_cpy[nBins];
  TH1D *r_hRD_upS_down_left_cpy[nBins], *r_hRD_upS_down_right_cpy[nBins];
  TH1D *r_hRD_downS_up_left_cpy[nBins], *r_hRD_downS_up_right_cpy[nBins];
  TH1D *r_hRD_downS_down_left_cpy[nBins], *r_hRD_downS_down_right_cpy[nBins];
  TCanvas* c3 = new TCanvas(); c3->Divide(4);
  Double_t cpy_minX =2.5, cpy_maxX =3.5;
  Double_t cpy_minY =0.7, cpy_maxY =1.3;
  for (Int_t bi=0, icolor=1; bi<nBins; bi++, icolor++) {
    c3->cd(1);
    r_hRD_upS_up_left_cpy[bi] = (TH1D*)r_hRD_upS_up_left[bi]->Clone();
    r_hRD_upS_up_left_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_upS_up_left_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_upS_up_left_cpy[bi]->Draw("e");
      DrawLine(r_hRD_upS_up_left_cpy[bi], 1.0);
    }
    else {
      r_hRD_upS_up_left_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_upS_up_left_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(1);
    r_hRD_upS_up_right_cpy[bi] = (TH1D*)r_hRD_upS_up_right[bi]->Clone();
    r_hRD_upS_up_right_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_upS_up_right_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_upS_up_right_cpy[bi]->Draw("e");
      DrawLine(r_hRD_upS_up_right_cpy[bi], 1.0);
    }
    else {
      r_hRD_upS_up_right_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_upS_up_right_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(2);
    r_hRD_upS_down_left_cpy[bi] = (TH1D*)r_hRD_upS_down_left[bi]->Clone();
    r_hRD_upS_down_left_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_upS_down_left_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_upS_down_left_cpy[bi]->Draw("e");
      DrawLine(r_hRD_upS_down_left_cpy[bi], 1.0);
    }
    else {
      r_hRD_upS_down_left_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_upS_down_left_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(2);
    r_hRD_upS_down_right_cpy[bi] = (TH1D*)r_hRD_upS_down_right[bi]->Clone();
    r_hRD_upS_down_right_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_upS_down_right_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_upS_down_right_cpy[bi]->Draw("e");
      DrawLine(r_hRD_upS_down_right_cpy[bi], 1.0);
    }
    else {
      r_hRD_upS_down_right_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_upS_down_right_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(3);
    r_hRD_downS_up_left_cpy[bi] = (TH1D*)r_hRD_downS_up_left[bi]->Clone();
    r_hRD_downS_up_left_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_downS_up_left_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_downS_up_left_cpy[bi]->Draw("e");
      DrawLine(r_hRD_downS_up_left_cpy[bi], 1.0);
    }
    else {
      r_hRD_downS_up_left_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_downS_up_left_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(3);
    r_hRD_downS_up_right_cpy[bi] = (TH1D*)r_hRD_downS_up_right[bi]->Clone();
    r_hRD_downS_up_right_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_downS_up_right_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_downS_up_right_cpy[bi]->Draw("e");
      DrawLine(r_hRD_downS_up_right_cpy[bi], 1.0);
    }
    else {
      r_hRD_downS_up_right_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_downS_up_right_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(4);
    r_hRD_downS_down_left_cpy[bi] = (TH1D*)r_hRD_downS_down_left[bi]->Clone();
    r_hRD_downS_down_left_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_downS_down_left_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_downS_down_left_cpy[bi]->Draw("e");
      DrawLine(r_hRD_downS_down_left_cpy[bi], 1.0);
    }
    else {
      r_hRD_downS_down_left_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_downS_down_left_cpy[bi]->SetLineColor(icolor);
    }

    c3->cd(4);
    r_hRD_downS_down_right_cpy[bi] = (TH1D*)r_hRD_downS_down_right[bi]->Clone();
    r_hRD_downS_down_right_cpy[bi]->GetXaxis()->SetRangeUser(cpy_minX, cpy_maxX);
    r_hRD_downS_down_right_cpy[bi]->GetYaxis()->SetRangeUser(cpy_minY, cpy_maxY);
    if (bi==0) {
      r_hRD_downS_down_right_cpy[bi]->Draw("e");
      DrawLine(r_hRD_downS_down_right_cpy[bi], 1.0);
    }
    else {
      r_hRD_downS_down_right_cpy[bi]->Draw("esame");
      if (icolor==5) icolor += 1;
      r_hRD_downS_down_right_cpy[bi]->SetLineColor(icolor);
    }
  }

  TString fOutput = Form("FitMass_%s_%.2f_%.2f_",
			 physBinned.Data(), Mmin, Mmax);
  if (PolCorr) fOutput += "corr.root";
  else fOutput += "noCorr.root";
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    /*g_AN_JPsi_upS_up->Write("JPsi_upS_up");
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
    g_JPsi_Right_downS_down->Write("JPsi_Right_downS_down");//*/
  }

  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Mass Range is: " << Mmin << "  -  " << Mmax << endl;
  cout << "Polarization was performed:  " << PolCorr << endl;
  cout << "Physics binning is:  " << physBinned << endl;
  cout << "Fit failed:  " << Failure << "  times" << endl;
  cout << "Number of histogram bins used:  "<< hbins << endl;
  cout << "Outputting asymmetry for    JPsi" << endl;
  cout << "Fit used was:  " << whichFit << endl;
  cout << "nPars used in fit:  "  << nPar << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
