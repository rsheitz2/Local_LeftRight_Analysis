#ifndef FIT_ONE_H
#define FIT_ONE_H


Double_t Get_1_85_one_Ratio(){
  Double_t ratioPsi = 1.145;
  
  return ratioPsi;
}


Double_t Fit_one(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi = Get_1_85_one_Ratio();
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[1]*ratioPsi)/(par[2]*ratioPsi);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t OC = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


void Paras_one(TF1 *fitFunc, Int_t nPar){
  Int_t ipar =0;
  fitFunc->SetParameter(ipar, 1e5); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, 3e3); ipar++;//psi
  fitFunc->SetParameter(ipar, 1e6); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.5); ipar++;
  fitFunc->SetParameter(ipar, 1e6); ipar++;//OC
  fitFunc->SetParameter(ipar, -2.5); ipar++;
  fitFunc->SetParameter(ipar, 8e3); ipar++;//DY
  fitFunc->SetParameter(ipar, -0.9); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }

  
  ipar=0;
  fitFunc->SetParLimits(ipar, 0, 1e6); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 0, 1e6); ipar++;//psi
  fitFunc->SetParLimits(ipar, 1e5, 1e7); ipar++;//Bg
  fitFunc->SetParLimits(ipar, -6.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 5e4, 1e7); ipar++;//OC
  fitFunc->SetParLimits(ipar, -3.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 1e3, 1e4); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void NicePars_one(TH1D *h, TF1 *fitFunc,
		  Double_t *nicePars, Double_t *eNicePars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *ePars = fitFunc->GetParErrors();

  Double_t ratioPsi = Get_1_85_one_Ratio();
  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0){
    nicePars[0] = pars[1];
    nicePars[2] = pars[2];
    nicePars[4] = pars[1]*ratioPsi;
    nicePars[6] = pars[2]*ratioPsi;

    eNicePars[0] = ePars[1];
    eNicePars[2] = ePars[2];
    eNicePars[4] = ePars[1]*ratioPsi;
    eNicePars[6] = ePars[2]*ratioPsi;
  }
  else {
    nicePars[1] = pars[1];
    nicePars[3] = pars[2];
    nicePars[5] = pars[1]*ratioPsi;
    nicePars[7] = pars[2]*ratioPsi;

    eNicePars[1] = ePars[1];
    eNicePars[3] = ePars[2];
    eNicePars[5] = ePars[1]*ratioPsi;
    eNicePars[7] = ePars[2]*ratioPsi;
  }
    
}

#endif
