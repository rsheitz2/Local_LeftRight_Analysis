#ifndef FIT_SIX_H
#define FIT_SIX_H

Double_t Fit_six(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[4])/(par[5]);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_six(TF1 *fitFunc, Int_t nPar){
  Int_t ipar =0;
  fitFunc->SetParameter(ipar, 9e4); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, 2e3); ipar++;//psi
  fitFunc->SetParameter(ipar, 3.6); ipar++;
  fitFunc->SetParameter(ipar, 0.25); ipar++;
  fitFunc->SetParameter(ipar, 8.8e5); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.8); ipar++;
  fitFunc->SetParameter(ipar, 8e3); ipar++;//DY
  fitFunc->SetParameter(ipar, -0.84); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }

  
  ipar=0;
  fitFunc->SetParLimits(ipar, 0, 1e6); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 1e3, 1e6); ipar++;//psi
  fitFunc->SetParLimits(ipar, 3.3, 4.0); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.4); ipar++;
  fitFunc->SetParLimits(ipar, 1e5, 1e7); ipar++;//Bg
  fitFunc->SetParLimits(ipar, -6.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 1e3, 2e4); ipar++;//DY
  fitFunc->SetParLimits(ipar, -0.9, 0); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void NicePars_six(TH1D *h, TF1 *fitFunc,
		  Double_t *nicePars, Double_t *eNicePars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *ePars = fitFunc->GetParErrors();

  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0){
    nicePars[0] = pars[1];
    nicePars[2] = pars[2];
    nicePars[4] = pars[4];
    nicePars[6] = pars[5];

    eNicePars[0] = ePars[1];
    eNicePars[2] = ePars[2];
    eNicePars[4] = ePars[4];
    eNicePars[6] = ePars[5];
  }
  else {
    nicePars[1] = pars[1];
    nicePars[3] = pars[2];
    nicePars[5] = pars[4];
    nicePars[7] = pars[5];

    eNicePars[1] = ePars[1];
    eNicePars[3] = ePars[2];
    eNicePars[5] = ePars[4];
    eNicePars[7] = ePars[5];
  }
}

#endif
