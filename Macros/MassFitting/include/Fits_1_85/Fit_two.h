#ifndef FIT_TWO_H
#define FIT_TWO_H

Double_t Get_1_85_two_Ratio(TString targ){
  Double_t ratioPsi_upS = 1.14;
  Double_t ratioPsi_downS = 1.13;
  
  if (targ=="UpS") return ratioPsi_upS;
  else if (targ=="DownS") return ratioPsi_downS;
  else {
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t Fit_two_upS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi_upS = Get_1_85_two_Ratio("UpS");
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[1]*ratioPsi_upS)/(par[2]*ratioPsi_upS);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t OC = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


Double_t Fit_two_downS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi_downS = Get_1_85_two_Ratio("DownS");
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi = ( x[0] - par[1]*ratioPsi_downS)/(par[2]*ratioPsi_downS);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t OC = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + OC + JPsi + psi + DY;
}


void Paras_two(TF1 *fitFunc, Int_t nPar){
  Int_t ipar =0;
  fitFunc->SetParameter(ipar, 1e4); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, 3e2); ipar++;//psi
  fitFunc->SetParameter(ipar, 1e5); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.5); ipar++;
  fitFunc->SetParameter(ipar, 1e5); ipar++;//OC
  fitFunc->SetParameter(ipar, -2.5); ipar++;
  fitFunc->SetParameter(ipar, 1e3); ipar++;//DY
  fitFunc->SetParameter(ipar, -0.9); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }

  
  ipar=0;
  fitFunc->SetParLimits(ipar, 0, 1e5); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 0, 1e5); ipar++;//psi
  fitFunc->SetParLimits(ipar, 1e4, 1e6); ipar++;//Bg
  fitFunc->SetParLimits(ipar, -6.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 1e4, 1e6); ipar++;//OC
  fitFunc->SetParLimits(ipar, -3.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 1e2, 1e4); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void NicePars_two(TH1D *h, TF1 *fitFunc,
		  Double_t *nicePars, Double_t *eNicePars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *ePars = fitFunc->GetParErrors();

  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0){
    Double_t ratioPsi_upS = Get_1_85_two_Ratio("UpS");
    nicePars[0] = pars[1];
    nicePars[2] = pars[2];
    nicePars[4] = pars[1]*ratioPsi_upS;
    nicePars[6] = pars[2]*ratioPsi_upS;

    eNicePars[0] = ePars[1];
    eNicePars[2] = ePars[2];
    eNicePars[4] = ePars[1]*ratioPsi_upS;
    eNicePars[6] = ePars[2]*ratioPsi_upS;
  }
  else {
    Double_t ratioPsi_downS = Get_1_85_two_Ratio("DownS");
    nicePars[1] = pars[1];
    nicePars[3] = pars[2];
    nicePars[5] = pars[1]*ratioPsi_downS;
    nicePars[7] = pars[2]*ratioPsi_downS;

    eNicePars[1] = ePars[1];
    eNicePars[3] = ePars[2];
    eNicePars[5] = ePars[1]*ratioPsi_downS;
    eNicePars[7] = ePars[2]*ratioPsi_downS;
  }
    
}

#endif
