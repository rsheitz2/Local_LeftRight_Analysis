#ifndef FIT_FIVE_H
#define FIT_FIVE_H

Double_t Get_1_85_five_Ratio(TString targ){
  Double_t ratioPsi_upS_2exp = 1.14;
  Double_t ratioPsi_downS_2exp = 1.13;
  
  if (targ=="UpS") return ratioPsi_upS_2exp;
  else if (targ=="DownS") return ratioPsi_downS_2exp;
  else {
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t Fit_five_upS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi_upS_2exp = Get_1_85_five_Ratio("UpS");
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi
    = ( x[0] - par[1]*ratioPsi_upS_2exp)/(par[2]*ratioPsi_upS_2exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t DY = par[6]*TMath::Exp( par[7]*xShift );
    
  return CombBg + JPsi + psi + DY;
}




Double_t Fit_five_downS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi_downS_2exp = Get_1_85_five_Ratio("DownS");
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi
    = ( x[0] - par[1]*ratioPsi_downS_2exp)/(par[2]*ratioPsi_downS_2exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t DY = par[6]*TMath::Exp( par[7]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_five(TF1 *fitFunc, Int_t nPar){
  Int_t ipar =0;
  fitFunc->SetParameter(ipar, 1e5); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, 3e3); ipar++;//psi
  fitFunc->SetParameter(ipar, 1e6); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.5); ipar++;
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
  fitFunc->SetParLimits(ipar, 1e3, 1e4); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void NicePars_five(TH1D *h, TF1 *fitFunc,
		  Double_t *nicePars, Double_t *eNicePars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *ePars = fitFunc->GetParErrors();

  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0){
    Double_t ratioPsi_upS_2exp = Get_1_85_five_Ratio("UpS");
    nicePars[0] = pars[1];
    nicePars[2] = pars[2];
    nicePars[4] = pars[1]*ratioPsi_upS_2exp;
    nicePars[6] = pars[2]*ratioPsi_upS_2exp;

    eNicePars[0] = ePars[1];
    eNicePars[2] = ePars[2];
    eNicePars[4] = ePars[1]*ratioPsi_upS_2exp;
    eNicePars[6] = ePars[2]*ratioPsi_upS_2exp;
  }
  else {
    Double_t ratioPsi_downS_2exp = Get_1_85_five_Ratio("DownS");
    nicePars[1] = pars[1];
    nicePars[3] = pars[2];
    nicePars[5] = pars[1]*ratioPsi_downS_2exp;
    nicePars[7] = pars[2]*ratioPsi_downS_2exp;

    eNicePars[1] = ePars[1];
    eNicePars[3] = ePars[2];
    eNicePars[5] = ePars[1]*ratioPsi_downS_2exp;
    eNicePars[7] = ePars[2]*ratioPsi_downS_2exp;
  }
    
}


void ProcessPars_five(TF1 *fitFunc, Double_t *processPars, Double_t *LR_cov,
		      TFitResultPtr &status, TString process,
		      Int_t nPars, Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();

  Double_t psiMW;
  if (hIsUpS) psiMW = Get_1_85_five_Ratio("UpS");
  else psiMW = Get_1_85_five_Ratio("DownS");

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[3];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[6];//DY
  processPars[7] = pars[7];

  Double_t *cov = status->GetCovarianceMatrix().GetMatrixArray();
  if (process=="JPsi"){
    LR_cov[0] = cov[0];
    LR_cov[1] = cov[1];
    LR_cov[2] = cov[2];
    
    LR_cov[3] = cov[nPars];
    LR_cov[4] = cov[nPars+1];
    LR_cov[5] = cov[nPars+2];
    
    LR_cov[6] = cov[nPars*2];
    LR_cov[7] = cov[nPars*2+1];
    LR_cov[8] = cov[nPars*2+2];
  }
  else if (process=="psi"){
    LR_cov[0] = cov[nPars*3+3];
    LR_cov[1] = cov[nPars*3+1]*psiMW;
    LR_cov[2] = cov[nPars*3+2]*psiMW;
    
    LR_cov[3] = cov[nPars+3]*psiMW;
    LR_cov[4] = cov[nPars+1]*psiMW*psiMW;
    LR_cov[5] = cov[nPars+2]*psiMW*psiMW;
    
    LR_cov[6] = cov[nPars*2+3]*psiMW;
    LR_cov[7] = cov[nPars*2+1]*psiMW*psiMW;
    LR_cov[8] = cov[nPars*2+2]*psiMW*psiMW;
  }
  else if (process=="DY"){
    LR_cov[0] = cov[nPars*6+6];
    LR_cov[1] = cov[nPars*6+7];
        
    LR_cov[3] = cov[nPars*7+6];
    LR_cov[4] = cov[nPars*7+7];
  }
  else{
    cout << "Process not defined well in Fit_five" << endl;
    exit(EXIT_FAILURE);
  }
}


#endif
