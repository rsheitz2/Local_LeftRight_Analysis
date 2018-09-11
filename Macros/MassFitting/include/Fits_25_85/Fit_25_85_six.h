#ifndef FIT_25_85_SIX_H
#define FIT_25_85_SIX_H


Double_t Get_25_85_six_Ratio(TString targ){
  Double_t ratioPsi_UpS_2exp = 1.15;//12.338
  Double_t ratioPsi_DownS_2exp = 1.135;//4.54284
  
  if (targ=="UpS") return ratioPsi_UpS_2exp;
  else if (targ=="DownS") return ratioPsi_DownS_2exp;
  else {
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t Fit_25_85_six_UpS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi_UpS_2exp = Get_25_85_six_Ratio("UpS");
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi
    = ( x[0] - par[1]*ratioPsi_UpS_2exp)/(par[2]*ratioPsi_UpS_2exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t DY = par[6]*TMath::Exp( par[7]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_25_85_six_DownS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-Mmin;

  Double_t ratioPsi_DownS_2exp = Get_25_85_six_Ratio("DownS");
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi
    = ( x[0] - par[1]*ratioPsi_DownS_2exp)/(par[2]*ratioPsi_DownS_2exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t CombBg = par[4]*TMath::Exp( par[5]*xShift );
  Double_t DY = par[6]*TMath::Exp( par[7]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_25_85_six(TF1 *fitFunc, Int_t nPar, Int_t nBinsPhys=5){
  Int_t ipar =0;
  if (nBinsPhys==5){
    fitFunc->SetParameter(ipar, 4e3); ipar++;//JPsi
    fitFunc->SetParameter(ipar, 3.1); ipar++;
    fitFunc->SetParameter(ipar, 0.15); ipar++;
    fitFunc->SetParameter(ipar, 1e2); ipar++;//psi
    fitFunc->SetParameter(ipar, 1e3); ipar++;//CombBg
    fitFunc->SetParameter(ipar, -3.5); ipar++;
    fitFunc->SetParameter(ipar, 1e2); ipar++;//DY
    fitFunc->SetParameter(ipar, -0.9); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }
  
    ipar=0;
    fitFunc->SetParLimits(ipar, 0, 1e4); ipar++;//JPsi
    fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
    fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
    fitFunc->SetParLimits(ipar, 0, 5e3); ipar++;//psi
    fitFunc->SetParLimits(ipar, 1e2, 1e4); ipar++;//Bg
    fitFunc->SetParLimits(ipar, -4.0, -1.0); ipar++;
    fitFunc->SetParLimits(ipar, 0, 1e4); ipar++;//DY
    fitFunc->SetParLimits(ipar, -2.0, 0); ipar++;
  }//nBinsPhys==5
  else if (nBinsPhys==1){
    fitFunc->SetParameter(ipar, 8e4); ipar++;//JPsi
    fitFunc->SetParameter(ipar, 3.1); ipar++;
    fitFunc->SetParameter(ipar, 0.15); ipar++;
    fitFunc->SetParameter(ipar, 2e3); ipar++;//psi
    fitFunc->SetParameter(ipar, 1e4); ipar++;//CombBg
    fitFunc->SetParameter(ipar, -3.5); ipar++;
    fitFunc->SetParameter(ipar, 1e3); ipar++;//DY
    fitFunc->SetParameter(ipar, -0.9); ipar++;
    if (ipar != nPar){
      cout << "ipar problem" << endl;
      exit(EXIT_FAILURE);
    }
  
    ipar=0;
    fitFunc->SetParLimits(ipar, 1e3, 1e7); ipar++;//JPsi
    fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
    fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
    fitFunc->SetParLimits(ipar, 1e2, 1e5); ipar++;//psi
    fitFunc->SetParLimits(ipar, 1e3, 1e5); ipar++;//Bg
    fitFunc->SetParLimits(ipar, -4.0, -1.0); ipar++;
    fitFunc->SetParLimits(ipar, 1e2, 1e4); ipar++;//DY
    fitFunc->SetParLimits(ipar, -2.0, 0); ipar++;
  }//nBinsPhys==1
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void Paras_25_85_six(TH1D *h, TF1 *fitFunc, Int_t nPar){
  Int_t ipar =0;
  Double_t ratioPsi = Get_25_85_six_Ratio("UpS");
  ratioPsi += Get_25_85_six_Ratio("DownS");
  ratioPsi /= 2.0;
  
  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(4.5) );
  A_DY /= (TMath::Exp(-0.84*(5.0-Mmin)));

  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_Bg); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.5); ipar++;
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, -0.84); ipar++;
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  
  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*5.0); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 2.5, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_Bg/factor, A_Bg*factor); ipar++;//Bg
  fitFunc->SetParLimits(ipar, -5.0, -1.5); ipar++;
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -1.5, -0.5); ipar++;
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void NicePars_25_85_six(TH1D *h, TF1 *fitFunc,
		  Double_t *nicePars, Double_t *eNicePars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *ePars = fitFunc->GetParErrors();

  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0){
    Double_t ratioPsi_UpS_2exp = Get_25_85_six_Ratio("UpS");
    
    nicePars[0] = pars[1];
    nicePars[2] = pars[2];
    nicePars[4] = pars[1]*ratioPsi_UpS_2exp;
    nicePars[6] = pars[2]*ratioPsi_UpS_2exp;

    eNicePars[0] = ePars[1];
    eNicePars[2] = ePars[2];
    eNicePars[4] = ePars[1]*ratioPsi_UpS_2exp;
    eNicePars[6] = ePars[2]*ratioPsi_UpS_2exp;
  }
  else {
    Double_t ratioPsi_DownS_2exp = Get_25_85_six_Ratio("DownS");
    
    nicePars[1] = pars[1];
    nicePars[3] = pars[2];
    nicePars[5] = pars[1]*ratioPsi_DownS_2exp;
    nicePars[7] = pars[2]*ratioPsi_DownS_2exp;

    eNicePars[1] = ePars[1];
    eNicePars[3] = ePars[2];
    eNicePars[5] = ePars[1]*ratioPsi_DownS_2exp;
    eNicePars[7] = ePars[2]*ratioPsi_DownS_2exp;
  }
    
}


void ProcessPars_25_85_six(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
			   TMatrixDSym &cov, TString process, Int_t nPars,
			   Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();

  Double_t psiMW;
  if (hIsUpS) psiMW = Get_25_85_six_Ratio("UpS");
  else psiMW = Get_25_85_six_Ratio("DownS");

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[3];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[6];//DY
  processPars[7] = pars[7];

  if (process=="JPsi"){
    LR_cov[0] = cov(0, 0);
    LR_cov[1] = cov(0, 1);
    LR_cov[2] = cov(0, 2);
    
    LR_cov[3] = cov(1, 0);
    LR_cov[4] = cov(1, 1);
    LR_cov[5] = cov(1, 2);
    
    LR_cov[6] = cov(2, 0);
    LR_cov[7] = cov(2, 1);
    LR_cov[8] = cov(2, 2);
  }
  else if (process=="psi"){
    LR_cov[0] = cov(3, 3);
    LR_cov[1] = cov(3, 1)*psiMW;
    LR_cov[2] = cov(3, 2)*psiMW;
    
    LR_cov[3] = cov(1, 3)*psiMW;
    LR_cov[4] = cov(1, 1)*psiMW*psiMW;
    LR_cov[5] = cov(1, 2)*psiMW*psiMW;
    
    LR_cov[6] = cov(2, 3)*psiMW;
    LR_cov[7] = cov(2, 1)*psiMW*psiMW;
    LR_cov[8] = cov(2, 2)*psiMW*psiMW;
  }
  else if (process=="DY"){
    LR_cov[0] = cov(6, 6);
    LR_cov[1] = cov(6, 7);
    LR_cov[2] = 0.0;
        
    LR_cov[3] = cov(7, 6);
    LR_cov[4] = cov(7, 7);
    LR_cov[5] = 0.0;

    LR_cov[6] = 0.0;
    LR_cov[7] = 0.0;
    LR_cov[8] = 0.0;
  }
  else{
    cout << "Process not defined well in Fit_25_85_six" << endl;
    exit(EXIT_FAILURE);
  }
}


#endif
