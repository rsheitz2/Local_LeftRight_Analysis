#ifndef FIT_NINE_H
#define FIT_NINE_H
//2 Gaussians for JPsi and psi'
//      psi' M/W = Get_nine_Ratio(targ)*JPsi M/W
//1 Exponentials for background and DY

Double_t Get_nine_Ratio(){
  Double_t ratioPsi = 1.15;
  
  return ratioPsi;
}


Double_t Fit_nine(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[8];
  Double_t ratioPsi = Get_nine_Ratio();
  //par[3] =alpha, par[4] =n
  Double_t A = TMath::Power(par[4]/par[3], par[4])*TMath::Exp(-par[3]*par[3]/2.0);
  Double_t B = par[4]/par[3] - par[3];
  Double_t C = (par[4]/par[3])*(1.0/(par[4]-1.0))*TMath::Exp(-par[3]*par[3]/2.0);
  Double_t D = TMath::Pi()*(1+TMath::Erf( par[3]/TMath::Sqrt(2)) );

  Double_t Norm = 1.0/(par[2]*(C+D) );

  Double_t arg_JPsi = (x[0] - par[1])/par[2];
  Double_t JPsi;
  if (arg_JPsi > -par[3] ) JPsi = par[0]*Norm*TMath::Exp(-0.5*arg_JPsi*arg_JPsi);
  else JPsi = par[0]*Norm*A*TMath::Power((B - arg_JPsi), -par[1]);
    
  Double_t arg_psi = ( x[0] - par[1]*ratioPsi)/(par[2]*ratioPsi);
  Double_t psi;
  if (arg_psi > -par[3] ) psi = par[5]*Norm*TMath::Exp(-0.5*arg_psi*arg_psi);
  else psi = par[5]*Norm*A*TMath::Power((B - arg_psi), -par[1]);

  Double_t DY = par[6]*TMath::Exp( par[7]*xShift );
    
  return JPsi + psi + DY;
}


void Paras_nine(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax){
  Int_t ipar =0;
  Double_t ratioPsi = Get_nine_Ratio();
  
  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_DY =  h->GetBinContent(Mmin);

  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, 1.1); ipar++;
  fitFunc->SetParameter(ipar, 3.0); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, -0.84); ipar++;
  fitFunc->SetParameter(ipar, Mmin); ipar++;//Mmin
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  
  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 1.0, 5.0); ipar++;
  fitFunc->SetParLimits(ipar, 1.1, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -4.0, -0.5); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* SetupFunc_nine(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc,
		   Double_t Mmin, Double_t Mmax, Int_t *nPar){
  //Setup fit function
  *nPar =9;

  fitFunc = new TF1("fitFunc", Fit_nine, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_nine(h, fitFunc, *nPar, Mmin, Mmax);

  return fitFunc;
}


void ProcessPars_nine(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
			   TMatrixDSym &cov, TString process, Int_t nPars,
			   Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  Double_t psiMW = Get_nine_Ratio();
  
  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[5];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[6];//DY
  processPars[7] = pars[7];
  
  if (process=="JPsi"){
    for (Int_t i=0; i<5; i++) {
      LR_cov[i] = cov(0, i);
      LR_cov[i+5] = cov(1, i);
      LR_cov[i+10] = cov(2, i);
      LR_cov[i+15] = cov(3, i);
      LR_cov[i+20] = cov(4, i);
    }
  }
  else if (process=="psi"){
    LR_cov[0] = cov(5, 5);
    LR_cov[5] = cov(1, 5)*psiMW;
    LR_cov[10] = cov(2, 5)*psiMW;
    LR_cov[15] = cov(3, 5);
    LR_cov[20] = cov(4, 5);
    
    for (Int_t i=1; i<5; i++) {
      Double_t multi = (i == 1) ? psiMW : 1.0;
      
      LR_cov[i] = cov(5, i)*multi;
      LR_cov[i+5] = cov(1, i)*psiMW*multi;
      LR_cov[i+10] = cov(2, i)*psiMW*multi;
      LR_cov[i+15] = cov(3, i)*multi;
      LR_cov[i+20] = cov(4, i)*multi;
    }
  }
  else if (process=="DY"){
    cout << "covariance matrix doesn't work yet" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "Process not defined well in Fit_nine" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_nine(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		       TMatrixDSym &cov, TString process, Int_t nPars,
		       Bool_t hIsUpS, Double_t *e_processPars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  
  Double_t psiMW = Get_nine_Ratio();

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[5];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[6];//DY
  processPars[7] = pars[7];

  e_processPars[0] = e_pars[0];//JPsi
  e_processPars[1] = e_pars[1];
  e_processPars[2] = e_pars[2];

  e_processPars[3] = e_pars[5];//psi
  e_processPars[4] = e_pars[1]*psiMW;
  e_processPars[5] = e_pars[2]*psiMW;

  e_processPars[6] = e_pars[6];//DY
  e_processPars[7] = e_pars[7];
  
  if (process=="JPsi"){
    for (Int_t i=0; i<5; i++) {
      LR_cov[i] = cov(0, i);
      LR_cov[i+5] = cov(1, i);
      LR_cov[i+10] = cov(2, i);
      LR_cov[i+15] = cov(3, i);
      LR_cov[i+20] = cov(4, i);
    }
  }
  else if (process=="psi"){
    LR_cov[0] = cov(5, 5);
    LR_cov[5] = cov(1, 5)*psiMW;
    LR_cov[10] = cov(2, 5)*psiMW;
    LR_cov[15] = cov(3, 5);
    LR_cov[20] = cov(4, 5);
    
    for (Int_t i=1; i<5; i++) {
      Double_t multi = (i == 1) ? psiMW : 1.0;
      
      LR_cov[i] = cov(5, i)*multi;
      LR_cov[i+5] = cov(1, i)*psiMW*multi;
      LR_cov[i+10] = cov(2, i)*psiMW*multi;
      LR_cov[i+15] = cov(3, i)*multi;
      LR_cov[i+20] = cov(4, i)*multi;
    }
  }
  else if (process=="DY"){
    cout << "covariance matrix doesn't work yet" << endl;
    exit(EXIT_FAILURE);
    /*LR_cov[0] = cov(4, 4);
    LR_cov[1] = cov(4, 5);
    LR_cov[2] = 0.0;
        
    LR_cov[3] = cov(5, 4);
    LR_cov[4] = cov(5, 5);
    LR_cov[5] = 0.0;

    LR_cov[6] = 0.0;
    LR_cov[7] = 0.0;
    LR_cov[8] = 0.0;*/
  }
  else{
    cout << "Process not defined well in Fit_nine" << endl;
    exit(EXIT_FAILURE);
  }
}



TF1* ComponentFuncts_nine(Double_t *pars, Double_t Mmin, Double_t Mmax,
			  TString process){
  Double_t psiMW = Get_nine_Ratio();
  
  TF1 *JPsi = new TF1("f_JPsi", f_CrystalBall, Mmin, Mmax);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);

  TF1 *psi = new TF1("f_psi", f_CrystalBall, Mmin, Mmax);
  psi->SetParameters(pars[5], pars[1]*psiMW, pars[2]*psiMW, pars[3], pars[4]);

  TF1 *DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  DY->SetParameters(pars[6], pars[7], Mmin);

  JPsi->SetLineColor(kGreen); 
  psi->SetLineColor(kGreen);
  DY->SetLineColor(kGreen); 

  JPsi->Draw("same");
  psi->Draw("same");
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
    cout << "Wrong input process Fit_nine" << endl;
    exit(EXIT_FAILURE);
  }
}


void IntegrateLR_nine(TF1 *f, Double_t *pars, Double_t *LR_cov,
		       Double_t LR_Mmin, Double_t LR_Mmax,
		       Double_t *LR, Double_t *e_LR,
		       Double_t Mmin, Double_t Mmax, TString process){
  (*LR) = f->Integral(LR_Mmin, LR_Mmax);
  
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, pars, LR_cov);
  }
  else if (process =="psi"){
    Double_t psi_pars[] = {pars[5], pars[1], pars[2], pars[3], pars[4]};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, psi_pars, LR_cov);
  }
  else if (process =="DY"){
    Double_t DY_pars[] = {pars[6], pars[7], Mmin};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov);
  }
}


#endif
