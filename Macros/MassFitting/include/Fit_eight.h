#ifndef FIT_EIGHT_H
#define FIT_EIGHT_H
//2 Gaussians for JPsi and psi'
//      psi' M/W = Get_eight_Ratio(targ)*JPsi M/W
//2 Exponentials for background and DY

Double_t Get_eight_Ratio(TString targ){
  if (targ =="UpS") return 1.155;
  else if (targ =="DownS") return 1.135;
  else{
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t Fit_eight_UpS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_eight_Ratio("UpS");
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

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_eight_DownS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_eight_Ratio("DownS");
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

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_eight_scan(Double_t *x, Double_t *par){
  //Feed one scan parameter as par[11]=psiMW
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = par[11];
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

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_eight(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax){
  Int_t ipar =0;
  Double_t ratioPsi = Get_eight_Ratio("UpS");
  ratioPsi += Get_eight_Ratio("DownS");
  ratioPsi /= 2.0;

  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(5.0) );
  Double_t DY_slope = -1.0, DY_Mmin = 5.0;
  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));

  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.12); ipar++;
  fitFunc->SetParameter(ipar, 0.16); ipar++;
  fitFunc->SetParameter(ipar, 1.4); ipar++;
  fitFunc->SetParameter(ipar, 3.2); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_Bg); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.5); ipar++;
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, DY_slope); ipar++;
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
  fitFunc->SetParLimits(ipar, 1.0, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, 2.5, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_Bg/factor, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(ipar, -7.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0.0); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void Paras_eight(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax, Double_t psiMW){
  Int_t ipar =0;
  Double_t ratioPsi =psiMW;

  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(5.0) );
  Double_t DY_slope = -1.0, DY_Mmin = 5.0;
  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));

  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.12); ipar++;
  fitFunc->SetParameter(ipar, 0.16); ipar++;
  fitFunc->SetParameter(ipar, 1.4); ipar++;
  fitFunc->SetParameter(ipar, 3.2); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_Bg); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.5); ipar++;
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, DY_slope); ipar++;
  fitFunc->SetParameter(ipar, Mmin); ipar++;//Mmin
  fitFunc->SetParameter(ipar, psiMW); ipar++;//psiMW
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  
  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 1.0, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, 2.5, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_Bg/factor, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(ipar, -7.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0.0); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  fitFunc->SetParLimits(ipar, psiMW, psiMW); ipar++;//psiMW
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* SetupFunc_eight(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc,
		     Double_t Mmin, Double_t Mmax, Int_t *nPar){
  //Setup fit function
  *nPar =11;

  if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_eight_UpS, Mmin, Mmax, *nPar);
  else fitFunc = new TF1("fitFunc", Fit_eight_DownS, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_eight(h, fitFunc, *nPar, Mmin, Mmax);

  return fitFunc;
}


TF1* SetupFunc_eight(TH1D *h, TF1 *fitFunc, Double_t Mmin, Double_t Mmax,
		     Int_t *nPar, Double_t psiMW){
  //Setup fit function for scanning psiMW
  *nPar =12;

  fitFunc = new TF1("fitFunc", Fit_eight_scan, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_eight(h, fitFunc, *nPar, Mmin, Mmax, psiMW);

  return fitFunc;
}


void ProcessPars_eight(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		       TMatrixDSym &cov, TString process, Int_t nPars,
		       Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");
  
  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[5];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[8];//DY
  processPars[7] = pars[9];
  
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
    cout << "Process not defined well in Fit_eight" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_eight(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		       TMatrixDSym &cov, TString process, Int_t nPars,
		       Bool_t hIsUpS, Double_t *e_processPars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[5];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[8];//DY
  processPars[7] = pars[9];

  e_processPars[0] = e_pars[0];//JPsi
  e_processPars[1] = e_pars[1];
  e_processPars[2] = e_pars[2];

  e_processPars[3] = e_pars[5];//psi
  e_processPars[4] = e_pars[1]*psiMW;
  e_processPars[5] = e_pars[2]*psiMW;

  e_processPars[6] = e_pars[8];//DY
  e_processPars[7] = e_pars[9];
  
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
    
    /*LR_cov[0] = cov(3, 3);
      LR_cov[1] = cov(3, 1)*psiMW;
      LR_cov[2] = cov(3, 2)*psiMW;
    
      LR_cov[3] = cov(1, 3)*psiMW;
      LR_cov[4] = cov(1, 1)*psiMW*psiMW;
      LR_cov[5] = cov(1, 2)*psiMW*psiMW;
    
      LR_cov[6] = cov(2, 3)*psiMW;
      LR_cov[7] = cov(2, 1)*psiMW*psiMW;
      LR_cov[8] = cov(2, 2)*psiMW*psiMW;*/
  }
  else{
    cout << "Process not defined well in Fit_eight" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* ComponentFuncts_eight(Double_t *pars, Double_t Mmin, Double_t Mmax,
			   TString process, Bool_t hIsUpS){
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");
  
  TF1 *JPsi = new TF1("f_JPsi", f_CrystalBall, Mmin, Mmax, 5);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
  
  TF1 *psi = new TF1("f_psi", f_CrystalBall, Mmin, Mmax, 5);
  psi->SetParameters(pars[5], pars[1]*psiMW, pars[2]*psiMW, pars[3], pars[4]);

  TF1 *Bg = new TF1("f_Bg", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  Bg->SetParameters(pars[6], pars[7], Mmin);
  
  TF1 *DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
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
    cout << "Wrong input process Fit_eight" << endl;
    exit(EXIT_FAILURE);
  }
}


void IntegrateLR_eight(TF1 *f, Double_t *pars, Double_t *LR_cov,
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
    Double_t DY_pars[] = {pars[8], pars[9], Mmin};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov);
  }
}


void IntegrateLR_eight(TF1 *f, Double_t *pars, Double_t *LR_cov,
		       Double_t LR_width, Double_t *LR, Double_t *e_LR,
		       Bool_t hIsUpS, TString process){
  //Integrate around mass value
  Double_t Mass = pars[1];
  if (process == "psi") {
    Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS"):Get_eight_Ratio("DownS");
    Mass *= psiMW;
  }
  Double_t LR_Mmin = Mass - LR_width;
  Double_t LR_Mmax = Mass + LR_width;
  
  (*LR) = f->Integral(LR_Mmin, LR_Mmax);
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, pars, LR_cov);
  }
  else if (process =="psi"){
    Double_t psi_pars[] = {pars[5], pars[1], pars[2], pars[3], pars[4]};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, psi_pars, LR_cov);
  }
  else {
    cout << "Process not defined for this IntegrateLR_eight" << endl;
  }
}


#endif
