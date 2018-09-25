#ifndef FIT_SEVEN_H
#define FIT_SEVEN_H
//2 Gaussians for JPsi and psi'
//      psi' M/W = Get_seven_Ratio(targ)*JPsi M/W(different constraint per target)
//1 Exponentials for background and DY

Double_t Get_seven_Ratio(TString targ){
  Double_t ratioPsi_UpS_2exp = 1.15;
  Double_t ratioPsi_DownS_2exp = 1.135;
  
  if (targ=="UpS") return ratioPsi_UpS_2exp;
  else if (targ=="DownS") return ratioPsi_DownS_2exp;
  else {
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t Fit_seven_UpS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[6];
  Double_t ratioPsi_exp = Get_seven_Ratio("UpS");
    
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi =( x[0] - par[1]*ratioPsi_exp)/(par[2]*ratioPsi_exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[4]*TMath::Exp( par[5]*xShift );
        
  return JPsi + psi + DY;
}


Double_t Fit_seven_DownS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[6];
  Double_t ratioPsi_exp = Get_seven_Ratio("DownS");
    
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi =( x[0] - par[1]*ratioPsi_exp)/(par[2]*ratioPsi_exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[4]*TMath::Exp( par[5]*xShift );
    
  return JPsi + psi + DY;
}


Double_t Fit_seven_scan(Double_t *x, Double_t *par){
  //Feed one scan parameter as par[7]=psiMW
  Double_t xShift = x[0]-par[6];
  Double_t ratioPsi_exp = par[7];
    
  Double_t arg_JPsi = ( x[0] - par[1] )/par[2];
  Double_t JPsi = par[0]*TMath::Exp( -0.5*arg_JPsi*arg_JPsi );
  Double_t arg_psi =( x[0] - par[1]*ratioPsi_exp)/(par[2]*ratioPsi_exp);
  Double_t psi = par[3]*TMath::Exp( -0.5*arg_psi*arg_psi );

  Double_t DY = par[4]*TMath::Exp( par[5]*xShift );
    
  return JPsi + psi + DY;
}


void Paras_seven(TH1D *h, TF1 *fitFunc, Int_t nPar, Double_t Mmin, Double_t Mmax){
  Int_t ipar =0;
  Double_t ratioPsi = Get_seven_Ratio("UpS");
  ratioPsi += Get_seven_Ratio("DownS");
  ratioPsi /= 2.0;
  
  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(Mmin) );

  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
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
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*5.0); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -3.0, -0.5); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void Paras_seven(TH1D *h, TF1 *fitFunc, Int_t nPar, Double_t Mmin,
		 Double_t Mmax, Double_t psiMW){
  Int_t ipar =0;
  Double_t ratioPsi =psiMW;
  
  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(Mmin) );

  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.1); ipar++;
  fitFunc->SetParameter(ipar, 0.15); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, -0.84); ipar++;
  fitFunc->SetParameter(ipar, Mmin); ipar++;//Mmin
  fitFunc->SetParameter(ipar, psiMW); ipar++;//psiMW
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  
  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*5.0); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -3.0, -0.5); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  fitFunc->SetParLimits(ipar, psiMW, psiMW); ipar++;//psiMW
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* SetupFunc_seven(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc,
		     Double_t Mmin, Double_t Mmax, Int_t *nPar){
  //Setup fit function
  *nPar =7;

  if (hIsUpS)
    fitFunc = new TF1("fitFunc", Fit_seven_UpS, Mmin, Mmax, *nPar);
  else
    fitFunc = new TF1("fitFunc", Fit_seven_DownS, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_seven(h, fitFunc, *nPar, Mmin, Mmax);

  return fitFunc;
}


TF1* SetupFunc_seven(TH1D *h, TF1 *fitFunc, Double_t Mmin, Double_t Mmax,
		     Int_t *nPar, Double_t psiMW){
  //Setup fit function for scanning psiMW
  *nPar =8;

  fitFunc = new TF1("fitFunc", Fit_seven_scan, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_seven(h, fitFunc, *nPar, Mmin, Mmax, psiMW);

  return fitFunc;
}


void ProcessPars_seven(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		       TMatrixDSym &cov, TString process, Int_t nPars,
		       Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  
  Double_t psiMW;
  if (hIsUpS) psiMW = Get_seven_Ratio("UpS");
  else psiMW = Get_seven_Ratio("DownS");

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[3];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[4];//DY
  processPars[7] = pars[5];
  
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
    LR_cov[0] = cov(4, 4);
    LR_cov[1] = cov(4, 5);
    LR_cov[2] = 0.0;
        
    LR_cov[3] = cov(5, 4);
    LR_cov[4] = cov(5, 5);
    LR_cov[5] = 0.0;

    LR_cov[6] = 0.0;
    LR_cov[7] = 0.0;
    LR_cov[8] = 0.0;
  }
  else{
    cout << "Process not defined well in Fit_seven" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_seven(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		       TMatrixDSym &cov, TString process, Int_t nPars,
		       Bool_t hIsUpS, Double_t *e_processPars){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  
  Double_t psiMW;
  if (hIsUpS) psiMW = Get_seven_Ratio("UpS");
  else psiMW = Get_seven_Ratio("DownS");

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[3];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
  processPars[6] = pars[4];//DY
  processPars[7] = pars[5];

  e_processPars[0] = e_pars[0];//JPsi
  e_processPars[1] = e_pars[1];
  e_processPars[2] = e_pars[2];

  e_processPars[3] = e_pars[3];//psi
  e_processPars[4] = e_pars[1]*psiMW;
  e_processPars[5] = e_pars[2]*psiMW;

  e_processPars[6] = e_pars[4];//DY
  e_processPars[7] = e_pars[5];
  
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
    LR_cov[0] = cov(4, 4);
    LR_cov[1] = cov(4, 5);
    LR_cov[2] = 0.0;
        
    LR_cov[3] = cov(5, 4);
    LR_cov[4] = cov(5, 5);
    LR_cov[5] = 0.0;

    LR_cov[6] = 0.0;
    LR_cov[7] = 0.0;
    LR_cov[8] = 0.0;
  }
  else{
    cout << "Process not defined well in Fit_seven" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* ComponentFuncts_seven(Double_t *processPars, Double_t Mmin, Double_t Mmax,
			   TString process){
  TF1 *JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*(([1]-x)*([1]-x))/([2]*[2]))",
	    Mmin, Mmax);
  JPsi->SetParameters(processPars[0], processPars[1], processPars[2]);

  TF1 *psi =
    new TF1("f_psi", "[0]*TMath::Exp(-0.5*( ([1]-x)*([1]-x) )/( [2]*[2]) )",
	    Mmin, Mmax);
  psi->SetParameters(processPars[3], processPars[4], processPars[5]);

  TF1 *DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  DY->SetParameters(processPars[6], processPars[7], Mmin);

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
    cout << "Wrong input process Fit_seven" << endl;
    exit(EXIT_FAILURE);
  }
}


void IntegrateLR_seven(TF1 *f, Double_t *processPars, Double_t *LR_cov,
		       Double_t LR_Mmin, Double_t LR_Mmax,
		       Double_t *LR, Double_t *e_LR, Double_t Mmin, Double_t Mmax,
		       TString process){

  (*LR) = f->Integral(LR_Mmin, LR_Mmax);
  
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, processPars, LR_cov);
  }
  else if (process =="psi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, &processPars[3], LR_cov);
  }
  else if (process =="DY"){
    Double_t DY_pars[] = {processPars[6], processPars[7], Mmin};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov);
  }
  
}


#endif
