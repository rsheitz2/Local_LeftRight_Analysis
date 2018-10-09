#ifndef FIT_TEN_H
#define FIT_TEN_H
//2 Crystal Balls for JPsi and psi'
//2 Exponentials for background and DY


Double_t Fit_ten(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[12];
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
    
  Double_t arg_psi = ( x[0] - par[6])/(par[7]);
  Double_t psi;
  if (arg_psi > -par[3] ) psi = par[5]*Norm*TMath::Exp(-0.5*arg_psi*arg_psi);
  else psi = par[5]*Norm*A*TMath::Power((B - arg_psi), -par[1]);

  Double_t CombBg = par[8]*TMath::Exp( par[9]*xShift );
  Double_t DY = par[10]*TMath::Exp( par[11]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


/*Double_t Fit_ten_scan(Double_t *x, Double_t *par){
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
}*/


void Paras_ten(TH1D *h, TF1 *fitFunc, Int_t nPar,
	       Double_t Mmin, Double_t Mmax){
  Int_t ipar =0;
  Double_t ratioPsi = 1.145;

  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(5.0) );
  Double_t DY_slope = -1.0, DY_Mmin = 5.0;
  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));
  A_Bg = A_Bg - A_DY;

  if (A_JPsi == 0) A_JPsi = 200.;
  if (A_psi == 0) A_psi = 200.;
  if (A_Bg <= 0) A_Bg = 100.;
  if (A_DY == 0) A_DY = 200.;
  
  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.12); ipar++;
  fitFunc->SetParameter(ipar, 0.16); ipar++;
  fitFunc->SetParameter(ipar, 1.4); ipar++;
  fitFunc->SetParameter(ipar, 3.2); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, 3.12*ratioPsi); ipar++;
  fitFunc->SetParameter(ipar, 0.16*ratioPsi); ipar++;
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
  Double_t factor =10.0;
  //Double_t factor =5.0;
  //fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  ipar++;
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 1.0, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, 2.5, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  //ipar++;
  fitFunc->SetParLimits(ipar, 3.4, 4.0); ipar++;
  fitFunc->SetParLimits(ipar, 0.1, 0.3); ipar++;
  fitFunc->SetParLimits(ipar, 0.0, A_Bg*factor); ipar++;//CombBg
  //ipar++;
  //fitFunc->SetParLimits(ipar, A_Bg/factor, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(ipar, -7.0, -2.0); ipar++;
  //ipar++;
  //fitFunc->SetParLimits(ipar, 0.0, A_DY*factor); ipar++;//DY
  //ipar++;
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0.0); ipar++;
  //ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
    }
}


/*void Paras_ten(TH1D *h, TF1 *fitFunc, Int_t nPar,
	       Double_t Mmin, Double_t Mmax){
  Int_t ipar =0;
  Double_t ratioPsi = 1.145;

  Double_t A_JPsi = h->GetBinContent(h->FindBin(3.12) );
  Double_t A_psi = h->GetBinContent(h->FindBin(3.12*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(5.0) );
  Double_t DY_slope = -1.0, DY_Mmin = 5.0;
  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));
  A_Bg = A_Bg - A_DY;

  if (A_JPsi == 0) A_JPsi = 200.;
  if (A_psi == 0) A_psi = 200.;
  if (A_Bg <= 0) A_Bg = 100.;
  if (A_DY == 0) A_DY = 200.;
  
  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, 3.12); ipar++;
  fitFunc->SetParameter(ipar, 0.16); ipar++;
  fitFunc->SetParameter(ipar, 1.4); ipar++;
  fitFunc->SetParameter(ipar, 3.2); ipar++;
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, 3.12*ratioPsi); ipar++;
  fitFunc->SetParameter(ipar, 0.16*ratioPsi); ipar++;
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
  Double_t factor =10.0;
  //Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 1.0, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, 2.5, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, 3.4, 4.0); ipar++;
  fitFunc->SetParLimits(ipar, 0.1, 0.3); ipar++;
  fitFunc->SetParLimits(ipar, 0.0, A_Bg*factor); ipar++;//CombBg
  //fitFunc->SetParLimits(ipar, A_Bg/factor, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(ipar, -7.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 0.0, A_DY*factor); ipar++;//DY
  //fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0.0); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  }//*/


/*void Paras_ten(TH1D *h, TF1 *fitFunc, Int_t nPar,
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
  }*/


TF1* SetupFunc_ten(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc,
		   Double_t Mmin, Double_t Mmax, Int_t *nPar){
  //Setup fit function
  *nPar =13;
  fitFunc = new TF1("fitFunc", Fit_ten, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_ten(h, fitFunc, *nPar, Mmin, Mmax);

  return fitFunc;
}


/*TF1* SetupFunc_ten(TH1D *h, TF1 *fitFunc, Double_t Mmin, Double_t Mmax,
  Int_t *nPar, Double_t psiMW){
  //Setup fit function for scanning psiMW
  *nPar =12;

  fitFunc = new TF1("fitFunc", Fit_ten_scan, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_ten(h, fitFunc, *nPar, Mmin, Mmax, psiMW);

  return fitFunc;
  }*/


void ProcessPars_ten(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		     TMatrixDSym &cov, TString process, Int_t nPars,
		     Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  
  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[5];//psi
  processPars[4] = pars[6];
  processPars[5] = pars[7];
  
  processPars[6] = pars[10];//DY
  processPars[7] = pars[11];
  
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
    cout << "covariance matrix doesn't work yet" << endl;
    exit(EXIT_FAILURE);
    /*LR_cov[0] = cov(5, 5);
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
      }*/
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
    cout << "Process not defined well in Fit_ten" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_ten(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		     TMatrixDSym &cov, TString process, Int_t nPars,
		     Bool_t hIsUpS, Double_t *e_processPars){
  ProcessPars_ten(fitFunc, processPars, LR_cov, cov, process, nPars, hIsUpS);
    
  const Double_t *e_pars = fitFunc->GetParErrors();

  e_processPars[0] = e_pars[0];//JPsi
  e_processPars[1] = e_pars[1];
  e_processPars[2] = e_pars[2];

  e_processPars[3] = e_pars[5];//psi
  e_processPars[4] = e_pars[6];
  e_processPars[5] = e_pars[7];

  e_processPars[6] = e_pars[10];//DY
  e_processPars[7] = e_pars[11];
}


void ProcessPars_ten(TF1 *fitFunc, Double_t *processPars,
		     Double_t *e_processPars, TString process, Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();

  if (process=="JPsi"){
    processPars[0] = pars[1]; //Mass
    processPars[1] = pars[2]; //Width

    e_processPars[0] = e_pars[1]; //Mass
    e_processPars[1] = e_pars[2]; //Width
  }
  else if (process=="psi"){
    processPars[0] = pars[5]; //Mass
    processPars[1] = pars[6]; //Width

    e_processPars[0] = e_pars[5]; //Mass
    e_processPars[1] = e_pars[6]; //Width
  }
  else if (process=="DY"){
    processPars[0] = pars[10]; //Amplitude
    processPars[1] = pars[11]; //Slope

    e_processPars[0] = e_pars[10]; //Amplitude
    e_processPars[1] = e_pars[11]; //Slope
  }
  else{
    cout << "Process not defined well in Fit_ten" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* ComponentFuncts_ten(Double_t *pars, Double_t Mmin, Double_t Mmax,
			 TString process, Bool_t hIsUpS){
  TF1 *JPsi = new TF1("f_JPsi", f_CrystalBall, Mmin, Mmax, 5);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
  
  TF1 *psi = new TF1("f_psi", f_CrystalBall, Mmin, Mmax, 5);
  psi->SetParameters(pars[5], pars[6], pars[7], pars[3], pars[4]);

  TF1 *Bg = new TF1("f_Bg", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  Bg->SetParameters(pars[8], pars[9], Mmin);
  
  TF1 *DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  DY->SetParameters(pars[10], pars[11], Mmin);

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
    cout << "Wrong input process Fit_ten" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* ComponentFuncts_ten(Double_t *pars, Double_t Mmin, Double_t Mmax,
			 TString process, Bool_t hIsUpS, TH1D *hRatio){
  //Draws component functions and
  //makes ratio of RD to fit based on component functions
  //Doesn't considered fit correlation errors in ratio calculation
  TF1 *JPsi = new TF1("f_JPsi", f_CrystalBall, Mmin, Mmax, 5);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
  
  TF1 *psi = new TF1("f_psi", f_CrystalBall, Mmin, Mmax, 5);
  psi->SetParameters(pars[5], pars[6], pars[7], pars[3], pars[4]);

  TF1 *Bg = new TF1("f_Bg", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  Bg->SetParameters(pars[8], pars[9], Mmin);
  
  TF1 *DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", Mmin, Mmax);
  DY->SetParameters(pars[10], pars[11], Mmin);

  JPsi->SetLineColor(kGreen); 
  psi->SetLineColor(kGreen);
  Bg->SetLineColor(kGreen); 
  DY->SetLineColor(kGreen); 

  JPsi->Draw("same");
  psi->Draw("same");
  Bg->Draw("same");
  DY->Draw("same");

  for (Int_t bi=1; bi<hRatio->GetNbinsX()+1; bi++) {
    Double_t RealData = hRatio->GetBinContent(bi);
    Double_t RDerror = hRatio->GetBinError(bi);

    Double_t nJPsi = JPsi->Eval(hRatio->GetBinCenter(bi) );
    Double_t npsi = psi->Eval(hRatio->GetBinCenter(bi) );
    Double_t nBg = Bg->Eval(hRatio->GetBinCenter(bi) );
    Double_t nDY = DY->Eval(hRatio->GetBinCenter(bi) );
    
    Double_t fitValue = nJPsi + npsi + nBg + nDY;
    Double_t fitError = TMath::Sqrt(fitValue);
    
    Double_t error = RatioError(RealData, fitValue, RDerror, fitError);
    Double_t ratio = RealData/fitValue;
    
    hRatio->SetBinContent(bi, ratio);
    hRatio->SetBinError(bi, error);
  }
  
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
    cout << "Wrong input process Fit_ten" << endl;
    exit(EXIT_FAILURE);
  }
}


void IntegrateLR_ten(TF1 *f, Double_t *pars, Double_t *LR_cov,
		     Double_t binWidth, Double_t LR_Mmin, Double_t LR_Mmax,
		     Double_t *LR, Double_t *e_LR, Double_t Mmin,
		     Double_t Mmax, TString process){
  (*LR) = f->Integral(LR_Mmin, LR_Mmax)/binWidth;
  
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, pars, LR_cov)/binWidth;
  }
  else if (process =="psi"){
    Double_t psi_pars[] = {pars[5], pars[6], pars[7], pars[3], pars[4]};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, psi_pars, LR_cov)/binWidth;
  }
  else if (process =="DY"){
    Double_t DY_pars[] = {pars[10], pars[11], Mmin};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov)/binWidth;
  }
}


void IntegrateLR_ten(TF1 *f, Double_t *pars, Double_t *LR_cov,
		     Double_t binWidth, Double_t LR_width, Double_t *LR,
		     Double_t *e_LR, Bool_t hIsUpS, Double_t Mmin,Double_t Mmax,
		     TString process){
  //Integrate around mass value
  Double_t Mass = (process=="JPsi") ? pars[1] : pars[6];
  
  Double_t LR_Mmin = Mass - LR_width;
  Double_t LR_Mmax = Mass + LR_width;

  if (process != "JPsi"  && process != "psi") {
    cout << "Process not defined for this IntegrateLR_ten" << endl;
    exit(EXIT_FAILURE);
  }

  IntegrateLR_ten(f, pars, LR_cov, binWidth, LR_Mmin, LR_Mmax, LR, e_LR,
		  Mmin, Mmax, process);
}


void IntegrateLRpercent_ten(TF1 *f, Double_t *pars, Double_t *LR_cov,
			    Double_t binWidth, Double_t widthFraction,
			    Double_t *LR, Double_t *e_LR, Bool_t hIsUpS,
			    Double_t Mmin,Double_t Mmax, TString process){
  //Integrate around mass value with an interval as a percentage of the width
  if (process != "JPsi"  && process != "psi") {
    cout << "Process not defined for this IntegrateLR_ten" << endl;
    exit(EXIT_FAILURE);
  }
  
  //Integrate around mass value
  Double_t Mass = (process=="JPsi") ? pars[1] : pars[6];
  Double_t Width = (process=="JPsi") ? pars[2] : pars[7];
  
  Double_t LR_Mmin = Mass - Width*widthFraction;
  Double_t LR_Mmax = Mass + Width*widthFraction;

  IntegrateLR_ten(f, pars, LR_cov, binWidth, LR_Mmin, LR_Mmax, LR, e_LR,
		  Mmin, Mmax, process);
}


void IntegrateLR_ten(TF1 *f, TH1D *h, Double_t *pars, Double_t LR_width,
		     Double_t *LR, Double_t *e_LR, Int_t bin, TString process){
  //Define process result histogram
  Double_t multiNbins =100.0;
  Int_t nBins = h->GetNbinsX(); nBins = nBins*multiNbins;
  Double_t Mmin = h->GetXaxis()->GetXmin();
  Double_t Mmax = h->GetXaxis()->GetXmax();
  TH1D *hLR = new TH1D("hLR", "hLR", nBins, Mmin, Mmax);
  for (Int_t bi=1; bi<nBins+1; bi++) {
    Double_t xval = hLR->GetBinCenter(bi);
    Double_t f_val =f->Eval(xval);
    Double_t eFval = TMath::Sqrt(f_val);
    
    hLR->SetBinContent(bi, f_val);
    hLR->SetBinError(bi, eFval);
  }
  
  //Integrate around center value of histogram
  Double_t center;
  if (process =="JPsi") center = pars[1];
  else if (process =="psi") center = pars[6];
  else {
    cout << "Process not defined well in Fit_ten" << endl;
    exit(EXIT_FAILURE);
  }
  Double_t LR_Mmin = center - LR_width;
  Double_t LR_Mmax = center + LR_width;
  Int_t lowerBin = hLR->FindBin(LR_Mmin);
  Int_t upperBin = hLR->FindBin(LR_Mmax);
  LR[bin] = hLR->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  LR[bin] /= multiNbins;
  e_LR[bin] /= TMath::Sqrt(multiNbins);

  delete hLR;
}


#endif
