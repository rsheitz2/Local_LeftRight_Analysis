#ifndef FIT_EIGHT_H
#define FIT_EIGHT_H
//2 Crystal Balls for JPsi and psi'
//      psi' M/W = Get_eight_Ratio(targ)*JPsi M/W
//      psi' alpha/n parameters are the same as JPsi
//2 Exponentials for background and DY


Double_t Get_eight_Ratio(TString targ){
  if (targ =="UpS") return 1.155;
  else if (targ =="DownS") return 1.135;
  else{
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t _CrystalBall(Double_t *x, Double_t *par){
  //f(x; alpha, n, xBar, sigma, Amp)
  //   Amp = par[0], xBar = par[1], sigma = par[2]
  //   alpha = par[3], n = par[4]
  Double_t A = TMath::Power(par[4]/par[3], par[4])*
    TMath::Exp(-par[3]*par[3]/2.0);
  
  Double_t B = par[4]/par[3] - par[3];
  
  Double_t C = (par[4]/par[3])*(1.0/(par[4]-1.0))*
    TMath::Exp(-par[3]*par[3]/2.0);
  
  Double_t D =TMath::Sqrt( TMath::Pi()/2.0 )*
  (1+TMath::Erf( par[3]/TMath::Sqrt(2) ) );

  Double_t Norm = 1.0/(par[2]*(C+D) );

  Double_t arg = (x[0] - par[1])/par[2];
  if (arg > -par[3] ) return par[0]*Norm*TMath::Exp(-0.5*arg*arg);
  else return par[0]*Norm*A*TMath::Power((B - arg), -par[4]);
}


Double_t Fit_eight_UpS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_eight_Ratio("UpS");

  Double_t par_JPsi[] = {par[0], par[1], par[2], par[3], par[4]};
  Double_t JPsi = _CrystalBall(x, par_JPsi);
  Double_t par_psi[] = {par[5], par[1]*ratioPsi, par[2]*ratioPsi,par[3],par[4]};
  Double_t psi = _CrystalBall(x, par_psi);

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_eight_DownS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_eight_Ratio("DownS");

  Double_t par_JPsi[] = {par[0], par[1], par[2], par[3], par[4]};
  Double_t JPsi = _CrystalBall(x, par_JPsi);
  Double_t par_psi[] = {par[5], par[1]*ratioPsi, par[2]*ratioPsi,par[3],par[4]};
  Double_t psi = _CrystalBall(x, par_psi);

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_eight_scan(Double_t *x, Double_t *par){
  //Feed one scan parameter as par[11]=psiMW
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = par[11];

  Double_t par_JPsi[] = {par[0], par[1], par[2], par[3], par[4]};
  Double_t JPsi = _CrystalBall(x, par_JPsi);
  Double_t par_psi[] = {par[5], par[1]*ratioPsi, par[2]*ratioPsi,par[3],par[4]};
  Double_t psi = _CrystalBall(x, par_psi);

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_eight(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax){
  //Older parameter setup (no longer used)
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
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 3.0, 3.6); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.2); ipar++;
  fitFunc->SetParLimits(ipar, 1.0, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, 2.5, 8.0); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, 0.0, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(ipar, -7.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 0.0, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0.0); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


void Paras_eight(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax, Bool_t hIsUpS){
  //Newer (currently used) parameter setup
  Int_t ipar =0;
  Double_t ratioPsi = (hIsUpS) ?
    Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");

  Double_t M_JPsi = (hIsUpS) ? 3.13 : 3.12;
  Double_t A_JPsi = h->GetBinContent(h->FindBin(M_JPsi) );
  Double_t A_psi = h->GetBinContent(h->FindBin(M_JPsi*ratioPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  //Double_t DY_slope = -0.95, DY_Mmin = 4.5;
  Double_t DY_slope = -0.90, DY_Mmin = 4.5;
  Double_t A_DY =  h->GetBinContent(h->FindBin(DY_Mmin) );
  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));
  A_Bg = A_Bg - A_DY;

  if (A_JPsi == 0) A_JPsi = 200.;
  if (A_psi == 0) A_psi = 200.;
  if (A_Bg <= 0) A_Bg = 100.;
  if (A_DY == 0) A_DY = 200.;
  
  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, M_JPsi); ipar++;
  fitFunc->SetParameter(ipar, 0.16); ipar++;
  fitFunc->SetParameter(ipar, 1.7); ipar++;
  //fitFunc->SetParameter(ipar, 3.2); ipar++;
  fitFunc->SetParameter(ipar, 8.0); ipar++;
  
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_Bg); ipar++;//CombBg
  fitFunc->SetParameter(ipar, -3.0); ipar++;
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  fitFunc->SetParameter(ipar, DY_slope); ipar++;
  fitFunc->SetParameter(ipar, Mmin); ipar++;//Mmin
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }

  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(0, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(1, 2.8, 3.5); ipar++;
  fitFunc->SetParLimits(2, 0, 0.25); ipar++;
  fitFunc->SetParLimits(3, 0.5, 2.0); ipar++;
  fitFunc->SetParLimits(4, 0.5, 20.0); ipar++;
  fitFunc->SetParLimits(5, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(6, 0.0, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(7, -7.0, -1.0); ipar++;
  fitFunc->SetParLimits(8, 0.0, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(9, -2.0, 0.0); ipar++;
  fitFunc->SetParLimits(10, Mmin, Mmin); ipar++;//Mmin
  
  /*fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, 2.8, 3.5); ipar++;
  fitFunc->SetParLimits(ipar, 0, 0.25); ipar++;
  fitFunc->SetParLimits(ipar, 1.0, 2.0); ipar++;
  fitFunc->SetParLimits(ipar, 2.5, 4.5); ipar++;
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, 0.0, A_Bg*factor); ipar++;//CombBg
  fitFunc->SetParLimits(ipar, -4.0, -2.0); ipar++;
  fitFunc->SetParLimits(ipar, 0.0, A_DY*factor); ipar++;//DY
  fitFunc->SetParLimits(ipar, -2.0, 0.0); ipar++;
  fitFunc->SetParLimits(ipar, Mmin, Mmin); ipar++;//Mmin
  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }//*/
}


void Paras_eight(TH1D *h, TF1 *fitFunc, Int_t nPar,
		 Double_t Mmin, Double_t Mmax, Double_t psiMW){
  //Nominally used for setting up to scan psiMW value
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
  Paras_eight(h, fitFunc, *nPar, Mmin, Mmax, hIsUpS); //Newer more accurate pars

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
  //Gets some physical parameters, sets up cov matrix
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
  }
  else{
    cout << "Process not defined well in Fit_eight" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_eight(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
		       TMatrixDSym &cov, TString process, Int_t nPars,
		       Bool_t hIsUpS, Double_t *e_processPars){
  //Gets some physical parameters, sets up cov matrix
  //    and gets errors of physical parameters
  ProcessPars_eight(fitFunc, processPars, LR_cov, cov, process,nPars,hIsUpS);
  const Double_t *e_pars = fitFunc->GetParErrors();
  
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");

  e_processPars[0] = e_pars[0];//JPsi
  e_processPars[1] = e_pars[1];
  e_processPars[2] = e_pars[2];

  e_processPars[3] = e_pars[5];//psi
  e_processPars[4] = e_pars[1]*psiMW;
  e_processPars[5] = e_pars[2]*psiMW;

  e_processPars[6] = e_pars[8];//DY
  e_processPars[7] = e_pars[9];
}


void ProcessPars_eight(TF1 *fitFunc, Double_t *processPars,
		       Double_t *e_processPars, TString process, Bool_t hIsUpS){
  //Used to monitor some physical parameters
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  
  Double_t psiMW;
  if (hIsUpS) psiMW = Get_eight_Ratio("UpS");
  else psiMW = Get_eight_Ratio("DownS");

  if (process=="JPsi"){
    processPars[0] = pars[1]; //Mass
    processPars[1] = pars[2]; //Width

    e_processPars[0] = e_pars[1]; //Mass
    e_processPars[1] = e_pars[2]; //Width
  }
  else if (process=="psi"){
    processPars[0] = pars[1]*psiMW; //Mass
    processPars[1] = pars[2]*psiMW; //Width

    e_processPars[0] = e_pars[1]*psiMW; //Mass
    e_processPars[1] = e_pars[2]*psiMW; //Width
  }
  else if (process=="DY"){
    processPars[0] = pars[8]; //Amplitude
    processPars[1] = pars[9]; //Slope

    e_processPars[0] = e_pars[8]; //Amplitude
    e_processPars[1] = e_pars[9]; //Slope
  }
  else{
    cout << "Process not defined well in Fit_eight" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_eight(TF1 *fitFunc, vector<vector<Double_t> > &processPars,
		       vector<vector<Double_t> > &e_processPars, Bool_t hIsUpS){
  //Used to monitor all output paramters
  //For macro functParas.C
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();

  for (Int_t i=0; i<10; i++) {
    processPars[i].push_back(pars[i]);
    e_processPars[i].push_back(e_pars[i]);
  }
}


void ProcessPhysicsPars_eight(TF1 *fitFunc, Double_t *processPars,
			      Double_t *e_processPars, TString process,
			      Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();

  Double_t psiMW;
  if (hIsUpS) psiMW = Get_eight_Ratio("UpS");
  else psiMW = Get_eight_Ratio("DownS");

  processPars[0] = pars[1]; //JPsi Mass
  processPars[1] = pars[2]; //JPsi Width
  processPars[2] = psiMW*pars[1]; //psi' Mass
  processPars[3] = psiMW*pars[2]; //psi' Width
  processPars[4] = pars[9]; //DY Slope
    
  e_processPars[0] = e_pars[1]; //JPsi Mass
  e_processPars[1] = e_pars[2]; //JPsi Width
  e_processPars[2] = psiMW*e_pars[1]; //psi' Mass
  e_processPars[3] = psiMW*e_pars[2]; //psi' Width
  e_processPars[4] = e_pars[9]; //DY Slope
}


TF1* ComponentFuncts_eight(Double_t *pars, Double_t Mmin, Double_t Mmax,
			   TString process, Bool_t hIsUpS){
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");
  
  TF1 *JPsi = new TF1("f_JPsi", _CrystalBall, Mmin, Mmax, 5);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
  
  TF1 *psi = new TF1("f_psi", _CrystalBall, Mmin, Mmax, 5);
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


void IntegrateLR_eight(TF1 *f, Double_t binWidth,
		       Double_t LR_Mmin, Double_t LR_Mmax,
		       Double_t *LR, Double_t *e_LR){
  //Integration without covariance matrix taken into account
  (*LR) = f->Integral(LR_Mmin, LR_Mmax)/binWidth;
  (*e_LR) = TMath::Sqrt(*LR);
}


void IntegrateLR_eight(TF1 *f, Double_t *pars, Double_t *LR_cov,
		       Double_t binWidth, Double_t LR_Mmin, Double_t LR_Mmax,
		       Double_t *LR, Double_t *e_LR, Double_t Mmin,
		       Double_t Mmax, Bool_t hIsUpS, TString process){
  //Integration considering correlation errors
  (*LR) = f->Integral(LR_Mmin, LR_Mmax)/binWidth;
  
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, pars, LR_cov)/binWidth;
  }
  else if (process =="psi"){
    Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS"):Get_eight_Ratio("DownS");
    Double_t psi_pars[] =
      {pars[5], psiMW*pars[1], psiMW*pars[2], pars[3], pars[4]};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, psi_pars, LR_cov)/binWidth;
  }
  else if (process =="DY"){
    Double_t DY_pars[] = {pars[8], pars[9], Mmin};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov)/binWidth;
  }
}


void IntegrateLR_eight(TF1 *f, Double_t *pars, Double_t *LR_cov,
		       Double_t binWidth, Double_t LR_width, Double_t *LR,
		       Double_t *e_LR, Bool_t hIsUpS, TString process){
  //Integration considering correlation errors
  //Integrate around mass value
  Double_t Mass = pars[1];
  Double_t psiMW;
  if (process == "psi") {
    psiMW = (hIsUpS) ? Get_eight_Ratio("UpS"):Get_eight_Ratio("DownS");
    Mass *= psiMW;
  }
  Double_t LR_Mmin = Mass - LR_width;
  Double_t LR_Mmax = Mass + LR_width;
  
  (*LR) = f->Integral(LR_Mmin, LR_Mmax)/binWidth;
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, pars, LR_cov)/binWidth;
  }
  else if (process =="psi"){
    Double_t psi_pars[] =
      {pars[5], psiMW*pars[1], psiMW*pars[2], pars[3], pars[4]};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, psi_pars, LR_cov)/binWidth;
  }
  else {
    cout << "Process not defined for this IntegrateLR_eight" << endl;
  }
}


void IntegrateLRpercent_eight(TF1 *f, Double_t *pars, Double_t *LR_cov,
			      Double_t binWidth, Double_t widthFraction,
			      Double_t *LR, Double_t *e_LR, Bool_t hIsUpS,
			      Double_t Mmin,Double_t Mmax, TString process){
  //Integrate around mass value with an interval as a percentage of the width
  if (process != "JPsi"  && process != "psi") {
    cout << "Process not defined for this IntegrateLR_ten" << endl;
    exit(EXIT_FAILURE);
  }

  //Integrate around mass value
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");
  Double_t Mass = (process=="JPsi") ? pars[1] : psiMW*pars[1];
  Double_t Width = (process=="JPsi") ? pars[2] : psiMW*pars[2];
  
  Double_t LR_Mmin = Mass - Width*widthFraction;
  Double_t LR_Mmax = Mass + Width*widthFraction;
  
  IntegrateLR_eight(f, pars, LR_cov, binWidth, LR_Mmin, LR_Mmax, LR, e_LR,
		    Mmin, Mmax, hIsUpS, process);
}


void IntegrateLRpercent_eight(TF1 *f, Double_t *pars, Double_t binWidth,
 			      Double_t widthFraction, Double_t *LR,
			      Double_t *e_LR, Bool_t hIsUpS, TString process){
  //Integrate around mass value with an interval as a percentage of the width
  //Correlation errors not taken into account
  if (process != "JPsi"  && process != "psi") {
    cout << "Process not defined for this IntegrateLR_ten" << endl;
    exit(EXIT_FAILURE);
  }

  //Integrate around mass value
  Double_t psiMW = (hIsUpS) ? Get_eight_Ratio("UpS") : Get_eight_Ratio("DownS");
  Double_t Mass = (process=="JPsi") ? pars[1] : psiMW*pars[1];
  Double_t Width = (process=="JPsi") ? pars[2] : psiMW*pars[2];
  
  Double_t LR_Mmin = Mass - Width*widthFraction;
  Double_t LR_Mmax = Mass + Width*widthFraction;
  
  IntegrateLR_eight(f, binWidth, LR_Mmin, LR_Mmax, LR, e_LR);
}


#endif
