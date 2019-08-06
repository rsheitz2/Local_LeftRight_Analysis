#ifndef FIT_THIRTEEN_H
#define FIT_THIRTEEN_H
//2 Crystal Balls for JPsi and psi'
//      psi' M/W/alpha = Get_thirteen_Ratio(targ)*JPsi M/W
//      psi' n parameter is the same as JPsi
//2 Exponentials for background and DY


Double_t Get_thirteen_Ratio(TString targ){
  if (targ =="UpS") return 1.1485;
  else if (targ =="DownS") return 1.1418;
  else{
    cout << "Wrong input target" << endl;
    exit(EXIT_FAILURE);
  }
}


Double_t thirteen_CrystalBall(Double_t *x, Double_t *par){
  //f(x; alpha, n, xBar, sigma, Amp)
  //   Amp = par[0], xBar = par[1], sigma = par[2]
  //   alpha = par[3], n = par[4]
  Double_t A = TMath::Power(par[4]/par[3], par[4])*
    TMath::Exp(-par[3]*par[3]/2.0);
  
  Double_t B = par[4]/par[3] - par[3];

  Double_t C
    =(par[4]/par[3])*(1.0/(par[4]-1.0))*TMath::Exp(-par[3]*par[3]/2.0);
  
  Double_t D
    =TMath::Sqrt( TMath::Pi()/2.0 )*(1+TMath::Erf( par[3]/TMath::Sqrt(2) ) );

  Double_t Norm = 1.0/(par[2]*(C+D) );

  Double_t arg = (x[0] - par[1])/par[2];
  if (arg > -par[3] ) return par[0]*Norm*TMath::Exp(-0.5*arg*arg);
  else return par[0]*Norm*A*TMath::Power((B - arg), -par[4]);
}


Double_t Fit_thirteen_UpS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_thirteen_Ratio("UpS");

  Double_t par_JPsi[] = {par[0], par[1], par[2], par[3], par[4]};
  Double_t JPsi = thirteen_CrystalBall(x, par_JPsi);
  
  Double_t par_psi[] =
    {par[5], par[1]*ratioPsi, par[2]*ratioPsi, par[3]*ratioPsi, par[4]};
  Double_t psi = thirteen_CrystalBall(x, par_psi);

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


Double_t Fit_thirteen_DownS(Double_t *x, Double_t *par){
  Double_t xShift = x[0]-par[10];
  Double_t ratioPsi = Get_thirteen_Ratio("DownS");

  Double_t par_JPsi[] = {par[0], par[1], par[2], par[3], par[4]};
  Double_t JPsi = thirteen_CrystalBall(x, par_JPsi);
  
  Double_t par_psi[] =
    {par[5], par[1]*ratioPsi, par[2]*ratioPsi, par[3]*ratioPsi, par[4]};
  Double_t psi = thirteen_CrystalBall(x, par_psi);

  Double_t CombBg = par[6]*TMath::Exp( par[7]*xShift );
  Double_t DY = par[8]*TMath::Exp( par[9]*xShift );
    
  return CombBg + JPsi + psi + DY;
}


void Paras_thirteen(TH1D *h, TF1 *fitFunc, Int_t nPar, Double_t Mmin,
		    Double_t Mmax, Bool_t hIsUpS, Int_t bin,TString physBinned){
  //Normal parameter setup
  //     Fixed parameter setup
  Double_t massJPsi =3.131, widthJPsi=0.17;
  Double_t alphaJPsi =1.5, nJPsi =1.5;
  Double_t Bg_slope =-2.0;
  Double_t DY_slope = -0.90, DY_Mmin = 4.5;

  Double_t ratioPsi = (hIsUpS) ?
    Get_thirteen_Ratio("UpS") : Get_thirteen_Ratio("DownS");
  
  Double_t A_JPsi = h->GetBinContent(h->FindBin(massJPsi) );
  Double_t A_Bg =  h->GetBinContent(h->FindBin(Mmin) );
  Double_t A_DY =  h->GetBinContent(h->FindBin(DY_Mmin) );
  
  A_DY /= (TMath::Exp(DY_slope*(DY_Mmin-Mmin)));
  A_Bg = A_Bg - A_DY;

  Double_t C 
    =(nJPsi/alphaJPsi)*(1.0/(nJPsi-1.0))*TMath::Exp(-alphaJPsi*alphaJPsi/2.0);
  Double_t D
    =TMath::Sqrt( TMath::Pi()/2.0 )*(1+TMath::Erf( alphaJPsi/TMath::Sqrt(2) ) );
  Double_t Norm = 1.0/(widthJPsi*(C+D) );
  A_JPsi /= Norm;

  Double_t A_psi = A_JPsi/30.0;

  if (A_JPsi == 0) A_JPsi = 200.;
  if (A_psi < 10) A_psi = 200.;
  if (A_Bg <= 0) A_Bg = 100.;
  if (A_DY == 0) A_DY = 200.;

  fitFunc->SetParameter(0, A_JPsi); //JPsi
  fitFunc->SetParameter(1, massJPsi); 
  fitFunc->SetParameter(2, widthJPsi); 
  fitFunc->SetParameter(3, alphaJPsi); 
  fitFunc->SetParameter(4, nJPsi); 
  fitFunc->SetParameter(5, A_psi); //psi
  fitFunc->SetParameter(6, A_Bg); //CombBg
  fitFunc->SetParameter(7, Bg_slope); 
  fitFunc->SetParameter(8, A_DY); //DY
  fitFunc->SetParameter(9, DY_slope); 
  fitFunc->SetParameter(10, Mmin); //Mmin

  //Constraints
  Double_t factor =100.0;//pT setup
  fitFunc->SetParLimits(0, 10, A_JPsi*factor); //A_JPsi
  fitFunc->SetParLimits(4, 0., 10.0); //nJPsi
  fitFunc->SetParLimits(5, 1, A_psi*factor); //A_psi
  fitFunc->SetParLimits(7, -7.0, -1.2); //Bg_slope
  fitFunc->SetParLimits(9, -2.0, 0.0); //DY_slope

  fitFunc->SetParLimits(10, Mmin, Mmin);//Mmin*/
}


TF1* SetupFunc_thirteen(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc, Double_t Mmin,
			Double_t Mmax, Int_t *nPar, Int_t bin,
			TString physBinned){
  //Setup fit function
  *nPar =11;

  if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_thirteen_UpS, Mmin, Mmax, *nPar);
  else fitFunc = new TF1("fitFunc", Fit_thirteen_DownS, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_thirteen(h, fitFunc, *nPar, Mmin, Mmax, hIsUpS, bin, physBinned);

  return fitFunc;
}


void Basic_thirteenChecks(Double_t *pars){
  Int_t failure = 0;
  if (pars[0] < 0 ) {
    cout << pars[0] << " A_JPsi" << endl;
    failure++;
  }
  if ( (pars[1] < 2.9) || (pars[1] > 3.3) ) {
    cout << pars[1] << " massJPsi" << endl;
    failure++;
  }
  if ( (pars[2] < 0) || (pars[2] > 0.4) ) {
    cout << pars[2] << " widthJPsi" << endl;
    failure++;
  }
  if (pars[5] < 0){
    cout << pars[5] << " A_psi" << endl;
    failure++;
  }
  if (pars[6] < 0){
    cout << pars[6] << " A_Bg" << endl;
    failure++;
  }
  if (pars[8] < 0){
    cout << pars[9] << " A_DY" << endl;
    failure++;
  }

  if (failure){
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_thirteen(TF1 *fitFunc, Double_t *processPars,Double_t *LR_cov,
			  TMatrixDSym &cov, TString process, Int_t nPars,
			  Bool_t hIsUpS){
  //Gets some physical parameters, sets up cov matrix
  Double_t *pars = fitFunc->GetParameters();
  Double_t psiMW = (hIsUpS) ? Get_thirteen_Ratio("UpS") : Get_thirteen_Ratio("DownS");
  Basic_thirteenChecks(pars);

  processPars[0] = pars[0];//JPsi
  processPars[1] = pars[1];
  processPars[2] = pars[2];

  processPars[3] = pars[5];//psi
  processPars[4] = pars[1]*psiMW;
  processPars[5] = pars[2]*psiMW;
  
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
    cout << "covariance matrix probably doesn't work yet" << endl;
    exit(EXIT_FAILURE);
    LR_cov[0] = cov(5, 5);
    LR_cov[5] = cov(1, 5)*psiMW;
    LR_cov[10] = cov(2, 5)*psiMW;
    LR_cov[15] = cov(6, 5);
    LR_cov[20] = cov(7, 5);
    
    for (Int_t i=1; i<5; i++) {
      Double_t multi = (i == 1) ? psiMW : 1.0;
      
      LR_cov[i] = cov(5, i)*multi;
      LR_cov[i+5] = cov(1, i)*psiMW*multi;
      LR_cov[i+10] = cov(2, i)*psiMW*multi;
      LR_cov[i+15] = cov(6, i)*multi;
      LR_cov[i+20] = cov(7, i)*multi;
    }
  }
  else if (process=="DY"){
    cout << "covariance matrix doesn't work yet" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "Process not defined well in Fit_thirteen" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPars_thirteen(TF1 *fitFunc, Double_t *processPars,
			  Double_t *e_processPars, TString process,
			  Bool_t hIsUpS){
  //Used to monitor some physical parameters
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  Basic_thirteenChecks(pars);
  
  Double_t psiMW;
  if (hIsUpS) psiMW = Get_thirteen_Ratio("UpS");
  else psiMW = Get_thirteen_Ratio("DownS");

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
    cout << "Process not defined well in Fit_thirteen" << endl;
    exit(EXIT_FAILURE);
  }
}


void ProcessPhysicsPars_thirteen(TF1 *fitFunc, Double_t *processPars,
  Double_t *e_processPars, TString process,
  Bool_t hIsUpS){
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  Basic_thirteenChecks(pars);

  Double_t psiMW;
  if (hIsUpS) psiMW = Get_thirteen_Ratio("UpS");
  else psiMW = Get_thirteen_Ratio("DownS");

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


TF1* ComponentFuncts_thirteen(Double_t *pars, Double_t Mmin, Double_t Mmax,
			      TString process, Bool_t hIsUpS){
  Double_t psiMW = (hIsUpS) ? Get_thirteen_Ratio("UpS") : Get_thirteen_Ratio("DownS");
  
  TF1 *JPsi = new TF1("f_JPsi", thirteen_CrystalBall, Mmin, Mmax, 5);
  JPsi->SetParameters(pars[0], pars[1], pars[2], pars[3], pars[4]);
  
  TF1 *psi = new TF1("f_psi", thirteen_CrystalBall, Mmin, Mmax, 5);
  psi->SetParameters(pars[5], pars[1]*psiMW, pars[2]*psiMW, pars[3]*psiMW, pars[4]);

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
    cout << "Wrong input process Fit_thirteen" << endl;
    exit(EXIT_FAILURE);
  }
}


void IntegrateLR_thirteen(TF1 *f, Double_t binWidth,
			  Double_t LR_Mmin, Double_t LR_Mmax,
			  Double_t *LR, Double_t *e_LR){
  //Integration without covariance matrix taken into account
  (*LR) = f->Integral(LR_Mmin, LR_Mmax)/binWidth;
  (*e_LR) = TMath::Sqrt(*LR);
}


void IntegrateLR_thirteen(TF1 *f, Double_t *pars, Double_t *LR_cov,
			  Double_t binWidth, Double_t LR_Mmin, Double_t LR_Mmax,
			  Double_t *LR, Double_t *e_LR, Double_t Mmin,
			  Double_t Mmax, Bool_t hIsUpS, TString process){
  //Integration considering correlation errors
  (*LR) = f->Integral(LR_Mmin, LR_Mmax)/binWidth;
  
  if (process =="JPsi"){
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, pars, LR_cov)/binWidth;
  }
  else if (process =="psi"){
    Double_t psiMW = (hIsUpS) ? Get_thirteen_Ratio("UpS"):Get_thirteen_Ratio("DownS");
    Double_t psi_pars[] =
      {pars[5], psiMW*pars[1], psiMW*pars[2], psiMW*pars[3], pars[4]};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, psi_pars, LR_cov)/binWidth;
  }
  else if (process =="DY"){
    Double_t DY_pars[] = {pars[8], pars[9], Mmin};
    (*e_LR) = f->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov)/binWidth;
  }
}

#endif
