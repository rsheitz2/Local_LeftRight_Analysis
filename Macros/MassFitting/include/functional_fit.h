#ifndef FUNCTIONAL_FIT_H
#define FUNCTIONAL_FIT_H
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

#endif
