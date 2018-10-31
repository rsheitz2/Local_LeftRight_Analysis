#ifndef FIT_MC_5_H
#define FIT_MC_5_H
//Reconstructed Monte Carlo with:  JPsi, Psi', OC, Drell-Yan
//{JPsi, psi, OC, LMDY, HMDY}
TH1D *hist_5[5] = {NULL, NULL, NULL, NULL, NULL}; 

Double_t Fit_MC_5(Double_t *x, Double_t *par){
  Double_t nJPsi = hist_5[0]->GetBinContent(hist_5[0]->FindBin(x[0] ) );
  Double_t npsi = hist_5[1]->GetBinContent(hist_5[1]->FindBin(x[0] ) );
  Double_t nOC = hist_5[2]->GetBinContent(hist_5[2]->FindBin(x[0] ) );
  Double_t nLMDY = hist_5[3]->GetBinContent(hist_5[3]->FindBin(x[0] ) );
  Double_t nHMDY = hist_5[4]->GetBinContent(hist_5[4]->FindBin(x[0] ) );
    
  //return par[0]*nJPsi +par[1]*npsi +par[2]*nOC + par[3]*nLMDY + par[4]*nHMDY;//cleanup
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC + par[4]*nHMDY;
}


void Paras_MC_5(TH1D *h, TF1 *fitFunc, Int_t nPar, Double_t Mmin,Double_t Mmax){
  Double_t binJPsi = hist_5[0]->GetMaximumBin();
  Double_t binPsi = hist_5[1]->GetMaximumBin();
  Double_t binOC = hist_5[2]->FindBin(Mmin);
  Double_t binDY = hist_5[3]->FindBin(4.0);
  
  Double_t A_JPsi = h->GetBinContent(binJPsi);
  Double_t A_psi = h->GetBinContent(binPsi);
  Double_t A_OC =  h->GetBinContent(binOC);
  Double_t A_HMDY =  h->GetBinContent(binDY);
  Double_t A_LMDY =  h->GetBinContent(binDY);
  
  if (!(hist_5[0]->GetBinContent(binJPsi) == 0) ){
    A_JPsi /= hist_5[0]->GetBinContent(binJPsi); }
  if (!(hist_5[1]->GetBinContent(binPsi) == 0) ){
    A_psi /= hist_5[1]->GetBinContent(binPsi); }
  if (!(hist_5[2]->GetBinContent(binOC) == 0) ){
    A_OC /= hist_5[2]->GetBinContent(binOC); }
  if (!(hist_5[3]->GetBinContent(binDY) == 0) ){
    A_HMDY /= hist_5[3]->GetBinContent(binDY); }
  if (!(hist_5[4]->GetBinContent(binDY) == 0) ){
    A_LMDY /= hist_5[4]->GetBinContent(binDY); }
    
  if (A_JPsi == 0) A_JPsi =200.0;
  if (A_psi == 0) A_psi =200.0;
  if (A_OC == 0) A_OC =200.0;
  if (A_HMDY == 0) A_HMDY =200.0;
  if (A_LMDY == 0) A_LMDY =200.0;

  Int_t ipar =0;
  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_OC); ipar++;//OC
  fitFunc->SetParameter(ipar, A_HMDY); ipar++;//HMDY
  fitFunc->SetParameter(ipar, A_LMDY); ipar++;//LMDY
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  
  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_OC/factor, A_OC*factor); ipar++;//OC
  fitFunc->SetParLimits(ipar, A_HMDY/factor, A_HMDY*factor); ipar++;//HMDY
  fitFunc->SetParLimits(ipar, A_LMDY/factor, A_LMDY*factor); ipar++;//LMDY  
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* SetupFunc_MC_5(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc,
		  Double_t Mmin, Double_t Mmax, Int_t *nPar, TH1D **hMC_5){
  //Setup fit MC_5 hist_5ograms
  for (Int_t i=0; i<5; i++) {
    hist_5[i] = hMC_5[i];}
  
  //Setup fit function
  *nPar =5;
  fitFunc = new TF1("fitFunc", Fit_MC_5, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_MC_5(h, fitFunc, *nPar, Mmin, Mmax);

  return fitFunc;
}


TH1D* ComponentFuncts_MC_5(Double_t *processPars, Double_t Mmin, Double_t Mmax,
			 TString process, Bool_t hIsUpS){
  TH1D *hResult;
  for (Int_t i=0; i<5; i++) {
    hist_5[i]->Scale(processPars[i] );
    //hist_5[i]->SetLineColor(kGreen);
    hist_5[i]->SetLineColor(i+1);
    hist_5[i]->GetXaxis()->SetRangeUser(Mmin, Mmax);
    hist_5[i]->Draw("same");

    if (i==0) hResult=(TH1D*) hist_5[0]->Clone();
    //else hResult->Add(hist_5[i]);//cleanup
  }
  hResult->Draw("same"); hResult->SetLineColor(kRed);
  hResult->Add(hist_5[1]); hResult->Add(hist_5[2]); hResult->Add(hist_5[4]);
  
  if (process == "JPsi") { hist_5[0]->SetLineColor(kBlack); }
  else if (process == "psi") { hist_5[1]->SetLineColor(kBlack); }
  else if (process == "DY") { hist_5[4]->SetLineColor(kBlack); }
  else {
    cout << "Wrong input process Fit_MC_5" << endl;
    exit(EXIT_FAILURE);
  }

  return hResult;
}


void IntegrateLR_MC_5(Double_t LR_Mrange, Double_t *LR, Double_t *e_LR,
		    TString process, Int_t bin){
  if (process =="JPsi"){
    Double_t massBin = hist_5[0]->GetMaximumBin();
    Double_t JPsi_M = hist_5[0]->GetBinCenter(massBin);
    Double_t lowerBin = hist_5[0]->FindBin(JPsi_M - LR_Mrange);
    Double_t upperBin = hist_5[0]->FindBin(JPsi_M + LR_Mrange);

    LR[bin] = hist_5[0]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else if (process =="psi"){
    Double_t massBin = hist_5[1]->GetMaximumBin();
    Double_t psi_M = hist_5[1]->GetBinCenter(massBin);
    Double_t lowerBin = hist_5[1]->FindBin(psi_M - LR_Mrange);
    Double_t upperBin = hist_5[1]->FindBin(psi_M + LR_Mrange);

    LR[bin] = hist_5[1]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else { cout << "Process not defined for this IntegrateLR_MC_5" << endl; }
}


void IntegrateLR_MC_5(Double_t LR_Mmin, Double_t LR_Mmax,
		    Double_t *LR, Double_t *e_LR, TString process, Int_t bin){
  if (process =="JPsi"){
    Double_t lowerBin = hist_5[0]->FindBin(LR_Mmin);
    Double_t upperBin = hist_5[0]->FindBin(LR_Mmax);

    LR[bin] = hist_5[0]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else if (process =="psi"){
    Double_t lowerBin = hist_5[1]->FindBin(LR_Mmin);
    Double_t upperBin = hist_5[1]->FindBin(LR_Mmax);

    LR[bin] = hist_5[1]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else if (process =="DY"){
    Double_t lowerBin = hist_5[4]->FindBin(LR_Mmin);
    Double_t upperBin = hist_5[4]->FindBin(LR_Mmax);

    LR[bin] = hist_5[4]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else { cout << "Process not defined for this IntegrateLR_MC_5" << endl; }
}

#endif
