#ifndef FIT_MC_H
#define FIT_MC_H
//Reconstructed Monte Carlo with:  JPsi, Psi', OC, Drell-Yan
TH1D *hist[4] = {NULL, NULL, NULL, NULL}; //{JPsi, psi, OC, AMDY}

Double_t Fit_MC(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );
    
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}


void Paras_MC(TH1D *h, TF1 *fitFunc, Int_t nPar, Double_t Mmin, Double_t Mmax){
  Double_t binJPsi = hist[0]->GetMaximumBin();
  Double_t binPsi = hist[1]->GetMaximumBin();
  Double_t binOC = hist[2]->FindBin(Mmin);
  Double_t binDY = hist[3]->FindBin(4.0);
  
  Double_t A_JPsi = h->GetBinContent(binJPsi);
  Double_t A_psi = h->GetBinContent(binPsi);
  Double_t A_OC =  h->GetBinContent(binOC);
  Double_t A_DY =  h->GetBinContent(binDY);
  
  if (!(hist[0]->GetBinContent(binJPsi) == 0) ){
    A_JPsi /= hist[0]->GetBinContent(binJPsi); }
  if (!(hist[1]->GetBinContent(binPsi) == 0) ){
    A_psi /= hist[1]->GetBinContent(binPsi); }
  if (!(hist[2]->GetBinContent(binOC) == 0) ){
    A_OC /= hist[2]->GetBinContent(binOC); }
  if (!(hist[3]->GetBinContent(binDY) == 0) ){
    A_DY /= hist[3]->GetBinContent(binDY); }
    
  if (A_JPsi == 0) A_JPsi =200.0;
  if (A_psi == 0) A_psi =200.0;
  if (A_OC == 0) A_OC =200.0;
  if (A_DY == 0) A_DY =200.0;

  Int_t ipar =0;
  fitFunc->SetParameter(ipar, A_JPsi); ipar++;//JPsi
  fitFunc->SetParameter(ipar, A_psi); ipar++;//psi
  fitFunc->SetParameter(ipar, A_OC); ipar++;//OC
  fitFunc->SetParameter(ipar, A_DY); ipar++;//DY
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
  
  ipar=0;
  Double_t factor =5.0;
  fitFunc->SetParLimits(ipar, A_JPsi/factor, A_JPsi*factor); ipar++;//JPsi
  fitFunc->SetParLimits(ipar, A_psi/factor, A_psi*factor); ipar++;//psi
  fitFunc->SetParLimits(ipar, A_OC/factor, A_OC*factor); ipar++;//OC
  fitFunc->SetParLimits(ipar, A_DY/factor, A_DY*factor); ipar++;//DY
  if (ipar != nPar){
    cout << "ipar problem" << endl;
    exit(EXIT_FAILURE);
  }
}


TF1* SetupFunc_MC(TH1D *h, Bool_t hIsUpS, TF1 *fitFunc,
		  Double_t Mmin, Double_t Mmax, Int_t *nPar, TH1D **hMC){
  //Setup fit MC histograms
  for (Int_t i=0; i<4; i++) {
    hist[i] = hMC[i];}
  
  //Setup fit function
  *nPar =4;
  fitFunc = new TF1("fitFunc", Fit_MC, Mmin, Mmax, *nPar);

  //Setup intial parameters of fit and parameter constraints
  Paras_MC(h, fitFunc, *nPar, Mmin, Mmax);

  return fitFunc;
}


TH1D* ComponentFuncts_MC(Double_t *processPars, Double_t Mmin, Double_t Mmax,
			 TString process, Bool_t hIsUpS){
  for (Int_t i=0; i<4; i++) {
    hist[i]->Scale(processPars[i] );
    hist[i]->SetLineColor(kGreen);
    hist[i]->GetXaxis()->SetRangeUser(Mmin, Mmax);
    hist[i]->Draw("same");
  }
  
  TH1D *hResult = (TH1D*)hist[0]->Clone();
  hResult->Add(hist[1]); hResult->Add(hist[2]); hResult->Add(hist[3]);
  hResult->SetLineColor(kRed);
  hResult->Draw("same");
  
  if (process == "JPsi") { hist[0]->SetLineColor(kBlack); }
  else if (process == "psi") { hist[1]->SetLineColor(kBlack); }
  else if (process == "DY") { hist[3]->SetLineColor(kBlack); }
  else {
    cout << "Wrong input process Fit_MC" << endl;
    exit(EXIT_FAILURE);
  }

  return hResult;
}


void IntegrateLR_MC(Double_t LR_Mrange, Double_t *LR, Double_t *e_LR,
		    TString process, Int_t bin){
  if (process =="JPsi"){
    Double_t massBin = hist[0]->GetMaximumBin();
    Double_t JPsi_M = hist[0]->GetBinCenter(massBin);
    Double_t lowerBin = hist[0]->FindBin(JPsi_M - LR_Mrange);
    Double_t upperBin = hist[0]->FindBin(JPsi_M + LR_Mrange);

    LR[bin] = hist[0]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else if (process =="psi"){
    Double_t massBin = hist[1]->GetMaximumBin();
    Double_t psi_M = hist[1]->GetBinCenter(massBin);
    Double_t lowerBin = hist[1]->FindBin(psi_M - LR_Mrange);
    Double_t upperBin = hist[1]->FindBin(psi_M + LR_Mrange);

    LR[bin] = hist[1]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else { cout << "Process not defined for this IntegrateLR_MC" << endl; }
}


void IntegrateLR_MC(Double_t LR_Mmin, Double_t LR_Mmax,
		    Double_t *LR, Double_t *e_LR, TString process, Int_t bin){
  if (process =="JPsi"){
    Double_t lowerBin = hist[0]->FindBin(LR_Mmin);
    Double_t upperBin = hist[0]->FindBin(LR_Mmax);

    LR[bin] = hist[0]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else if (process =="psi"){
    Double_t lowerBin = hist[1]->FindBin(LR_Mmin);
    Double_t upperBin = hist[1]->FindBin(LR_Mmax);

    LR[bin] = hist[1]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else if (process =="DY"){
    Double_t lowerBin = hist[3]->FindBin(LR_Mmin);
    Double_t upperBin = hist[3]->FindBin(LR_Mmax);

    LR[bin] = hist[3]->IntegralAndError(lowerBin, upperBin, e_LR[bin]);
  }
  else { cout << "Process not defined for this IntegrateLR_MC" << endl; }
}

#endif
