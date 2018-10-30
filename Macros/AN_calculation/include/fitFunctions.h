#ifndef FITFUNCTIONS_H
#define FITFUNCTIONS_H

//2 Gaussian w/ psi' M/W = A*JPsi M/W by target
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_six.h"

//2 Gaussian w/ psi' M/W = A*JPsi M/W by target
//1 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_seven.h"

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_eight.h"

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//1 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_nine.h"

//2 Crystal Ball 
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_ten.h"


void DoFit(TH1D *h, TMatrixDSym &cov, Double_t Mmin, Double_t Mmax){
  h->Sumw2(); //Try different fits till it works....
  TFitResultPtr status;
  for (Int_t i=0; i<20; i++) {
    status = h->Fit("fitFunc", "RLSQ", "", Mmin-0.05*i, Mmax+0.05*i);
    if (!status->Status() ) break;
    status = h->Fit("fitFunc", "RLSQ", "", Mmin-0.05*i, Mmax);
    if (!status->Status() ) break;
    status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax+0.05*i);
    
    if (status->Status() && (i+1)==20){
      cout << h->GetTitle() << "  Fit failed:  " << i+1 << endl;
      exit(EXIT_FAILURE);
    }
  }

  cov = status->GetCovarianceMatrix();
}


TF1* FitGetLR(TH1D **h, Int_t bin, Double_t *LR, Double_t *e_LR,
		Double_t LR_Mmin, Double_t LR_Mmax, TString process,
		TString whichFit, Double_t Mmin, Double_t Mmax){
  
  Double_t processPars[8] = {0.0};//{aJPsi,mJPsi,wJPsi,Apsi,Mpsi,Wpsi,aDY,cDY}
  Double_t LR_cov[25] = {0.0};
  Double_t binWidth =h[bin]->GetBinWidth(1);
  
  Bool_t hIsUpS;
  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_up", 12) == 0)
    hIsUpS =true;
  else if (strncmp(Form("%s", h[bin]->GetTitle() ),"MuMu_right_up",13)==0)
    hIsUpS =true;
  else if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_down", 14) == 0)
    hIsUpS =false;
  else if (strncmp(Form("%s", h[bin]->GetTitle() ),"MuMu_right_down", 15)==0)
    hIsUpS =false;
  else{
    cout << "Histogram target not defined well" << endl;
    exit(EXIT_FAILURE); }

  //Fit Setups_____
  TF1 *fitFunc = NULL; Int_t nPar;
  if (whichFit =="six"){
    fitFunc = SetupFunc_six(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_six(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    TF1 *f_LR =ComponentFuncts_six(processPars, Mmin, Mmax, process);
    IntegrateLR_six(f_LR, processPars, LR_cov, binWidth, LR_Mmin, LR_Mmax,
		    &(LR[bin]), &(e_LR[bin]), Mmin, Mmax, process);
  }
  else if (whichFit =="seven"){
    fitFunc = SetupFunc_seven(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_seven(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    TF1 *f_LR =ComponentFuncts_seven(processPars, Mmin, Mmax, process);
    IntegrateLR_seven(f_LR, processPars, LR_cov, binWidth, LR_Mmin, LR_Mmax,
		      &(LR[bin]), &(e_LR[bin]), Mmin, Mmax, process);
  }
  else if (whichFit =="eight"){
    fitFunc = SetupFunc_eight(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_eight(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    Double_t *pars = fitFunc->GetParameters();
    
    TF1 *f_LR =ComponentFuncts_eight(pars, Mmin, Mmax, process, hIsUpS);
    IntegrateLR_eight(f_LR, binWidth, LR_Mmin, LR_Mmax,
		      &(LR[bin]), &(e_LR[bin]) );//*/
    /*IntegrateLR_eight(f_LR, processPars, LR_cov, binWidth,
		      (LR_Mmax - LR_Mmin)/2.0, &(LR[bin]), &(e_LR[bin]),
		      hIsUpS, process);//*/
    /*IntegrateLRpercent_eight(f_LR, processPars, LR_cov, binWidth,
			     1.1, &(LR[bin]), &(e_LR[bin]),
			     hIsUpS, Mmin, Mmax, process);//*/
    /*IntegrateLRpercent_eight(f_LR, processPars, binWidth, 6.0,
      &(LR[bin]), &(e_LR[bin]), hIsUpS, process);*/
  }
  else if (whichFit =="nine"){
    fitFunc = SetupFunc_nine(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_nine(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    
    Double_t *pars = fitFunc->GetParameters();
    TF1 *f_LR =ComponentFuncts_nine(pars, Mmin, Mmax, process);
    IntegrateLR_nine(f_LR, pars, LR_cov, binWidth, LR_Mmin, LR_Mmax,
		    &(LR[bin]), &(e_LR[bin]), Mmin, Mmax, process);    
  }
  else if (whichFit =="ten"){
    fitFunc = SetupFunc_ten(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);

    ProcessPars_ten(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    
    Double_t *pars = fitFunc->GetParameters();
    TF1 *f_LR =ComponentFuncts_ten(pars, Mmin, Mmax, process, hIsUpS);
    /*IntegrateLR_ten(f_LR, processPars, LR_cov, binWidth,
		    (LR_Mmax - LR_Mmin)/2.0, &(LR[bin]), &(e_LR[bin]),
		    hIsUpS, Mmin, Mmax, process);//*/
    /*TF1 *f_LR =ComponentFuncts_ten(pars, Mmin, Mmax, process, hIsUpS,
				   hRatio[bin]);
    IntegrateLR_ten(f_LR, h[bin], pars, (LR_Mmax - LR_Mmin)/2.0, LR, e_LR,
    bin, process);//*/
    /*IntegrateLR_ten(f_LR, processPars, LR_cov, binWidth,
		    0.2, &(LR[bin]), &(e_LR[bin]),
		    hIsUpS, Mmin, Mmax, process);//*/
    IntegrateLRpercent_ten(f_LR, processPars, LR_cov, binWidth,
			   1.20, &(LR[bin]), &(e_LR[bin]),
			   hIsUpS, Mmin, Mmax, process);//*/
  }
  else{
    cout << "Int valid fit type:   " << whichFit << endl;
    exit(EXIT_FAILURE);
  }

  return fitFunc;
}


TF1* FitGetLR(TH1D **h, TH1D **hRatio, Int_t bin, Double_t *LR, Double_t *e_LR,
		Double_t LR_Mmin, Double_t LR_Mmax, TString process,
		TString whichFit, Double_t Mmin, Double_t Mmax){
  //Wrapper function to also get ratio histogram
  TF1 *fitFunc =
    FitGetLR(h, bin, LR, e_LR, LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax);
  
  //Determine ratio of fit to data
  TGraphErrors *gError = new TGraphErrors(hRatio[bin]->GetNbinsX() );
  for (Int_t i=1; i<hRatio[bin]->GetNbinsX()+1; i++) 
    gError->SetPoint(i, h[bin]->GetBinCenter(i), 0);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gError);

  for (Int_t bi=1; bi<h[bin]->GetNbinsX()+1; bi++) {
    Double_t RealData = hRatio[bin]->GetBinContent(bi);
    Double_t RDerror = hRatio[bin]->GetBinError(bi);
    
    Double_t fitValue = gError->Eval(h[bin]->GetBinCenter(bi) );
    Double_t fitError = gError->GetErrorY(bi);
    
    Double_t error = RatioError(RealData, fitValue, RDerror, fitError);
    Double_t ratio = RealData/fitValue;
    
    hRatio[bin]->SetBinContent(bi, ratio);
    hRatio[bin]->SetBinError(bi, error);
  }

  return fitFunc;
}


void SelectFitPars(TF1 *fitFunc, Double_t *pars, Double_t *e_pars,
		   TString process, TString whichFit, Bool_t hIsUpS){
  
  if (whichFit =="six") ProcessPars_six(fitFunc, pars, e_pars, process, hIsUpS);
  else if (whichFit =="seven") {
    ProcessPars_seven(fitFunc, pars, e_pars, process, hIsUpS);}
  else if (whichFit =="eight") {
    ProcessPars_eight(fitFunc, pars, e_pars, process, hIsUpS);}
  else if (whichFit =="nine") {
    ProcessPars_nine(fitFunc, pars, e_pars, process, hIsUpS);}
  else if (whichFit =="ten") {
    ProcessPars_ten(fitFunc, pars, e_pars, process, hIsUpS);}
  else{
    cout << "Int valid fit type:   " << whichFit << endl;
    exit(EXIT_FAILURE); }
}


void SelectFitPhysicsPars(TF1 *fitFunc, Double_t *pars, Double_t *e_pars,
		       TString process, TString whichFit, Bool_t hIsUpS){
  
  if (whichFit =="eight") {
    ProcessPhysicsPars_eight(fitFunc, pars, e_pars, process, hIsUpS);}
  else if (whichFit =="ten") {
    ProcessPhysicsPars_ten(fitFunc, pars, e_pars, process, hIsUpS);}
  else{
    cout << "Int valid fit type:   " << whichFit << endl;
    exit(EXIT_FAILURE); }
}

#endif
