#include "include/helperFunctions.h"

Double_t Mmin;
Double_t Mmax;

#include "include/Fits_1_85/Fit_one.h"//ratioPsi
#include "include/Fits_1_85/Fit_two.h"//ratioPsi_upS, ratioPsi_downS

#include "include/Fits_1_85/Fit_five.h"//ratioPsi_upS_2exp, ratioPsi_downS_2exp


Double_t DoFit(TH1D *h, Int_t iFit){
  Mmin = h->GetXaxis()->GetXmin();
  Mmax = h->GetXaxis()->GetXmax();

  Bool_t hIsUpS =false;
  if (strncmp(Form("%s", h->GetTitle() ), "h_mumu_up", 9) == 0)
    hIsUpS=true;
  
  
  TF1* fitFunc;
  if (iFit==1){
    Int_t nPar =10;
    fitFunc = new TF1("fitFunc", Fit_one, Mmin, Mmax, nPar);
    Paras_one(fitFunc, nPar);
  }
  else if (iFit==2){
    Int_t nPar =10;
    if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_two_upS, Mmin, Mmax, nPar);
    else fitFunc = new TF1("fitFunc", Fit_two_downS, Mmin, Mmax, nPar);
    Paras_two(fitFunc, nPar);
  }
  else if (iFit==5){
    Int_t nPar =8;
    if (hIsUpS) fitFunc = new TF1("fitFunc", Fit_five_upS, Mmin, Mmax, nPar);
    else fitFunc = new TF1("fitFunc", Fit_five_downS, Mmin, Mmax, nPar);
    Paras_five(fitFunc, nPar);
  }
  else{
    cout << "Incorrect iFit value:   " << iFit << " given to DoFit 1\n";
    exit(EXIT_FAILURE);
  }

  
  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ) return -1.0;

  Int_t ndf = fitFunc->GetNDF();
  Double_t chi2 = fitFunc->GetChisquare();
  
  return chi2/(1.0*ndf);
}


void ScanFit_One_Pars1_85(TString start=""){
  if (start=="") {
    cout << "Script scans a specific constrain parameter in a fit";
    cout << " to get the mimimum chi2 value" << endl;
    cout << "This script is specifically for mass range of 1.0-8.5 GeV" << endl;
    cout <<"\nUsage:" << endl;
    cout << "root \'ScanFit_One_Pars1_85(1)\'\n";
    exit(EXIT_FAILURE);
  }
  //void ScanFitPars1_85(Int_t iFit){
  //Setup
  Int_t iFit =1;
  Double_t start_one=1.05, end_one=1.19, res_one=0.005; //ifits 1/2/6
  //Double_t start_one=1.13, end_one=1.16, res_one=0.001; //ifit 5
  Bool_t chiFromBothTarg;
  Bool_t isUpS =false;

  Double_t *one;
  switch(iFit){
  case 1 : chiFromBothTarg = true;
    one = &ratioPsi;
    break;
  case 2 : chiFromBothTarg = false;
    cout << "Remember to check which target Chi2 is coming from!!!" << endl;
    (isUpS) ? one = &ratioPsi_upS : one = &ratioPsi_downS;
    break;
  case 5 : chiFromBothTarg = false;
    cout << "Remember to check which target Chi2 is coming from!!!" << endl;
    (isUpS) ? one = &ratioPsi_upS_2exp : one = &ratioPsi_downS_2exp;
    break;
  default :
    cout << "Wrong iFit entered";
    exit(EXIT_FAILURE);
  }
  
  vector <Double_t> val_one;
  for (Double_t val=start_one; val<end_one+res_one; val+=res_one)
    val_one.push_back(val);

  TString inputPath = "InputData/";
  TString inputData = "MuMu_1_85_notScaled.root";
  TFile *f = TFile::Open(inputPath+inputData);
  if (!f){
    cout << "File   " << inputData << "  does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  TH1D *h_mumu_up = (TH1D*)f->Get("h_mumu_up");
  TH1D *h_mumu_down = (TH1D*)f->Get("h_mumu_down");
  TCanvas* cFit = new TCanvas();

  Double_t minChi2 = 1000.0;
  Double_t minOne, minTwo;
  vector <Double_t> val_chi;
  vector <Double_t> xval;
  cout << "\nFit failure with input value(s):   " << endl;
  for(vector<Double_t>::iterator iOne=val_one.begin(); iOne!=val_one.end();
      iOne++){
    *one = *iOne;
      
    Double_t Chi2_upS = DoFit(h_mumu_up, iFit);
    Double_t Chi2_downS = DoFit(h_mumu_down, iFit);

    Double_t Chi2;
    if (Chi2_upS > 0 && Chi2_downS > 0){
      if (chiFromBothTarg) {
	Chi2 = Chi2_upS + Chi2_downS;
	Chi2 /=2;
      }
      else if (isUpS) Chi2 = Chi2_upS;
      else Chi2 = Chi2_downS;
    }
    else {
      cout << *iOne << endl;//Fit failed input value
      continue;
    }
    
    
    val_chi.push_back(Chi2);

    if (Chi2 < minChi2) {
      minChi2 = Chi2;
      minOne = *iOne;
    }

    xval.push_back(*iOne);
  }


  TCanvas* cChi = new TCanvas();
  TGraph *gChi = new TGraph(xval.size(), &xval[0], &val_chi[0]);
  SetUp(gChi);
  gChi->Draw("AP");
  
  //End of Macro Output
  cout << " " << endl;
  cout << "Input data is:  " << inputData << endl;
  cout << "iFit:  " << iFit << endl;
  switch(iFit){
  case 1 : cout << "2 Gaussians Constrained (psi'=a*JPsi) && 3 Exponentials\n";
    break;
  case 2 : cout << "2 Gaussians Constrained by target (psi'=a*JPsi) ";
    cout << "&& 3 Exponentials\n";
    break;
  case 4 : cout<< "2 Gaussians Constrained (psi'=a*JPsi) && 2 Exponentials\n";
    break;
  case 5 : cout << "2 Gaussians Constrained by target (psi'=a*JPsi) ";
    cout << "&& 2 Exponentials\n";
    break;
  }
  cout << "Size of val_one:  " << val_one.size() << endl;
  cout << " " << endl;
  cout << "Minimum Chi2:  " << minChi2 << "  iOne: " << minOne << endl;
  cout << "Chi2 from both targets:   " << chiFromBothTarg << endl;
  if (!chiFromBothTarg) {
    cout << "    Chi2 from    ";
    (isUpS) ? cout << "UpS" << endl : cout << "DownS" << endl;
  }
  cout << "Note:" << endl;
  cout << "-b/(2c) = min x val   for a + b*x + c*x^2" << endl;
}
