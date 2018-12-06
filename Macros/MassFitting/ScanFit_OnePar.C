#include "include/helperFunctions.h"

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

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//   psi' alpha/n parameters free
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_eleven.h"


Double_t DoFit(TH1D *h, Double_t Mmin, Double_t Mmax, TString whichFit,
	       Double_t onePar){

  TF1 *fitFunc = NULL;
  Int_t nPar;
  if (whichFit =="six"){
    fitFunc = SetupFunc_six(h, fitFunc, Mmin, Mmax, &nPar, onePar);
  }
  else if (whichFit =="seven"){
    fitFunc = SetupFunc_seven(h, fitFunc, Mmin, Mmax, &nPar, onePar);
  }
  else if (whichFit =="eight"){
    fitFunc = SetupFunc_eight(h, fitFunc, Mmin, Mmax, &nPar, onePar);
  }
  else if (whichFit =="eleven"){
    fitFunc = SetupFunc_eleven(h, fitFunc, Mmin, Mmax, &nPar, onePar);
  }
  else{
    cout << "Incorrect fit input:   " << whichFit << " given to DoFit 1\n";
    exit(EXIT_FAILURE);
  }

  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ) return -1.0;

  Int_t ndf = fitFunc->GetNDF();
  Double_t chi2 = fitFunc->GetChisquare();
  
  return chi2/(1.0*ndf);
}


void ScanFit_OnePar(TString start=""){
  //Setup_______________
  TString period_Mtype ="WAll1_85";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =1;
  TString physBinned ="";//"", "xN", "xPi", "xF", "pT"
  TString whichFit ="eleven";
  Double_t Mmin =2.0;//Fit Mass minimum
  Double_t Mmax =8.50;//Fit Mass maximum
  Bool_t hIsUpS =true;
  //Double_t start_par=1.05, end_par=1.19, res_par=0.005; //ifits 1/2/6
  Double_t start_par=1.13, end_par=1.16, res_par=0.001; //ifit 5
  //Double_t start_par=1.152, end_par=start_par, res_par=0.001; //fit once
  
  //Setup_______________

  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Macros/MassFitting/Data/ScanFit_OnePar/";
  TString RDfile =Form("MuMu_%s_%s%ibins.root", period_Mtype.Data(),
		       physBinned.Data(), nBins);
  
  if (start==""){
    cout << "Script scans a specific constrain parameter in a fit";
    cout << " to get the mimimum chi2 value" << endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C" << endl;
    cout << "Usage:" << endl;
    cout << "root \'ScanFit_OnePar.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Mass data path:             " << pathRD << endl;
    cout << "Mass data file considered:  " << RDfile << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
    cout << "Scanning upstream target:   " << hIsUpS << endl;
    cout << "Starting/endding parameter value:   " << start_par << " - "
	 << end_par << endl;
    cout << "Parameter increment:   "  << res_par << endl;
    exit(EXIT_FAILURE);
  }

  //Get Input File and input histogram
  TFile *fRD  = TFile::Open(pathRD + RDfile);
  if ( !fRD ){
    cout << "RD file does not exist " << endl;
    exit(EXIT_FAILURE);
  }
  TH1D *h_mumu = (hIsUpS) ? (TH1D*)fRD->Get("h_mumu_up") :
    (TH1D*)fRD->Get("h_mumu_down");

  //Scan over parameter values
  TCanvas* cFit = new TCanvas();
  Double_t minChi2 =1000.0, minPar =1000.0;
  vector <Double_t> val_chi, xval;
  cout << "\nFit failure with input value(s):   " << endl;
  for (Double_t val=start_par; val<end_par+res_par; val+=res_par){
    Double_t Chi2 = DoFit(h_mumu, Mmin, Mmax, whichFit, val);

    if (Chi2 < 0){
      cout << val << endl;//Fit failed input value
      continue;
    }
    
    val_chi.push_back(Chi2);

    if (Chi2 < minChi2) {
      minChi2 = Chi2;
      minPar = val;
    }

    xval.push_back(val);
  }


  TCanvas* cChi = new TCanvas();
  TGraph *gChi = new TGraph(xval.size(), &xval[0], &val_chi[0]);
  SetUp(gChi);
  gChi->Draw("AP");
  
  //End of Macro Output
  cout << " " << endl;
  cout << "Minimum Chi2:  " << minChi2 << "  iPar: " << minPar << endl;
  cout << "Note:" << endl;
  cout << "-b/(2c) = min x val   for a + b*x + c*x^2" << endl;
}
