#include "include/helperFunctions.h"

Double_t CalSysErrorPeriod(Double_t mean, Double_t sigma){
  if (mean < 0) mean *= -1.0;
  Double_t sigma2_m1 = sigma*sigma - 1;
  if (sigma2_m1 < 0) sigma2_m1 *= -1.0;
  
  return TMath::Sqrt(sigma2_m1) + mean/2.0;
}


void sysErrorPeriod(TString start =""){
  //Setup_______________
  const Int_t nBins =3;
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit ="true"; //not used right now

  Bool_t toWrite =true;
  //Setup_______________

  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Data/\
pullDist/";
  
  if (start==""){//Basic info
    cout <<"\nScript calculates systematic error due to period compatability\n";
    cout << "\n\nUtilization:" << endl;
    cout << "root \'sysErrorPeriod(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "False asymmetry data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Period and Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    if (whichFit == "true"){
      if (lrMrange != fitMrange){
	cout << "L/R mass range does not equal fit mass range for true fit\n";
	exit(EXIT_FAILURE);
      }
      cout << "L/R mass range:     " << lrMrange << endl;
    }
    else{
      cout << "LR integral mass range:     " << lrMrange << endl;
      cout << "Fit mass range:     " << fitMrange << endl;
    }
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //Needs to be cleaned up
  //File name setup && get file 
  TString inFile =Form("pullDist_%s_%s%s_%ibins.root", Mtype.Data(),
		       process.Data(), lrMrange.Data(), nBins);
  
  TFile *fPeriod = TFile::Open(path+inFile);
  if (!fPeriod ){
    cout << "Pull distribution file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data from files/Fit with gaussian
  TH1D *hPull = (TH1D*)fPeriod->Get("hPull");
  hPull->Fit("gaus");
  TF1 *f_gaus = hPull->GetFunction("gaus");
  Double_t mean = f_gaus->GetParameter(1);
  Double_t sigma = f_gaus->GetParameter(2);

  Double_t error = CalSysErrorPeriod(mean, sigma);
  Double_t sysError[nBins], xvals[nBins];
  for (Int_t i=0; i<nBins; i++) {
    sysError[i] = error;
    xvals[i] = 1 + i;
  }
  
  //Draw systematic error
  TGraph *gSys = new TGraphErrors(nBins, xvals, sysError);
  SetUp(gSys); gSys->GetYaxis()->SetRangeUser(0.0, sysError[0]*1.3);

  TCanvas* c1 = new TCanvas();
  gSys->Draw("AP");

  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Data/\
sysError";
  TString fOutput;//cleanup
  fOutput = Form("%s/sysErrorPeriod_WAll_%s_%s%s.root", thisDirPath.Data(),
		 Mtype.Data(), process.Data(), lrMrange.Data());
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    gSys->Write("gSys");
  }

  cout << "\nSettings________" << endl;
  cout << "False asymmetry data from:            " << path << endl;
  cout << "False asymmetry file used:                   " << inFile << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Period and Mass type considered:   " << Mtype << endl;
  cout << "AN physical process:        " << process << endl;
  if (whichFit == "true"){
    cout << "L/R mass range:     " << lrMrange << endl;
  }
  else{
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
  }
  cout << "Which fit considered:       " << whichFit << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
