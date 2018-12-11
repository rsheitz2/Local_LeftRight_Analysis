#include "include/helperFunctions.h"

void wAvg(){
  //Setup_______________
  const Int_t nBins =1;//# of physBinned bins
  const Int_t nHbins =16;
  TString Mtype ="HMDY";
  TString physBinned ="xN";//"xN", "xPi", "xF", "pT", "M"
  TString production ="slot1";//"t3", "slot1"
  
  Bool_t toWrite =false;
  //Setup_______________
  
  //Basic Setup
  TString pathDR="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio/Data";
  const Int_t nPer =9;
  TString period[nPer] ={"07", "08", "09", "10", "11", "12", "13", "14", "15"};
  
  Double_t yVal_wAvg[nBins] ={0.0}, e_yVal_wAvg[nBins] ={0.0};
  Double_t ex[nBins] ={0.0};
  Double_t *xvals;

  //Period Amps
  Double_t yPeriod[nPer], e_yPeriod[nPer];
  Double_t xPer[nPer] = {1, 2, 3, 4, 5, 6, 7, 8, 9}, ex_per[nPer] ={0.0};

  //Get Data and add for wAvg
  for (Int_t p=0; p<nPer; p++) {
    TString n_Wper =
      Form("%s/doubleRatio_W%s_%s_%s%i_%ihbins_%s.root", pathDR.Data(),
	   period[p].Data(), Mtype.Data(), physBinned.Data(), nBins, nHbins,
	   production.Data());
    TFile *f_Wper = OpenFile(n_Wper);
    TGraphErrors *g_Wper = (TGraphErrors*) f_Wper->Get("Amp");

    Double_t *yval = g_Wper->GetY();
    Double_t *e_yval = g_Wper->GetEY();
    for (Int_t i=0; i<nBins; i++) {
      yVal_wAvg[i] += yval[i]/(e_yval[i]*e_yval[i]);
      e_yVal_wAvg[i] += 1.0/(e_yval[i]*e_yval[i]);
    }

    //per period wAvg
    Double_t perWavg =0.0, e_perWavg =0.0;
    for (Int_t i=0; i<nBins; i++) {
      perWavg += yval[i]/(e_yval[i]*e_yval[i]);
      e_perWavg += 1.0/(e_yval[i]*e_yval[i]);
    }

    perWavg /= e_perWavg;
    e_perWavg = 1.0/TMath::Sqrt(e_perWavg);
    yPeriod[p] = perWavg;
    e_yPeriod[p] =e_perWavg;
    
    if (p==0) xvals = g_Wper->GetX();
  }

  //Final wAvg calculation
  for (Int_t i=0; i<nBins; i++) {
    yVal_wAvg[i] /= e_yVal_wAvg[i];
    e_yVal_wAvg[i] = TMath::Sqrt(1.0/e_yVal_wAvg[i]);
  }

  TGraphErrors* g_WAvg =
    new TGraphErrors(nBins, xvals, yVal_wAvg, ex, e_yVal_wAvg);
  SetUp(g_WAvg); 

  //Draw data
  TCanvas* c1 = new TCanvas();
  g_WAvg->Draw("AP");
  if (Mtype=="HMDY"){ g_WAvg->GetYaxis()->SetRangeUser(-0.5, 0.5); }
  else { g_WAvg->GetYaxis()->SetRangeUser(-0.08, 0.08); }
  g_WAvg->SetTitle("Final Amplitude from wAveraging");
  g_WAvg->GetXaxis()->SetTitle(physBinned);
  DrawLine(g_WAvg, 0.0);

  //per Period Wavg
  TGraphErrors* g_perWavg =
    new TGraphErrors(nPer, xPer, yPeriod, ex_per, e_yPeriod);
  SetUp(g_perWavg); g_perWavg->SetTitle("By Period");
  g_perWavg->GetXaxis()->SetTitle("period");
  g_perWavg->GetXaxis()->SetNdivisions(510);

  TCanvas* c2 = new TCanvas();
  g_perWavg->Draw("AP"); DrawLine(g_perWavg, 0.0);

  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio/Data/WAvg";
  TString fOutput =
    Form("%s/wAvg_%s_%s%i_%ihbins_%s.root", thisDirPath.Data(), Mtype.Data(),
	 physBinned.Data(), nBins, nHbins, production.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_WAvg->Write("Amp");
  }

  cout << "\nSettings______" << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Fit in nHbins:              " << nHbins << endl;
  cout << "Mass range considered:      " << Mtype << endl;
  cout << "Binned in which DY physics: " << physBinned << endl;
  cout << "Production considered:      " << production << endl;
  cout << "\nTo write output file:       " << toWrite << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
