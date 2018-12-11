#include "include/helperFunctions.h"

void wAvg(){
  //Setup_______________
  const Int_t nBins =1;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  /*const Int_t nBins =1;//JPsi
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =true;
   //Setup_______________

  
  //Basic Setup
  TString pathAN="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";
  const Int_t nPer =9;
  TString period[nPer] ={"07", "08", "09", "10", "11", "12", "13", "14", "15"};
  
  Double_t yVal_wAvg[nBins] ={0.0}, e_yVal_wAvg[nBins] ={0.0};
  Double_t ex[nBins] ={0.0};
  Double_t *xvals;

  //Get Data and add for wAvg
  for (Int_t p=0; p<nPer; p++) {
    TString n_Wper;

    if (whichFit=="true"){
      n_Wper =
	Form("%s/GeoMean4Targ_%s_W%s_%s_%s%s_%s%s%i_%s_%s.root",pathAN.Data(),
	     whichFit.Data(), period[p].Data(), Mtype.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),physBinned.Data(),
	     nBins, production.Data(), additionalCuts.Data());
    }
    else{
      n_Wper =
	Form("%s/GeoMean4Targ_%s%s_W%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(), period[p].Data(),
	     Mtype.Data(), process.Data(), lrMrange.Data(), binRange.Data(),
	     physBinned.Data(), nBins, hbins, production.Data(),
	     additionalCuts.Data());
    }
    TFile *f_Wper = OpenFile(n_Wper);
    TGraphErrors *g_Wper = (TGraphErrors*) f_Wper->Get("AN");

    Double_t *yval = g_Wper->GetY();
    Double_t *e_yval = g_Wper->GetEY();
    for (Int_t i=0; i<nBins; i++) {
      yVal_wAvg[i] += yval[i]/(e_yval[i]*e_yval[i]);
      e_yVal_wAvg[i] += 1.0/(e_yval[i]*e_yval[i]);
    }

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
  if (Mtype=="HMDY"){ g_WAvg->GetYaxis()->SetRangeUser(-0.25, 0.25); }
  else { g_WAvg->GetYaxis()->SetRangeUser(-0.08, 0.08); }
  DrawLine(g_WAvg, 0.0);

  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data";
  TString fOutput;

  if (whichFit=="true"){
    fOutput = Form("%s/WAvg/wAvg_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), Mtype.Data(),
		   process.Data(), lrMrange.Data(), binRange.Data(),
		   physBinned.Data(), nBins, hbins, production.Data(),
		   additionalCuts.Data());
  }
  else{
    fOutput = Form("%s/WAvg/wAvg_%s%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		   Mtype.Data(), process.Data(), lrMrange.Data(),
		   physBinned.Data(), binRange.Data(), nBins, hbins,
		   production.Data(), additionalCuts.Data());
  }
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_WAvg->Write("AN");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathAN << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Mass type considered:   " << Mtype << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
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
