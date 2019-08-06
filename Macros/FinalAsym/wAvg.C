#include "include/helperFunctions.h"

void AddData(TGraphErrors *g, Double_t *wAvg, Double_t *eWavg){
  Double_t *yval = g->GetY();
  Double_t *e_yval = g->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    wAvg[i] += yval[i]/(e_yval[i]*e_yval[i]);
    eWavg[i] += 1.0/(e_yval[i]*e_yval[i]);
  }
}


void TableFormat(Int_t p, TString physBinned, TGraphErrors *g);

void TableFormat(Int_t p, TString physBinned, TGraphErrors *g,
		 TGraphErrors *g_uncorr);

void TableFormat(TString p, TString physBinned, TGraphErrors *g,
		 TGraphErrors *g_uncorr);

void FinalWavg(Double_t *wAvg, Double_t *eWavg, Int_t nBins){
  for (Int_t i=0; i<nBins; i++) {
    wAvg[i] /= eWavg[i];
    eWavg[i] = TMath::Sqrt(1.0/eWavg[i]);
  }
}


void wAvg(){
  //Setup_______________
  /*const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString physBinned ="M";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  const Int_t nBins =4;//JPsi
  TString Mtype ="LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="M";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.87_3.38";
  TString fitMrange ="2.87_3.38";
  TString binRange ="29_34";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
   //Setup_______________

  
  //Basic Setup
  TString pathAN="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";
  const Int_t nPer =9;
  TString period[nPer] ={"07", "08", "09", "10", "11", "12", "13", "14", "15"};
  
  Double_t yVal_wAvg[nBins] ={0.0}, e_yVal_wAvg[nBins] ={0.0};
  Double_t yVal_wAvg_uncorr[nBins] ={0.0}, e_yVal_wAvg_uncorr[nBins] ={0.0};
  Double_t yVal_wAvg_upS[nBins] ={0.0}, e_yVal_wAvg_upS[nBins] ={0.0};
  Double_t yVal_wAvg_downS[nBins] ={0.0}, e_yVal_wAvg_downS[nBins] ={0.0};
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
    TGraphErrors *g_Wper_uncorr = (TGraphErrors*) f_Wper->Get("ANuncorr");
    TGraphErrors *g_Wper_upS = (TGraphErrors*) f_Wper->Get("AN_ups");
    TGraphErrors *g_Wper_downS = (TGraphErrors*) f_Wper->Get("AN_downs");

    AddData(g_Wper, yVal_wAvg, e_yVal_wAvg);
    AddData(g_Wper_uncorr, yVal_wAvg_uncorr, e_yVal_wAvg_uncorr);
    AddData(g_Wper_upS, yVal_wAvg_upS, e_yVal_wAvg_upS);
    AddData(g_Wper_downS, yVal_wAvg_downS, e_yVal_wAvg_downS);

    TableFormat(p, physBinned, g_Wper, g_Wper_uncorr);

    if (p==0) xvals = g_Wper->GetX();
  }

  //Final wAvg calculation
  FinalWavg(yVal_wAvg, e_yVal_wAvg, nBins);
  FinalWavg(yVal_wAvg_uncorr, e_yVal_wAvg_uncorr, nBins);
  FinalWavg(yVal_wAvg_upS, e_yVal_wAvg_upS, nBins);
  FinalWavg(yVal_wAvg_downS, e_yVal_wAvg_downS, nBins);

  TGraphErrors* g_WAvg =
    new TGraphErrors(nBins, xvals, yVal_wAvg, ex, e_yVal_wAvg);
  TGraphErrors* g_WAvg_uncorr =
    new TGraphErrors(nBins, xvals, yVal_wAvg_uncorr, ex, e_yVal_wAvg_uncorr);
  TGraphErrors* g_WAvg_upS =
    new TGraphErrors(nBins, xvals, yVal_wAvg_upS, ex, e_yVal_wAvg_upS);
  TGraphErrors* g_WAvg_downS =
    new TGraphErrors(nBins, xvals, yVal_wAvg_downS, ex, e_yVal_wAvg_downS);
  SetUp(g_WAvg); SetUp(g_WAvg_upS); SetUp(g_WAvg_downS);
  if (physBinned=="xN"){//Specific setups
    g_WAvg->GetXaxis()->SetNdivisions(503);
    g_WAvg_upS->GetXaxis()->SetNdivisions(503);
    g_WAvg_downS->GetXaxis()->SetNdivisions(503);
  }
  else if (physBinned=="xF"){
    g_WAvg->GetXaxis()->SetLimits(0.05, 0.65);
    g_WAvg_upS->GetXaxis()->SetLimits(0.05, 0.65);
    g_WAvg_downS->GetXaxis()->SetLimits(0.05, 0.65);
  }

  TableFormat("avg", physBinned, g_WAvg, g_WAvg_uncorr);

  //Draw data
  TCanvas* c1 = new TCanvas();
  g_WAvg->Draw("AP");
  if (Mtype=="HMDY"){ g_WAvg->GetYaxis()->SetRangeUser(-0.25, 0.25); }
  else { g_WAvg->GetYaxis()->SetRangeUser(-0.08, 0.08); }
  DrawLine(g_WAvg, 0.0);

  TAxis *xaxis = g_WAvg_upS->GetXaxis();
  Double_t range = xaxis->GetXmax() - xaxis->GetXmin();
  OffSet(g_WAvg_downS, range/30.0);
  g_WAvg_upS->Draw("Psame"); g_WAvg_upS->SetMarkerColor(kBlue);
  g_WAvg_downS->Draw("Psame"); g_WAvg_downS->SetMarkerColor(kRed);

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
    g_WAvg_upS->Write("AN_ups");
    g_WAvg_downS->Write("AN_downs");
  }

  /*cout << " " << endl;
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
  else cout << "File: " << fOutput << " was NOT written" << endl;*/
}


void TableFormat(Int_t p, TString physBinned, TGraphErrors *g){
  Double_t *yval = g->GetY();
  Double_t *e_yval = g->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    cout << "W" << p+7 << ", " << physBinned << ", "
	 << yval[i] << ", " << e_yval[i] << endl;
  }
}


void TableFormat(Int_t p, TString physBinned, TGraphErrors *g,
		 TGraphErrors *g_uncorr){
  Double_t *yval = g->GetY();
  Double_t *e_yval = g->GetEY();
  Double_t *yval_uncorr = g_uncorr->GetY();
  Double_t *e_yval_uncorr = g_uncorr->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    cout << "W" << p+7 << ", " << physBinned << ", "
	 << yval[i] << ", " << e_yval[i] << ", "
	 << yval_uncorr[i] << ", " << e_yval_uncorr[i] << ", "
	 <<  yval_uncorr[i]/yval[i] << endl;
  }
}


void TableFormat(TString p, TString physBinned, TGraphErrors *g,
		 TGraphErrors *g_uncorr){
  Double_t *yval = g->GetY();
  Double_t *e_yval = g->GetEY();
  Double_t *yval_uncorr = g_uncorr->GetY();
  Double_t *e_yval_uncorr = g_uncorr->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    cout << p << ", " << physBinned << ", "
	 << yval[i] << ", " << e_yval[i] << ", "
	 << yval_uncorr[i] << ", " << e_yval_uncorr[i] << ", "
	 <<  yval_uncorr[i]/yval[i] << endl;
  }
}
