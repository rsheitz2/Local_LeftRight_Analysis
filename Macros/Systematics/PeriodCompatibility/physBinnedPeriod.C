#include "include/helperFunctions.h"

void physBinnedPeriod(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="M";//xN, xPi, xF, pT, M
  TString process ="DY";
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/
  
  /*const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="t3";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN ="/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/\
GeoMean4Targ";
  
  if (start==""){
    cout << "Script takes the AN data for a given DY kinematic binning/n";
    cout << "    made per period and makes a plot of it" << endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C  ->  GeoMean4Targ.C  ->";
    cout << "  physBinnedPeriod.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'physBinnedPeriod.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << fitMrangeType << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  //Period name setups
  const Int_t nPeriods =9;
  TString periods[nPeriods] =
    {"07", "08", "09", "10", "11", "12", "13", "14", "15"};
  TGraphErrors *g_periods[nPeriods];

  //Aesthetic setups
  TCanvas* cAllper = new TCanvas();
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  Double_t yMax = (fitMrangeType=="HMDY") ? 0.75 : 0.25;
  Int_t icolor[nPeriods] = {kBlue+2, kRed+2, kGreen+2, kMagenta+2, kCyan+2,
			    kBlue, kRed, kGreen, kMagenta};
  Int_t imarker[nPeriods] = {20, 21, 22, 23, 24, 25, 26, 27, 28};
  Double_t offset;
  if (physBinned =="xF") offset = 0.003;
  else if (physBinned =="pT") offset = 0.01;
  else if (physBinned =="xN") offset = 0.0009;
  else if (physBinned =="xPi") offset = 0.002;
  else if (physBinned =="xM") offset = 0.035;

  //Compute Weighed Average per Period 
  Double_t xVal[nPeriods], perWavg[nPeriods];
  Double_t ex[nPeriods] ={0.0}, e_perWavg[nPeriods];

  //Get Data files/TGraphs and Draw
  for (Int_t p=0; p<nPeriods; p++) {
    TString fname;
    if  (whichFit=="true"){
      fname =
	Form("%s/GeoMean4Targ_true_W%s_%s_%s%s_%s%s%i_%s_%s.root",pathAN.Data(),
	     periods[p].Data(), fitMrangeType.Data(), process.Data(),
	     fitMrange.Data(), binRange.Data(), physBinned.Data(), nBins,
	     production.Data(), additionalCuts.Data());
    }
    else{
      fname =
	Form("%s/GeoMean4Targ_%s%s_W%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(),periods[p].Data(),
	     fitMrangeType.Data(), process.Data(), lrMrange.Data(),
	     binRange.Data(), physBinned.Data(), nBins, hbins,production.Data(),
	     additionalCuts.Data());
    }
    cout << "Using data from:  " << fname << endl;
    TFile *f = TFile::Open(fname);

    if (!f){
      cout << "RD or RD_noCorr file does not exist " << endl;
      exit(EXIT_FAILURE);
    }

    g_periods[p] = (TGraphErrors*) f->Get("AN");
    SetUp(g_periods[p], icolor[p], imarker[p], p*offset);

    Double_t wAvg = WeightedAvg(g_periods[p]);
    Double_t e_Wavg = WeightedErr(g_periods[p]);
    legend->AddEntry(g_periods[p],
		     Form("W%s= %.3f +/- %.3f",periods[p].Data(), wAvg, e_Wavg),
		     "p");
    
    if (p==0) {
      g_periods[p]->Draw("AP");
      g_periods[p]->GetYaxis()->SetRangeUser(-yMax, yMax);
      g_periods[p]->SetTitle("AN");
      g_periods[p]->GetXaxis()->SetTitle(physBinned);
      if (physBinned=="M")
	g_periods[p]->GetXaxis()->SetLimits(4.3, 7.0);
      DrawLine(g_periods[p], 0.0);
    }
    else g_periods[p]->Draw("Psame");

    //Weighted average per period
    perWavg[p] = wAvg;
    e_perWavg[p] = e_Wavg;
    xVal[p] = p+1;
    
  }//end p loop
  legend->Draw("same");

  //Graph aesthetic setups/Draw Graph
  TGraphErrors *g_Wavg =
    new TGraphErrors(nPeriods, xVal, perWavg, ex, e_perWavg);
  SetUp(g_Wavg);
  
  TCanvas* cWavg = new TCanvas();
  g_Wavg->Draw("AP");
  g_Wavg->Fit("pol0");
  g_Wavg->SetTitle("Weighted Average");
  g_Wavg->GetXaxis()->SetTitle("Period");
  DrawLine(g_Wavg, 0.0);

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility";
  TString fOutput;
  if  (whichFit=="true"){
    fOutput
      =Form("%s/Data/physBinned/physBinnedPeriod_true_%s_%s%s_%s%s%i_%s_%s.root",
	    thisDirPath.Data(), fitMrangeType.Data(), process.Data(),
	    lrMrange.Data(), binRange.Data(), physBinned.Data(), nBins,
	    production.Data(), additionalCuts.Data());
  }
  else {
    cout << "Need to fix this file naming..." << endl;
    exit(EXIT_FAILURE);
    fOutput =
      Form("%s/Data/physBinned/\
physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
	   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
	   fitMrangeType.Data(), process.Data(), lrMrange.Data(),
	   physBinned.Data(), nBins, hbins, production.Data(),
	   additionalCuts.Data());
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    for (Int_t p=0; p<nPeriods; p++)
      g_periods[p]->Write(Form("AN_W%s", periods[p].Data() ));

    cAllper->Write("cAllper");

    g_Wavg->Write("wAvg");
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << fitMrangeType << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
  }
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
