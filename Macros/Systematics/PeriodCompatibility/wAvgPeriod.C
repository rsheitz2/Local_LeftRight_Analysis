#include "include/helperFunctions.h"


void wAvgPeriod(TString start=""){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/\
Data/physBinned";
    
  //Setup_______________
  const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xPi";//xN, xPi, xF, pT
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit ="true";

  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout << "Script takes the AN data for a given DY kinematic binning/n";
    cout << "    made per period and makes a plot of it" << endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C  ->  GeoMean4Targ.C  ->";
    cout << "  physBinnedPeriod.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'physBinnedPeriod.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Data coming from:            " << path << endl;
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
  
  //Get Data file
  TString fname;
  if (whichFit == "true"){
    if (fitMrange != lrMrange){
      cout << "fit Mass range != left/right mass range with true fit" << endl;
      exit(EXIT_FAILURE);
    }
    
    fname =Form("%s/physBinnedPeriod_%s_%s_%s%s_%s%i.root", path.Data(),
		whichFit.Data(), fitMrangeType.Data(), process.Data(),
		lrMrange.Data(), physBinned.Data(), nBins);
  }
  else {
    fname =Form("%s/physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin.root",path.Data(),
		whichFit.Data(), fitMrange.Data(), fitMrangeType.Data(),
		process.Data(),lrMrange.Data(),physBinned.Data(), nBins, hbins);
  }
  TFile *f_in = TFile::Open(fname);

  if (!f_in){
    cout << "RD or RD_noCorr file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Period name setups
  const Int_t nPeriods =9;
  TString periods[nPeriods]
    = {"07", "08", "09", "10", "11", "12", "13", "14", "15"};

  //Compute Weighed per Period in period loop
  Double_t xVal[nPeriods], perWavg[nPeriods];
  Double_t ex[nPeriods] ={0.0}, e_perWavg[nPeriods];
  
  for (Int_t p=0; p<nPeriods; p++) {
    TGraphErrors *g_AN
      = (TGraphErrors*) f_in->Get(Form("AN_%s", periods[p].Data() ));
    
    Double_t wAvg = WeightedAvg(g_AN);
    Double_t e_Wavg = WeightedErr(g_AN);

    perWavg[p] = wAvg;
    e_perWavg[p] = e_Wavg;
    xVal[p] = p+1;
  }//end p loop

  //Graph aesthetic setups/Draw Graph
  TGraphErrors *g_period_AN =
    new TGraphErrors(nPeriods, xVal, perWavg, ex, e_perWavg);
  g_period_AN->SetMarkerStyle(21);
  g_period_AN->GetYaxis()->SetLabelFont(22);
  g_period_AN->GetYaxis()->SetLabelSize(0.08);
  g_period_AN->GetXaxis()->SetLabelFont(22);
  g_period_AN->GetXaxis()->SetLabelSize(0.08);

  TCanvas* c1 = new TCanvas();
  g_period_AN->Draw("AP");
  g_period_AN->Fit("pol0");
  DrawLine(g_period_AN, 0.0);
  
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility";
  TString fOutput;
  if (whichFit == "true"){
    fOutput =Form("%s/Data/wAvgPeriod/wAvgPeriod_%s_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrangeType.Data(),
		  process.Data(), lrMrange.Data(), physBinned.Data(), nBins);
  }
  else {
    fOutput =Form("%s/Data/wAvgPeriod/wAvgPeriod_%s%s_%s_%s%s_%s%i_%ihbin.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  fitMrangeType.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TList *doc = new TList();
    TString docName = fname;
    doc->Add((TObject*)(new TObjString(docName)) );
    doc->Write("InputData");
    
    g_period_AN->Write("wAvg_AN");
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << path << endl;
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
