#include "include/helperFunctions.h"


void physBinnedPeriod(TString start=""){
  TString path4TargAN ="/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/\
GeoMean4Targ";
  
  //Setup_______________
  const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xPi";//xN, xPi, xF, pT, M
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
    cout << "Data coming from:            " << path4TargAN << endl;
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
  TCanvas* c1 = new TCanvas();
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  Double_t yMax = 0.25;
  Int_t icolor[nPeriods] = {4, 6, 15, 2, 3, 7, 5, 9, 30};
  Int_t imarker[nPeriods] = {20, 21, 22, 23, 24, 25, 26, 27, 28};
  Double_t offset;
  if (physBinned =="xF") offset = 0.003;
  else if (physBinned =="pT") offset = 0.008;
  else if (physBinned =="xN") offset = 0.0009;
  else if (physBinned =="xPi") offset = 0.002;
  TList *doc = new TList(); TString docName="";

  //Get Data files/TGraphs and Draw
  for (Int_t p=0; p<nPeriods; p++) {
    TString fname;
    if  (whichFit=="true"){
      fname =Form("%s/GeoMean4Targ_true_W%s_%s_%s%s_%s%i.root",
		  path4TargAN.Data(), periods[p].Data(),
		  fitMrangeType.Data(),
		  process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins);
    }
    else{
      fname =Form("%s/GeoMean4Targ_%s%s_W%s_%s_%s%s_%s%i.root",
		  path4TargAN.Data(), whichFit.Data(),
		  fitMrange.Data(), periods[p].Data(),
		  fitMrangeType.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins);
    }
    TFile *f = TFile::Open(fname);
    docName+= fname+"\n";

    if (!f){
      cout << "RD or RD_noCorr file does not exist " << endl;
      exit(EXIT_FAILURE);
    }

    g_periods[p] = (TGraphErrors*) f->Get("AN");
    SetUp(g_periods[p], icolor[p], imarker[p], p*offset);

    Double_t wAvg = WeightedAvg(g_periods[p]);
    Double_t e_Wavg = WeightedErr(g_periods[p]);
    legend->AddEntry(g_periods[p], Form("W%s= %.3f +/- %.3f",
					periods[p].Data(), wAvg, e_Wavg),
		     "p");
    
    if (p==0) {
      g_periods[p]->Draw("AP");
      g_periods[p]->GetYaxis()->SetRangeUser(-yMax, yMax);
      g_periods[p]->SetTitle(physBinned);
      DrawLine(g_periods[p], 0.0);
    }
    else g_periods[p]->Draw("Psame");
    
  }//end p loop
  legend->Draw("same");

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility";
  TString fOutput;
  if  (whichFit=="true"){
    fOutput
      =Form("%s/Data/physBinned/physBinnedPeriod_true_%s_%s%s_%s%i.root",
	    thisDirPath.Data(), fitMrangeType.Data(),
	    process.Data(), lrMrange.Data(), physBinned.Data(), nBins);
  }
  else {
    fOutput
      =Form("%s/Data/physBinned/physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin.root",
	    thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
	    fitMrangeType.Data(), process.Data(), lrMrange.Data(),
	    physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    doc->Add((TObject*)(new TObjString(docName)) );
    doc->Write("InputData");
    
    for (Int_t p=0; p<nPeriods; p++)
      g_periods[p]->Write(Form("AN_%s", periods[p].Data() ));

    c1->Write("c1");
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << path4TargAN << endl;
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
