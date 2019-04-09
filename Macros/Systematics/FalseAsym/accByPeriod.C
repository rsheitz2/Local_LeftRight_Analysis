#include "include/helperFunctions.h"

void accByPeriod(TString start=""){
  //Setup_______________
  /*const Int_t nBins =1;//HMDY
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/
  
  const Int_t nBins =1;
  TString fitMrangeType ="LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.80_3.44";
  TString fitMrange ="2.80_3.44";
  TString binRange ="25_43";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN ="/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/\
GeoMean4Targ";
  TString pathFA ="/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data";
  
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

  //Aesthetic setups
  Double_t yMax = (fitMrangeType=="HMDY") ? 1.2 : 1.05;
  TCanvas* cAlpha = new TCanvas(); cAlpha->Divide(1,2);

  //Get Data files/TGraphs and Draw
  for (Int_t p=0; p<nPeriods; p++) {
    TString fname_FA;
    if  (whichFit=="true"){
      fname_FA = Form("%s/acceptanceFourTargRatio/\
acceptanceFourTargRatio_true_W%s_%s_%s%s_%s%s%i_%s_%s.root",
		   pathFA.Data(), periods[p].Data(), fitMrangeType.Data(),
		   process.Data(), fitMrange.Data(), binRange.Data(),
		   physBinned.Data(), nBins, production.Data(),
		   additionalCuts.Data());
    }
    else{
      cout <<"not setup" <<endl;
      exit(EXIT_FAILURE);
      fname_FA =
	Form("%s/GeoMean4Targ_%s%s_W%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(),periods[p].Data(),
	     fitMrangeType.Data(), process.Data(), lrMrange.Data(),
	     binRange.Data(), physBinned.Data(), nBins, hbins,production.Data(),
	     additionalCuts.Data());
    }
    cout << "Using data from:  " << fname_FA << endl;
    TFile *f = OpenFile(fname_FA);

    TGraphErrors *g_alpha_subper = (TGraphErrors*) f->Get("alpha_subper");
    TGraphErrors *g_alpha_pol = (TGraphErrors*) f->Get("alpha_pol");
    
    TGraphErrors *g_alpha_upS = (TGraphErrors*) f->Get("alpha_upS");
    TGraphErrors *g_alpha_downS = (TGraphErrors*) f->Get("alpha_downS");

    OffSet(g_alpha_subper, p+1.0); OffSet(g_alpha_pol, p+1.0);
    OffSet(g_alpha_upS, p+1.0); OffSet(g_alpha_downS, p+1.0); 

    //Draw acceptance
    cAlpha->cd(1);
    if (p==0) {
      g_alpha_subper->Draw("AP");
      g_alpha_subper->GetYaxis()->SetRangeUser(2-yMax , yMax);
      g_alpha_subper->GetXaxis()->SetLimits(0, nPeriods+1);
      DrawLine(g_alpha_subper, 1.0);
    }
    else g_alpha_subper->Draw("Psame");

    g_alpha_pol->Draw("Psame");

    cAlpha->cd(2);
    if (p==0) {
      g_alpha_upS->Draw("AP");
      g_alpha_upS->GetYaxis()->SetRangeUser(2-yMax, yMax);
      g_alpha_upS->GetXaxis()->SetLimits(0, nPeriods+1);
      DrawLine(g_alpha_upS, 1.0);
    }
    else g_alpha_upS->Draw("Psame");

    g_alpha_downS->Draw("Psame");
    
  }//end p loop

  /*//Write Output/Final Settings
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
  else cout << "File: " << fOutput << " was NOT written" << endl;//*/
}
