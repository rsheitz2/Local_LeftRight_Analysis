#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void DrawLegend(TGraphErrors *g, TString whichName){
  g->SetName("gr");
  Double_t sigma;
  Double_t avg = WeightedAvgAndError(g, &sigma);

  TLegend *leg = new TLegend(0.1,0.9,0.7,0.99);
  leg->SetBorderSize(0); leg->SetTextFont(132); leg->SetTextSize(0.08);
  
  if(whichName=="Acc"){
    leg->AddEntry("gr", Form("#bar{Acc} = %0.2f #pm %0.2f", avg, sigma), "p");
  }
  else if(whichName=="FA"){
    leg->AddEntry("gr", Form("#bar{A}_{N,false} = %0.2f #pm %0.2f", avg, sigma),
		  "p");
  }
  
  leg->Draw("same");
}


void DrawLegend(TGraph *g, TString whichName){
  g->SetName("gr");
  Double_t avg = Avg(g);

  TLegend *leg = new TLegend(0.1,0.9,0.7,0.99);
  leg->SetBorderSize(0); leg->SetTextFont(132); leg->SetTextSize(0.08);
  
  if(whichName=="SysStat"){
    leg->AddEntry("gr",
		  Form("#bar{#sigma_{systemaic}/#sigma_{statistical}} = %0.2f",
		       avg), "p");
  }
  
  leg->Draw("same");
}


void all_FA_Acc(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";
  
  //const Int_t nBins =5;//JPsi
  //TString Mtype ="W12_LowM_AMDY";
  //Int_t hbins =150;
  //TString process ="JPsi";//JPsi, psi, DY
  //TString lrMrange ="2.00_5.00";
  //TString fitMrange ="2.00_7.50";
  //TString binRange ="25_43";
  //TString whichFit ="thirteen";
  //TString production ="slot1";//"t3", "slot1"
  //TString additionalCuts ="phiS0.195";

  Bool_t toWrite =false;
  //Setup_______________
  
  TString pathFA="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data";
  if (start==""){//Basic info
    cout << "\nScript draws the input files with the false asymmetry " << endl;
    cout << "and the acceptance ratio together" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'all_FA_Acc.C(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << pathFA << endl;
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
  
  const Int_t nPhysBinned =6;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "M", "xN"};
  //2nd xN is for integrated
  
  //Aesthetics setup
  gStyle->SetLineWidth(2);
  TCanvas* cwFA = new TCanvas(); cwFA->Divide(nPhysBinned, 1, 0, 0.01);
  TCanvas* cwAcc = new TCanvas(); cwAcc->Divide(nPhysBinned, 1, 0, 0.01);
  TCanvas* cwSysStat = new TCanvas();cwSysStat->Divide(nPhysBinned, 1, 0, 0.01);

  //Loop over physics binning
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    
    Int_t nBinsName =nBins;
    if (phys == nPhysBinned-1) {//used for integrated
      nBinsName =1; }
    TString fname;
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      fname=Form("%s/WAvg/wAvg_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
		 pathFA.Data(), whichFit.Data(), Mtype.Data(), process.Data(),
		 lrMrange.Data(), binRange.Data(), physBinned[phys].Data(),
		 nBinsName, hbins, production.Data(), additionalCuts.Data());
    }
    else {
      cout << "Only works for true fit for now" << endl;
      exit(EXIT_FAILURE);
      fname =Form("%s/WAvg/wAvg_%s%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
		  pathFA.Data(), whichFit.Data(), Mtype.Data(),
		  Mtype.Data(), process.Data(),lrMrange.Data(),
		  physBinned[phys].Data(), nBinsName, hbins, production.Data(),
		  additionalCuts.Data() );
    }

    //Get Data
    TFile *f_in = OpenFile(fname);

    //FA
    cwFA->cd(phys+1); gPad->SetFrameLineWidth(2);
    TGraphErrors *g_WAvg_FAsubper = (TGraphErrors*) f_in->Get("falseAN_subper");
    FinalSetup(g_WAvg_FAsubper); FinalClearTitles(g_WAvg_FAsubper);
    g_WAvg_FAsubper->SetMarkerColor(kRed); g_WAvg_FAsubper->SetMarkerStyle(20);
    g_WAvg_FAsubper->Draw("AP");
    DrawLine(g_WAvg_FAsubper, 0.0); DrawLegend(g_WAvg_FAsubper, "FA");

    if (phys==nPhysBinned-1)
      g_WAvg_FAsubper->GetXaxis()->SetLimits(0.09, 0.26);

    //Acceptance
    cwAcc->cd(phys+1); gPad->SetFrameLineWidth(2);
    TGraphErrors *g_WAvg_Acc = (TGraphErrors*) f_in->Get("acc_subper");
    FinalSetup(g_WAvg_Acc); FinalClearTitles(g_WAvg_Acc);
    g_WAvg_Acc->SetMarkerColor(kRed); g_WAvg_Acc->SetMarkerStyle(20);
    g_WAvg_Acc->Draw("AP"); g_WAvg_Acc->GetYaxis()->SetRangeUser(0.95, 1.05);
    DrawLine(g_WAvg_Acc, 1.0); DrawLegend(g_WAvg_Acc, "Acc");

    if (phys==nPhysBinned-1)
      g_WAvg_Acc->GetXaxis()->SetLimits(0.09, 0.26);

    //Sys_Stat
    cwSysStat->cd(phys+1); gPad->SetFrameLineWidth(2);
    TGraph *g_WAvg_Sys_Stat = (TGraph*) f_in->Get("gSys_Stat");
    FinalSetup(g_WAvg_Sys_Stat); FinalClearTitles(g_WAvg_Sys_Stat);
    g_WAvg_Sys_Stat->SetMarkerColor(kRed); g_WAvg_Sys_Stat->SetMarkerStyle(20);
    g_WAvg_Sys_Stat->Draw("AP");
    g_WAvg_Sys_Stat->GetYaxis()->SetRangeUser(0, 0.6);
    DrawLegend(g_WAvg_Sys_Stat, "SysStat");

    if (phys==nPhysBinned-1)
      g_WAvg_Sys_Stat->GetXaxis()->SetLimits(0.09, 0.26);
  }//phys binned loop

  //Write output/final settings
  TString fOutput =
    Form("%s/AllFAAcc/all_FA_Acc_%s_%s_%s%s_%s%ibin_%ihbin_%s_%s.root",
	 pathFA.Data(), whichFit.Data(), Mtype.Data(), process.Data(),
	 lrMrange.Data(), binRange.Data(), nBins, hbins, production.Data(),
	 additionalCuts.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    cwFA->Write("cwFA");
    cwAcc->Write("cwAcc");
    cwSysStat->Write("cwSysStat");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathFA << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Mass type considered:   " << Mtype << endl;
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
