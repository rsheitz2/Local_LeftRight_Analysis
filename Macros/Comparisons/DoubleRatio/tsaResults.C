#include "include/helperFunctions.h"

void DrawLegend(TGraphErrors *g);

void tsaResults(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//# of physBinned bins
  const Int_t nHbins =8;
  TString Mtype ="HMDY";
  TString production ="slot1";//"t3", "slot1"
  
  Bool_t toWrite =false;
  //Setup_______________
  
  //Basic Setup
  TString pathDR="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio/Data";
  const Int_t nPhys =6;//1st xN=integrated
  TString physBinned[nPhys] ={"xN", "xN", "xPi", "xF", "pT", "M"};
  const Int_t nTSA =3;
  TString whichTSA[nTSA] = {"Siv", "Pretz", "Trans"};

  TCanvas* cTSA = new TCanvas();
  cTSA->SetLeftMargin(0.2); cTSA->Divide(nPhys, nTSA, 0, 0);
  for (Int_t tsa=0, ipad=1; tsa<nTSA; tsa++) {
    for (Int_t phys=0; phys<nPhys; phys++, ipad++) {

      Bool_t integrated =false;
      if (phys == 0) 
	integrated =true;
      Int_t nBinsName = (integrated) ? 1 : nBins;
      /*cout << "Using WAll, not weighted average" << endl;
      TString fname =
	Form("%s/DoubleRatio/doubleRatio_WAll_%s_%s%i_%ihbins_%s_%s.root",
	     pathDR.Data(), Mtype.Data(), physBinned[phys].Data(), nBinsName,
	     nHbins, production.Data(), whichTSA[tsa].Data());//*/
      TString fname = Form("%s/WAvg/wAvg_%s_%s%i_%ihbins_%s_%s.root",
			   pathDR.Data(), Mtype.Data(), physBinned[phys].Data(),
			   nBinsName, nHbins, production.Data(),
			   whichTSA[tsa].Data());//*/
      TFile *fdata = OpenFile(fname);
      TGraphErrors *gAmp = (TGraphErrors*) fdata->Get("Amp");

      cTSA->cd(ipad);
      gAmp->Draw("AP"); gAmp->SetTitle("");
      gAmp->SetMarkerColor(kBlue);
      DrawLegend(gAmp);
      if (integrated){
	gAmp->GetXaxis()->SetLimits(0.0, 0.35);
      }
      gAmp->SetMarkerStyle(20);
      DrawLine(gAmp, 0.0);

      if(whichTSA[tsa] == "Siv")
	gAmp->GetYaxis()->SetRangeUser(-0.5, 0.5);
      else if(whichTSA[tsa] == "Pretz")
	gAmp->GetYaxis()->SetRangeUser(-0.8, 0.8);
      else if(whichTSA[tsa] == "Trans")
	gAmp->GetYaxis()->SetRangeUser(-0.8, 0.8);
    }  
  }

  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio/Data/TSA";
  TString fOutput =
    Form("%s/tsa_%s_%ibins_%ihbins_%s.root", thisDirPath.Data(),
	 Mtype.Data(), nBins, nHbins, production.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    cTSA->Write("cTSA");
  }

  cout << "\nSettings______" << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Fit in nHbins:              " << nHbins << endl;
  cout << "Mass range considered:      " << Mtype << endl;
  cout << "Production considered:      " << production << endl;
  cout << "\nTo write output file:       " << toWrite << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}

void DrawLegend(TGraphErrors *g){
  g->SetName("Atsa");
  Double_t sigma;
  Double_t avg = WeightedAvgAndError(g, &sigma);

  TLegend *leg = new TLegend(0.1,0.9,0.7,0.99);
  leg->AddEntry("Atsa", Form("#bar{A} = %0.2f #pm %0.2f", avg, sigma), "p");
  leg->SetBorderSize(0); leg->SetTextFont(132); leg->SetTextSize(0.08);
  leg->Draw("same");
}
