#include "include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void DrawLegend(TGraphErrors *g);

void FinalSetupLocal(TGraphErrors *g, TString xName, Double_t yMax,
		     Bool_t same=false);

void FinalLocalIntegrated(TGraphErrors *g);

void physBinnedData(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  /*const Int_t nBins =4;//JPsi
  TString Mtype ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.87_3.38";//"3.08_3.17"; //"2.87_3.38";
  TString fitMrange =lrMrange;
  TString binRange ="29_34";//"31_32"; //"29_34"; //"25_43";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/WAvg";
  const Int_t nPhysBinned =6;
  TString physBinned[nPhysBinned] ={"xN", "xN", "xPi", "xF", "pT", "M"};
  //First xN used for integrated values
    
  if (start==""){
    cout <<"Script draws Final AN asymmetry with systematics binned in physics";
    cout << "\nSystematics are input as a TGraph for each physics bin";
    cout << "\n\n\nUsage:" << endl;
    cout << "root \'binned4TargGeoMean.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Asym data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " <<  whichFit << endl;
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  
  //Aesthetics setup
  TCanvas* cAsym = new TCanvas();
  cAsym->SetLeftMargin(0.2); cAsym->Divide(nPhysBinned, 1, 0, 0);
  Double_t yMax =(process=="DY") ? 0.25 : 0.06;
  TCanvas* cAsym_1targ= new TCanvas();
  cAsym_1targ->SetLeftMargin(0.2); cAsym_1targ->Divide(nPhysBinned, 1, 0, 0);
  TString xNames[nPhysBinned] =
    {"", "x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphErrors *g_AN[nPhysBinned];
  TGraphErrors *g_AN_upS[nPhysBinned], *g_AN_downS[nPhysBinned];;
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName;

    Bool_t integrated =false;
    if (phys == 0) 
      integrated =true;
    Int_t nBinsName = (integrated) ? 1 : nBins;
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      AsymName =
	Form("%s/wAvg_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), Mtype.Data(), process.Data(),
	     lrMrange.Data(), binRange.Data(), physBinned[phys].Data(),
	     nBinsName, hbins, production.Data(), additionalCuts.Data());
    }
    else {
      AsymName =
	Form("%s/wAvg_%s%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(), Mtype.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),
	     physBinned[phys].Data(), nBinsName, hbins, production.Data(),
	     additionalCuts.Data());
    }
        
    TFile *f_AN = OpenFile(AsymName);

    cAsym->cd(phys+1); gPad->SetFrameLineWidth(2);
    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");
    g_AN[phys]->SetName("AN");
    FinalSetupLocal(g_AN[phys], xNames[phys], yMax);
    if (integrated) { FinalLocalIntegrated(g_AN[phys]); }

    cAsym_1targ->cd(phys+1); gPad->SetFrameLineWidth(2);
    g_AN_upS[phys] =(TGraphErrors*)f_AN->Get("AN_ups");
    g_AN_upS[phys]->SetName("AN_upS");
    FinalSetupLocal(g_AN_upS[phys], xNames[phys], 2*yMax);
    g_AN_upS[phys]->GetYaxis()->SetNdivisions(505);
    if (integrated) { FinalLocalIntegrated(g_AN_upS[phys]); }
    g_AN_upS[phys]->SetMarkerStyle(22);

    g_AN_downS[phys] =(TGraphErrors*)f_AN->Get("AN_downs");
    g_AN_downS[phys]->SetName("AN_downs");
    FinalSetupLocal(g_AN_downS[phys], xNames[phys], 2*yMax, true);
    g_AN_downS[phys]->GetYaxis()->SetNdivisions(505);
    if (integrated) { FinalLocalIntegrated(g_AN_downS[phys]); }
    g_AN_downS[phys]->SetMarkerColor(kBlue);
    g_AN_downS[phys]->SetMarkerStyle(23);
    DrawLegend(g_AN_downS[phys]);
  }//phys binned loop
  
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/PhysBinnedData";
  TString fOutput;
  if (whichFit=="true"){
    fOutput = Form("%s/physBinnedData_%s_%s_%s%s_%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), Mtype.Data(),
		   process.Data(), lrMrange.Data(), nBins,
		   hbins, production.Data(), additionalCuts.Data());
  }
  else{
    fOutput = Form("%s/physBinnedData_%s%s_%s_%s%s_%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		   Mtype.Data(), process.Data(), lrMrange.Data(),
		   nBins, hbins, production.Data(), 
		   additionalCuts.Data());
  }

  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    for (Int_t i=0; i<nPhysBinned; i++) {
      g_AN[i]->Write(Form("AN_%s", physBinned[i].Data() ));
    }

    cAsym->Write("Asym2Targ");
    cAsym_1targ->Write("Asym1Targs");
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Asymmetry data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " <<  whichFit << endl;
  }
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
}


void DrawLegend(TGraphErrors *g){
  TString name = g->GetName();
  Double_t sigma;
  Double_t avg = WeightedAvgAndError(g, &sigma);

  TLegend *leg = new TLegend(0.25,0.9,0.7,0.99);
  leg->AddEntry(name, Form("#bar{A} = %0.2f #pm %0.2f", avg, sigma),
		"p");
  leg->SetBorderSize(0); leg->SetTextFont(133); leg->SetTextSize(25);
  leg->Draw("same");
}


void FinalSetupLocal(TGraphErrors *g, TString xName, Double_t yMax,
		     Bool_t same=false){
  (same) ? g->Draw("Psame") : g->Draw("AP");
  g->GetYaxis()->SetRangeUser(-yMax, yMax);
  
  FinalSetup(g);

  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.8);
  g->SetFillStyle(0);
  g->SetLineWidth(2);
  g->SetTitle("");

  SetTitleName(g, xName, "x");

  DrawLine(g, 0.0);
  if (!same) DrawLegend(g);
}


void FinalLocalIntegrated(TGraphErrors *g){
  SetTitleName(g, "A_{lr}", "y");
  g->GetXaxis()->SetLimits(0.0, 0.35);
  g->GetXaxis()->SetLabelSize(0.0);
  g->GetXaxis()->SetTickSize(0.0);

  DrawLine(g, 0.0);
}
