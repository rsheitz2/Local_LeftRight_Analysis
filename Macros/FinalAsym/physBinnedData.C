#include "include/helperFunctions.h"

void physBinnedData(TString start=""){
  //Setup_______________
  /*const Int_t nBins =3;//HMDY
    TString Mtype ="HMDY";
    Int_t hbins =150;
    TString process ="DY";//JPsi, psi, DY
    TString lrMrange ="4.30_8.50";
    TString fitMrange ="4.30_8.50";
    TString binRange ="43_85";
    TString whichFit ="true";
    TString production ="slot1";
    TString additionalCuts ="phiS0.53";//*/

  const Int_t nBins =5;//JPsi
  TString Mtype ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.00_5.00";
  TString fitMrange ="2.00_8.50";
  TString binRange ="25_43";
  TString whichFit ="thirteen";
  TString production ="slot1";
  TString additionalCuts ="phiS0.53";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/WAvg";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  //TString physBinned[nPhysBinned] ={"xN"};
  
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
  TCanvas* cAsym = new TCanvas(); cAsym->Divide(nPhysBinned+1, 1, 0, 0.01);
  Double_t yMax =(process=="DY") ? 0.25 : 0.1;
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphErrors *g_AN[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName;
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      AsymName =
	Form("%s/wAvg_%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(),
	     nBins, hbins, production.Data(), additionalCuts.Data());
    }
    else {
      AsymName =
	Form("%s/wAvg_%s%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(), Mtype.Data(),
	     process.Data(), lrMrange.Data(), physBinned[phys].Data(),
	     nBins, hbins, production.Data(), additionalCuts.Data());
    }
        
    TFile *f_AN = OpenFile(AsymName);

    cAsym->cd(phys+1);
    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");
    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    g_AN[phys]->SetTitle("");
    DrawLine(g_AN[phys], 0.0);
  }//phys binned loop

  //Integrated value
  cout << "\nIntegrated valued determined from fit...\n" << endl;
  cAsym->cd(nPhysBinned+1);
  g_AN[0]->Fit("pol0", "0");
  TF1 *f_pol0 = (TF1*)g_AN[0]->GetFunction("pol0");
  Double_t yInt[] ={ f_pol0->GetParameters()[0] };
  Double_t e_yInt[] = { f_pol0->GetParErrors()[0] };
  Double_t xvals[] ={1.0};
  Double_t ex[] = {0.0};
  TGraphErrors *g_AN_int = new TGraphErrors(1, xvals, yInt, ex, e_yInt);
  SetUp(g_AN_int);
  g_AN_int->GetYaxis()->SetRangeUser(-yMax, yMax);
  g_AN_int->GetXaxis()->SetLimits(0.0, 2.0);
  g_AN_int->Draw("AP"); DrawLine(g_AN_int, 0.0);

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

    cAsym->Write();
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

