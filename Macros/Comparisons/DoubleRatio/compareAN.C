#include "include/helperFunctions.h"

void compareAN(TString start=""){
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

  /*const Int_t nBins =5;//JPsi
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
  TString pathDR = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio/Data/WAvg";
  const Int_t nPhysBinned =5;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "xN"};
  
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
  TCanvas* cAsym = new TCanvas(); cAsym->Divide(nPhysBinned, 1, 0, 0.01);
  Double_t yMax =(process=="DY") ? 0.5 : 0.1;
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphErrors *g_AN[nPhysBinned], *g_DR[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName, DRName;
    Int_t binsName =nBins;
    if (phys == nPhysBinned-1){ binsName =1; }//Integrated
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      AsymName =
	Form("%s/wAvg_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), Mtype.Data(), process.Data(),
	     lrMrange.Data(), binRange.Data(), physBinned[phys].Data(),
	     binsName, hbins, production.Data(), additionalCuts.Data());

      Int_t nHbins =8;
      if (phys ==1){
	cout << "\n" << nHbins << " bins used for doubleratio\n" << endl;
      }
      DRName =
	Form("%s/wAvg_%s_%s%i_%ihbins_%s.root", pathDR.Data(), Mtype.Data(),
	     physBinned[phys].Data(), binsName, nHbins, production.Data());
    }
    else {
      cout << "Doesn't work" << endl;
      exit(EXIT_FAILURE);
    }
        
    TFile *f_AN = OpenFile(AsymName);
    TFile *f_DR = OpenFile(DRName);

    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");
    g_DR[phys] =(TGraphErrors*)f_DR->Get("Amp");

    cAsym->cd(phys+1);
    Double_t *yAN =g_AN[phys]->GetY();
    Double_t *e_yAN =g_AN[phys]->GetEY();
    for (Int_t bi=0; bi<g_AN[phys]->GetN(); bi++) {
      yAN[bi] = yAN[bi]*TMath::Pi()/2.0;
      e_yAN[bi] = e_yAN[bi]*TMath::Pi()/2.0;
    }

    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    g_AN[phys]->SetTitle(""); g_AN[phys]->SetMarkerColor(kBlue);
    DrawLine(g_AN[phys], 0.0);

    g_DR[phys]->SetMarkerColor(kRed);
    g_DR[phys]->Draw("Psame");
    
    Double_t offset;
    if (physBinned[phys] =="xF") offset = 0.003;
    else if (physBinned[phys] =="pT") offset = 0.01;
    else if (physBinned[phys] =="xN") offset = 0.0009;
    else if (physBinned[phys] =="xPi") offset = 0.002;
    offset *= 5.0;
    OffSet(g_DR[phys], offset);
  }//phys binned loop
  
  //Write Output/Final Settings
  /*TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis \
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
  else cout << "File: " << fOutput << " was NOT written" << endl;*/
}

