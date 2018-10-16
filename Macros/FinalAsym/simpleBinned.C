#include "include/helperFunctions.h"

void simpleBinned(TString start=""){
  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  //const Int_t nPhysBinned =2;
  //TString physBinned[nPhysBinned] ={"xF", "pT"};
  
  //Setup_______________
  const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit[nPhysBinned] = {"true", "true", "true", "true"};//*/

  /*const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.90_3.30";
  TString fitMrange ="2.00_7.50";
  TString whichFit[nPhysBinned] = {"eight", "eight", "eight", "eight"};//*/

  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout <<"Script draws Final AN asymmetry WITHOUT systematics binned in phys";
    cout << "\n\n\nUsage:" << endl;
    cout << "root \'simpleBinned.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Asym data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << period_Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "\nWhich fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << physBinned[i] << "   ";}
    cout << " " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << whichFit[i] << " ";
    }
    cout << "\n\nOutput is to be written:     " << toWrite << "\n" << endl;
    exit(EXIT_FAILURE);
  }
  
  //Aesthetics setup
  TCanvas* cAsym = new TCanvas(); cAsym->Divide(nPhysBinned, 1, 0, 0.01);
  Double_t yMax = (nBins==3) ? 0.25 : 0.1;

  //Systematic error setup
  TGraphAsymmErrors *g_sys[nPhysBinned];
  Double_t yvals[nBins], ey[nBins] = {0.0};
  Double_t ysys = (nBins==3) ? -0.15 : -0.075;
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphErrors *g_AN[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName;
    if (whichFit[phys] == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      AsymName =
	Form("%s/GeoMean4Targ_%s_%s_%s%s_%s%i.root", pathAN.Data(),
	     whichFit[phys].Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);
    }
    else {
      AsymName =
	Form("%s/GeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbin.root",
	     pathAN.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);
    }

    TFile *f_AN = TFile::Open(AsymName);
    if ( !f_AN ){
      cout << "Asymmetries file does not exist"<<endl;
      exit(EXIT_FAILURE);
    }
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit[phys]+" ";

    cAsym->cd(phys+1);
    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");
    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    g_AN[phys]->SetTitle("");
    DrawLine(g_AN[phys], 0.0);

    //tmp
    /*Double_t *e_yAN = g_AN[phys]->GetEY();
    for (Int_t i=0; i<nBins; i++) {
      Double_t divide = (nBins == 3) ? 0.01+0.01*3 : 0.001 + 0.001*3;
      cout << divide/e_yAN[i] << "  ";
    }
    cout << "\n";//tmp*/

    //Systematic error
    Double_t *xvals = g_AN[phys]->GetX();
    Double_t exR[nBins], exL[nBins];
    Double_t sysErr[nBins];
    for (Int_t i=0; i<nBins; i++) {
      yvals[i] = ysys;
      if (i==0) {exL[i] = 0.0; }
      else{ exL[i] = (xvals[i] - xvals[i-1])/2.0; }
      
      if (i==nBins-1){exR[i] = 0.0;}
      else{ exR[i] = (xvals[i+1] - xvals[i])/2.0; }

      if (nBins ==3) {sysErr[i] = (0.01+0.01*3)/2.0; }
      else {sysErr[i] = (0.001+0.001*3)/2.0; }
    }
    g_sys[phys] = new TGraphAsymmErrors(nBins, xvals, yvals, exL, exR,
					ey, sysErr);
    SetUp(g_sys[phys]); g_sys[phys]->Draw("2same");
    g_sys[phys]->SetFillColor(kRed); g_sys[phys]->SetFillStyle(3002);
  }//phys binned loop

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/binned4TargGeoMean";
  TString fOutput;
  if (whichFit[0] == "true"){
    fOutput =
      Form("%s/binned4TargGeoMean_true_%s_%s%s_%ibins.root", thisDirPath.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(), nBins);
  }
  else {
    fOutput =
      Form("%s/binned4TargGeoMean_%s_%s_%s%s_%ibins_%ihbin.root",
	   thisDirPath.Data(), fitMrange.Data(), period_Mtype.Data(),
	   process.Data(), lrMrange.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TNamed pBinNam ("physBinned", physBinnedNames.Data());
    TNamed fitNam ("fitNames", fitNames.Data());
    pBinNam.Write();
    fitNam.Write();

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
    cout << "Mass type considered:   " << period_Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << physBinned[i] << ":   " << whichFit[i] << endl; }
  }
  if (toWrite) cout << "\nFile:  " << fOutput << "   was written" << endl;
  else cout << "\nFile: " << fOutput << " was NOT written" << endl;
}

