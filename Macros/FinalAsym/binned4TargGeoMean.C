#include "include/helperFunctions.h"


void CalTotalSysErrors(Double_t *sysErr, TGraphErrors *g_AN,
		       vector<TGraphErrors*> &vec, Int_t nBins){
  for (vector<TGraphErrors*>::iterator it=vec.begin(); it!=vec.end(); it++){
    
    Double_t *yvals = (*it)->GetY();
    Double_t *e_AN = g_AN->GetEY();
    for (Int_t bi=0; bi<nBins; bi++) {
      sysErr[bi] += e_AN[bi]*e_AN[bi]*yvals[bi]*yvals[bi];
    }
    
  }

  for (Int_t bi=0; bi<nBins; bi++) {
    sysErr[bi] = TMath::Sqrt(sysErr[bi]);
  }
}


void binned4TargGeoMean(TString start=""){
  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";
  TString pathSys = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  
  //Setup_______________
  const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit[nPhysBinned] = {"true", "true", "true", "true"};
  //TString whichFit[nPhysBinned] = {"six", "six", "seven", "seven"};

  Bool_t falseAsym =true;

  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout <<"Script draws Final AN asymmetry with systematics binned in physics";
    cout << "\n\n\nUsage:" << endl;
    cout << "root \'binned4TargGeoMean.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Asym data coming from:            " << pathAN << endl;
    cout << "Systematic data coming from:            " << pathSys << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << period_Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "\nWhich fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << physBinned[i] << "   ";
    }
    cout << " " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << whichFit[i] << " ";
    }
    cout << "\n\nSytemtaics considered:" << endl;
    cout << "    False asymmetries   " << falseAsym << endl;
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  
  //Aesthetics setup
  TCanvas* cAsym = new TCanvas(); cAsym->Divide(4, 1, 0, 0.01);
  Double_t yMax =0.3;
  Double_t ysys =-0.25;
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphAsymmErrors *g_sys[nPhysBinned];
  TGraphErrors *g_AN[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName, FAname;
    if (whichFit[phys] == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      AsymName =
	Form("%s/GeoMean4Targ_%s_%s_%s%s_%s%i.root", pathAN.Data(),
	     whichFit[phys].Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);

      FAname =
	Form("%s/FalseAsym/Data/sysError/sysErrorFA_%s_%s_%s%s_%s%i.root",
	     pathSys.Data(),
	     whichFit[phys].Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);

    }
    else {
      AsymName =
	Form("%s/GeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbin.root",
	     pathAN.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);

      FAname =
	Form("%s/FalseAsym/Data/sysError/sysErrorFA_%s%s_%s_%s%s_%s%i_%ihbin.root",
	     pathSys.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);
    }
    
    TFile *f_AN = TFile::Open(AsymName);
    TFile *f_sysFA = TFile::Open(FAname);
    if ( !f_AN || !f_sysFA ){
      cout << "Asymmetries file or systematic file does not exist"<<endl;
      exit(EXIT_FAILURE);
    }
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit[phys]+" ";

    cAsym->cd(phys+1);
    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");
    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    g_AN[phys]->SetTitle("");
    DrawLine(g_AN[phys], 0.0);

    //Calculate and Draw systematic error bars
    TGraphErrors *g_FA =(TGraphErrors*)f_sysFA->Get("gSys");
    
    vector<TGraphErrors*> vec_sys;
    if (falseAsym) vec_sys.push_back(g_FA);

    Double_t *xvals = g_FA->GetX();
    Double_t yvals[nBins], ex[nBins] = {0.0};
    for (Int_t i=0; i<nBins; i++) {
      yvals[i] = ysys;
    }

    Double_t sysErr[nBins] = {0.0};
    CalTotalSysErrors(sysErr, g_AN[phys], vec_sys, nBins);

    g_sys[phys] = new TGraphAsymmErrors(nBins, xvals, yvals, ex, ex, ex,sysErr);
    SetUp(g_sys[phys]);
    g_sys[phys]->Draw("3Same");
    g_sys[phys]->SetFillColor(kRed);
    g_sys[phys]->SetFillStyle(3002);
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
      g_sys[i]->Write(Form("sys_%s", physBinned[i].Data() ));
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
      cout << physBinned[i] << ":   " << whichFit[i] << endl;
    }
  }
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
}

