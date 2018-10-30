#include "include/helperFunctions.h"


void CalTotalSysErrors(Double_t *sysErr, vector<TGraph*> &vec,
		       Int_t nBins){
  for (vector<TGraph*>::iterator it=vec.begin(); it!=vec.end(); it++){
    Double_t *yvals = (*it)->GetY();
    for (Int_t bi=0; bi<nBins; bi++) {
      sysErr[bi] += yvals[bi]*yvals[bi];
    }
  }

  for (Int_t bi=0; bi<nBins; bi++) {
    sysErr[bi] = TMath::Sqrt(sysErr[bi]);
  }
}


void FourTarg_inputSystematics(TString start=""){
  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";
  TString pathSys = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  
  //Setup_______________
  /*const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit[nPhysBinned] = {"true", "true", "true", "true"};//*/

  const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.00_5.00";
  TString fitMrange ="2.00_7.50";
  TString whichFit[nPhysBinned] = {"eight", "eight", "eight", "eight"};//*/

  //Which systematics to include
  Bool_t alphaAcc =true;
  
  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout <<"Script draws Final AN asymmetry with systematics binned in physics";
    cout << "\nSystematics are input as a TGraph for each physics bin";
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
    cout << "    Acceptance alpha      " << alphaAcc << endl;
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  
  //Aesthetics setup
  TCanvas* cAsym = new TCanvas(); cAsym->Divide(nPhysBinned, 1, 0, 0.01);
  Double_t yMax =(process=="DY") ? 0.25 : 0.1;
  Double_t ysys =(process=="DY") ? -0.15 : -0.075;
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphAsymmErrors *g_sys[nPhysBinned];
  TGraphErrors *g_AN[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName;
    TString alphaAccName;//Systematic file names
    if (whichFit[phys] == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      AsymName =
	Form("%s/GeoMean4Targ_%s_%s_%s%s_%s%i.root", pathAN.Data(),
	     whichFit[phys].Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);

      TString specficPath =
	"FalseAsym/Data/acceptanceFourTargRatio/SystematicError";
      alphaAccName =
	Form("%s/%s/accSys4TargRatio_true_%s_%s%s_%s%i.root", pathSys.Data(),
	     specficPath.Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);
    }
    else {
      AsymName =
	Form("%s/GeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbin.root",
	     pathAN.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);

      TString specficPath =
	"FalseAsym/Data/acceptanceFourTargRatio/SystematicError";
      alphaAccName =
	Form("%s/%s/accSys4TargRatio_%s%s_%s_%s%s_%s%i_%ihbins.root",
	     pathSys.Data(), specficPath.Data(), whichFit[phys].Data(),
	     fitMrange.Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins, hbins);
    }
        
    TFile *f_AN = OpenFile(AsymName);
    TFile *f_alphaAcc = (alphaAcc) ? OpenFile(alphaAccName) : NULL;

    cAsym->cd(phys+1);
    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");
    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    g_AN[phys]->SetTitle("");
    DrawLine(g_AN[phys], 0.0);

    //Calculate and Draw systematic error bars
    TGraph *g_alphaAcc =
      (alphaAcc) ? (TGraph*)f_alphaAcc->Get("gSys") : NULL;
    
    vector<TGraph*> vec_sys;
    if (alphaAcc) vec_sys.push_back(g_alphaAcc);

    Double_t *xvals = g_AN[phys]->GetX();
    Double_t yvals[nBins], eyb[nBins] = {0.0};
    Double_t exR[nBins], exL[nBins];
    for (Int_t i=0; i<nBins; i++) {
      yvals[i] = ysys;
      if (i==0) {exL[i] = 0.0; }
      else{ exL[i] = (xvals[i] - xvals[i-1])/2.0; }
      
      if (i==nBins-1){exR[i] = 0.0;}
      else{ exR[i] = (xvals[i+1] - xvals[i])/2.0; }
    }

    Double_t sysErr[nBins] = {0.0};
    CalTotalSysErrors(sysErr, vec_sys, nBins);

    g_sys[phys] = new TGraphAsymmErrors(nBins, xvals, yvals, exL, exR,
					eyb, sysErr);
    SetUp(g_sys[phys]);
    //g_sys[phys]->Draw("3Same");
    g_sys[phys]->Draw("2same");
    g_sys[phys]->SetFillColor(kRed);
    g_sys[phys]->SetFillStyle(3002);
  }//phys binned loop

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/FourTarg_inputSystematics";
  TString fOutput;
  if (whichFit[0] == "true"){
    fOutput =
      Form("%s/FourTarg_inputSystematics_true_%s_%s%s_%ibins.root",
	   thisDirPath.Data(), period_Mtype.Data(), process.Data(),
	   lrMrange.Data(), nBins);
  }
  else {
    fOutput =
      Form("%s/FourTarg_inputSystematics_%s%s_%s_%s%s_%ibins_%ihbin.root",
	   thisDirPath.Data(), whichFit[0].Data(), fitMrange.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    //TNamed pBinNam ("physBinned", physBinnedNames.Data());
    //TNamed fitNam ("fitNames", fitNames.Data());
    //pBinNam.Write();
    //fitNam.Write();

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

