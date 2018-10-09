#include "include/helperFunctions.h"

void CalAbsSys(TGraphErrors *gFA_1, TGraphErrors *gFA_2, Double_t *systematics){
  Double_t *yFA_1 = gFA_1->GetY();
  Double_t *yFA_2 = gFA_2->GetY();

  Double_t *e_yFA_1 = gFA_1->GetEY();
  Double_t *e_yFA_2 = gFA_2->GetEY();

  if (gFA_1->GetN() != gFA_2->GetN() ) {
    cout << "Binning number inconsistency" << endl;
    exit(EXIT_FAILURE);
  }

  for (Int_t i=0; i<gFA_1->GetN(); i++) {
    Double_t sys_1 = yFA_1[i]/e_yFA_1[i];
    if (sys_1 < 0.0) sys_1 *= -1.0;
    if (sys_1 < 0.68) sys_1 = 0;
    else sys_1 = TMath::Sqrt(sys_1*sys_1 - 0.68*0.68);

    Double_t sys_2 = yFA_2[i]/e_yFA_2[i];
    if (sys_2 < 0.0) sys_2 *= -1.0;
    if (sys_2 < 0.68) sys_2 = 0;
    else sys_2 = TMath::Sqrt(sys_2*sys_2 - 0.68*0.68);

    systematics[i] = (sys_1 + sys_2)/2.0;
  }
}//CalAbsSys


void CalAbsSys(vector<TGraphErrors*> &v_FA, Double_t *systematics){
  Int_t nBins=0;
  Int_t nFA=0;
  for(vector<TGraphErrors*>::iterator it=v_FA.begin(); it!=v_FA.end(); it++){
    Double_t *yFA = (*it)->GetY();
    Double_t *e_yFA = (*it)->GetEY();

    if (nBins==0) nBins = (*it)->GetN();
    else if ((*it)->GetN() != nBins ) {
      cout << "Binning number inconsistency" << endl;
      exit(EXIT_FAILURE);
    }

    for (Int_t bi=0; bi<nBins; bi++) {
      Double_t sys = yFA[bi]/e_yFA[bi];
      
      if (sys < 0.0) sys *= -1.0;
      if (sys < 0.68) sys = 0.0;
      else sys = TMath::Sqrt(sys*sys - 0.68*0.68);

      systematics[bi] += sys;
    }

    nFA++;
  }
  
  for (Int_t bi=0; bi<nBins; bi++) systematics[bi] /= nFA;
}//CalAbsSys


void SysError(TGraphErrors* g_FA, Double_t *sysErr){
  Double_t *yFA = g_FA->GetY();
  Double_t *e_yFA = g_FA->GetEY();
  Int_t nBins = g_FA->GetN();
  for (Int_t bi=0; bi<nBins; bi++) {
    Double_t sys = yFA[bi]/e_yFA[bi];
      
    if (sys < 0.0) sys *= -1.0;
    if (sys < 0.68) sys = 0.0;
    else sys = TMath::Sqrt(sys*sys - 0.68*0.68);

    sysErr[bi] += sys;
  }
}//SysError


void OffSet(TGraphErrors *g, Double_t xMin, Double_t xMax, Int_t i){
  Double_t offset = i*(xMax - xMin)/100.0;

  Double_t *xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xval[i] += offset;
}


void OffSet(TGraph *g, Int_t i){
  Double_t xMin = g->GetXaxis()->GetXmin();
  Double_t xMax = g->GetXaxis()->GetXmax();
  Double_t offset = i*(xMax - xMin)/100.0;

  Double_t *xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xval[i] += offset;
}


void DrawFA(vector<TGraphErrors*> &v_FA){
  Int_t i =0, nBins;
  Double_t xMin, xMax;
  for(vector<TGraphErrors*>::iterator it=v_FA.begin(); it!=v_FA.end(); it++){
    if (i==0) {
      (*it)->Draw("AP");
      //(*it)->GetYaxis()->SetRangeUser(-0.5, 0.5);
      (*it)->GetYaxis()->SetRangeUser(-0.25, 0.25);
      (*it)->SetTitle("False Asymmetries");
      xMin = (*it)->GetXaxis()->GetXmin();
      xMax = (*it)->GetXaxis()->GetXmax();
      nBins = (*it)->GetN();
    }
    else {
      OffSet((*it), xMin, xMax, i);
      (*it)->Draw("Psame");
    }
    
    i++;
  }
}


void sysErrorFA(TString start =""){
  //Setup_______________
  const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.90_3.30";
  TString fitMrange ="2.00_7.50";
  TString whichFit ="ten";

  Bool_t toWrite =false;
  //Setup_______________

  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/";
  
  if (start==""){//Basic info
    cout << "\nScript calculates systematic error due false AN asymmetries\n";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C or trueCount.C ";
    cout <<"-> fasleGeoMean4Targ_targFlips.C  -> sysErrorFA.C" << endl;
    cout << "systematic_leftRight  ->  (for mass fitting only) functMFit.C ";
    cout <<"-> falseGeoMean4Targ_splitTarg.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'sysErrorFA(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "False asymmetry data coming from:            " << pathFA << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Period and Mass type considered:   " << period_Mtype << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
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

  //File name setup && get file
  TString inFileFA, inSplitTarg, inRunNum;
  if (whichFit=="true"){
    inFileFA =
      Form("TargFlip/falseGeoMean4Targ_true_%s_%s%s_%s%i.root",
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   physBinned.Data(), nBins);
    
    inSplitTarg =
      Form("SplitTarg/falseGeoMeanSplitTarg_true_%s_%s%s_%s%i.root",
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   physBinned.Data(), nBins);

    inRunNum =
      Form("RunNums/falseGeoMeanRunNums_true_%s_%s%s_%s%i.root",
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   physBinned.Data(), nBins);
  }
  else{
    inFileFA =
      Form("TargFlip/falseGeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbins.root",
	   whichFit.Data(), fitMrange.Data(),period_Mtype.Data(),
	   process.Data(),lrMrange.Data(), physBinned.Data(), nBins, hbins);

    inSplitTarg =
      Form("SplitTarg/falseGeoMeanSplitTarg_%s%s_%s_%s%s_%s%i_%ihbins.root",
	   whichFit.Data(), fitMrange.Data(),period_Mtype.Data(), process.Data(),
	   lrMrange.Data(), physBinned.Data(), nBins, hbins);

    inRunNum =
      Form("RunNums/falseGeoMeanRunNums_%s%s_%s_%s%s_%s%i_%ihbins.root",
	   whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
	   process.Data(), lrMrange.Data(), physBinned.Data(), nBins, hbins);
  }
  
  TFile *fFA = TFile::Open(pathFA+inFileFA);
  TFile *fFA_sT = TFile::Open(pathFA+inSplitTarg);
  TFile *fFA_runNum = TFile::Open(pathFA+inRunNum);
  if (!fFA || !fFA_sT || !fFA_runNum ){
    cout << "One or more false asym files does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data from files
  //     TargFlip
  TGraphErrors *gFA_pol = (TGraphErrors*)fFA->Get("falseAN_pol");
  TGraphErrors *gFA_subper = (TGraphErrors*)fFA->Get("falseAN_subper");
  //     SplitTarg
  TGraphErrors *gsT_sb1_center= (TGraphErrors*)fFA_sT->Get("sb1_center");
  TGraphErrors *gsT_sb2_center= (TGraphErrors*)fFA_sT->Get("sb2_center");
  TGraphErrors *gsT_sb1_Sup= (TGraphErrors*)fFA_sT->Get("sb1_Sup");
  TGraphErrors *gsT_sb2_Sup= (TGraphErrors*)fFA_sT->Get("sb2_Sup");
  //     RunNum
  TGraphErrors *gFA_runNum= (TGraphErrors*)fFA_runNum->Get("FA_RunNums");
  
  //Which False asymmetries to consider
  vector<TGraphErrors*> vecFA{gFA_pol, gFA_subper,
      gsT_sb1_center, gsT_sb2_center, gsT_sb1_Sup, gsT_sb2_Sup,
      gFA_runNum};
    
  //Make FA/statisics average
  Double_t FA_sys[nBins] ={0.0};
  Double_t *xvals = gFA_pol->GetX();
  CalAbsSys(vecFA, FA_sys);

  //Individual sys errors
  Double_t sys_pol[nBins], sys_subper[nBins];
  Double_t sys_sb1_center[nBins], sys_sb2_center[nBins];
  Double_t sys_sb1_Sup[nBins], sys_sb2_Sup[nBins];
  Double_t sys_runNum[nBins];
  SysError(gFA_pol, sys_pol);
  SysError(gFA_subper, sys_subper);
  SysError(gsT_sb1_center, sys_sb1_center);
  SysError(gsT_sb2_center, sys_sb2_center);
  SysError(gsT_sb1_Sup, sys_sb1_Sup);
  SysError(gsT_sb2_Sup, sys_sb2_Sup);
  SysError(gFA_runNum, sys_runNum);
  TGraph *gSys_pol = new TGraph(nBins, xvals, sys_pol);
  TGraph *gSys_subper = new TGraph(nBins, xvals, sys_subper);
  TGraph *gSys_sb1_center = new TGraph(nBins, xvals, sys_sb1_center);
  TGraph *gSys_sb2_center = new TGraph(nBins, xvals, sys_sb2_center);
  TGraph *gSys_sb1_Sup = new TGraph(nBins, xvals, sys_sb1_Sup);
  TGraph *gSys_sb2_Sup = new TGraph(nBins, xvals, sys_sb2_Sup);
  TGraph *gSys_runNum = new TGraph(nBins, xvals, sys_runNum);
  SetUp(gSys_pol); SetUp(gSys_subper); 
  SetUp(gSys_sb1_center); SetUp(gSys_sb2_center);
  SetUp(gSys_sb1_Sup); SetUp(gSys_sb2_Sup);
  SetUp(gSys_runNum);
  OffSet(gSys_pol, 1); OffSet(gSys_subper, 2);
  OffSet(gSys_sb1_center, 3); OffSet(gSys_sb2_center, 4);
  OffSet(gSys_sb1_Sup, 5); OffSet(gSys_sb2_Sup, 6);
  OffSet(gSys_runNum, 7);
  
  //Draw systematic error
  TGraph *gSys = new TGraphErrors(nBins, xvals, FA_sys);
  SetUp(gSys);

  TCanvas* c1 = new TCanvas();
  //gSys->Draw("AP");  gSys->GetYaxis()->SetRangeUser(0.0, 1.5);
  gSys->Draw("AP");  gSys->GetYaxis()->SetRangeUser(0.0, 10);
  gSys->SetTitle("Systematic Errors/Stat Errors");
  gSys_pol->Draw("PSame"); gSys_pol->SetMarkerColor(kBlue);
  gSys_subper->Draw("PSame"); gSys_subper->SetMarkerColor(kRed);
  gSys_sb1_center->Draw("PSame"); gSys_sb1_center->SetMarkerColor(36);
  gSys_sb2_center->Draw("PSame"); gSys_sb2_center->SetMarkerColor(kRed);
  gSys_sb1_Sup->Draw("PSame"); gSys_sb1_Sup->SetMarkerColor(kBlue);
  gSys_sb2_Sup->Draw("PSame"); gSys_sb2_Sup->SetMarkerColor(38);
  gSys_runNum->Draw("PSame"); gSys_runNum->SetMarkerColor(kGreen);
  
  TCanvas* cAll = new TCanvas();
  DrawFA(vecFA);  DrawLine(gFA_pol, 0.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/sysError";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =
      Form("%s/sysErrorFA_true_%s_%s%s_%s%i.root", thisDirPath.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),
	   physBinned.Data(), nBins);
  }
  else{
    fOutput =
      Form("%s/sysErrorFA_%s%s_%s_%s%s_%s%i_%ihbin.root", thisDirPath.Data(),
	   whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
	   process.Data(), lrMrange.Data(), physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    gSys->Write("gSys");
    c1->Write("SysErrors");
    cAll->Write("All_FAs");
  }

  cout << "\nSettings________" << endl;
  cout << "False asymmetry data from:            " << pathFA << endl;
  cout << "False asymmetry file used:                   " << inFileFA << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Period and Mass type considered:   " << period_Mtype << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
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
  cout << "Warning not using pol flip FA" << endl;
}
