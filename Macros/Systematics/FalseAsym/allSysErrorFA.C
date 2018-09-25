#include "include/helperFunctions.h"


void OffSet(TGraphErrors *g, Double_t offset){
  Double_t *xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xval[i] += offset;
}


void FillValues(TH1D *h, TGraphErrors *g){
  Double_t *yVals = g->GetY();
  
  for (Int_t i=0; i<g->GetN(); i++){
    h->Fill(yVals[i]);
  }
}


void CalWavg(vector<TGraphErrors*> &v_FA, Double_t *wAvg, Double_t *wErr){
  Int_t nBins=0, nFA=0;
  for(vector<TGraphErrors*>::iterator it=v_FA.begin(); it!=v_FA.end(); it++){
    Double_t *yFA = (*it)->GetY();
    Double_t *e_yFA = (*it)->GetEY();

    if (nBins==0) nBins = (*it)->GetN();
    else if ((*it)->GetN() != nBins ) {
      cout << "Binning number inconsistency" << endl;
      exit(EXIT_FAILURE);
    }
    
    for (Int_t bi=0; bi<nBins; bi++) {
      wAvg[bi] += yFA[bi]/(e_yFA[bi]*e_yFA[bi]);
      wErr[bi] += 1.0/(e_yFA[bi]*e_yFA[bi]);
    }
  }

  for (Int_t bi=0; bi<nBins; bi++) {
    wAvg[bi] /= wErr[bi];
    wErr[bi] = TMath::Sqrt( 1.0/wErr[bi] );
  }
}


void allSysErrorFA(TString start=""){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/\
Data";
  //const Int_t nPhysBinned =4;
  //TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  const Int_t nPhysBinned =2;
  TString physBinned[nPhysBinned] ={"xF", "pT"};
  
  //Setup_______________
  const Int_t nBins =5;
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.90_3.30";
  TString fitMrange ="1.00_8.50";
  //TString whichFit[] ={"true", "true", "true", "true"};
  //TString whichFit[] ={"seven", "seven", "six", "six"};
  TString whichFit[] ={"eight", "eight"};

  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout << "Script draws false asymmetries and systematic error from false ";
    cout << " asymmetries per period and physics binning on a nice plot";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C  ->  GeoMean4Targ.C  ->";
    cout << "falseGeoMean4Targ_targFlip.C/falseGeoMean4Targ_splitTarg.C  -> ";
    cout << " sysErrorFA.C  ->  allSysError.C";
    cout << "\n\nUsage:" << endl;
    cout << "root \'allSysError.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << period_Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) 
      cout << whichFit[i] << " ";
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  
  //Aesthetics setup
  //TCanvas* cFA = new TCanvas(); cFA->Divide(4, 1, 0, 0.01);
  TCanvas* cFA = new TCanvas(); cFA->Divide(2, 1, 0.01, 0.01);
  //TCanvas* cSysFA = new TCanvas(); cSysFA->Divide(4, 1, 0, 0.01);
  TCanvas* cSysFA = new TCanvas(); cSysFA->Divide(2, 1, 0.01, 0.01);
  //TCanvas* cWavg = new TCanvas(); cWavg->Divide(4, 1, 0, 0.01);
  TCanvas* cWavg = new TCanvas(); cWavg->Divide(2, 1, 0.01, 0.01);
  //Double_t offsets[nPhysBinned] = {0.006, 0.01, 0.01, 0.05};
  Double_t offsets[nPhysBinned] = {0.01, 0.05};
  Double_t yMax =0.2;
  
  //Get Data file/Get graphs and plot
  TH1D* hSys = new TH1D("hSys", "hSys", 5, 0, 0.60);
  TString physBinnedNames ="", fitNames="";
  Double_t FA_wAvg[nPhysBinned][nBins], eFA_wAvg[nPhysBinned][nBins];
  Double_t ex[nBins] = {0.0};
  TGraphErrors *gFA_wAvg[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString inFileFA, inSplitTarg, sysName;
    if (whichFit[phys] == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      inFileFA =
	Form("%s/TargFlip/falseGeoMean4Targ_true_%s_%s%s_%s%i.root",path.Data(),
	     period_Mtype.Data(), process.Data(), lrMrange.Data(),
	     physBinned[phys].Data(), nBins);
      inSplitTarg =
	Form("%s/SplitTarg/falseGeoMeanSplitTarg_true_%s_%s%s_%s%i.root",
	     path.Data(), period_Mtype.Data(), process.Data(), lrMrange.Data(),
	     physBinned[phys].Data(), nBins);
      sysName =
	Form("%s/sysError/sysErrorFA_true_%s_%s%s_%s%i.root", path.Data(),
	     period_Mtype.Data(), process.Data(), lrMrange.Data(),
	     physBinned[phys].Data(), nBins);
    }
    else {
      inFileFA =
	Form("%s/TargFlip/falseGeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbins.root",
	     path.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);
      inSplitTarg =
	Form("%s/SplitTarg/falseGeoMeanSplitTarg_%s%s_%s_%s%s_%s%i_%ihbins.root",
	     path.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(), lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);
      sysName =
	Form("%s/sysError/sysErrorFA_%s%s_%s_%s%s_%s%i_%ihbin.root", 
	     path.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);
    }
    
    TFile *f_FA = TFile::Open(inFileFA);
    TFile *fFA_sT = TFile::Open(inSplitTarg);
    TFile *f_sys = TFile::Open(sysName);
    if (!f_FA || !fFA_sT || !f_sys){
      cout << "False asymmetries or systematic error file does not exist"<<endl;
      exit(EXIT_FAILURE);
    }
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit[phys]+" ";

    cFA->cd(phys+1);
    //     TargFlip
    TGraphErrors *g_FA_pol =(TGraphErrors*)f_FA->Get("falseAN_pol");
    TGraphErrors *g_FA_subper =(TGraphErrors*)f_FA->Get("falseAN_subper");
    OffSet(g_FA_subper, offsets[phys]);
    //     SplitTarg
    TGraphErrors *gsT_sb1_center= (TGraphErrors*)fFA_sT->Get("sb1_center");
    TGraphErrors *gsT_sb2_center= (TGraphErrors*)fFA_sT->Get("sb2_center");
    TGraphErrors *gsT_sb1_Sup= (TGraphErrors*)fFA_sT->Get("sb1_Sup");
    TGraphErrors *gsT_sb2_Sup= (TGraphErrors*)fFA_sT->Get("sb2_Sup");
    vector<TGraphErrors*> vecFA{g_FA_subper,
      gsT_sb1_center, gsT_sb2_center, gsT_sb1_Sup, gsT_sb2_Sup};
    //vector<TGraphErrors*> vecFA{g_FA_pol, g_FA_subper};
    
    OffSet(gsT_sb1_center, offsets[phys]);
    OffSet(gsT_sb2_center, offsets[phys]);
    OffSet(gsT_sb1_Sup, offsets[phys]);
    OffSet(gsT_sb2_Sup, offsets[phys]);
    gsT_sb2_center->SetMarkerColor(7);
    gsT_sb2_Sup->SetMarkerColor(28);
    
    //g_FA_pol->Draw("AP");
    //g_FA_subper->Draw("Psame");
    g_FA_subper->Draw("AP"); g_FA_subper->GetYaxis()->SetRangeUser(-0.3, 0.3);
    gsT_sb1_center->Draw("Psame"); gsT_sb2_center->Draw("Psame");
    gsT_sb1_Sup->Draw("Psame"); gsT_sb2_Sup->Draw("Psame");
    g_FA_pol->SetTitle("False Asym");
    DrawLine(g_FA_pol, 0.0);

    cSysFA->cd(phys+1);
    TGraphErrors *g_sys =(TGraphErrors*)f_sys->Get("gSys");
    g_sys->Draw("AP"); g_sys->SetTitle("Sys Error/Stat Error");
    g_sys->GetYaxis()->SetRangeUser(0, 0.6);

    FillValues(hSys, g_sys);

    //Weighted average calculation
    for (Int_t bi=0; bi<nBins; bi++) {
      FA_wAvg[phys][bi] = 0.0; eFA_wAvg[phys][bi] = 0.0;
    }
    CalWavg(vecFA, FA_wAvg[phys], eFA_wAvg[phys]);

    Double_t *xvals = g_FA_pol->GetX();
    gFA_wAvg[phys] = new TGraphErrors(nBins, xvals, FA_wAvg[phys],
				      ex, eFA_wAvg[phys]);
    SetUp(gFA_wAvg[phys]); gFA_wAvg[phys]->GetYaxis()->SetRangeUser(-0.2, 0.2);
    cWavg->cd(phys+1);
    gFA_wAvg[phys]->Draw("AP"); gFA_wAvg[phys]->SetTitle("WAvg");
  }//phys binned loop

  TCanvas* cDist = new TCanvas();
  hSys->Draw("E");
  SetUp(hSys);

    
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/allSysErrorFA";
  TString fOutput;
  if (whichFit[0] == "true"){
    fOutput =
      Form("%s/allSysErrorFA_true_%s_%s%s_%ibins.root", thisDirPath.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(), nBins);
  }
  else {
    fOutput =
      Form("%s/allSysErrorFA_%s_%s_%s%s_%ibins_%ihbin.root", thisDirPath.Data(),
	   fitMrange.Data(), period_Mtype.Data(), process.Data(),
	   lrMrange.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TNamed pBinNam ("physBinned", physBinnedNames.Data());
    TNamed fitNam ("fitNames", fitNames.Data());
    pBinNam.Write();
    fitNam.Write();
    
    cFA->Write();
    cSysFA->Write();
    hSys->Write();
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << period_Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) 
      cout << whichFit[i] << " ";
  }
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
  cout << "\nTargFlip pol false asym was not used!!!" << endl;
}

