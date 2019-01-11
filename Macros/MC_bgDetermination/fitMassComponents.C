#include "Include/helperFunctions.h"

const Int_t nFiles = 5;
TH1D *hist[nFiles]; //{JPsi, psi, OC, AMDY, Real}

Double_t FitMCs(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );
    
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}


Double_t FitMCHMDY(Double_t *x, Double_t *par){
  Double_t nHMDY = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
    
  return par[0]*nHMDY;
}


void fitMassComponents(TString start=""){
  //Setup__________
  TString whichMC ="Charles"; //"Yu", "Charles"
  const Int_t nBins =200;
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
  TString whichRD ="slot1"; //"t3", "slot1"
  Bool_t fromHist =true;
  Bool_t HMDYonly =true;
  
  Bool_t toWrite =true;
  //Setup__________
  
  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);

  if (start=="") {
    cout << "Script determines the amount of each MC component in real data";
    cout << "\nThis is done fitting MC distributions as a function of Mass to";
    cout << " to get the contribution from each physical process" << endl;
    cout << "Script only works for AMDY already made for now..." << endl;
    cout << "\nCurrent Settings:" << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Which MC considered:          " << whichMC << endl;
    cout << "Which real data considered:  " << whichRD << endl;
    cout << "Minimum mass considered:      " << Mmin << endl;
    cout << "Maximum mass considered:      " << Mmax << endl;
    cout << "Number bins in mass hist:     " << nBins << endl;
    cout << "Mass distributions from premade histogram:   " << fromHist << endl;
    cout << "HM Drell-Yan only:            " << HMDYonly << endl;
    cout << "\nTo Write:                   " << toWrite << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'fitMassComponenets.C(1)\'" << endl;
    exit(EXIT_FAILURE);
  }
  
  //Aesthetic setup
  gStyle->SetOptFit(1111);
  Int_t icolor[nFiles] = {6,3,9,4,1};//{JPsi, psi, OC, AMDY, RealData}
  
  TString localData = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Presents/DATA/";
  TString thisDirPath = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/MC_bgDetermination/";
  TString dataPathNames[nFiles];
  if (fromHist) {
    if (whichMC=="Yu"){
      cout << "Does not work for Yu with hist for now..." << endl;
      exit(EXIT_FAILURE);
    }
    else if (HMDYonly){
      TString pathMC =thisDirPath+"Data/MassDist/";
      dataPathNames[0] = pathMC+Form("M_%s_HMDY_%ibins%0.1f_%0.1f.root",
				     whichMC.Data(), nBins, Mmin, Mmax);
      dataPathNames[1] = pathMC+Form("M_%s_%ibins%0.1f_%0.1f.root",
				     whichRD.Data(), nBins, Mmin, Mmax);
    }
    else if (whichMC=="Charles") {
      TString pathMC =thisDirPath+"Data/MassDist/";
      dataPathNames[0] = pathMC+Form("M_%s_Jpsi_%ibins%0.1f_%0.1f.root",
				     whichMC.Data(), nBins, Mmin, Mmax);
      dataPathNames[1] = pathMC+Form("M_%s_Psi_%ibins%0.1f_%0.1f.root",
				     whichMC.Data(), nBins, Mmin, Mmax);
      dataPathNames[2] = pathMC+Form("M_%s_OC_%ibins%0.1f_%0.1f.root",
				     whichMC.Data(), nBins, Mmin, Mmax);
      dataPathNames[3] = pathMC+Form("M_%s_AMDY_%ibins%0.1f_%0.1f.root",
				     whichMC.Data(), nBins, Mmin, Mmax);
      dataPathNames[4] = pathMC+Form("M_%s_%ibins%0.1f_%0.1f.root",
				     whichRD.Data(), nBins, Mmin, Mmax);
    }
    else {
      cout << "MC specified does not exist" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else{
    if (HMDYonly){
      cout << "Not setup for without hist MC fitting for now..." << endl;
      exit(EXIT_FAILURE);
    }
	
    if (whichMC=="Yu"){
      TString pathMC =localData+"MC_Data/YuShiangMC/";
      dataPathNames[0] = pathMC+"JPsi/Yu_Wall_full_main_JPsi_20bins.root";
      dataPathNames[1] = pathMC+"psi/Yu_Wall_full_main_psi_20bins.root";
      dataPathNames[2] = pathMC+"OC/Yu_Wall_full_main_OC_20bins.root";
      dataPathNames[3] = pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root";
    }
    else if (whichMC=="Charles") {
      cout << "Does not work for Charles without hist for now..." << endl;
      exit(EXIT_FAILURE);
    }
    else {
      cout << "MC specified does not exist" << endl;
      exit(EXIT_FAILURE);
    }

    if (whichRD=="t3"){
      dataPathNames[4] = localData+"RealData/LowM_AMDY/WAll_LowM_AMDY.root";    
    }
    else if (whichRD=="slot1"){
      dataPathNames[4] = localData+"RealData/LowM_AMDY/slot1WAll_LowM_AMDY.root";
    }
    else {
      cout << "RD specifed does not exist" << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Open data files && get MC mass histograms
  TFile *fJPsi, *fpsi, *fOC, *fAMDY, *fReal, *fHMDY;
  TCanvas* cResults = new TCanvas(); gPad->SetLogy();
  TRatioPlot *rp;
  TF1 *fitFunc;
  TH1D *histResult[nFiles];//{JPsi, psi, OC, AMDY, Sum}  
  if (HMDYonly){//HM only
    fHMDY = OpenFile(dataPathNames[0]);
    fReal = OpenFile(dataPathNames[1]);
    
    if (fromHist){
      hist[0] = (TH1D*)fHMDY->Get("h_HMDY");
      hist[0]->Sumw2();
      hist[0]->Scale(1/(hist[0]->Integral()) );
      hist[0]->SetLineColor(icolor[0]);
      Setup(hist[0]);

      hist[1] = (TH1D*)fReal->Get(Form("h_%s", whichRD.Data()) );
      Setup(hist[1]);
    }

    //Do fit
    fitFunc = new TF1("fitFunc", FitMCHMDY, Mmin, Mmax, 1);
    hist[1]->Fit("fitFunc", "RWL0", "", Mmin, Mmax);

    histResult[0] = (TH1D*)hist[0]->Clone();
    histResult[0]->Scale(fitFunc->GetParameter(0) );
    histResult[0]->Draw("same");

    rp =new TRatioPlot(histResult[0], hist[1]);
  }
  else{//All physics processes
    fJPsi = OpenFile(dataPathNames[0]);
    fpsi = OpenFile(dataPathNames[1]);
    fOC = OpenFile(dataPathNames[2]);
    fAMDY = OpenFile(dataPathNames[3]);
    fReal = OpenFile(dataPathNames[4]);
    
    TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fReal};
    TString type[nFiles] = {"Jpsi", "Psi", "OC", "AMDY", whichRD};
    
    if (fromHist){//From Hist
      for (Int_t i=0; i<nFiles; i++) {
	hist[i] = (TH1D*)Files[i]->Get(Form("h_%s", type[i].Data()));
	Setup(hist[i]);
      }
    }
    else{//From MC fitting
      for (Int_t i=0; i<nFiles; i++) {
	hist[i] = new TH1D(Form("h_%s", type[i].Data()),
			   Form("h_%s", type[i].Data()), nBins, Mmin, Mmax);
	Setup(hist[i]);
      }
      TTree *tJPsi = (TTree*) fJPsi->Get("pT_Weighted");
      TTree *tpsi = (TTree*) fpsi->Get("pT_Weighted");
      TTree *tOC = (TTree*) fOC->Get("pT_Weighted");
      TTree *tAMDY = (TTree*) fAMDY->Get("pT_Weighted");
      TTree *tReal = (TTree*) fReal->Get("pT_Weighted");

      cResults->cd();
      tJPsi->Draw("Mmumu>>h_JPsi", MassCut, "0");
      tpsi->Draw("Mmumu>>h_psi", MassCut, "0");
      tOC->Draw("Mmumu>>h_OC", MassCut, "0");
      tAMDY->Draw("Mmumu>>h_AMDY", MassCut, "0");
      tReal->Draw("Mmumu>>h_Real", MassCut, "0");
    }

    //Clone MC histograms to draw on added up graph
    hist[nFiles-1]->Draw();
    histResult[nFiles-1] = new TH1D("hResult_Sum", "hResult_Sum",nBins,Mmin,Mmax);
    for (Int_t i=0; i<nFiles-1; i++) {
      hist[i]->Sumw2();
      hist[i]->Scale(1/(hist[i]->Integral()) );
      hist[i]->SetLineColor(icolor[i]);

      histResult[i] = (TH1D*)hist[i]->Clone();
    }

    //Do fit
    fitFunc = new TF1("fitFunc", FitMCs, Mmin, Mmax, nFiles-1);
    hist[nFiles-1]->Fit("fitFunc", "RWL0", "", Mmin, Mmax);

    for (Int_t i=0; i<nFiles-1; i++) {
      histResult[i]->Scale(fitFunc->GetParameter(i) );
      histResult[i]->Draw("same");

      histResult[nFiles-1]->Add(histResult[i]);
    }

    rp =new TRatioPlot(histResult[nFiles-1], hist[nFiles-1]);
  }//All physics processes  

  Double_t chi =fitFunc->GetChisquare();
  Int_t ndf =fitFunc->GetNDF();
  cout << "\n\nChi2:  " << chi << "  / ndf:  " << ndf << endl;
  cout << "reduced chi2:  " << chi/(1.0*ndf) << endl;
  
  //Compare Fit with real data
  TCanvas* cRatio = new TCanvas();
  cRatio->SetTicks(0, 1);
  rp->Draw();
  cRatio->Update();
  cRatio->cd(1);

  //Make parameters tgraph
  Double_t *pars =fitFunc->GetParameters();
  const Double_t *e_pars =fitFunc->GetParErrors();
  Double_t xvals[nFiles-1], ex[nFiles-1] = {0.0};
  for (Int_t i=0; i<nFiles-1; i++) { xvals[i] = i+1; }
  TGraphErrors* gPar;
  if (HMDYonly) gPar = new TGraphErrors(1, xvals, pars, ex, e_pars);
  else gPar = new TGraphErrors(nFiles-1, xvals, pars, ex, e_pars);
  Setup(gPar);
  TCanvas* cPars = new TCanvas();
  gPar->Draw("AP"); gPar->SetTitle("Components Coefficents");

  //Write output file
  TString outName;
  if (HMDYonly){
    outName=Form("fitComponents_%s_%s_%0.2f_%0.2f_histFit%i_HMonly.root",
		 whichMC.Data(), whichRD.Data(), Mmin, Mmax, fromHist);
  }
  else{
    outName=Form("fitComponents_%s_%s_%0.2f_%0.2f_histFit%i.root",
		 whichMC.Data(), whichRD.Data(), Mmin, Mmax, fromHist);
  }
  if (toWrite){
    TFile* fOut = new TFile("Data/FitMassComponents/"+outName, "RECREATE");
    gPar->Write("gPar");

    if (HMDYonly){
      hist[0]->Write();
    }
    else{
      for (Int_t i=0; i<nFiles-1; i++) {
	hist[i]->Write();
      }
    }

    fOut->Close();
  }

  //Final Output
  cout << "\nSettings !!!!" << endl;
  cout << "Mass Range is		:  " << MassCut << endl;
  cout << "Which MC considered:          " << whichMC << endl;
  cout << "Which reald data considered:  " << whichRD << endl;
  cout << "Minimum mass considered:      " << Mmin << endl;
  cout << "Maximum mass considered:      " << Mmax << endl;
  cout << "Number bins in mass hist:     " << nBins << endl;
  cout << "Mass distributions from premade histogram:   " << fromHist << endl;
  cout << "HM Drell-Yan only:            " << HMDYonly << endl;
  if (toWrite) { cout << outName << "   was written" << endl; }
  else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
