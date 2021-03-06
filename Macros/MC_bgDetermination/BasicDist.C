#include "Include/helperFunctions.h"

void basicDist(TString start=""){
  //Setup__________
  TString physType = "PhiS";
  Double_t xMin=-TMath::Pi(), xMax=TMath::Pi(); Int_t nBins=100;
  TString additionalCut = "";
  
  TString whichMC ="Charles"; //"Yu", "Charles"
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
  TString whichRD ="slot1"; //"t3", "slot1"

  Bool_t fromHist =true;
  Bool_t HMDYonly =false;

  Bool_t toWrite =false;
  //Setup__________

  //  Other frequently used setup options
  /*TString physType = "Mmumu"; Double_t xMin=2., xMax=8.5; Int_t nBins=200;
  //TString physType = "x_beam"; Double_t xMin=0.1, xMax=0.9; Int_t nBins=100;
  //TString physType = "x_target"; Double_t xMin=0.04, xMax=0.5; Int_t nBins=200;
  TString additionalCut = "&&x_beam>0.2&&x_beam<0.8";//*/
  
  const Int_t nFiles=6;
  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);
    
  if (start=="") {
    cout << "Script plots the components that add up to a basic distribution";
    cout << "\nThis is done fitting MC distributions as a function of Mass to";
    cout << " to get the contribution from each physical process" << endl;
    cout << "\nCurrent Settings:" << endl;
    cout << "Basic distribution plotted:   " << physType << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Additional cuts considered:   " << additionalCut << endl;
    cout << "Which MC considered:          " << whichMC << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'BasicDist(1)\'" << endl;
    exit(EXIT_FAILURE);
  }
  
  //{JPsi, psi, OC, AMDY}//Fit results from FitMassComponents.C
  TString fitFileName =Form("fitComponents_%s_%s_%0.2f_%0.2f_histFit%i.root",
			    whichMC.Data(), whichRD.Data(), Mmin, Mmax,
			    fromHist);
  TFile *fFits = OpenFile("Data/FitMassComponents/"+fitFileName);
  TGraphErrors *gPars = (TGraphErrors*)fFits->Get("gPar");
  Double_t *Counts = gPars->GetY();

  //Aesthetic setup
  Int_t icolor[nFiles+1] = {6,3,9,4,2,1, 1};//{JPsi, psi, OC, AMDY,Sum,RealData}

  //Data file names and paths
  TString localData = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/";
  TString dataPathNames[nFiles];
  if (whichMC=="Yu"){
    TString pathMC =localData+"MC_Data/YuShiangMC/";
    dataPathNames[0] = pathMC+"JPsi/Yu_Wall_full_main_JPsi_20bins.root";
    dataPathNames[1] = pathMC+"psi/Yu_Wall_full_main_psi_20bins.root";
    dataPathNames[2] = pathMC+"OC/Yu_Wall_full_main_OC_20bins.root";
    dataPathNames[3] = pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root";
  }
  else if (whichMC=="Charles") {
    TString pathMC ="/Volumes/Seagate/DrellYan/Charles_Official/";
    dataPathNames[0] = pathMC+"Charles_W12_Jpsi.root";
    dataPathNames[1] = pathMC+"Charles_W12_Psi.root";
    dataPathNames[2] = pathMC+"Charles_W12_OC.root";
    dataPathNames[3] = pathMC+"Charles_W12_LMDY.root";
    dataPathNames[4] = pathMC+"Charles_W12_HMDY.root";
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
  
  TFile *fJPsi = TFile::Open(dataPathNames[0]);
  TFile *fpsi = TFile::Open(dataPathNames[1]);
  TFile *fOC = TFile::Open(dataPathNames[2]);
  TFile *fAMDY = TFile::Open(dataPathNames[3]);
  TFile *fReal = TFile::Open(dataPathNames[4]);
  TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fReal};
  TH1D *hist[nFiles];
  
  TString type[nFiles] = {"JPsi", "psi", "OC", "AMDY", "Real"};
  for (Int_t i=0; i<nFiles; i++) {
    if (!Files[i]) {
      cout << "File:   " << i << "   does not exist " << endl;
      exit(EXIT_FAILURE); }

    hist[i] = new TH1D(Form("h_%s", type[i].Data()),
		       Form("h_%s", type[i].Data()), nBins, xMin, xMax);
    Setup(hist[i]);
  }
  TTree *tJPsi = (TTree*) fJPsi->Get("pT_Weighted");
  TTree *tpsi = (TTree*) fpsi->Get("pT_Weighted");
  TTree *tOC = (TTree*) fOC->Get("pT_Weighted");
  TTree *tAMDY = (TTree*) fAMDY->Get("pT_Weighted");
  TTree *tReal = (TTree*) fReal->Get("pT_Weighted");

  //Draw distributions from tree
  TCanvas* c1 = new TCanvas();
  tJPsi->Draw(Form("%s>>h_JPsi",physType.Data() ), MassCut+additionalCut,
	      "0");
  tpsi->Draw(Form("%s>>h_psi"  ,physType.Data() ), MassCut+additionalCut,
	     "0");
  tOC->Draw(Form("%s>>h_OC"    ,physType.Data() ), MassCut+additionalCut,
	    "0");
  tAMDY->Draw(Form("%s>>h_AMDY",physType.Data() ), MassCut+additionalCut,
	      "0");
  tReal->Draw(Form("%s>>h_Real",physType.Data() ), MassCut+additionalCut,
	      "0");
  hist[nFiles-1]->Draw(); hist[nFiles-1]->SetLineColor(icolor[nFiles]);

  
  for (Int_t i=0; i<nFiles-1; i++) {
    hist[i]->Sumw2();
    hist[i]->Scale(Counts[i]/(hist[i]->Integral()) );
    hist[i]->SetLineColor(icolor[i]);

    hist[i]->Draw("sames");
  }

  TH1D* histTotal = (TH1D*)hist[0]->Clone();
  for (Int_t i=1; i<nFiles-1; i++) histTotal->Add(hist[i]);
    
  histTotal->SetLineColor(icolor[nFiles-1]);
  histTotal->Draw("same");
  
  
  TCanvas* c2 = new TCanvas();
  TRatioPlot *rp = new TRatioPlot(hist[nFiles-1], histTotal);
  rp->Draw(); Setup(rp); 
  rp->GetLowerRefGraph()->SetMinimum(0.8);
  rp->GetLowerRefGraph()->SetMaximum(1.2);
  c2->Update();

  //Write output file
  TString outName =Form("basicDist_%s%s_%s_%s_%0.2f_%0.2f.root",
			physType.Data(), additionalCut.Data(),
			whichMC.Data(), whichRD.Data(), Mmin, Mmax);
  if (toWrite){
    TFile* fOut = new TFile("Data/BasicDist/"+outName, "RECREATE");
    histTotal->Write("SumOfFits");
    hist[nFiles-1]->Write();

    fOut->Close();
  }

  //Final Setup output
  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Physics type considered	:  " << physType << endl;
  cout << "Mass Range is		:  " << MassCut << endl;
  cout << "Additional cuts		:  " << additionalCut << endl;
  cout << "Which MC considered:          " << whichMC << endl;
  if (toWrite) { cout << outName << "   was written" << endl; }
  else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
