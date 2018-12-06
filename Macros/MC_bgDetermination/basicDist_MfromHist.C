#include "Include/helperFunctions.h"

void basicDist_MfromHist(TString start=""){
  //Setup__________
  TString physType = "PhiS";
  Double_t xMin=-TMath::Pi(), xMax=TMath::Pi(); Int_t nBins=100;
  TString additionalCut = "";
  
  TString whichMC ="Charles"; //"Yu", "Charles"
  const Double_t Mmin=2.0, Mmax=8.5; //Mass cut
  TString whichRD ="slot1"; //"t3", "slot1"

  Bool_t toWrite =false;
  //Setup__________

  const Int_t nFiles=6;
  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);
    
  if (start=="") {
    cout << "Script plots the components that add up to a basic distribution";
    cout << "\nThis is done fitting MC distributions as a function of Mass to";
    cout << " to get the contribution from each physical process" << endl;
    cout << "\nTo be used specifically with Charles MC";
    cout << "  with combined LMDY & HMDY\n";
    cout << "\nCurrent Settings:" << endl;
    cout << "Basic distribution plotted:   " << physType << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Additional cuts considered:   " << additionalCut << endl;
    cout << "Which MC considered:          " << whichMC << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'BasicDist(1)\'" << endl;
    cout << "\nWarning!!! not setup for highM DY" << endl;
    exit(EXIT_FAILURE);
  }

      if (Mmin > 2.1) {
	cout << "\nWarning!!! not setup for highM DY" << endl;
	exit(EXIT_FAILURE);
      }
  
  //{JPsi, psi, OC, LMDY, HMDY}//Fit results from FitMassComponents_fromHist.C
  TString fitFileName =Form("fitComponents_%s_%s_%0.2f_%0.2f.root",
			    whichMC.Data(), whichRD.Data(), Mmin, Mmax);
  TFile *fFits = TFile::Open("Data/FitMassComponents/"+fitFileName);
  if (!fFits){
    cout << "fit file does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  TGraphErrors *gPars = (TGraphErrors*)fFits->Get("gPar");
  Double_t *Counts = gPars->GetY();

  //Aesthetic setup
  Int_t icolor[nFiles+1] ={6,3,9,4,2,1,1};//{JPsi, psi, OC, LMDY, HMDY, Sum, RD}
  //Data file names and paths
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData/LowM_AMDY/";
  TString pathMC ="/Volumes/Seagate/DrellYan/Charles_Official/";
  TString dataPathNames[nFiles];
  if (whichMC=="Charles"){
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
    dataPathNames[5] = pathRD+"WAll_LowM_AMDY.root";    
  }
  else if (whichRD=="slot1"){
    dataPathNames[5] = pathRD+"slot1WAll_LowM_AMDY.root";
  }
  else {
    cout << "RD specifed does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  TFile *fJPsi = TFile::Open(dataPathNames[0]);
  TFile *fpsi = TFile::Open(dataPathNames[1]);
  TFile *fOC = TFile::Open(dataPathNames[2]);
  TFile *fLMDY = TFile::Open(dataPathNames[3]);
  TFile *fHMDY = TFile::Open(dataPathNames[4]);
  TFile *fReal = TFile::Open(dataPathNames[5]);
  TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fLMDY, fHMDY, fReal};
  TH1D *hist[nFiles];
  
  TString type[nFiles] = {"JPsi", "psi", "OC", "LMDY","HMDY", "Real"};
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
  TTree *tLMDY = (TTree*) fLMDY->Get("pT_Weighted");
  TTree *tHMDY = (TTree*) fHMDY->Get("pT_Weighted");
  TTree *tReal = (TTree*) fReal->Get("pT_Weighted");

  //Draw distributions from tree
  TCanvas *cTmp = new TCanvas("cTmp", "cTmp");
  tJPsi->Draw(Form("%s>>h_JPsi",physType.Data() ), MassCut+additionalCut,
	      "0");
  tpsi->Draw(Form("%s>>h_psi"  ,physType.Data() ), MassCut+additionalCut,
	     "0");
  tOC->Draw(Form("%s>>h_OC"    ,physType.Data() ), MassCut+additionalCut,
	    "0");
  tLMDY->Draw(Form("%s>>h_LMDY",physType.Data() ), MassCut+additionalCut,
	      "0");
  tHMDY->Draw(Form("%s>>h_HMDY",physType.Data() ), MassCut+additionalCut,
	      "0");
  tReal->Draw(Form("%s>>h_Real",physType.Data() ), MassCut+additionalCut,
	      "0");
  
  //Scaling all but     Not LMDY or HMDY
  TCanvas* c1 = new TCanvas("c1", "c1");
  hist[nFiles-1]->Draw(); hist[nFiles-1]->SetLineColor(icolor[nFiles-1]);
  hist[nFiles-1]->GetYaxis()->SetRangeUser(0, hist[nFiles-1]->GetMaximum()*1.1);
  for (Int_t i=0; i<nFiles-3; i++) {
    cout << "Drawing " << type[i] << endl;
    hist[i]->Sumw2();
    hist[i]->Scale(Counts[i]/(hist[i]->Integral()) );
    hist[i]->SetLineColor(icolor[i]);

    hist[i]->Draw("sames");
  }
  
  //Add LMDY/HMDY
  TFile *f_M_AMDY = TFile::Open("Data/MassDist/M_Charles_AMDY.root");
  TGraph *g_LM_HM_scale = (TGraph*) f_M_AMDY->Get("LMDY_HMDY_scale");
  Double_t scale_LM_HM = g_LM_HM_scale->GetY()[0];
  TH1D *h_AMDY = (TH1D*)hist[3]->Clone();//LMDY
  h_AMDY->Scale(1.0/h_AMDY->Integral() );
  hist[4]->Scale(scale_LM_HM/(hist[4]->Integral()) );//HMDY
  h_AMDY->Add(hist[4]);//HMDY
  h_AMDY->Scale(Counts[3]/(h_AMDY->Integral()) );
  h_AMDY->SetLineColor(icolor[3]); h_AMDY->Draw("sames");//cleanup
      
  TH1D* histTotal = (TH1D*)hist[0]->Clone();
  histTotal->Add(hist[1]); histTotal->Add(hist[2]); histTotal->Add(h_AMDY);
  histTotal->SetLineColor(kRed);
  histTotal->Draw("sames");

  TCanvas* c2 = new TCanvas("c2", "c2");
  TRatioPlot *rp = new TRatioPlot(hist[nFiles-1], histTotal);
  rp->Draw(); Setup(rp); 
  rp->GetLowerRefGraph()->SetMinimum(0.9);
  rp->GetLowerRefGraph()->SetMaximum(1.1);
  c2->Update();

  //Make Ratio
  TH1D *hRatio =(TH1D*)histTotal->Clone();
  hRatio->Divide(hist[nFiles-1]); hRatio->SetLineColor(kBlack);
  TCanvas* cRatio = new TCanvas();
  hRatio->Draw();
  DrawLine(hRatio, 1.0);
  
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
