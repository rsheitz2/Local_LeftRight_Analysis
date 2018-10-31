const Int_t nFiles = 5;
TH1D *hist[nFiles]; //{JPsi, psi, OC, AMDY, Real}


Double_t FitMCs(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );
    
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}


void fitCheck_mc(TString start=""){
  if (start==""){
    cout << "Script fits data using a sum of MC distributions and ";
    cout << "compares real data with the fit" << endl;
    cout << "Usage:" << endl;
    cout << "root \'MCFitComparison(1)\'" << endl;
    exit(EXIT_FAILURE);
  }
  const Int_t nBins=200;
  //const Double_t minX=2.5, maxX=8.5;
  const Double_t minX=4.3, maxX=8.5;
  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/";
  TString pathMC = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/MC_Data/YuShiangMC/";
  
  TFile *fJPsi = TFile::Open(pathMC+"JPsi/Yu_Wall_full_main_JPsi_20bins.root");
  TFile *fpsi = TFile::Open(pathMC+"psi/Yu_Wall_full_main_psi_20bins.root");
  TFile *fOC = TFile::Open(pathMC+"OC/Yu_Wall_full_main_OC_20bins.root");
  TFile *fAMDY = TFile::Open(pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root");
  TFile *fReal = TFile::Open(path+"RealData/AMDY/WAll_AllMass.root");

  TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fReal};
  for (Int_t i=0; i<nFiles; i++) {
    if (!Files[i]) {
      cout << "File:   " << i << "   does not exist " << endl;
      exit(EXIT_FAILURE); }
  }
  TTree *tJPsi = (TTree*) fJPsi->Get("pT_Weighted");
  TTree *tpsi = (TTree*) fpsi->Get("pT_Weighted");
  TTree *tOC = (TTree*) fOC->Get("pT_Weighted");
  TTree *tAMDY = (TTree*) fAMDY->Get("pT_Weighted");
  TTree *tReal = (TTree*) fReal->Get("pT_Weighted");

  TString type[nFiles] = {"JPsi", "psi", "OC", "AMDY", "Real"};
  for (Int_t i=0; i<nFiles; i++) {
    hist[i] = new TH1D(Form("h_%s", type[i].Data()),
		       Form("h_%s", type[i].Data()), nBins, minX, maxX);
  }

  TCanvas* c1 = new TCanvas();
  TString MassCut = "Mmumu<8.5&&Mmumu>2.5";
  tJPsi->Draw("Mmumu>>h_JPsi", MassCut, "0");
  tpsi->Draw("Mmumu>>h_psi", MassCut, "0");
  tOC->Draw("Mmumu>>h_OC", MassCut, "0");
  tAMDY->Draw("Mmumu>>h_AMDY", MassCut, "0");
  tReal->Draw("Mmumu>>h_Real", MassCut, "0");

  
  hist[nFiles-1]->Draw();
  TH1D *histResult[nFiles];//{JPsi, psi, OC, AMDY, Sum}
  Int_t icolor[nFiles+1] = {6,3,9,4,2,1};//{JPsi, psi, OC, AMDY, Sum, RealData}
  histResult[nFiles-1] = new TH1D("hResult_Sum", "hResult_Sum",nBins,minX,maxX);
  for (Int_t i=0; i<nFiles-1; i++) {
    hist[i]->Sumw2();
    hist[i]->Scale(1/(hist[i]->Integral()) );
    hist[i]->SetLineColor(icolor[i]);

    histResult[i] = (TH1D*)hist[i]->Clone();
  }

  TF1 *fitFunc = new TF1("fitFunc", FitMCs, minX, maxX, nFiles-1);
  hist[nFiles-1]->Fit("fitFunc", "RWL0", "", minX, maxX);

  for (Int_t i=0; i<nFiles-1; i++) {
    histResult[i]->Scale(fitFunc->GetParameter(i) );
    histResult[i]->Draw("same");

    histResult[nFiles-1]->Add(histResult[i]);
  }
  histResult[nFiles-1]->SetLineColor(icolor[nFiles-1]);
  histResult[nFiles-1]->Draw("same");

  
  TCanvas* c2 = new TCanvas();
  auto rp = new TRatioPlot(histResult[nFiles-1], hist[nFiles-1]);
  c2->SetTicks(0, 1);
  rp->Draw();
  c2->Update();
  c2->cd(1);


  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Mass Range is: " << minX << "  -  " << maxX << endl;
  cout << " " << endl;
}
