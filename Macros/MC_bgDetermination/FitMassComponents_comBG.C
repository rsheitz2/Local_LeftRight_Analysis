//const Int_t nFiles = 6;
const Int_t nFiles = 5;
TH1D *hist[nFiles]; //{JPsi, psi, OC, AMDY, CombBg, Real}


Double_t FitMCs(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );
  //Double_t nCombBg = hist[4]->GetBinContent(hist[4]->FindBin(x[0] ) );
    
  //return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY +par[4]*nCombBg;
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}


void FitMassComponents(){
  const Int_t nBins=200;
  const Double_t minX=2.5, maxX=8.5;
  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/";
  TString pathMC = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/MC_Data/YuShiangMC/";
  
  TFile *fJPsi = TFile::Open(pathMC+"JPsi/Yu_Wall_full_main_JPsi_20bins.root");
  TTree *tJPsi = (TTree*) fJPsi->Get("pT_Weighted");
  TFile *fpsi = TFile::Open(pathMC+"psi/Yu_Wall_full_main_psi_20bins.root");
  TTree *tpsi = (TTree*) fpsi->Get("pT_Weighted");
  TFile *fOC = TFile::Open(pathMC+"OC/Yu_Wall_full_main_OC_20bins.root");
  TTree *tOC = (TTree*) fOC->Get("pT_Weighted");
  TFile *fAMDY = TFile::Open(pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root");
  TTree *tAMDY = (TTree*) fAMDY->Get("pT_Weighted");
  //TFile *fCombBg = TFile::Open(path+"RealData/CombinatorialBg.root");
  TFile *fReal = TFile::Open(path+"RealData/WAll_AllMass.root");
  TTree *tReal = (TTree*) fReal->Get("pT_Weighted");


  //TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fCombBg, fReal};
  TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fReal};
  for (Int_t i=0; i<nFiles; i++) {
    if (!Files[i]) {
      cout << "File:   " << i << "   does not exist " << endl;
      exit(EXIT_FAILURE); }
  }

  //TString type[nFiles] = {"JPsi", "psi", "OC", "AMDY", "CombBg", "Real"};
  TString type[nFiles] = {"JPsi", "psi", "OC", "AMDY", "Real"};
  for (Int_t i=0; i<nFiles; i++) {
    //if (i==4) continue;//CombBg is already made
    
    hist[i] = new TH1D(Form("h_%s", type[i].Data()),
		       Form("h_%s", type[i].Data()), nBins, minX, maxX);
  }

  TCanvas* c1 = new TCanvas();
  TString MassCut = "Mmumu<8.5&&Mmumu>2.5";
  tJPsi->Draw("Mmumu>>h_JPsi", MassCut, "0");
  tpsi->Draw("Mmumu>>h_psi", MassCut, "0");
  tOC->Draw("Mmumu>>h_OC", MassCut, "0");
  tAMDY->Draw("Mmumu>>h_AMDY", MassCut, "0");
  //hist[4] = (TH1D*)fCombBg->Get("h_Bg");
  tReal->Draw("Mmumu>>h_Real", MassCut, "0");

  
  hist[nFiles-1]->Draw();
  TCanvas* c3 = new TCanvas();//cleanup
  TH1D *histResult[nFiles];//{JPsi, psi, OC, AMDY, CombBg, Sum}
  histResult[nFiles-1] = new TH1D("hResult_Sum", "hResult_Sum",nBins,minX,maxX);
  for (Int_t i=0; i<nFiles-1; i++) {
    hist[i]->Sumw2();
    //hist[i]->Scale(1/(hist[i]->GetEntries()) );
    //cout << hist[i]->
    //hist[i]->Scale(1/(hist[i]->ComputeIntegral()) );//cleanup
    hist[i]->Scale(1/(hist[i]->Integral()) );//cleanup
    hist[i]->SetLineColor(5+i);

    c3->cd();//cleanup
    if (i==0) hist[i]->Draw(); else hist[i]->Draw("same");//clean up

    histResult[i] = (TH1D*)hist[i]->Clone();
  }

  TF1 *fitFunc = new TF1("fitFunc", FitMCs, minX, maxX, nFiles-1);
  hist[nFiles-1]->Fit("fitFunc", "RWL0", "", minX, maxX);

  c1->cd();//cleanup
  for (Int_t i=0; i<nFiles-1; i++) {
    histResult[i]->Scale(fitFunc->GetParameter(i) );
    histResult[i]->Draw("same");

    histResult[nFiles-1]->Add(histResult[i]);
  }

  TCanvas* c2 = new TCanvas();
  //hist[nFiles-1]->GetXaxis()->SetTitle("test1");
  //hist[nFiles-1]->GetYaxis()->SetTitle("test2");
  //histTotal->GetYaxis()->SetTitle("test3");
  //histTotal->SetFillColor(kBlue);
  auto rp = new TRatioPlot(histResult[nFiles-1], hist[nFiles-1]);
  //auto rp = new TRatioPlot(histTotal, hist[nFiles-1]);
  c2->SetTicks(0, 1);
  rp->Draw();
  c2->Update();
  c2->cd(1);

}
