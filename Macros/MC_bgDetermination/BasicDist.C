const Int_t nFiles=5;
const Int_t nMrange=4;

void BasicDist(){
  //TString physType = "Mmumu"; Double_t xMin=2., xMax=8.5; Int_t nBins=200;
  TString physType = "x_beam"; Double_t xMin=0.1, xMax=0.9; Int_t nBins=100;
  //TString physType = "x_target"; Double_t xMin=0.04, xMax=0.5; Int_t nBins=200;

  Int_t iMcut = 2;
  TString MassCut[nMrange] = {"Mmumu>2.5&&Mmumu<8.5", "Mmumu>2.5&&Mmumu<4.3",
			      "Mmumu>4.3&&Mmumu<8.5", ""};
  //TString additionalCut = "&&x_beam>0.2&&x_beam<0.8";
  TString additionalCut = "";
  
  
  //{JPsi, psi, OC, AMDY}//Fit results from FitMassComponents.C
  Double_t Counts[nMrange][nFiles-1] = {
    {1.53946e+06, 2.98861e+04, 2.13914e+05, 1.62672e+05},//2.5-8.5GeV
    {2.73913e+06, 4.22147e+04, 3.33131e+05, 2.96746e+05},//2.5-4.3GeV
    {2.58264e+02, 7.16396e+02, 8.86266e+01, 3.38407e+04},//4.3-8.5GeV
    {1.53946e+06, 2.98861e+04, 2.13914e+05, 1.62672e+05},//""
  };

  Int_t icolor[nFiles+1] = {6,3,9,4,2,1};//{JPsi, psi, OC, AMDY, Sum, RealData}
  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/";
  TString pathMC = path+"MC_Data/YuShiangMC/";
  TFile *fJPsi = TFile::Open(pathMC+"JPsi/Yu_Wall_full_main_JPsi_20bins.root");
  TTree *tJPsi = (TTree*) fJPsi->Get("pT_Weighted");
  TFile *fpsi = TFile::Open(pathMC+"psi/Yu_Wall_full_main_psi_20bins.root");
  TTree *tpsi = (TTree*) fpsi->Get("pT_Weighted");
  TFile *fOC = TFile::Open(pathMC+"OC/Yu_Wall_full_main_OC_20bins.root");
  TTree *tOC = (TTree*) fOC->Get("pT_Weighted");
  TFile *fAMDY = TFile::Open(pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root");
  TTree *tAMDY = (TTree*) fAMDY->Get("pT_Weighted");
  TFile *fReal = TFile::Open(path+"RealData/WAll_AllMass.root");
  TTree *tReal = (TTree*) fReal->Get("pT_Weighted");

  TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fReal};
  TH1D *hist[nFiles];
  TString type[nFiles] = {"JPsi", "psi", "OC", "AMDY", "Real"};
  for (Int_t i=0; i<nFiles; i++) {
    if (!Files[i]) {
      cout << "File:   " << i << "   does not exist " << endl;
      exit(EXIT_FAILURE); }

    hist[i] = new TH1D(Form("h_%s", type[i].Data()),
		       Form("h_%s", type[i].Data()), nBins, xMin, xMax);
  }

  
  TCanvas* c1 = new TCanvas();
  tJPsi->Draw(Form("%s>>h_JPsi",physType.Data() ), MassCut[iMcut]+additionalCut,
	      "0");
  tpsi->Draw(Form("%s>>h_psi"  ,physType.Data() ), MassCut[iMcut]+additionalCut,
	     "0");
  tOC->Draw(Form("%s>>h_OC"    ,physType.Data() ), MassCut[iMcut]+additionalCut,
	    "0");
  tAMDY->Draw(Form("%s>>h_AMDY",physType.Data() ), MassCut[iMcut]+additionalCut,
	      "0");
  tReal->Draw(Form("%s>>h_Real",physType.Data() ), MassCut[iMcut]+additionalCut,
	      "0");
  hist[nFiles-1]->Draw(); hist[nFiles-1]->SetLineColor(icolor[nFiles]);

  
  for (Int_t i=0; i<nFiles-1; i++) {
    hist[i]->Sumw2();
    hist[i]->Scale(Counts[iMcut][i]/(hist[i]->Integral()) );
    hist[i]->SetLineColor(icolor[i]);

    hist[i]->Draw("sames");
  }

  TH1D* histTotal = (TH1D*)hist[0]->Clone();
  for (Int_t i=1; i<nFiles-1; i++) histTotal->Add(hist[i]);
    
  histTotal->SetLineColor(icolor[nFiles-1]);
  histTotal->Draw("same");
  
  
  TCanvas* c2 = new TCanvas();
  auto rp = new TRatioPlot(hist[nFiles-1], histTotal);
  rp->GetLowYaxis()->SetNdivisions(505);
  rp->Draw();
  c2->Update();


  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Physics type considered	:  " << physType << endl;
  cout << "Mass Range is		:  " << MassCut[iMcut] << endl;
  cout << "Additional cuts		:  " << additionalCut << endl;
  cout << " " << endl;
}
