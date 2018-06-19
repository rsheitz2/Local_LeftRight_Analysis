void MakePretty(TH1D *h1){
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetNdivisions(9);

  h1->GetYaxis()->SetLabelSize(0.06);
  h1->Draw();
}

void DrawAngleBasic(){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/MC_Data/YuShiangMC_W11_run1.1/";
  TFile* f_AMDY =
    TFile::Open(path+"AMDY/Yu_Wall_full_main_AMDY_20bins.root");
    
  TFile* f_JPsi =
    TFile::Open(path+"JPsi/Yu_Wall_full_main_JPsi_20bins.root");
    
  TFile* f_psi =
    TFile::Open(path+"psi/Yu_Wall_full_main_psi_20bins.root");
  
  TFile* f_OC =
    TFile::Open(path+"OC/Yu_Wall_full_main_OC_20bins.root");

  //Basic file checks
  if (!(f_AMDY) || !(f_JPsi) || !(f_psi) || !(f_OC) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  
  const Int_t nProcess = 4;
  TString processName[nProcess] = {"psi", "OC", "AMDY", "JPsi"};//axes drawn on [0]
  
  Int_t iColor[nProcess];
  TTree *tree[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    if (processName[i] == "AMDY" ) {
      iColor[i] = 1;
      tree[i] = (TTree*)f_AMDY->Get("pT_Weighted");
    }
    else if (processName[i] == "psi" ) {
      iColor[i] = 2;
      tree[i] = (TTree*)f_psi->Get("pT_Weighted");
    }
    else if (processName[i] == "OC" ) {
      iColor[i] = 3;
      tree[i] = (TTree*)f_OC->Get("pT_Weighted");
    }
    else if (processName[i] == "JPsi" ) {
      iColor[i] = 4;
      tree[i] = (TTree*)f_JPsi->Get("pT_Weighted");
    }
  }

  
  TH1D *h_PhiS[nProcess], *h_Phi_CS[nProcess], *h_CosTheta_CS[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    h_PhiS[i] = new TH1D(Form("h_PhiS_%s", processName[i].Data() ),
			 Form("h_PhiS_%s", processName[i].Data() ),
			 150, -TMath::Pi(), TMath::Pi() );

    h_Phi_CS[i] = new TH1D(Form("h_Phi_CS_%s", processName[i].Data() ),
			   Form("h_Phi_CS_%s", processName[i].Data() ),
			   150, -TMath::Pi(), TMath::Pi() );

    h_CosTheta_CS[i] = new TH1D(Form("h_CosTheta_CS_%s",processName[i].Data()),
				Form("h_CosTheta_CS_%s",processName[i].Data()),
				150, -1, 1);
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3);
  Double_t PhiS, Phi_CS, Theta_CS;
  for (Int_t i=0; i<nProcess; i++) {

    //Tree stuff
    tree[i]->SetBranchAddress("PhiS", &PhiS);
    tree[i]->SetBranchAddress("Phi_CS", &Phi_CS);
    tree[i]->SetBranchAddress("Theta_CS", &Theta_CS);

    for (Int_t ev=0; ev<tree[i]->GetEntries(); ev++) {
      tree[i]->GetEntry(ev);
      
      h_PhiS[i]->Fill(PhiS);
      h_Phi_CS[i]->Fill(Phi_CS);
      h_CosTheta_CS[i]->Fill(TMath::Cos(Theta_CS) );
    }//tree loop


    //Draw stuff
    h_PhiS[i]->SetLineColor(iColor[i]);
    h_Phi_CS[i]->SetLineColor(iColor[i]);
    h_CosTheta_CS[i]->SetLineColor(iColor[i]);

    if (i==0){
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      MakePretty(h_PhiS[i]);
      h_PhiS[i]->GetYaxis()->SetRangeUser(1, 8200);

      c1->cd(ipad); ipad++;
      MakePretty(h_Phi_CS[i]);
      h_Phi_CS[i]->GetYaxis()->SetRangeUser(1, 8500);

      c1->cd(ipad); ipad++;
      MakePretty(h_CosTheta_CS[i]);
    }
    else{
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      h_PhiS[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_Phi_CS[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_CosTheta_CS[i]->Draw("Same");
    }
  }

}
