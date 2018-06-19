void MakePretty(TH1D *h1, Bool_t normalized=false){
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetNdivisions(9);

  h1->GetYaxis()->SetLabelSize(0.06);
  //if (normalized) ? h1->DrawNormalized() : h1->Draw();
}

void DrawMuBasic(){
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

  
  TH1D *h_theta_traj1[nProcess], *h_qP_traj1[nProcess], *h_phi_traj1[nProcess];
  TH1D *h_theta_traj2[nProcess], *h_qP_traj2[nProcess], *h_phi_traj2[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    h_theta_traj1[i] = new TH1D(Form("h_theta_traj1_%s", processName[i].Data() ),
				Form("h_theta_traj1_%s", processName[i].Data() ),
				150, 0, 0.2);

    h_qP_traj1[i] = new TH1D(Form("h_qP_traj1_%s", processName[i].Data() ),
			     Form("h_qP_traj1_%s", processName[i].Data() ),
			     150, 0, 200);

    h_phi_traj1[i] = new TH1D(Form("h_phi_traj1_%s",processName[i].Data()),
			      Form("h_phi_traj1_%s",processName[i].Data()),
			      150, -TMath::Pi(), TMath::Pi() );

    h_theta_traj2[i] = new TH1D(Form("h_theta_traj2_%s", processName[i].Data() ),
				Form("h_theta_traj2_%s", processName[i].Data() ),
				150, 0, 0.2);
			      
    h_qP_traj2[i] = new TH1D(Form("h_qP_traj2_%s", processName[i].Data() ),
			     Form("h_qP_traj2_%s", processName[i].Data() ),
			     150, 0, 200);
    
    h_phi_traj2[i] = new TH1D(Form("h_phi_traj2_%s",processName[i].Data()),
			      Form("h_phi_traj2_%s",processName[i].Data()),
			      150, -TMath::Pi(), TMath::Pi() );
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3,2);
  Double_t theta_traj1, qP_traj1, phi_traj1;
  Double_t theta_traj2, qP_traj2, phi_traj2;
  for (Int_t i=0; i<nProcess; i++) {
    //Tree stuff
    tree[i]->SetBranchAddress("theta_traj1", &theta_traj1);
    tree[i]->SetBranchAddress("qP_traj1", &qP_traj1);
    tree[i]->SetBranchAddress("phi_traj1", &phi_traj1);
    tree[i]->SetBranchAddress("theta_traj2", &theta_traj2);
    tree[i]->SetBranchAddress("qP_traj2", &qP_traj2);
    tree[i]->SetBranchAddress("phi_traj2", &phi_traj2);

    for (Int_t ev=0; ev<tree[i]->GetEntries(); ev++) {
      //for (Int_t ev=0; ev<1000; ev++) {
      tree[i]->GetEntry(ev);
      
      h_theta_traj1[i]->Fill(theta_traj1);
      h_qP_traj1[i]->Fill(qP_traj1);
      h_phi_traj1[i]->Fill(phi_traj1);
      h_theta_traj2[i]->Fill(theta_traj2);
      h_qP_traj2[i]->Fill(-1.0*qP_traj2);
      h_phi_traj2[i]->Fill(phi_traj2);
    }//tree loop


    //Draw stuff
    h_theta_traj1[i]->SetLineColor(iColor[i]);
    h_qP_traj1[i]->SetLineColor(iColor[i]);
    h_phi_traj1[i]->SetLineColor(iColor[i]);
    h_theta_traj2[i]->SetLineColor(iColor[i]);
    h_qP_traj2[i]->SetLineColor(iColor[i]);
    h_phi_traj2[i]->SetLineColor(iColor[i]);

    if (i==0){
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_theta_traj1[i]);
      h_theta_traj1[i]->Draw();

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_qP_traj1[i]);
      h_qP_traj1[i]->Draw();

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_phi_traj1[i]);
      h_phi_traj1[i]->DrawNormalized();

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_theta_traj2[i]);
      h_theta_traj2[i]->Draw();

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_qP_traj2[i]);
      h_qP_traj2[i]->Draw();
      
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_phi_traj2[i]);
      h_phi_traj2[i]->DrawNormalized();
    }
    else{
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      h_theta_traj1[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_qP_traj1[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_phi_traj1[i]->DrawNormalized("Same");

      c1->cd(ipad); ipad++;
      h_theta_traj2[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_qP_traj2[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_phi_traj2[i]->DrawNormalized("Same");
    }
  }

}
