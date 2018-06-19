void MakePretty(TH1D *h1, Bool_t normalized=false){
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetNdivisions(9);

  h1->GetYaxis()->SetLabelSize(0.06);
  
  if (normalized) {
    h1->DrawNormalized();
  }
  else h1->Draw();
}

void DrawMu_genBasic(){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/Gen_MC_Data/Yu_BW/";
  TFile* f_AMDY =
    TFile::Open(path+"AMDY.root");
    
  TFile* f_JPsi =
    TFile::Open(path+"JPsi.root");
    
  TFile* f_psi =
    TFile::Open(path+"psi_minusW10_SP2.root");
  
  TFile* f_OC =
    TFile::Open(path+"OC.root");

  //Basic file checks
  if (!(f_AMDY) || !(f_JPsi) || !(f_psi) || !(f_OC) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  
  const Int_t nProcess = 4;
  TString processName[nProcess] = {"AMDY", "psi", "OC", "JPsi"};//axes drawn on [0]
  
  Int_t iColor[nProcess];
  TTree *tree[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    if (processName[i] == "AMDY" ) {
      iColor[i] = 1;
      tree[i] = (TTree*)f_AMDY->Get("Event");
    }
    else if (processName[i] == "psi" ) {
      iColor[i] = 2;
      tree[i] = (TTree*)f_psi->Get("Event");
    }
    else if (processName[i] == "OC" ) {
      iColor[i] = 3;
      tree[i] = (TTree*)f_OC->Get("Event");
    }
    else if (processName[i] == "JPsi" ) {
      iColor[i] = 4;
      tree[i] = (TTree*)f_JPsi->Get("Event");
    }
  }

  
  TH1D *h_theta_muP[nProcess], *h_qP_muP[nProcess], *h_phi_muP[nProcess];
  TH1D *h_theta_muM[nProcess], *h_qP_muM[nProcess], *h_phi_muM[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    h_theta_muP[i] = new TH1D(Form("h_theta_muP_%s", processName[i].Data() ),
				Form("h_theta_muP_%s", processName[i].Data() ),
				150, 0, 0.2);

    h_qP_muP[i] = new TH1D(Form("h_qP_muP_%s", processName[i].Data() ),
			     Form("h_qP_muP_%s", processName[i].Data() ),
			     150, 0, 200);

    h_phi_muP[i] = new TH1D(Form("h_phi_muP_%s",processName[i].Data()),
			      Form("h_phi_muP_%s",processName[i].Data()),
			      150, -TMath::Pi(), TMath::Pi() );

    h_theta_muM[i] = new TH1D(Form("h_theta_muM_%s", processName[i].Data() ),
				Form("h_theta_muM_%s", processName[i].Data() ),
				150, 0, 0.2);
			      
    h_qP_muM[i] = new TH1D(Form("h_qP_muM_%s", processName[i].Data() ),
			     Form("h_qP_muM_%s", processName[i].Data() ),
			     150, 0, 200);
    
    h_phi_muM[i] = new TH1D(Form("h_phi_muM_%s",processName[i].Data()),
			      Form("h_phi_muM_%s",processName[i].Data()),
			      150, -TMath::Pi(), TMath::Pi() );
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3,2);
  Double_t theta_muP, qP_muP, phi_muP;
  Double_t theta_muM, qP_muM, phi_muM;
  for (Int_t i=0; i<nProcess; i++) {
    //Tree stuff
    tree[i]->SetBranchAddress("theta_muP", &theta_muP);
    tree[i]->SetBranchAddress("qP_muP", &qP_muP);
    tree[i]->SetBranchAddress("phi_muP", &phi_muP);
    tree[i]->SetBranchAddress("theta_muM", &theta_muM);
    tree[i]->SetBranchAddress("qP_muM", &qP_muM);
    tree[i]->SetBranchAddress("phi_muM", &phi_muM);

    for (Int_t ev=0; ev<tree[i]->GetEntries(); ev++) {
      //for (Int_t ev=0; ev<1000; ev++) { cout << "debugging" << endl;
      tree[i]->GetEntry(ev);
      
      h_theta_muP[i]->Fill(theta_muP);
      h_qP_muP[i]->Fill(qP_muP);
      h_phi_muP[i]->Fill(phi_muP);
      h_theta_muM[i]->Fill(theta_muM);
      h_qP_muM[i]->Fill(qP_muM);
      h_phi_muM[i]->Fill(phi_muM);
    }//tree loop


    //Draw stuff
    h_theta_muP[i]->SetLineColor(iColor[i]);
    h_qP_muP[i]->SetLineColor(iColor[i]);
    h_phi_muP[i]->SetLineColor(iColor[i]);
    h_theta_muM[i]->SetLineColor(iColor[i]);
    h_qP_muM[i]->SetLineColor(iColor[i]);
    h_phi_muM[i]->SetLineColor(iColor[i]);

    if (i==0){
      Int_t ipad=1; Bool_t mormalized = true;
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_theta_muP[i]);
      
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_qP_muP[i]);
            
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_phi_muP[i], mormalized);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_theta_muM[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_qP_muM[i]);
      
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_phi_muM[i], mormalized);
    }
    else{
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      h_theta_muP[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_qP_muP[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_phi_muP[i]->DrawNormalized("Same");

      c1->cd(ipad); ipad++;
      h_theta_muM[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_qP_muM[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_phi_muM[i]->DrawNormalized("Same");
    }
  }

}
