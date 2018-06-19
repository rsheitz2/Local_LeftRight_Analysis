void MakePretty(TH1D *h1, Bool_t normalized=false){
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetNdivisions(9);

  h1->GetYaxis()->SetLabelSize(0.06);

  if (normalized) {
    h1->DrawNormalized();
  }
  else h1->Draw();
}

void DrawPhiS_byTarg(){
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
  TString processName[nProcess] = {"OC", "psi", "AMDY", "JPsi"};//axes drawn on [0]
  
  Int_t iColor[nProcess];
  Double_t iScale[nProcess];
  TTree *tree[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    if (processName[i] == "AMDY" ) {
      iColor[i] = 1; iScale[i] = 1;
      tree[i] = (TTree*)f_AMDY->Get("Event");
    }
    else if (processName[i] == "psi" ) {
      iColor[i] = 2; iScale[i] = 1;
      tree[i] = (TTree*)f_psi->Get("Event");
    }
    else if (processName[i] == "OC" ) {
      iColor[i] = 3; iScale[i] = 1.1;
      tree[i] = (TTree*)f_OC->Get("Event");
    }
    else if (processName[i] == "JPsi" ) {
      iColor[i] = 4; iScale[i] = 1.2;
      tree[i] = (TTree*)f_JPsi->Get("Event");
    }
  }

  
  TH1D *h_PhiS[nProcess], *h_PhiS_upS[nProcess], *h_PhiS_downS[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    h_PhiS[i] = new TH1D(Form("h_PhiS_%s", processName[i].Data() ),
			 Form("h_PhiS_%s", processName[i].Data() ),
			 70, -TMath::Pi(), TMath::Pi() );

    h_PhiS_upS[i] = new TH1D(Form("h_PhiS_upS_%s", processName[i].Data() ),
			     Form("h_PhiS_upS_%s", processName[i].Data() ),
			     70, -TMath::Pi(), TMath::Pi() );

    h_PhiS_downS[i] = new TH1D(Form("h_PhiS_downS_%s", processName[i].Data() ),
			       Form("h_PhiS_downS_%s", processName[i].Data() ),
			       70, -TMath::Pi(), TMath::Pi() );
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3);
  Double_t PhiS, vx_z;
  for (Int_t i=0; i<nProcess; i++) {

    //Tree stuff
    tree[i]->SetBranchAddress("PhiS", &PhiS);
    tree[i]->SetBranchAddress("vx_z", &vx_z);
    
    for (Int_t ev=0; ev<tree[i]->GetEntries(); ev++) {
      tree[i]->GetEntry(ev);

      h_PhiS[i]->Fill(PhiS);
      if (vx_z < -229.4) h_PhiS_upS[i]->Fill(PhiS);
      else h_PhiS_downS[i]->Fill(PhiS);

      if (vx_z < -320.0 || vx_z > -135.0) {
	cout << "Target vx_z problems: vx_z= " << vx_z << endl;
	exit(EXIT_FAILURE);
      }
    }//tree loop


    //Draw stuff
    h_PhiS[i]->SetLineColor(iColor[i]);
    h_PhiS_upS[i]->SetLineColor(iColor[i]);
    h_PhiS_downS[i]->SetLineColor(iColor[i]);
    h_PhiS[i]->Scale(iScale[i]);
    h_PhiS_upS[i]->Scale(iScale[i]);
    h_PhiS_downS[i]->Scale(iScale[i]);
    
    if (i==0){
      Bool_t normalized=true;
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      MakePretty(h_PhiS[i]);
      
      c1->cd(ipad); ipad++;
      MakePretty(h_PhiS_upS[i]);
      
      c1->cd(ipad); ipad++;
      MakePretty(h_PhiS_downS[i]);
    }
    else{
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      h_PhiS[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_PhiS_upS[i]->Draw("Same");      

      c1->cd(ipad); ipad++;
      h_PhiS_downS[i]->Draw("Same");
    }
  }

}
