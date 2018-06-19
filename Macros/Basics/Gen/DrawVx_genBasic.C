void MakePretty(TH1D *h1){
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetNdivisions(9);

  h1->GetYaxis()->SetLabelSize(0.06);
  h1->Draw();
}

void DrawVx_genBasic(){
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
  TString processName[nProcess] = {"psi", "OC", "AMDY", "JPsi"};//axes drawn on [0]
  
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

  
  TH1D *h_vx_z[nProcess], *h_vx_x[nProcess], *h_vx_y[nProcess];
  TH1D *h_vx_zVar[nProcess], *h_vx_xVar[nProcess], *h_vx_yVar[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    h_vx_z[i] = new TH1D(Form("h_vx_z_%s", processName[i].Data() ),
			 Form("h_vx_z_%s", processName[i].Data() ),
			 150, -310, -150);

    h_vx_x[i] = new TH1D(Form("h_vx_x_%s", processName[i].Data() ),
			 Form("h_vx_x_%s", processName[i].Data() ),
			 150, -3.5, 3.5);

    h_vx_y[i] = new TH1D(Form("h_vx_y_%s",processName[i].Data()),
			 Form("h_vx_y_%s",processName[i].Data()),
			 150, -3.5, 3.5);
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3,2);
  Double_t vx_z, vx_x, vx_y;
  for (Int_t i=0; i<nProcess; i++) {
    //Tree stuff
    tree[i]->SetBranchAddress("vx_z", &vx_z);
    tree[i]->SetBranchAddress("vx_x", &vx_x);
    tree[i]->SetBranchAddress("vx_y", &vx_y);
    
    for (Int_t ev=0; ev<tree[i]->GetEntries(); ev++) {
    //for (Int_t ev=0; ev<1000; ev++) {
      tree[i]->GetEntry(ev);
      
      h_vx_z[i]->Fill(vx_z);
      h_vx_x[i]->Fill(vx_x);
      h_vx_y[i]->Fill(vx_y);
    }//tree loop


    //Draw stuff
    h_vx_z[i]->SetLineColor(iColor[i]);
    h_vx_x[i]->SetLineColor(iColor[i]);
    h_vx_y[i]->SetLineColor(iColor[i]);

    if (i==0){
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_vx_z[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_vx_x[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_vx_y[i]);
    }
    else{
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      h_vx_z[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_vx_x[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      h_vx_y[i]->Draw("Same");
    }
  }

}
