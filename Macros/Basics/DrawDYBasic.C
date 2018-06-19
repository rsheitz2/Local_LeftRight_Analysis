void MakePretty(TH1D *h1){
  h1->GetXaxis()->SetLabelSize(0.06);
  h1->GetXaxis()->SetNdivisions(9);

  h1->GetYaxis()->SetLabelSize(0.06);
  h1->Draw();
}

void DrawDYBasic(){

  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/MC_Data/YuShiangMC_W11_run1.1/";
  TFile* f_AMDY =
    TFile::Open(path+"AMDY/Yu_Wall_full_main_AMDY_20bins.root");
  TTree *t_AMDY = (TTree*)f_AMDY->Get("pT_Weighted");
  
  TFile* f_JPsi =
    TFile::Open(path+"JPsi/Yu_Wall_full_main_JPsi_20bins.root");
  TTree *t_JPsi = (TTree*)f_JPsi->Get("pT_Weighted");
    
  TFile* f_psi =
    TFile::Open(path+"psi/Yu_Wall_full_main_psi_20bins.root");
  TTree *t_psi = (TTree*)f_psi->Get("pT_Weighted");
  
  TFile* f_OC =
    TFile::Open(path+"OC/Yu_Wall_full_main_OC_20bins.root");
  TTree *t_OC = (TTree*)f_OC->Get("pT_Weighted");
  
  //Basic file checks
  if (!(f_AMDY) || !(f_JPsi) || !(f_psi) || !(f_OC) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nProcess = 4;
  TString processName[nProcess] = {"JPsi", "AMDY", "psi", "OC"};//axes drawn on [0]
  Int_t iColor[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    if (processName[i] == "AMDY" ) iColor[i] = 1;
    else if (processName[i] == "psi" ) iColor[i] = 2;
    else if (processName[i] == "OC" ) iColor[i] = 3;
    else if (processName[i] == "JPsi" ) iColor[i] = 4;
  }

  TH1D *h_x_target[nProcess], *h_x_beam[nProcess], *h_xF[nProcess];
  TH1D *h_pT[nProcess], *h_M[nProcess];
  for (Int_t i=0; i<nProcess; i++) {
    h_x_target[i] = new TH1D(Form("h_x_target_%s", processName[i].Data() ),
			     Form("h_x_target_%s", processName[i].Data() ),
			     150, 0, 0.5);

    h_x_beam[i] = new TH1D(Form("h_x_beam_%s", processName[i].Data() ),
			   Form("h_x_beam_%s", processName[i].Data() ),
			   150, 0, 1);

    h_xF[i] = new TH1D(Form("h_xF_%s", processName[i].Data() ),
		       Form("h_xF_%s", processName[i].Data() ),
		       150, -0.4, 1);

    h_pT[i] = new TH1D(Form("h_pT_%s", processName[i].Data() ),
		       Form("h_pT_%s", processName[i].Data() ),
		       150, 0, 6);

    h_M[i] = new TH1D(Form("h_M_%s", processName[i].Data() ),
		      Form("h_M_%s", processName[i].Data() ),
		      150, 0, 12);
  }

  t_AMDY->Draw("x_target>>h_x_target_AMDY", "", "0");
  t_JPsi->Draw("x_target>>h_x_target_JPsi", "", "0");
  t_psi->Draw("x_target>>h_x_target_psi", "", "0");
  t_OC->Draw("x_target>>h_x_target_OC", "", "0");

  t_AMDY->Draw("x_beam>>h_x_beam_AMDY", "", "0");
  t_JPsi->Draw("x_beam>>h_x_beam_JPsi", "", "0");
  t_psi->Draw("x_beam>>h_x_beam_psi", "", "0");
  t_OC->Draw("x_beam>>h_x_beam_OC", "", "0");

  t_AMDY->Draw("x_feynman>>h_xF_AMDY", "", "0");
  t_JPsi->Draw("x_feynman>>h_xF_JPsi", "", "0");
  t_psi->Draw("x_feynman>>h_xF_psi", "", "0");
  t_OC->Draw("x_feynman>>h_xF_OC", "", "0");

  t_AMDY->Draw("q_transverse>>h_pT_AMDY", "", "0");
  t_JPsi->Draw("q_transverse>>h_pT_JPsi", "", "0");
  t_psi->Draw("q_transverse>>h_pT_psi", "", "0");
  t_OC->Draw("q_transverse>>h_pT_OC", "", "0");

  t_AMDY->Draw("Mmumu>>h_M_AMDY", "", "0");
  t_JPsi->Draw("Mmumu>>h_M_JPsi", "", "0");
  t_psi->Draw("Mmumu>>h_M_psi", "", "0");
  t_OC->Draw("Mmumu>>h_M_OC", "", "0");
  
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3,2);
  for (Int_t i=0; i<nProcess; i++) {
    h_x_target[i]->SetLineColor(iColor[i]);
    h_x_beam[i]->SetLineColor(iColor[i]);
    h_xF[i]->SetLineColor(iColor[i]);
    h_pT[i]->SetLineColor(iColor[i]);
    h_M[i]->SetLineColor(iColor[i]);

    if (i==0){
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_x_target[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_x_beam[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_xF[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_pT[i]);

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      MakePretty(h_M[i]);
    }
    else{
      Int_t ipad=1;
      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      h_x_target[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      h_x_beam[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      h_xF[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      h_pT[i]->Draw("Same");

      c1->cd(ipad); ipad++;
      gPad->SetLogy();
      h_M[i]->Draw("Same");
    }
  }

}
