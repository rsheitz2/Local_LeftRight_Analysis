void AcceptanceDraw(){
  //TString phys = "x_beam"; Double_t xMin =0., xMax =1.;
  //TString phys = "x_target"; Double_t xMin =0., xMax =1.;
  TString phys = "x_feynman"; Double_t xMin =-0.2, xMax =0.8;
  //TString phys = "q_transverse"; Double_t xMin =0., xMax =5.;
  //TString phys = "Mmumu"; Double_t xMin =4., xMax =9.;

  Int_t nBins = 100;
  
  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/";
  TFile* f_phast = TFile::Open(path+"MC_Data/YuShiangMC/HMDY/\
Merged_W07run0_HMDY.root");
  TFile* f_gen = TFile::Open(path+"Gen_MC_Data/Yu_BW/HMDY_W07_run0.root");

  TTree *T_phast = (TTree*)f_phast->Get("pT_Weighted");
  TTree *T_gen = (TTree*)f_gen->Get("Event");
  
  if ( !(f_phast) || !(f_gen) ){//Basic file checks
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  TString zCut_up = "vx_z>-294.5&&vx_z<-239.3";
  TString zCut_down = "vx_z>-219.5&&vx_z<-164.3";
  TString radCut = "vx_x*vx_x+vx_y*vx_y<=1.9*1.9";
  TString qTCut = "q_transverse>0.4&&q_transverse<5.0";
  TString MassCut = "Mmumu>4.3&&Mmumu<8.5";
  TString leftCut = "PhiS_simple>0&&PhiS_simple<=TMath::Pi()";
  TString rightCut = "PhiS_simple>-1.0*TMath::Pi()&&PhiS_simple<=0.0";
  
  TString zCut = zCut_up + " || " + zCut_down;
  TString genCut = radCut +" && "+ qTCut +" && "+ MassCut +" && ("+ zCut +")";
  TString genUpSCut = radCut +" && "+ qTCut +" && "+ MassCut +" && "+ zCut_up;
  TString genDownSCut = radCut +" && "+ qTCut +" && "+ MassCut +" && "+ zCut_down;
  TString genLeftCut = genCut +" && "+ leftCut;
  TString genRightCut = genCut +" && "+ rightCut;

  
  TH1D* h_phast = new TH1D("h_phast", "h_phast", nBins, xMin, xMax);
  TH1D* h_phastUpS = new TH1D("h_phastUpS", "h_phastUpS", nBins, xMin, xMax);
  TH1D* h_phastDownS = new TH1D("h_phastDownS", "h_phastDownS", nBins, xMin, xMax);
  TH1D* h_phastLeft = new TH1D("h_phastLeft", "h_phastLeft", nBins, xMin, xMax);
  TH1D* h_phastRight = new TH1D("h_phastRight", "h_phastRight",nBins,xMin,xMax);
  TH1D* h_gen = new TH1D("h_gen", "h_gen", nBins, xMin, xMax);
  TH1D* h_genUpS = new TH1D("h_genUpS", "h_genUpS", nBins, xMin, xMax);
  TH1D* h_genDownS = new TH1D("h_genDownS", "h_genDownS", nBins, xMin, xMax);
  TH1D* h_genLeft = new TH1D("h_genLeft", "h_genLeft", nBins, xMin, xMax);
  TH1D* h_genRight = new TH1D("h_genRight", "h_genRight",nBins,xMin,xMax);
  
  T_phast->Draw(phys+">>h_phast", "", "0");
  T_phast->Draw(phys+">>h_phastUpS", zCut_up, "0");
  T_phast->Draw(phys+">>h_phastDownS", zCut_down, "0");
  T_phast->Draw(phys+">>h_phastLeft", leftCut, "0");
  T_phast->Draw(phys+">>h_phastRight", rightCut, "0");
  
  T_gen->Draw(phys+">>h_gen", genCut, "0");
  T_gen->Draw(phys+">>h_genUpS", genUpSCut, "0");
  T_gen->Draw(phys+">>h_genDownS", genDownSCut, "0");
  T_gen->Draw(phys+">>h_genLeft", genLeftCut, "0");
  T_gen->Draw(phys+">>h_genRight", genRightCut, "0");
  h_phast->Sumw2(); h_phastUpS->Sumw2(); h_phastDownS->Sumw2();
  h_phastLeft->Sumw2(); h_phastRight->Sumw2();
  h_gen->Sumw2(); h_genUpS->Sumw2(); h_genDownS->Sumw2(); h_genLeft->Sumw2();
  h_genRight->Sumw2();
    
  TH1D *h_acc = (TH1D*)h_phast->Clone();
  h_acc->Divide(h_gen);
  h_acc->GetYaxis()->SetRangeUser(0, 1);
  
  
  gStyle->SetOptStat(11111111);
  TCanvas* c1 = new TCanvas();
  c1->Divide(2);

  c1->cd(1);
  h_gen->Draw();
  h_gen->SetLineColor(kRed);
  h_phast->Draw("sames");
  gPad->SetLogy();
  
  c1->cd(2);
  h_acc->Draw();


  TH1D *h_accUpS = (TH1D*)h_phastUpS->Clone();
  h_accUpS->Divide(h_genUpS);
  h_accUpS->GetYaxis()->SetRangeUser(0, 1);
  TH1D *h_accDownS = (TH1D*)h_phastDownS->Clone();
  h_accDownS->Divide(h_genDownS);
  h_accDownS->GetYaxis()->SetRangeUser(0, 1);
  TCanvas* cUpDown = new TCanvas();
  cUpDown->Divide(2, 2);

  cUpDown->cd(1);
  h_genUpS->Draw();
  h_genUpS->SetLineColor(kRed);
  h_phastUpS->Draw("sames");
  gPad->SetLogy();
  
  cUpDown->cd(2);
  h_accUpS->Draw();
  
  cUpDown->cd(3);
  h_genDownS->Draw();
  h_genDownS->SetLineColor(kRed);
  h_phastDownS->Draw("sames");
  gPad->SetLogy();
  
  cUpDown->cd(4);
  h_accDownS->Draw();

  
  TCanvas* cR_udS = new TCanvas();
  TRatioPlot *R_updown = new TRatioPlot(h_accUpS, h_accDownS);
  R_updown->Draw();
  cR_udS->Update();


  TH1D *h_accLeft = (TH1D*)h_phastLeft->Clone();
  h_accLeft->Divide(h_genLeft);
  h_accLeft->GetYaxis()->SetRangeUser(0, 1);
  TH1D *h_accRight = (TH1D*)h_phastRight->Clone();
  h_accRight->Divide(h_genRight);
  h_accRight->GetYaxis()->SetRangeUser(0, 1);
  TCanvas* cLR = new TCanvas();
  cLR->Divide(2, 2);

  cLR->cd(1);
  h_genLeft->Draw();
  h_genLeft->SetLineColor(kRed);
  h_phastLeft->Draw("sames");
  gPad->SetLogy();
  
  cLR->cd(2);
  h_accLeft->Draw();
  
  cLR->cd(3);
  h_genRight->Draw();
  h_genRight->SetLineColor(kRed);
  h_phastRight->Draw("sames");
  gPad->SetLogy();
  
  cLR->cd(4);
  h_accRight->Draw();


  TCanvas* cR = new TCanvas();
  TRatioPlot *R_lr = new TRatioPlot(h_accLeft, h_accRight);
  R_lr->Draw();
  cR->Update();

  
  cout << " " << endl;
  cout << "_____Setup_____" << endl;
  cout << "Ploting    " << phys << endl;
  cout << "Range    min:  " << xMin << "   max: " << xMax << endl;
  cout << " " << endl;
  cout << " " << endl;
  cout << "Generated Cuts:" << endl;
  cout << genCut << endl;
  cout << " " << endl;  
  cout << "Generated Left Cuts:" << endl;
  cout << "genCut && " << leftCut << endl;
  cout << " " << endl;
  cout << "Generated Right Cuts:" << endl;
  cout << "genCut && " << rightCut << endl;
}
