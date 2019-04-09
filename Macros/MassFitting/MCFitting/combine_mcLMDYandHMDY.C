void SetUp(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);
  
  h->GetXaxis()->SetNdivisions(508);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);
}

void combine_mcLMDYandHMDY(){
  //Setup_____
  const Int_t nBins =5;
  TString physBinned ="xN";
  //Setup_____

  cout << "Macro combines low and high mass drell yan by fitting with and ";
  cout << "exponential near the boundary and extrapolating in between" << endl;

  TString lrPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/";
  /*TString nameLMDY ="Data/leftRight_byTarget_Charles_LMDY1.00_8.50_5bins25_43_150hbin.root";
    TString nameHMDY ="Data/leftRight_byTarget_Charles_HMDY1.00_8.50_5bins25_43_150hbin.root";//*/
  TString nameLMDY ="Data/leftRight_byTarget_Charles_LMDY1.00_8.50_4bins29_34_150hbin.root";
  TString nameHMDY ="Data/leftRight_byTarget_Charles_HMDY1.00_8.50_4bins29_34_150hbin.root";//*/
  TFile *fHMDY = TFile::Open(lrPath+nameHMDY);
  TFile *fLMDY = TFile::Open(lrPath+nameLMDY);

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas(); c1->Divide(2, 3);
  TCanvas* c2 = new TCanvas(); c2->Divide(2, 3);
  TH1D *hHMDY[nBins], *hLMDY[nBins];
  TH1D *hSum[nBins];
  for (Int_t i=0; i<nBins; i++) {
    c1->cd(i+1);
    gPad->SetLogy();

    //HMDY
    hHMDY[i] =
      (TH1D*)fHMDY->Get(Form("MuMu_left_upstream_up_%s%i", physBinned.Data(), i) );
    hHMDY[i]->Sumw2();
    hHMDY[i]->Scale(1.0/(hHMDY[i]->Integral() ) );
    hHMDY[i]->SetLineColor(kGreen);
    if (physBinned == "pT") hHMDY[i]->Fit("expo", "LQ0", "", 4.0, 4.3);//pT
    else hHMDY[i]->Fit("expo", "LQ0", "", 4.0, 4.2);//xN, xPi, xF
        
    TF1 *fexp = (TF1*)hHMDY[i]->GetFunction("expo");
    Double_t HMDY_val = fexp->Eval(3.5);

    //LMDY
    hLMDY[i] =
      (TH1D*)fLMDY->Get(Form("MuMu_left_upstream_up_%s%i", physBinned.Data(), i) );
    hLMDY[i]->Sumw2();
    hLMDY[i]->Scale(1.0/(hLMDY[i]->Integral() ) );
    if (physBinned == "pT") hLMDY[i]->Fit("expo", "LQ0", "" , 3.0, 3.2);//pT
    else hLMDY[i]->Fit("expo", "LQ0", "" , 2.9, 3.1);//xN, xPi, xF
        
    fexp = (TF1*)hLMDY[i]->GetFunction("expo");
    Double_t LMDY_val = fexp->Eval(3.5);

    //Combine and Draw
    hSum[i] = (TH1D*) hHMDY[i]->Clone(); SetUp(hSum[i]);
    hSum[i]->Scale(LMDY_val/HMDY_val);
    hSum[i]->Add(hLMDY[i]);
    hSum[i]->Draw(); hSum[i]->SetLineColor(kRed);
    
    hHMDY[i]->Scale(LMDY_val/HMDY_val);
    hHMDY[i]->Draw("sames");
    hLMDY[i]->Draw("sames");

    c2->cd(i+1);
    gPad->SetLogy();
    hSum[i]->Draw();
  }
    
}
