void SetUp(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);
  
  h->GetXaxis()->SetNdivisions(508);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);
}

void writeFile_mcLMDYandHMDY(){
  //Setup_____
  const Int_t nBins =4;
  TString binRange="29_34";
  TString physBinned[4] ={"xN", "xPi", "xF", "pT"};
  TString targets[4] ={"upstream_up", "upstream_down",
		       "downstream_up", "downstream_down"};
  //Setup_____

  cout << "Macro combines low and high mass drell yan by fitting with and ";
  cout << "exponential near the boundary and extrapolating in between" << endl;
  cout << "Output is written to a file" << endl;
  
  TString lrPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/";
  TString nameLMDY =Form("Data/leftRight_byTarget_Charles_LMDY1.00_8.50_%ibins%s_150hbin.root", nBins, binRange.Data());
  TString nameHMDY =Form("Data/leftRight_byTarget_Charles_HMDY1.00_8.50_%ibins%s_150hbin.root", nBins, binRange.Data());
  TFile *fHMDY = TFile::Open(lrPath+nameHMDY);
  TFile *fLMDY = TFile::Open(lrPath+nameLMDY);
  
  TFile *fOutput = new TFile(Form("%sData/leftRight_byTarget_Charles_AMDY1.00_8.50_%ibins%s_150hbin.root", lrPath.Data(), nBins, binRange.Data()), "RECREATE");
  for (Int_t tr=0; tr<4; tr++) {
    for (Int_t phys=0; phys<4; phys++) {
      for (Int_t bi=0; bi<nBins; bi++) {
	//HMDY
	TH1D *hHMDY_left =
	  (TH1D*)fHMDY->Get(Form("MuMu_left_%s_%s%i", targets[tr].Data(), physBinned[phys].Data(), bi) );
	hHMDY_left->Sumw2();
	hHMDY_left->Scale(1.0/(hHMDY_left->Integral() ) );
	if (physBinned[phys] == "pT") hHMDY_left->Fit("expo", "LQ0", "", 4.0, 4.3);//pT
	else hHMDY_left->Fit("expo", "LQ0", "", 4.0, 4.2);//xN, xPi, xF
    
	TF1 *fexp = (TF1*)hHMDY_left->GetFunction("expo");
	Double_t HMDY_val = fexp->Eval(3.5);

	//LMDY
	TH1D *hLMDY_left =
	  (TH1D*)fLMDY->Get(Form("MuMu_left_%s_%s%i", targets[tr].Data(), physBinned[phys].Data(), bi) );
	hLMDY_left->Sumw2();
	hLMDY_left->Scale(1.0/(hLMDY_left->Integral() ) );
	if (physBinned[phys] == "pT") hLMDY_left->Fit("expo", "LQ0", "" , 3.0, 3.2);//pT
	else hLMDY_left->Fit("expo", "LQ0", "" , 2.9, 3.1);//xN, xPi, xF
    
	fexp = (TF1*)hLMDY_left->GetFunction("expo");
	Double_t LMDY_val = fexp->Eval(3.5);

	//Combine and Draw
	TH1D *hSum_left = (TH1D*) hHMDY_left->Clone(); SetUp(hSum_left);
	hSum_left->Scale(LMDY_val/HMDY_val);
	hSum_left->Add(hLMDY_left);
	hSum_left->Scale(hSum_left->Integral() );
	hSum_left->Write(Form("MuMu_left_%s_%s%i", targets[tr].Data(), physBinned[phys].Data(), bi) );

	/////Right////
	//HMDY
	TH1D *hHMDY_right =
	  (TH1D*)fHMDY->Get(Form("MuMu_right_%s_%s%i", targets[tr].Data(), physBinned[phys].Data(), bi) );
	hHMDY_right->Sumw2();
	hHMDY_right->Scale(1.0/(hHMDY_right->Integral() ) );
	if (physBinned[phys] == "pT") hHMDY_right->Fit("expo", "LQ0", "", 4.0, 4.3);//pT
	else hHMDY_right->Fit("expo", "LQ0", "", 4.0, 4.2);//xN, xPi, xF
    
	fexp = (TF1*)hHMDY_right->GetFunction("expo");
	HMDY_val = fexp->Eval(3.5);

	//LMDY
	TH1D *hLMDY_right =
	  (TH1D*)fLMDY->Get(Form("MuMu_right_%s_%s%i", targets[tr].Data(), physBinned[phys].Data(), bi) );
	hLMDY_right->Sumw2();
	hLMDY_right->Scale(1.0/(hLMDY_right->Integral() ) );
	if (physBinned[phys] == "pT") hLMDY_right->Fit("expo", "LQ0", "" , 3.0, 3.2);//pT
	else hLMDY_right->Fit("expo", "LQ0", "" , 2.9, 3.1);//xN, xPi, xF
    
	fexp = (TF1*)hLMDY_right->GetFunction("expo");
	LMDY_val = fexp->Eval(3.5);

	//Combine and Draw
	TH1D *hSum_right = (TH1D*) hHMDY_right->Clone(); SetUp(hSum_right);
	hSum_right->Scale(LMDY_val/HMDY_val);
	hSum_right->Add(hLMDY_right);
	hSum_right->Scale(hSum_right->Integral() );
	hSum_right->Write(Form("MuMu_right_%s_%s%i", targets[tr].Data(), physBinned[phys].Data(), bi) );
      }
    }
  }
  fOutput->Close();
    
}
