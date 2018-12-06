#include "Include/helperFunctions.h"

void writeFile_mcLMDYandHMDY(){
  //Setup_____
	Bool_t toWrite =false;
  //Setup_____

  cout << "\nMacro combines low and high mass drell yan by fitting with an ";
  cout << "exponential near the boundary and extrapolating in between" << endl;
  cout << "\nOutput is written to a file:   " << toWrite << endl;
	cout << "" << endl;
  
	TString path =
		"/u/sciteam/heitz/Analysis/TGeant/Main/Data/MainOutput/Charles_Official/";
  TString nameLMDY ="Charles_W12_LMDY.root";
  TString nameHMDY ="Charles_W12_HMDY.root";
  TFile *fHMDY = TFile::Open(path+nameHMDY);
  TFile *fLMDY = TFile::Open(path+nameLMDY);

  TTree *tLMDY = (TTree*) fLMDY->Get("pT_Weighted");
  TTree *tHMDY = (TTree*) fHMDY->Get("pT_Weighted");

	TH1D *hLMDY = new TH1D("h_LMDY", "h_LMDY", 200, 2.0, 8.5);
	TH1D *hHMDY = new TH1D("h_HMDY", "h_HMDY", 200, 2.0, 8.5);
	
  tLMDY->Draw("Mmumu>>h_LMDY", "Mmumu>2.0&&Mmumu<8.5", "0");
  tHMDY->Draw("Mmumu>>h_HMDY", "Mmumu>2.0&&Mmumu<8.5", "0");
  
	//HMDY
	hHMDY->Sumw2();
	hHMDY->Scale(1.0/(hHMDY->Integral() ) );
	hHMDY->Fit("expo", "LQ0", "", 4.0, 4.3);
    
	TF1 *fexp = (TF1*)hHMDY->GetFunction("expo");
	Double_t HMDY_val = fexp->Eval(3.5);

	//LMDY
	hLMDY->Sumw2();
	hLMDY->Scale(1.0/(hLMDY->Integral() ) );
	hLMDY->Fit("expo", "LQ0", "" , 3.0, 3.2);
    
	fexp = (TF1*)hLMDY->GetFunction("expo");
	Double_t LMDY_val = fexp->Eval(3.5);

	//Combine and Draw
	TH1D *hSum = (TH1D*) hHMDY->Clone(); Setup(hSum);
	hSum->Scale(LMDY_val/HMDY_val);
	hSum->Add(hLMDY);

	TCanvas *c1 = new TCanvas(); c1->Divide(2);
	c1->cd(1); gPad->SetLogy();
	hSum->Draw();
	hLMDY->Draw("same"); hHMDY->Draw("same");
	hLMDY->SetLineColor(kBlack);
	hHMDY->SetLineColor(kRed);
	hHMDY->Scale(LMDY_val/HMDY_val);
	c1->cd(2); gPad->SetLogy();
	hSum->Draw();

	Double_t scale[] = {LMDY_val/HMDY_val};
	Double_t xval[] = {1};
	TGraph *gScale = new TGraph(1, xval, scale);
	Setup(gScale);

	TString outName ="M_AMDY.root";
	if (toWrite){
  TFile *fOutput = new TFile(path+"Macros/Data/MassDist/"+outName, "RECREATE");
	hSum->Write("h_AMDY");
	gScale->Write("LMDY_HMDY_scale");
  fOutput->Close();
	}

	//Final output
	cout << "\nFile  " << outName << "   was written:   " << toWrite << endl;
	cout << "" << endl;
}
