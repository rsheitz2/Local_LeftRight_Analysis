#include "include/helperFunctions.h"

Double_t Mod (Double_t *x, Double_t *p){
  return p[0]*(1+p[1]*TMath::Sin(x[0]));
}


void GenLoop(Int_t N_gen, Double_t A_siv, Double_t phiScut,
	     Double_t alpha, TH1D *h, Double_t *Left, Double_t *Right){
  //Define generation function
  TF1* f = new TF1("f_upS_up", Mod, -TMath::Pi(), TMath::Pi(), 2);
  f->SetParameters(1.0/(2*TMath::Pi()), A_siv);

  //Random number generators
  TRandom2 *rand = new TRandom2();
  TRandom2 *phiRand = new TRandom2();
  rand->SetSeed();
  phiRand->SetSeed();
  
  for (Int_t i=0; i<N_gen*alpha; i++) {
    //Random generation
    Double_t phiS = phiRand->Uniform(-TMath::Pi(), TMath::Pi() );
    Double_t dist = f->Eval(phiS);
    if (rand->Uniform()*(1+A_siv)/(2.0*TMath::Pi()) > dist ) continue;
    
    //Generated
    h->Fill(phiS);
    
    //PhiS cut
    if ( (phiS < phiScut) && (phiS > -phiScut) ) continue;
    if ( (phiS < -TMath::Pi() + phiScut) || (phiS > TMath::Pi() -phiScut) ) {
      continue; 
    }

    if ( (phiS < TMath::Pi() ) && (phiS > 0) ) Left[0]++;
    else if ( (phiS > -TMath::Pi() ) && (phiS < 0) ) Right[0]++;
    else {
      cout << "PhiS generatored wrong" << endl;
      exit(EXIT_FAILURE);
    }
  }//N_gen loop
}//GenLoop


void genRej(){
  //Setup_______________
  Double_t phiScut =1.44;
  TString additionalCuts ="phiS1.44";
  Double_t A_siv =0.8;
  Int_t N_gen =1000;
  Bool_t alphaScale =true;

  TString production ="slot1";
  TString physBinned ="xN";//"xF", "pT"
  const Int_t nBins =1;//# of physBinned bins

  Bool_t toWrite =false;
  //Setup_______________
  if (nBins > 1){
    cout << "This only works for nbins == 1 for now" << endl;
    exit(EXIT_FAILURE);
  }
  if (A_siv > 1){
    cout << "Amplitude generation is too larger for generatio" << endl;
    cout << A_siv << "  > 1  " << endl;
    exit(EXIT_FAILURE);
  }

  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =
    Form("leftRight_byTarget_WAll_HMDY4.30_8.50_%ibins43_85_150hbin_%s_%s.root",
	 nBins, production.Data(), additionalCuts.Data());
  TFile *fRD  = OpenFile(pathRD + RDfile);

  TGraph *g_alpha_upS_up =
    (TGraph*)fRD->Get(Form("%s_left_upstream_up", physBinned.Data()));
  Double_t y_alpha_upS_up = g_alpha_upS_up->GetY()[0];
  TGraph *g_alpha_upS_down =
    (TGraph*)fRD->Get(Form("%s_left_upstream_down", physBinned.Data()));
  Double_t y_alpha_upS_down = g_alpha_upS_down->GetY()[0];
  TGraph *g_alpha_downS_up =
    (TGraph*)fRD->Get(Form("%s_left_downstream_up", physBinned.Data()));
  Double_t y_alpha_downS_up = g_alpha_downS_up->GetY()[0];
  TGraph *g_alpha_downS_down =
    (TGraph*)fRD->Get(Form("%s_left_downstream_down", physBinned.Data()));
  Double_t y_alpha_downS_down = g_alpha_downS_down->GetY()[0];

  Double_t norm = y_alpha_upS_up + y_alpha_upS_down +
    y_alpha_downS_up + y_alpha_downS_down;

  //Amount of data to generate
  Double_t alpha_upS_up, alpha_upS_down;
  Double_t alpha_downS_up, alpha_downS_down;
  if (alphaScale) {
    N_gen *= 4.0;
    alpha_upS_up =y_alpha_upS_up/norm;
    alpha_upS_down =y_alpha_upS_down/norm;
    alpha_downS_up =y_alpha_downS_up/norm;
    alpha_downS_down =y_alpha_downS_down/norm;
    cout << alpha_upS_up << " " << alpha_upS_down;
    cout << " " << alpha_downS_up << " " << alpha_downS_down << endl;
  }
  else{
    alpha_upS_up =1.0; alpha_upS_down =1.0;
    alpha_downS_up =1.0; alpha_downS_down =1.0;
  }
  
  //Aesthetic setup
  Double_t ex[nBins] ={0.0}, xvals[nBins] ={1.0};

  //upS_up
  TH1D* hGen_upS_up =
    new TH1D("hGen_upS_up", "hGen_upS_up", 100, -TMath::Pi(), TMath::Pi());

  Double_t Left_upS_up[nBins] ={0.0}, Right_upS_up[nBins] ={0.0};
  GenLoop(N_gen, A_siv, phiScut,
	  alpha_upS_up, hGen_upS_up, Left_upS_up, Right_upS_up);
  
  Double_t e_Left_upS_up[nBins] = {TMath::Sqrt(Left_upS_up[0])};
  Double_t e_Right_upS_up[nBins] = {TMath::Sqrt(Right_upS_up[0])};
  TGraphErrors* g_Left_upS_up =
    new TGraphErrors(nBins, xvals, Left_upS_up, ex, e_Left_upS_up);
  TGraphErrors* g_Right_upS_up =
    new TGraphErrors(nBins, xvals, Right_upS_up, ex, e_Right_upS_up);
  SetUp(g_Left_upS_up); SetUp(g_Right_upS_up);

  //upS_down
  TH1D* hGen_upS_down =
    new TH1D("hGen_upS_down", "hGen_upS_down", 100, -TMath::Pi(), TMath::Pi());

  Double_t Left_upS_down[nBins] ={0.0}, Right_upS_down[nBins] ={0.0};
  GenLoop(N_gen, A_siv, phiScut,
	  alpha_upS_down, hGen_upS_down, Left_upS_down, Right_upS_down);
  
  Double_t e_Left_upS_down[nBins] = {TMath::Sqrt(Left_upS_down[0])};
  Double_t e_Right_upS_down[nBins] = {TMath::Sqrt(Right_upS_down[0])};
  TGraphErrors* g_Left_upS_down =
    new TGraphErrors(nBins, xvals, Left_upS_down, ex, e_Left_upS_down);
  TGraphErrors* g_Right_upS_down =
    new TGraphErrors(nBins, xvals, Right_upS_down, ex, e_Right_upS_down);
  SetUp(g_Left_upS_down); SetUp(g_Right_upS_down);

  //downS_up
  TH1D* hGen_downS_up =
    new TH1D("hGen_downS_up", "hGen_downS_up", 100, -TMath::Pi(), TMath::Pi());

  Double_t Left_downS_up[nBins] ={0.0}, Right_downS_up[nBins] ={0.0};
  GenLoop(N_gen, A_siv, phiScut,
	  alpha_downS_up, hGen_downS_up, Left_downS_up, Right_downS_up);
  
  Double_t e_Left_downS_up[nBins] = {TMath::Sqrt(Left_downS_up[0])};
  Double_t e_Right_downS_up[nBins] = {TMath::Sqrt(Right_downS_up[0])};
  TGraphErrors* g_Left_downS_up =
    new TGraphErrors(nBins, xvals, Left_downS_up, ex, e_Left_downS_up);
  TGraphErrors* g_Right_downS_up =
    new TGraphErrors(nBins, xvals, Right_downS_up, ex, e_Right_downS_up);
  SetUp(g_Left_downS_up); SetUp(g_Right_downS_up);

  //downS_down
  TH1D* hGen_downS_down =
    new TH1D("hGen_downS_down", "hGen_downS_down",100,-TMath::Pi(),TMath::Pi());

  Double_t Left_downS_down[nBins] ={0.0}, Right_downS_down[nBins] ={0.0};
  GenLoop(N_gen, A_siv, phiScut,
	  alpha_downS_down, hGen_downS_down, Left_downS_down, Right_downS_down);
  
  Double_t e_Left_downS_down[nBins] = {TMath::Sqrt(Left_downS_down[0])};
  Double_t e_Right_downS_down[nBins] = {TMath::Sqrt(Right_downS_down[0])};
  TGraphErrors* g_Left_downS_down =
    new TGraphErrors(nBins, xvals, Left_downS_down, ex, e_Left_downS_down);
  TGraphErrors* g_Right_downS_down =
    new TGraphErrors(nBins, xvals, Right_downS_down, ex, e_Right_downS_down);
  SetUp(g_Left_downS_down); SetUp(g_Right_downS_down);

  //Draw LR results
  TCanvas* cLR = new TCanvas(); cLR->Divide(2, 4);
  cLR->cd(1); g_Left_upS_up->Draw("AP");
  cLR->cd(2); g_Right_upS_up->Draw("AP");
  cLR->cd(3); g_Left_upS_down->Draw("AP");
  cLR->cd(4); g_Right_upS_down->Draw("AP");
  cLR->cd(5); g_Left_downS_up->Draw("AP");
  cLR->cd(6); g_Right_downS_up->Draw("AP");
  cLR->cd(7); g_Left_downS_down->Draw("AP");
  cLR->cd(8); g_Right_downS_down->Draw("AP");
  
  //Draw Gen
  TCanvas* cGen = new TCanvas();
  hGen_upS_up->Draw();
  hGen_upS_down->Draw("sames"); hGen_upS_down->SetLineColor(kRed);
  hGen_downS_up->Draw("sames"); hGen_downS_up->SetLineColor(kGreen);
  hGen_downS_down->Draw("sames"); hGen_downS_down->SetLineColor(kBlue);

  //Final output
  if (alphaScale) N_gen /=4.0;//N_gen was scaled up by 4 bc alpha scaled down
  TString thisDirPath = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/ToyMC/Data";
  TString fOutput = Form("%s/GenRej/genRej_%0.3f_%s_gen%i_alpha%i.root",
			 thisDirPath.Data(), A_siv, additionalCuts.Data(),
			 N_gen, alphaScale);
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");

    g_Left_upS_up->Write("Left_upS_up");
    g_Left_upS_down->Write("Left_upS_down");
    g_Left_downS_up->Write("Left_downS_up");
    g_Left_downS_down->Write("Left_downS_down");
    
    g_Right_upS_up->Write("Right_upS_up");
    g_Right_upS_down->Write("Right_upS_down");
    g_Right_downS_up->Write("Right_downS_up");
    g_Right_downS_down->Write("Right_downS_down");
    
    Double_t Pol[1] = {1.0};
    TGraph *g_Pol = new TGraph(1, xvals, Pol);
    g_Pol->Write(Form("%s_Pol", physBinned.Data()) );

    hGen_upS_up->Write();
    hGen_upS_down->Write();
    hGen_downS_up->Write();
    hGen_downS_down->Write();
    
    fResults->Close();
  }

  cout << "\nSettings________" << endl;
  cout << "Data coming from:            " << pathRD << endl;
  cout << "Input P corrected data:        " << RDfile << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
  cout << "PhiS cut used:         " << phiScut << endl;
  cout << "additionalCuts used:   " << additionalCuts << endl;
  cout << "Sivers amplitude:      " << A_siv << endl;
  cout << "Number of events generated:   " << N_gen << endl;
  cout << "Events per target scaled like real data:   " << alphaScale << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
  cout << " " << endl;
}
