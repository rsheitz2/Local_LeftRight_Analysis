#include "include/helperFunctions.h"

Double_t Mod (Double_t *x, Double_t *p){
  return p[0]*(1+p[1]*TMath::Sin(x[0]));
}


void generator(){
  //Setup_______________
  Double_t phiScut =0.0;
  TString additionalCuts ="phiS0.0";
  Double_t A_siv =0.0;
  Int_t N_gen =10000;
  Bool_t alphaScale =false;

  TString production ="slot1";
  TString physBinned ="xN";//"xF", "pT"
  const Int_t nBins =1;//# of physBinned bins

  Bool_t toWrite =true;
  //Setup_______________
  if (nBins > 1){
    cout << "This only works for nbins == 1 for now" << endl;
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
  }
  else{
    alpha_upS_up =1.0; alpha_upS_down =1.0;
    alpha_downS_up =1.0; alpha_downS_down =1.0;
  }
    
  //Aesthetic setup
  Double_t ex[nBins] ={0.0}, xvals[nBins] ={1.0};

  //Generator functions
  TF1* f_upS_up = new TF1("f_upS_up", Mod, -TMath::Pi(), TMath::Pi(), 2);
  f_upS_up->SetParameters(1.0/(2*TMath::Pi()), A_siv);
  TH1D* hGen_upS_up =
    new TH1D("hGen_upS_up", "hGen_upS_up", 100, -TMath::Pi(), TMath::Pi());

  TF1* f_upS_down = new TF1("f_upS_down", Mod, -TMath::Pi(), TMath::Pi(), 2);
  f_upS_down->SetParameters(1.0/(2*TMath::Pi()), A_siv);
  TH1D* hGen_upS_down =
    new TH1D("hGen_upS_down", "hGen_upS_down", 100, -TMath::Pi(), TMath::Pi());

  TF1* f_downS_up = new TF1("f_downS_up", Mod, -TMath::Pi(), TMath::Pi(), 2);
  f_downS_up->SetParameters(1.0/(2*TMath::Pi()), A_siv);
  TH1D* hGen_downS_up =
    new TH1D("hGen_downS_up", "hGen_downS_up", 100, -TMath::Pi(), TMath::Pi());

  TF1* f_downS_down = new TF1("f_downS_down", Mod, -TMath::Pi(), TMath::Pi(),2);
  f_downS_down->SetParameters(1.0/(2*TMath::Pi()), A_siv);
  TH1D* hGen_downS_down =
    new TH1D("hGen_downS_down", "hGen_downS_down",100,-TMath::Pi(),TMath::Pi());

  //upS_up
  Double_t left_upS_up[nBins] ={0.0}, right_upS_up[nBins] ={0.0};
  for (Int_t i=0; i<N_gen*alpha_upS_up; i++) {
    Double_t phiS =f_upS_up->GetRandom();

    //Generated
    hGen_upS_up->Fill(phiS);
    
    //PhiS cut
    if ( (phiS < phiScut) && (phiS > -phiScut) ) continue;
    if ( (phiS < -TMath::Pi() + phiScut) || (phiS > TMath::Pi() -phiScut) ) {
      continue; 
    }

    if ( (phiS < TMath::Pi() ) && (phiS > 0) ) left_upS_up[0]++;
    else if ( (phiS > -TMath::Pi() ) && (phiS < 0) ) right_upS_up[0]++;
    else {
      cout << "PhiS generatored wrong" << endl;
      exit(EXIT_FAILURE);
    }
  }
  Double_t e_left_upS_up[nBins] = {TMath::Sqrt(left_upS_up[0])};
  Double_t e_right_upS_up[nBins] = {TMath::Sqrt(right_upS_up[0])};
  TGraphErrors* g_left_upS_up =
    new TGraphErrors(nBins, xvals, left_upS_up, ex, e_left_upS_up);
  SetUp(g_left_upS_up);
  TGraphErrors* g_right_upS_up =
    new TGraphErrors(nBins, xvals, right_upS_up, ex, e_right_upS_up);
  SetUp(g_right_upS_up);

  //upS_down
  Double_t left_upS_down[nBins] ={0.0}, right_upS_down[nBins] ={0.0};
  for (Int_t i=0; i<N_gen*alpha_upS_down; i++) {
    Double_t phiS =f_upS_down->GetRandom();

    //Generated
    hGen_upS_down->Fill(phiS);
    
    //PhiS cut
    if ( (phiS < phiScut) && (phiS > -phiScut) ) continue;
    if ( (phiS < -TMath::Pi() + phiScut) || (phiS > TMath::Pi() -phiScut) ) {
      continue; 
    }

    if ( (phiS < TMath::Pi() ) && (phiS > 0) ) left_upS_down[0]++;
    else if ( (phiS > -TMath::Pi() ) && (phiS < 0) ) right_upS_down[0]++;
    else {
      cout << "PhiS generatored wrong" << endl;
      exit(EXIT_FAILURE);
    }
  }
  Double_t e_left_upS_down[nBins] = {TMath::Sqrt(left_upS_down[0])};
  Double_t e_right_upS_down[nBins] = {TMath::Sqrt(right_upS_down[0])};
  TGraphErrors* g_left_upS_down =
    new TGraphErrors(nBins, xvals, left_upS_down, ex, e_left_upS_down);
  SetUp(g_left_upS_down);
  TGraphErrors* g_right_upS_down =
    new TGraphErrors(nBins, xvals, right_upS_down, ex, e_right_upS_down);
  SetUp(g_right_upS_down);

    //downS_up
  Double_t left_downS_up[nBins] ={0.0}, right_downS_up[nBins] ={0.0};
  for (Int_t i=0; i<N_gen*alpha_downS_up; i++) {
    Double_t phiS =f_downS_up->GetRandom();

    //Generated
    hGen_downS_up->Fill(phiS);
    
    //PhiS cut
    if ( (phiS < phiScut) && (phiS > -phiScut) ) continue;
    if ( (phiS < -TMath::Pi() + phiScut) || (phiS > TMath::Pi() -phiScut) ) {
      continue; 
    }

    if ( (phiS < TMath::Pi() ) && (phiS > 0) ) left_downS_up[0]++;
    else if ( (phiS > -TMath::Pi() ) && (phiS < 0) ) right_downS_up[0]++;
    else {
      cout << "PhiS generatored wrong" << endl;
      exit(EXIT_FAILURE);
    }
  }
  Double_t e_left_downS_up[nBins] = {TMath::Sqrt(left_downS_up[0])};
  Double_t e_right_downS_up[nBins] = {TMath::Sqrt(right_downS_up[0])};
  TGraphErrors* g_left_downS_up =
    new TGraphErrors(nBins, xvals, left_downS_up, ex, e_left_downS_up);
  SetUp(g_left_downS_up);
  TGraphErrors* g_right_downS_up =
    new TGraphErrors(nBins, xvals, right_downS_up, ex, e_right_downS_up);
  SetUp(g_right_downS_up);

  //downS_down
  Double_t left_downS_down[nBins] ={0.0}, right_downS_down[nBins] ={0.0};
  for (Int_t i=0; i<N_gen*alpha_downS_down; i++) {
    Double_t phiS =f_downS_down->GetRandom();

    //Generated
    hGen_downS_down->Fill(phiS);
    
    //PhiS cut
    if ( (phiS < phiScut) && (phiS > -phiScut) ) continue;
    if ( (phiS < -TMath::Pi() + phiScut) || (phiS > TMath::Pi() -phiScut) ) {
      continue; 
    }

    if ( (phiS < TMath::Pi() ) && (phiS > 0) ) left_downS_down[0]++;
    else if ( (phiS > -TMath::Pi() ) && (phiS < 0) ) right_downS_down[0]++;
    else {
      cout << "PhiS generatored wrong" << endl;
      exit(EXIT_FAILURE);
    }
  }
  Double_t e_left_downS_down[nBins] = {TMath::Sqrt(left_downS_down[0])};
  Double_t e_right_downS_down[nBins] = {TMath::Sqrt(right_downS_down[0])};
  TGraphErrors* g_left_downS_down =
    new TGraphErrors(nBins, xvals, left_downS_down, ex, e_left_downS_down);
  SetUp(g_left_downS_down);
  TGraphErrors* g_right_downS_down =
    new TGraphErrors(nBins, xvals, right_downS_down, ex, e_right_downS_down);
  SetUp(g_right_downS_down);

  //Draw LR results
  TCanvas* cLR = new TCanvas(); cLR->Divide(2, 4);
  cLR->cd(1); g_left_upS_up->Draw("AP");
  cLR->cd(2); g_right_upS_up->Draw("AP");
  cLR->cd(3); g_left_upS_down->Draw("AP");
  cLR->cd(4); g_right_upS_down->Draw("AP");
  cLR->cd(5); g_left_downS_up->Draw("AP");
  cLR->cd(6); g_right_downS_up->Draw("AP");
  cLR->cd(7); g_left_downS_down->Draw("AP");
  cLR->cd(8); g_right_downS_down->Draw("AP");

  //Draw Gen
  TCanvas* cGen = new TCanvas();
  hGen_upS_up->Draw();
  hGen_upS_down->Draw("same"); hGen_upS_down->SetLineColor(kRed);
  hGen_downS_up->Draw("same"); hGen_downS_up->SetLineColor(kGreen);
  hGen_downS_down->Draw("same"); hGen_downS_down->SetLineColor(kBlue);

  //Final output
  if (alphaScale) N_gen /= 4.0;
  TString thisDirPath = "/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/ToyMC/Data";
  TString fOutput = Form("%s/Gen/generator_%0.3f_%s_gen%i_alpha%i.root",
			 thisDirPath.Data(), A_siv, additionalCuts.Data(),
			 N_gen, alphaScale);
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");

    g_left_upS_up->Write("Left_upS_up");
    g_left_upS_down->Write("Left_upS_down");
    g_left_downS_up->Write("Left_downS_up");
    g_left_downS_down->Write("Left_downS_down");

    g_right_upS_up->Write("Right_upS_up");
    g_right_upS_down->Write("Right_upS_down");
    g_right_downS_up->Write("Right_downS_up");
    g_right_downS_down->Write("Right_downS_down");

    Double_t Pol[1] = {1.0};
    TGraph *g_Pol = new TGraph(1, xvals, Pol);
    g_Pol->Write(Form("%s_Pol", physBinned.Data()) );

    hGen_upS_up->Write();
    hGen_upS_down->Write();
    hGen_downS_up->Write();
    hGen_downS_down->Write();

    fResults->Close();
  }
}
