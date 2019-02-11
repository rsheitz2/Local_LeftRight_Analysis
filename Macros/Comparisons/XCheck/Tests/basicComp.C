#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void basicComp(TString start=""){
  //Setup_______________
  TString Mtype ="HMDY";
  TString production ="slot1";

  Bool_t toWrite =false;
  //Setup_______________

  //Basic Setup
  TString localPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/XCheck";
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData";
  Int_t nhbin =150;

  //Basic setup
  TCanvas* cComp = new TCanvas(); cComp->Divide(4);

  //My distributions
  TString myData =Form("%s/%s/%sW07_%s.root", pathData.Data(), Mtype.Data(),
		       production.Data(), Mtype.Data());
  TFile *fMe = OpenFile(myData);
  TTree *myT = (TTree*)fMe->Get("pT_Weighted");
  TH1D *hMe[4];
  hMe[0] = new TH1D ("hMe_xpi", "hMe_xpi", nhbin, 0, 1);
  hMe[1] = new TH1D ("hMe_xp", "hMe_xp", nhbin, 0, 1);
  hMe[2] = new TH1D ("hMe_qt", "hMe_qt", nhbin, 0, 7);
  hMe[3] = new TH1D ("hMe_phis", "hMe_phis", nhbin, -TMath::Pi(), TMath::Pi());
  
  cComp->cd(1);
  myT->Draw("x_beam>>hMe_xpi", "", "E1");

  cComp->cd(2);
  myT->Draw("x_target>>hMe_xp", "", "E1");

  cComp->cd(3);
  myT->Draw("q_transverse>>hMe_qt", "", "E1");

  cComp->cd(4);
  //myT->Draw("PhiS>>hMe_phis", "", "E1");
  myT->Draw("PhiS_simple>>hMe_phis", "", "E1");

  //XCheck distributions
  TString xCheckData =
    Form("%s/Data/BasicDist/W07_An_xcheck.root", localPath.Data() );
  TFile *fxCheck = OpenFile(xCheckData);
  TTree *xCheckT = (TTree*)fxCheck->Get("An");
  TH1D *hxCheck[4];
  hxCheck[0] = new TH1D ("hxCheck_xpi", "hxCheck_xpi", nhbin, 0, 1);
  hxCheck[1] = new TH1D ("hxCheck_xp", "hxCheck_xp", nhbin, 0, 1);
  hxCheck[2] = new TH1D ("hxCheck_qt", "hxCheck_qt", nhbin, 0, 7);
  hxCheck[3] = new TH1D ("hxCheck_phis", "hxCheck_phis",
			 nhbin, -TMath::Pi(), TMath::Pi());
  xCheckT->SetLineColor(kRed);
  
  cComp->cd(1);

  xCheckT->Draw("xpi>>hxCheck_xpi", "", "E1 sames");

  cComp->cd(2);
  xCheckT->Draw("xp>>hxCheck_xp", "", "E1 sames");

  cComp->cd(3);
  xCheckT->Draw("qt>>hxCheck_qt", "", "E1 sames");

  cComp->cd(4);
  xCheckT->Draw("PhiS>>hxCheck_phis", "", "E1 sames");

  //Divison comparison
  TCanvas* c2 = new TCanvas(); c2->Divide(4);
  TH1D *hDivide[4];
  for (Int_t i=0; i<4; i++) {
    hMe[i]->Sumw2(); hxCheck[i]->Sumw2();
    hxCheck[i]->SetLineColor(kRed);
    hDivide[i] = (TH1D*)hMe[i]->Clone();
    hDivide[i]->Divide(hxCheck[i]);
    
    c2->cd(i+1);
    hDivide[i]->Draw("E1");
  }

}
