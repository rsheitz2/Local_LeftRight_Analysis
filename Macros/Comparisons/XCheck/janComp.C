#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

Double_t ShiftAngle(Double_t angle, Double_t phast=TMath::Pi());

void janComp(TString start=""){
  //Setup_______________
  TString Mtype ="HMDY";
  TString production ="t3"; //"slot1", "t3"

  Bool_t toWrite =false;
  //Setup_______________

  //Basic Setup
  TString localPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/XCheck";
  TString pTPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
pT_weight/Data";
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData";
  Int_t nhbin =50;

  //Get Data/setup trees
  if (production == "t3") production = "";
  /*TString myData =Form("%s/%s/%sW07_%s.root", pathData.Data(), Mtype.Data(),
    production.Data(), Mtype.Data());//*/
  TString myData =Form("%s/Main/mainW07.root", pTPath.Data() );
  TFile *fMe = OpenFile(myData);
  TTree *myT = (TTree*)fMe->Get("pT_Weighted");
  if (production == "") production = "t3";

  TString xCheckData =
    Form("%s/JanData/dy15w07_hmass_asym.root", pTPath.Data() );
  TFile *fxCheck = OpenFile(xCheckData);
  TH1D *hxCheckPhiS = (TH1D*)
    fxCheck->Get("bins_integrated/kinematic_distr/hphiS_DR_integrated_0");

  Double_t x_beam_me, q_transverse_me, Spin_0, PhiS_me, PhiS_simple_me;
  Int_t targetPosition;
  myT->SetBranchAddress("x_beam", &x_beam_me);
  myT->SetBranchAddress("Spin_0", &Spin_0);
  myT->SetBranchAddress("PhiS", &PhiS_me);
  myT->SetBranchAddress("PhiS_simple", &PhiS_simple_me);
  myT->SetBranchAddress("targetPosition", &targetPosition);

  //Distributions
  TString whichTarget[] = {"upS_pUp", "upS_pDown", "downS_pUp", "downS_pDown",
			   "all"};
  TH1D *hMe[5], *hMe_simple[5], *hxCheck[5];
  for (Int_t i=0; i<5; i++) {
    hMe[i] = new TH1D (Form("hMe_phis_%s", whichTarget[i].Data()),
		       Form("hMe_phis_%s", whichTarget[i].Data()),
		       nhbin, -TMath::Pi(), TMath::Pi());
    hMe_simple[i] = new TH1D (Form("hMe_simple_phis_%s", whichTarget[i].Data()),
			      Form("hMe_simple_phis_%s", whichTarget[i].Data()),
			      nhbin, -TMath::Pi(), TMath::Pi());
    hxCheck[i] = new TH1D (Form("hxCheck_phis_%s", whichTarget[i].Data()),
			   Form("hxCheck_phis_%s", whichTarget[i].Data()),
			   nhbin, -TMath::Pi(), TMath::Pi());    
  }

  //Tree loop
  for (Int_t ev=0; ev<myT->GetEntries(); ev++) {
    myT->GetEntry(ev);

    Double_t PhiS_xc;
    if ((targetPosition==0) && (Spin_0 > 0)){//upstream Pup
      //PhiS_xc = PhiS_me;
      //PhiS_xc = ShiftAngle(PhiS_me);
      /*if (PhiS_me > 0)
	PhiS_xc = PhiS_me - TMath::Pi();
      else
	PhiS_xc = PhiS_me + TMath::Pi();//*/
      //PhiS_xc = PhiS_simple_me;
      PhiS_xc = ShiftAngle(PhiS_simple_me);
      
      hMe[0]->Fill(PhiS_me);
      hMe_simple[0]->Fill(PhiS_simple_me);
      hxCheck[0]->Fill(PhiS_xc);
    }
    else if ((targetPosition==0) && (Spin_0 < 0)){//upstream Pdown
      //PhiS_xc = PhiS_me;
      //PhiS_xc = PhiS_simple_me;
      /*if (PhiS_me > 0)
	PhiS_xc = PhiS_me - TMath::Pi();
      else
      PhiS_xc = PhiS_me + TMath::Pi();//*/
      PhiS_xc = ShiftAngle(PhiS_simple_me);

      hMe[1]->Fill(PhiS_me);
      hMe_simple[1]->Fill(PhiS_simple_me);
      hxCheck[1]->Fill(PhiS_xc);
    }
    else if ((targetPosition==1) && (Spin_0 > 0)){//downstream Pup
      //PhiS_xc = PhiS_me;
      //PhiS_xc = ShiftAngle(PhiS_me);
      /*if (PhiS_me > 0)
	PhiS_xc = PhiS_me - TMath::Pi();
      else
	PhiS_xc = PhiS_me + TMath::Pi();//*/
      //PhiS_xc = PhiS_simple_me;
      PhiS_xc = ShiftAngle(PhiS_simple_me);
      
      hMe[2]->Fill(PhiS_me);
      hMe_simple[2]->Fill(PhiS_simple_me);
      hxCheck[2]->Fill(PhiS_xc);
    }
    else if ((targetPosition==1) && (Spin_0 < 0)){//downstream Pdown
      //PhiS_xc = PhiS_me;
      //PhiS_xc = PhiS_simple_me;
      /*if (PhiS_me > 0)
	PhiS_xc = PhiS_me - TMath::Pi();
      else
      PhiS_xc = PhiS_me + TMath::Pi();//*/
      PhiS_xc = ShiftAngle(PhiS_simple_me);
	    
      hMe[3]->Fill(PhiS_me);
      hMe_simple[3]->Fill(PhiS_simple_me);
      hxCheck[3]->Fill(PhiS_xc);
    }

    hMe[4]->Fill(PhiS_me);
    hMe_simple[4]->Fill(PhiS_simple_me);
    hxCheck[4]->Fill(PhiS_xc);
  }//end Me tree loop

  /*//Draw distributions
  TCanvas* cMe = new TCanvas(); cMe->Divide(2,2);
  TCanvas* cMe_simple = new TCanvas(); cMe_simple->Divide(2,2);
  TCanvas* cXC = new TCanvas(); cXC->Divide(2,2);
  for (Int_t i=0; i<4; i++) {
    cMe->cd(i+1);
    hMe[i]->Draw();

    cMe_simple->cd(i+1);
    hMe_simple[i]->Draw();

    cXC->cd(i+1);
    hxCheck[i]->Draw();
    }//*/

  TCanvas* cTot = new TCanvas(); cTot->Divide(2,2);
  cTot->cd(1);
  hMe[4]->Draw("E1");

  cTot->cd(2);
  hMe_simple[4]->Draw("E1");

  cTot->cd(3);
  hxCheck[4]->Draw("E1");

  cTot->cd(4);
  hxCheckPhiS->Draw();

  TCanvas* cComp = new TCanvas();
  hxCheck[4]->Draw("E1");
  hxCheckPhiS->Draw("same"); hxCheckPhiS->SetLineColor(kRed);
}


Double_t ShiftAngle(Double_t angle, Double_t phase=TMath::Pi()){
  if (angle > 0)
    return angle - phase;
  else
    return angle + phase;
}
