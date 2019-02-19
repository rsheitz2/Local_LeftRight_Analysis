#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

Double_t ShiftAngle(Double_t angle, Double_t phast=TMath::Pi());

void mePhiSComp(TString start=""){
  //Setup_______________
  TString Mtype ="HMDY";
  TString production ="t3";

  Bool_t toWrite =false;
  //Setup_______________

  //Basic Setup
  TString localPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/XCheck";
  TString pTPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
pT_weight/Data";
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData";
  Int_t nhbin =150;

  //Get Data/setup trees
  if (production == "t3") production = "";
  TString newData =Form("%s/%s/%sW07_%s.root", pathData.Data(), Mtype.Data(),
		       production.Data(), Mtype.Data());
  TFile *fnew = OpenFile(newData);
  TTree *newT = (TTree*)fnew->Get("pT_Weighted");
  if (production == "") production = "t3";

  TString oldData =Form("%s/Main/mainW07.root", pTPath.Data() );
  TFile *fold = OpenFile(oldData);
  TTree *oldT = (TTree*)fold->Get("pT_Weighted");

  Double_t x_beam_new, x_feynman_new, Spin_0_new, PhiS_new, PhiS_simple_new;
  Int_t targetPosition_new;
  newT->SetBranchAddress("x_beam", &x_beam_new);
  newT->SetBranchAddress("x_feynman", &x_feynman_new);
  newT->SetBranchAddress("Spin_0", &Spin_0_new);
  newT->SetBranchAddress("PhiS", &PhiS_new);
  newT->SetBranchAddress("PhiS_simple", &PhiS_simple_new);
  newT->SetBranchAddress("targetPosition", &targetPosition_new);

  Double_t x_beam_old, x_feynman_old, Spin_0_old, PhiS_old, PhiS_simple_old;
  Int_t targetPosition_old;
  oldT->SetBranchAddress("x_beam", &x_beam_old);
  oldT->SetBranchAddress("x_feynman", &x_feynman_old);
  oldT->SetBranchAddress("Spin_0", &Spin_0_old);
  oldT->SetBranchAddress("PhiS", &PhiS_old);
  oldT->SetBranchAddress("PhiS_simple", &PhiS_simple_old);
  oldT->SetBranchAddress("targetPosition", &targetPosition_old);

  //Distributions
  TString whichTarget[] = {"upS_pUp", "upS_pDown", "downS_pUp", "downS_pDown"};
  TH1D *hnew[4], *hnew_simple[4];
  TH1D *hold[4], *hold_simple[4];
  TH1D *hnew_xF[4], *hnew_xPi[4];
  TH1D *hold_xF[4], *hold_xPi[4];
  for (Int_t i=0; i<4; i++) {
    hnew[i] = new TH1D (Form("hnew_phis_%s", whichTarget[i].Data()),
		       Form("hnew_phis_%s", whichTarget[i].Data()),
		       nhbin, -TMath::Pi(), TMath::Pi());
    hnew_simple[i] = new TH1D (Form("hnew_simple_phis_%s",
				    whichTarget[i].Data()),
			      Form("hnew_simple_phis_%s",
				   whichTarget[i].Data()),
			      nhbin, -TMath::Pi(), TMath::Pi());
    hnew_xF[i] = new TH1D (Form("hnew_xF_%s", whichTarget[i].Data()),
			   Form("hnew_xF_%s", whichTarget[i].Data()),
			   nhbin, -1, 1);
    hnew_xPi[i] = new TH1D (Form("hnew_xPi_%s", whichTarget[i].Data()),
			      Form("hnew_xPi_%s", whichTarget[i].Data()),
			      nhbin, 0, 1);

    hold[i] = new TH1D (Form("hold_phis_%s", whichTarget[i].Data()),
		       Form("hold_phis_%s", whichTarget[i].Data()),
		       nhbin, -TMath::Pi(), TMath::Pi());
    hold_simple[i] = new TH1D (Form("hold_simple_phis_%s",
				    whichTarget[i].Data()),
			      Form("hold_simple_phis_%s",
				   whichTarget[i].Data()),
			      nhbin, -TMath::Pi(), TMath::Pi());
    hold_xF[i] = new TH1D (Form("hold_xF_%s", whichTarget[i].Data()),
			   Form("hold_xF_%s", whichTarget[i].Data()),
			   nhbin, -1, 1);
    hold_xPi[i] = new TH1D (Form("hold_xPi_%s", whichTarget[i].Data()),
			      Form("hold_xPi_%s", whichTarget[i].Data()),
			      nhbin, 0, 1);
  }
  
  //Tree loop
  Double_t epsilon = 0.000001;
  for (Int_t ev_new=0; ev_new<newT->GetEntries(); ev_new++) {
    newT->GetEntry(ev_new);
    
    for (Int_t ev_old=0; ev_old<oldT->GetEntries(); ev_old++) {
      oldT->GetEntry(ev_old);

      if ((x_beam_new>=x_beam_old-epsilon) && (x_beam_new<=x_beam_old+epsilon)){
	if ((x_feynman_new >= x_feynman_old-epsilon) &&
	    (x_feynman_new <= x_feynman_old+epsilon)){

	  if ((targetPosition_new==0) && (Spin_0_new > 0)){//upstream Pup
	    hnew[0]->Fill(PhiS_new);
	    hnew_simple[0]->Fill(PhiS_simple_new);
	    hold[0]->Fill(PhiS_old);
	    hold_simple[0]->Fill(PhiS_simple_old);

	    hnew_xF[0]->Fill(x_feynman_new);
	    hnew_xPi[0]->Fill(x_beam_new);
	    hold_xF[0]->Fill(x_feynman_old);
	    hold_xPi[0]->Fill(x_beam_old);
	  }
	  else if ((targetPosition_new==0) && (Spin_0_new < 0)){//upstream Pdown
	    hnew[1]->Fill(PhiS_new);
	    hnew_simple[1]->Fill(PhiS_simple_new);
	    hold[1]->Fill(PhiS_old);
	    hold_simple[1]->Fill(PhiS_simple_old);

	    hnew_xF[1]->Fill(x_feynman_new);
	    hnew_xPi[1]->Fill(x_beam_new);
	    hold_xF[1]->Fill(x_feynman_old);
	    hold_xPi[1]->Fill(x_beam_old);
	  }
	  else if ((targetPosition_new==1) && (Spin_0_new > 0)){//downstream Pup
	    hnew[2]->Fill(PhiS_new);
	    hnew_simple[2]->Fill(PhiS_simple_new);
	    hold[2]->Fill(PhiS_old);
	    hold_simple[2]->Fill(PhiS_simple_old);

	    hnew_xF[2]->Fill(x_feynman_new);
	    hnew_xPi[2]->Fill(x_beam_new);
	    hold_xF[2]->Fill(x_feynman_old);
	    hold_xPi[2]->Fill(x_beam_old);
	  }
	  else if ((targetPosition_new==1) && (Spin_0_new < 0)){//downstream Pdown
	    hnew[3]->Fill(PhiS_new);
	    hnew_simple[3]->Fill(PhiS_simple_new);
	    hold[3]->Fill(PhiS_old);
	    hold_simple[3]->Fill(PhiS_simple_old);

	    hnew_xF[3]->Fill(x_feynman_new);
	    hnew_xPi[3]->Fill(x_beam_new);
	    hold_xF[3]->Fill(x_feynman_old);
	    hold_xPi[3]->Fill(x_beam_old);
	  }
	  break;
	}
      }

    }//end old tree loop
  }//end new tree loop

  //Draw distributions
  TCanvas* cnew = new TCanvas(); cnew->Divide(2,2);
  TCanvas* cnew_simple = new TCanvas(); cnew_simple->Divide(2,2);
  TCanvas* cold = new TCanvas(); cold->Divide(2,2);
  TCanvas* cold_simple = new TCanvas(); cold_simple->Divide(2,2);

  TCanvas* cxF = new TCanvas(); cxF->Divide(2,2);
  TCanvas* cxPi = new TCanvas(); cxPi->Divide(2,2);
  for (Int_t i=0; i<4; i++) {
    cnew->cd(i+1);
    hnew[i]->Draw();

    cnew_simple->cd(i+1);
    hnew_simple[i]->Draw();

    cold->cd(i+1);
    hold[i]->Draw();

    cold_simple->cd(i+1);
    hold_simple[i]->Draw();

    cxF->cd(i+1);
    hnew_xF[i]->Draw("E1");
    hold_xF[i]->Draw("E1 sames"); hold_xF[i]->SetLineColor(kRed);

    cxPi->cd(i+1);
    hnew_xPi[i]->Draw("E1");
    hold_xPi[i]->Draw("E1 sames"); hold_xPi[i]->SetLineColor(kRed);
  }
}


Double_t ShiftAngle(Double_t angle, Double_t phase=TMath::Pi()){
  if (angle > 0)
    return angle - phase;
  else
    return angle + phase;
}
