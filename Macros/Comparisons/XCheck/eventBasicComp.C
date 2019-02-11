#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void eventBasicComp(TString start=""){
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

  //Get Data/setup trees
  TString myData =Form("%s/%s/%sW07_%s.root", pathData.Data(), Mtype.Data(),
		       production.Data(), Mtype.Data());
  TFile *fMe = OpenFile(myData);
  TTree *myT = (TTree*)fMe->Get("pT_Weighted");

  TString xCheckData =
    Form("%s/Data/BasicDist/W07_An_xcheck.root", localPath.Data() );
  TFile *fxCheck = OpenFile(xCheckData);
  TTree *xCheckT = (TTree*)fxCheck->Get("An");

  Double_t x_beam_me, q_transverse_me, Spin_0, PhiS_me, PhiS_simple_me;
  Int_t targetPosition;
  myT->SetBranchAddress("x_beam", &x_beam_me);
  myT->SetBranchAddress("q_transverse", &q_transverse_me);
  myT->SetBranchAddress("Spin_0", &Spin_0);
  myT->SetBranchAddress("PhiS", &PhiS_me);
  myT->SetBranchAddress("PhiS_simple", &PhiS_simple_me);
  myT->SetBranchAddress("targetPosition", &targetPosition);

  Float_t xpi_xc, qt_xc, PhiS_xc;
  xCheckT->SetBranchAddress("xpi", &xpi_xc);
  xCheckT->SetBranchAddress("qt", &qt_xc);
  xCheckT->SetBranchAddress("PhiS", &PhiS_xc);

  //Distributions
  TString whichTarget[] = {"upS_pUp", "upS_pDown", "downS_pUp", "downS_pDown"};
  TH1D *hMe[4], *hMe_simple[4], *hxCheck[4];
  for (Int_t i=0; i<4; i++) {
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

  TH1D* hDiff = new TH1D("hDiff", "hDiff", nhbin, -2*TMath::Pi(), 2*TMath::Pi());
  
  //Tree loop
  Double_t epsilon = 0.000001;
  for (Int_t ev=0; ev<myT->GetEntries(); ev++) {
    //for (Int_t ev=0; ev<3; ev++) {
    myT->GetEntry(ev);
    
    for (Int_t j=0; j<xCheckT->GetEntries(); j++) {
      xCheckT->GetEntry(j);

      if ((x_beam_me >= xpi_xc-epsilon) && (x_beam_me <= xpi_xc+epsilon)){
	if ((q_transverse_me >= qt_xc-epsilon) &&
	    (q_transverse_me <= qt_xc+epsilon)){

	  if ((targetPosition==0) && (Spin_0 > 0)){//upstream Pup
	    hMe[0]->Fill(PhiS_me);
	    hMe_simple[0]->Fill(PhiS_simple_me);
	    hxCheck[0]->Fill(PhiS_xc);

	    //Double_t Diff = PhiS_me - PhiS_xc;
	    //cout << PhiS_me << " " << PhiS_simple_me << " " << PhiS_xc << " " << Diff << endl;
	  }
	  else if ((targetPosition==0) && (Spin_0 < 0)){//upstream Pdown
	    hMe[1]->Fill(PhiS_me);
	    hMe_simple[1]->Fill(PhiS_simple_me);
	    hxCheck[1]->Fill(PhiS_xc);

	    //Double_t Diff = PhiS_me - PhiS_xc;
	    Double_t Diff = PhiS_simple_me - PhiS_xc;
	    if (Diff < -6.1){
	      cout << PhiS_me << " " << PhiS_simple_me << " " << PhiS_xc
		   << " " << Diff << endl;  
	    }
	    hDiff->Fill(Diff);
	  }
	  else if ((targetPosition==1) && (Spin_0 > 0)){//downstream Pup
	    hMe[2]->Fill(PhiS_me);
	    hMe_simple[2]->Fill(PhiS_simple_me);
	    hxCheck[2]->Fill(PhiS_xc);
	  }
	  else if ((targetPosition==1) && (Spin_0 < 0)){//downstream Pdown
	    hMe[3]->Fill(PhiS_me);
	    hMe_simple[3]->Fill(PhiS_simple_me);
	    hxCheck[3]->Fill(PhiS_xc);
	  }
	  break;
	}
      }

    }//end xCheck tree loop

  }//end Me tree loop


  //Draw distributions
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
  }

  TCanvas* cDiff = new TCanvas();
  hDiff->Draw();
}
