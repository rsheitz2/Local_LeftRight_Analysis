#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

Double_t ShiftAngle(Double_t angle, Double_t phast=TMath::Pi());

void eventBasicComp(TString start=""){
  //Setup_______________
  TString Mtype ="HMDY";
  TString production ="slot1"; //"slot1", "t3"

  Bool_t toWrite =false;
  //Setup_______________

  //Basic Setup
  TString localPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/XCheck";
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData";
  Int_t nhbin =150;

  //Get Data/setup trees
  if (production == "t3") production = "";
  TString myData =Form("%s/%s/%sW07_%s.root", pathData.Data(), Mtype.Data(),
		       production.Data(), Mtype.Data());
  TFile *fMe = OpenFile(myData);
  TTree *myT = (TTree*)fMe->Get("pT_Weighted");
  if (production == "") production = "t3";

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
  TH1D *hMe_xPi[4], *hMe_qT[4];
  TH1D *hXC_xPi[4], *hXC_qT[4];
  TH1D *hDiff[4];
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

    hMe_xPi[i] = new TH1D (Form("hMe_xPi_%s", whichTarget[i].Data()),
			   Form("hMe_xPi_%s", whichTarget[i].Data()),
			   nhbin, 0, 1);
    hMe_qT[i] = new TH1D (Form("hMe_qT_%s", whichTarget[i].Data()),
			      Form("hMe_qT_%s", whichTarget[i].Data()),
			      nhbin, 0, 6);
    hXC_xPi[i] = new TH1D (Form("hXC_xPi_%s", whichTarget[i].Data()),
			   Form("hXC_xPi_%s", whichTarget[i].Data()),
			   nhbin, 0, 1);
    hXC_qT[i] = new TH1D (Form("hXC_qT_%s", whichTarget[i].Data()),
			      Form("hXC_qT_%s", whichTarget[i].Data()),
			      nhbin, 0, 6);

    hDiff[i] = new TH1D (Form("hDiff_%s", whichTarget[i].Data()),
			 Form("hDiff_%s", whichTarget[i].Data()),
			 nhbin, -2*TMath::Pi(), 2*TMath::Pi());
  }

  //Tree loop
  Double_t epsilon = 0.000001;
  Double_t minDiff = 10;
  for (Int_t ev=0; ev<myT->GetEntries(); ev++) {
    //for (Int_t ev=0; ev<3; ev++) {
    myT->GetEntry(ev);
    
    for (Int_t j=0; j<xCheckT->GetEntries(); j++) {
      xCheckT->GetEntry(j);

      if ((x_beam_me >= xpi_xc-epsilon) && (x_beam_me <= xpi_xc+epsilon)){
	if ((q_transverse_me >= qt_xc-epsilon) &&
	    (q_transverse_me <= qt_xc+epsilon)){

	  Double_t Diff;
	  if ((targetPosition==0) && (Spin_0 > 0)){//upstream Pup
	    //PhiS_me = ShiftAngle(PhiS_me);
	    PhiS_simple_me = ShiftAngle(PhiS_simple_me);
	    
	    hMe[0]->Fill(PhiS_me);
	    hMe_simple[0]->Fill(PhiS_simple_me);
	    hxCheck[0]->Fill(PhiS_xc);

	    hMe_xPi[0]->Fill(x_beam_me);
	    hMe_qT[0]->Fill(q_transverse_me);
	    hXC_xPi[0]->Fill(xpi_xc);
	    hXC_qT[0]->Fill(qt_xc);

	    Diff = PhiS_xc - PhiS_simple_me;
	    hDiff[0]->Fill(Diff);
	  }
	  else if ((targetPosition==0) && (Spin_0 < 0)){//upstream Pdown
	    PhiS_me = ShiftAngle(PhiS_me);
	    PhiS_simple_me = ShiftAngle(PhiS_simple_me);
	    
	    hMe[1]->Fill(PhiS_me);
	    hMe_simple[1]->Fill(PhiS_simple_me);
	    hxCheck[1]->Fill(PhiS_xc);

	    hMe_xPi[1]->Fill(x_beam_me);
	    hMe_qT[1]->Fill(q_transverse_me);
	    hXC_xPi[1]->Fill(xpi_xc);
	    hXC_qT[1]->Fill(qt_xc);

	    Diff = PhiS_xc - PhiS_simple_me;
	    hDiff[1]->Fill(Diff);
	  }
	  else if ((targetPosition==1) && (Spin_0 > 0)){//downstream Pup
	    //PhiS_me = ShiftAngle(PhiS_me);
	    PhiS_simple_me = ShiftAngle(PhiS_simple_me);
	    
	    hMe[2]->Fill(PhiS_me);
	    hMe_simple[2]->Fill(PhiS_simple_me);
	    hxCheck[2]->Fill(PhiS_xc);

	    hMe_xPi[2]->Fill(x_beam_me);
	    hMe_qT[2]->Fill(q_transverse_me);
	    hXC_xPi[2]->Fill(xpi_xc);
	    hXC_qT[2]->Fill(qt_xc);

	    Diff = PhiS_xc - PhiS_simple_me;
	    hDiff[2]->Fill(Diff);
	  }
	  else if ((targetPosition==1) && (Spin_0 < 0)){//downstream Pdown
	    PhiS_me = ShiftAngle(PhiS_me);
	    PhiS_simple_me = ShiftAngle(PhiS_simple_me);
	    
	    hMe[3]->Fill(PhiS_me);
	    hMe_simple[3]->Fill(PhiS_simple_me);
	    hxCheck[3]->Fill(PhiS_xc);

	    hMe_xPi[3]->Fill(x_beam_me);
	    hMe_qT[3]->Fill(q_transverse_me);
	    hXC_xPi[3]->Fill(xpi_xc);
	    hXC_qT[3]->Fill(qt_xc);

	    Diff = PhiS_xc - PhiS_simple_me;
	    hDiff[3]->Fill(Diff);
	  }

	  if (TMath::Abs(Diff) < minDiff) minDiff = TMath::Abs(Diff);
	  if ((Diff < epsilon) && (Diff > -epsilon)){
	    cout << PhiS_me << " " << PhiS_simple_me << " " << PhiS_xc << endl;
	  }
	  break;
	}
      }

    }//end xCheck tree loop

  }//end Me tree loop

  cout << "\n\nMinimum PhiS difference btw me and XC" << endl;
  cout << "   " << minDiff << "\n\n" << endl;

  //Draw distributions
  TCanvas* cMe = new TCanvas(); cMe->Divide(2,2);
  TCanvas* cMe_simple = new TCanvas(); cMe_simple->Divide(2,2);
  TCanvas* cXC = new TCanvas(); cXC->Divide(2,2);

  TCanvas* cxPi = new TCanvas(); cxPi->Divide(2,2);
  TCanvas* cqT = new TCanvas(); cqT->Divide(2,2);
  
  TCanvas* cDiff = new TCanvas(); cDiff->Divide(2,2);
  for (Int_t i=0; i<4; i++) {
    cMe->cd(i+1);
    hMe[i]->Draw();

    cMe_simple->cd(i+1);
    hMe_simple[i]->Draw();

    cXC->cd(i+1);
    hxCheck[i]->Draw();

    cxPi->cd(i+1);
    hMe_xPi[i]->Draw("E1");
    hXC_xPi[i]->Draw("E1 sames"); hXC_xPi[i]->SetLineColor(kRed);

    cqT->cd(i+1);
    hMe_qT[i]->Draw("E1");
    hXC_qT[i]->Draw("E1 sames"); hXC_qT[i]->SetLineColor(kRed);

    cDiff->cd(i+1);
    hDiff[i]->Draw();
  }
}


Double_t ShiftAngle(Double_t angle, Double_t phase=TMath::Pi()){
  if (angle > 0)
    return angle - phase;
  else
    return angle + phase;
}
