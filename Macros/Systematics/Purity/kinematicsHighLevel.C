#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Double_t DeterminePhi(Double_t vX, Double_t vY);

Double_t  ShiftPhiSimple (Double_t PhiS_simple);

void kinematicsHighLevel(){
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/";

  Double_t PhiS, PhiS_simple, Spin_0, Polarization, dilutionFactor, Mmumu;
  Double_t vPhoton_X, vPhoton_Y;
  Double_t vOpenAngle, x_beam, x_target, x_feynman, q_transverse;
  Double_t rapidity;
  Int_t targetPosition;

  TString fname ="slot1WAll_LowM_AMDY.root";
  TFile *f = OpenFile(pathData+fname);
  TTree *tree = (TTree*)f->Get("pT_Weighted");

  tree->SetBranchAddress("PhiS", &PhiS);
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("Spin_0", &Spin_0);
  tree->SetBranchAddress("Polarization", &Polarization);
  tree->SetBranchAddress("dilutionFactor", &dilutionFactor);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  tree->SetBranchAddress("Mmumu", &Mmumu);
  tree->SetBranchAddress("vPhoton_X", &vPhoton_X);
  tree->SetBranchAddress("vPhoton_Y", &vPhoton_Y);
  tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  tree->SetBranchAddress("x_beam", &x_beam);
  tree->SetBranchAddress("x_target", &x_target);
  tree->SetBranchAddress("x_feynman", &x_feynman);
  tree->SetBranchAddress("q_transverse", &q_transverse);
  tree->SetBranchAddress("rapidity", &rapidity);

  const Int_t nTarg =4;
  TString targNames[nTarg] = {"ups_up", "ups_down", "downs_up", "downs_down"};
  TH1D *h_vOpenAngle[nTarg], *h_xBeam[nTarg], *h_xTarg[nTarg], *h_xFey[nTarg];
  TH1D *h_qT[nTarg], *h_rapidity[nTarg];
  for (Int_t i=0; i<nTarg; i++) {
    h_vOpenAngle[i] = new TH1D(Form("vOpenAngle_%s", targNames[i].Data()),
			       Form("vOpenAngle_%s", targNames[i].Data()),
			       100, 0, 0.32);

    h_xBeam[i] = new TH1D(Form("xBeam_%s", targNames[i].Data()),
			  Form("xBeam_%s", targNames[i].Data()),
			  100, 0, 1);

    h_xTarg[i] = new TH1D(Form("xTarg_%s", targNames[i].Data()),
			  Form("xTarg_%s", targNames[i].Data()),
			  100, 0, 1);

    h_xFey[i] = new TH1D(Form("xFey_%s", targNames[i].Data()),
			 Form("xFey_%s", targNames[i].Data()),
			 100, -1, 1);

    h_qT[i] = new TH1D(Form("qT_%s", targNames[i].Data()),
		       Form("qT_%s", targNames[i].Data()),
		       100, 0, 6);

    h_rapidity[i] = new TH1D(Form("rapidity_%s", targNames[i].Data()),
			     Form("rapidity_%s", targNames[i].Data()),
			     100, -3, 3);
  }


  for (Int_t ev=0; ev<tree->GetEntries(); ev++) {
    //for (Int_t ev=0; ev<1000; ev++) { if (ev==0) cout << "debug" << endl;
    tree->GetEntry(ev);

    /*if (Mmumu > 8.5) continue;
    if (Mmumu < 4.3) continue;//*/
    
    if (Mmumu > 3.44) continue;
    if (Mmumu < 2.80) continue;//*/

    /*//forced acceptance cut
    Double_t phiPhoton = DeterminePhi(vPhoton_X, vPhoton_Y);
    Double_t phiLab = ShiftPhiSimple(PhiS_simple);
    Bool_t labLeft =false, labRight =false;
    if ( (phiLab > -TMath::Pi()/2) && (phiLab < TMath::Pi()/2) )
      labLeft = true;
    else labRight = true;
    Bool_t tfLeft =false, tfRight =false;
    if ( (phiPhoton > -TMath::Pi()/2) && (phiPhoton < TMath::Pi()/2) )
      tfLeft = true;
    else tfRight = true;

    if (labLeft != tfLeft) continue;
    if (labRight != tfRight) continue;//*/

    if (targetPosition==0){
      if (Spin_0 > 0 ){
	h_vOpenAngle[0]->Fill(vOpenAngle);
	h_xBeam[0]->Fill(x_beam);
	h_xTarg[0]->Fill(x_target);
	h_xFey[0]->Fill(x_feynman);
	h_qT[0]->Fill(q_transverse);
	h_rapidity[0]->Fill(rapidity);
	
      }
      else{
	h_vOpenAngle[1]->Fill(vOpenAngle);
	h_xBeam[1]->Fill(x_beam);
	h_xTarg[1]->Fill(x_target);
	h_xFey[1]->Fill(x_feynman);
	h_qT[1]->Fill(q_transverse);
	h_rapidity[1]->Fill(rapidity);
      }
    }
    else{
      if (Spin_0 > 0 ){
	h_vOpenAngle[2]->Fill(vOpenAngle);
	h_xBeam[2]->Fill(x_beam);
	h_xTarg[2]->Fill(x_target);
	h_xFey[2]->Fill(x_feynman);
	h_qT[2]->Fill(q_transverse);
	h_rapidity[2]->Fill(rapidity);
	
      }
      else{
	h_vOpenAngle[3]->Fill(vOpenAngle);
	h_xBeam[3]->Fill(x_beam);
	h_xTarg[3]->Fill(x_target);
	h_xFey[3]->Fill(x_feynman);
	h_qT[3]->Fill(q_transverse);
	h_rapidity[3]->Fill(rapidity);
	
      }
    }
  }//event loop


  TCanvas* cOA = new TCanvas(); cOA->Divide(2);
  for (Int_t i=0; i<nTarg; i++) {
    cOA->cd(1);
    if (i==0) h_vOpenAngle[i]->Draw();
    else h_vOpenAngle[i]->Draw("sames");

    h_vOpenAngle[i]->SetLineColor(i+1);

    cOA->cd(2);
    if (i==0) h_rapidity[i]->Draw();
    else h_rapidity[i]->Draw("sames");

    h_rapidity[i]->SetLineColor(i+1);
  }

  TCanvas* cKin = new TCanvas(); cKin->Divide(2, 2);
  for (Int_t i=0; i<nTarg; i++) {
    cKin->cd(1);
    if (i==0) h_xBeam[i]->Draw();
    else h_xBeam[i]->Draw("sames");

    h_xBeam[i]->SetLineColor(i+1);

    cKin->cd(2);
    if (i==0) h_xTarg[i]->Draw();
    else h_xTarg[i]->Draw("sames");

    h_xTarg[i]->SetLineColor(i+1);

    cKin->cd(3);
    if (i==0) h_xFey[i]->Draw();
    else h_xFey[i]->Draw("sames");

    h_xFey[i]->SetLineColor(i+1);

    cKin->cd(4);
    if (i==0) h_qT[i]->Draw();
    else h_qT[i]->Draw("sames");

    h_qT[i]->SetLineColor(i+1);
  }
  
}


Double_t DeterminePhi(Double_t vX, Double_t vY){

  Double_t phi;
  if ( (vX > 0) && (vY > 0) ) phi = TMath::ATan(vY/vX);
  else if ( (vX > 0) && (vY < 0) ){
    phi = TMath::ATan(-1.0*vY/vX);
    phi *= -1.0;
  }
  else if ( (vX < 0) && (vY > 0) ){
    phi = TMath::ATan(-1.0*vY/vX);
    phi = TMath::Pi() - phi;
  }
  else{
    phi = TMath::ATan(vY/vX);
    phi = phi - TMath::Pi();
  }

  return phi;
}


Double_t  ShiftPhiSimple (Double_t PhiS_simple) {
  //Lab reference frame:
  //    -Pi/2 ->   P/i/2 = left
  //    Pi/2  -> 3*P/i/2 = right
  Double_t phi = TMath::Pi()/2 - PhiS_simple;
  if (phi > TMath::Pi()) phi = phi - 2*TMath::Pi();
  
  return phi;
}//ShiftPhiSimple
