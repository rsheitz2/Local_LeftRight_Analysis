#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Double_t DeterminePhi(Double_t vX, Double_t vY);

Double_t  ShiftPhiSimple (Double_t PhiS_simple);

void kinematicsLowLevel(){
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/";
  
  Double_t PhiS, PhiS_simple, Spin_0, Polarization, dilutionFactor, Mmumu;
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E, vOpenAngle;
  Int_t MasterTrigMask;
  Double_t theta_traj1, phi_traj1, qP_traj1;
  Double_t theta_traj2, phi_traj2, qP_traj2;
  Double_t theta_trajPIn, phi_trajPIn, qP_trajPIn;
  Double_t vx_z, vx_x, vx_y;
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
  tree->SetBranchAddress("vPhoton_Z", &vPhoton_Z);
  tree->SetBranchAddress("vPhoton_E", &vPhoton_E);
  tree->SetBranchAddress("MasterTrigMask", &MasterTrigMask);
  tree->SetBranchAddress("theta_traj1", &theta_traj1);
  tree->SetBranchAddress("phi_traj1", &phi_traj1);
  tree->SetBranchAddress("qP_traj1", &qP_traj1);
  tree->SetBranchAddress("theta_traj2", &theta_traj2);
  tree->SetBranchAddress("phi_traj2", &phi_traj2);
  tree->SetBranchAddress("qP_traj2", &qP_traj2);
  tree->SetBranchAddress("theta_trajPIn", &theta_trajPIn);
  tree->SetBranchAddress("phi_trajPIn", &phi_trajPIn);
  tree->SetBranchAddress("qP_trajPIn", &qP_trajPIn);
  tree->SetBranchAddress("vx_z", &vx_z);
  tree->SetBranchAddress("vx_y", &vx_y);
  tree->SetBranchAddress("vx_x", &vx_x);
    
  const Int_t nTarg =4;
  TString targNames[nTarg] = {"ups_up", "ups_down", "downs_up", "downs_down"};
  TH1D *h_MasterTrig[nTarg];
  TH1D *h_theta_traj1[nTarg], *h_phi_traj1[nTarg], *h_qP_traj1[nTarg];
  TH1D *h_theta_traj2[nTarg], *h_phi_traj2[nTarg], *h_qP_traj2[nTarg];
  TH1D *h_theta_trajPIn[nTarg], *h_phi_trajPIn[nTarg], *h_qP_trajPIn[nTarg];
  TH1D *h_vPhotonZ[nTarg], *h_vPhotonE[nTarg];
  TH1D *h_vx_z[nTarg], *h_vx_y[nTarg], *h_vx_x[nTarg];
  TH2D *h_vx_r[nTarg], *h_qP_theta2[nTarg];
  TH2D *h_theta1_phi1[nTarg], *h_theta2_phi2[nTarg];
  for (Int_t i=0; i<nTarg; i++) {
    h_MasterTrig[i] = new TH1D(Form("MasterTrig_%s", targNames[i].Data()),
			       Form("MasterTrig_%s", targNames[i].Data()),
			       300, -10, 1030);

    h_theta_traj1[i] = new TH1D(Form("theta_traj1_%s", targNames[i].Data()),
				Form("theta_traj1_%s", targNames[i].Data()),
				100, 0, 0.5);
    h_phi_traj1[i] = new TH1D(Form("phi_traj1_%s", targNames[i].Data()),
			      Form("phi_traj1_%s", targNames[i].Data()),
			      100, -3.2, 3.2);
    h_qP_traj1[i] = new TH1D(Form("qP_traj1_%s", targNames[i].Data()),
			     Form("qP_traj1_%s", targNames[i].Data()),
			     100, 0, 200);

    h_theta_traj2[i] = new TH1D(Form("theta_traj2_%s", targNames[i].Data()),
				Form("theta_traj2_%s", targNames[i].Data()),
				100, 0, 0.5);
    h_phi_traj2[i] = new TH1D(Form("phi_traj2_%s", targNames[i].Data()),
			      Form("phi_traj2_%s", targNames[i].Data()),
			      100, -3.2, 3.2);
    h_qP_traj2[i] = new TH1D(Form("qP_traj2_%s", targNames[i].Data()),
			     Form("qP_traj2_%s", targNames[i].Data()),
			     100, -200, 0);

    h_theta_trajPIn[i] = new TH1D(Form("theta_trajPIn_%s", targNames[i].Data()),
				  Form("theta_trajPIn_%s", targNames[i].Data()),
				  20, -0.1, 0.009);
    h_phi_trajPIn[i] = new TH1D(Form("phi_trajPIn_%s", targNames[i].Data()),
				Form("phi_trajPIn_%s", targNames[i].Data()),
				100, -3.2, 3.2);
    h_qP_trajPIn[i] = new TH1D(Form("qP_trajPIn_%s", targNames[i].Data()),
			       Form("qP_trajPIn_%s", targNames[i].Data()),
			       10, -195, -185);

    h_vPhotonZ[i] = new TH1D(Form("vPhotonZ_%s", targNames[i].Data()),
			     Form("vPhotonZ_%s", targNames[i].Data()),
			     100, 0, 190);

    h_vPhotonE[i] = new TH1D(Form("vPhotonE_%s", targNames[i].Data()),
			     Form("vPhotonE_%s", targNames[i].Data()),
			     100, 0, 190);

    h_vx_z[i] = new TH1D(Form("vx_z_%s", targNames[i].Data()),
			 Form("vx_z_%s", targNames[i].Data()),
			 500, -350, -150);
    h_vx_x[i] = new TH1D(Form("vx_x_%s", targNames[i].Data()),
			 Form("vx_x_%s", targNames[i].Data()),
			 100, -2, 2);
    h_vx_y[i] = new TH1D(Form("vx_y_%s", targNames[i].Data()),
			 Form("vx_y_%s", targNames[i].Data()),
			 100, -2, 2);
    
    h_vx_r[i] = new TH2D(Form("vx_r_%s", targNames[i].Data()),
			 Form("vx_r_%s", targNames[i].Data()),
			 100, -2, 2, 100, -2, 2);
    h_qP_theta2[i] = new TH2D(Form("qP_theta2_%s", targNames[i].Data()),
			      Form("qP_theta2_%s", targNames[i].Data()),
			      100, -200, 0, 100, 0, 0.35);
    h_theta1_phi1[i] = new TH2D(Form("theta1_phi1_%s", targNames[i].Data()),
				Form("theta1_phi1_%s", targNames[i].Data()),
				100, 0, 0.5, 100, -3.2, 3.2);
    h_theta2_phi2[i] = new TH2D(Form("theta2_phi2_%s", targNames[i].Data()),
				Form("theta2_phi2_%s", targNames[i].Data()),
				100, 0, 0.5, 100, -3.2, 3.2);
    
  }


  for (Int_t ev=0; ev<tree->GetEntries(); ev++) {
  //for (Int_t ev=0; ev<1000; ev++) { if (ev==0) cout << "debug" << endl;
    tree->GetEntry(ev);

    ////////////
    //Cuts
    //Mass cut
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

    //open angle cut
    /*if (vOpenAngle > 0.17) continue;
    if (vOpenAngle < 0.05) continue;//*/

    //Target radius cut
    Double_t vx_rad2 = vx_x*vx_x + vx_y*vx_y;
    if (vx_rad2 > 0.35*0.35) continue;//*/

    if (targetPosition==0){
      if (Spin_0 > 0 ){
	h_MasterTrig[0]->Fill(MasterTrigMask);
	h_theta_traj1[0]->Fill(theta_traj1);
	h_phi_traj1[0]->Fill(phi_traj1);
	h_qP_traj1[0]->Fill(qP_traj1);
	h_theta_traj2[0]->Fill(theta_traj2);
	h_phi_traj2[0]->Fill(phi_traj2);
	h_qP_traj2[0]->Fill(qP_traj2);
	h_theta_trajPIn[0]->Fill(theta_trajPIn);
	h_phi_trajPIn[0]->Fill(phi_trajPIn);
	h_qP_trajPIn[0]->Fill(qP_trajPIn);
	h_vx_z[0]->Fill(vx_z);
	h_vx_y[0]->Fill(vx_y);
	h_vx_x[0]->Fill(vx_x);
	h_vPhotonZ[0]->Fill(vPhoton_Z);
	h_vPhotonE[0]->Fill(vPhoton_E);

	h_vx_r[0]->Fill(vx_x, vx_y);
	h_qP_theta2[0]->Fill(qP_traj2, theta_traj2);
	h_theta1_phi1[0]->Fill(theta_traj1, phi_traj1);
	h_theta2_phi2[0]->Fill(theta_traj2, phi_traj2);
		
      }
      else{
	h_MasterTrig[1]->Fill(MasterTrigMask);
	h_theta_traj1[1]->Fill(theta_traj1);
	h_phi_traj1[1]->Fill(phi_traj1);
	h_qP_traj1[1]->Fill(qP_traj1);
	h_theta_traj2[1]->Fill(theta_traj2);
	h_phi_traj2[1]->Fill(phi_traj2);
	h_qP_traj2[1]->Fill(qP_traj2);
	h_theta_trajPIn[1]->Fill(theta_trajPIn);
	h_phi_trajPIn[1]->Fill(phi_trajPIn);
	h_qP_trajPIn[1]->Fill(qP_trajPIn);
	h_vx_z[1]->Fill(vx_z);
	h_vx_y[1]->Fill(vx_y);
	h_vx_x[1]->Fill(vx_x);
	h_vPhotonZ[1]->Fill(vPhoton_Z);
	h_vPhotonE[1]->Fill(vPhoton_E);

	h_vx_r[1]->Fill(vx_x, vx_y);
	h_qP_theta2[1]->Fill(qP_traj2, theta_traj2);
	h_theta1_phi1[1]->Fill(theta_traj1, phi_traj1);
	h_theta2_phi2[1]->Fill(theta_traj2, phi_traj2);
      }
    }
    else{
      if (Spin_0 > 0 ){
	h_MasterTrig[2]->Fill(MasterTrigMask);
	h_theta_traj1[2]->Fill(theta_traj1);
	h_phi_traj1[2]->Fill(phi_traj1);
	h_qP_traj1[2]->Fill(qP_traj1);
	h_theta_traj2[2]->Fill(theta_traj2);
	h_phi_traj2[2]->Fill(phi_traj2);
	h_qP_traj2[2]->Fill(qP_traj2);
	h_theta_trajPIn[2]->Fill(theta_trajPIn);
	h_phi_trajPIn[2]->Fill(phi_trajPIn);
	h_qP_trajPIn[2]->Fill(qP_trajPIn);
	h_vx_z[2]->Fill(vx_z);
	h_vx_y[2]->Fill(vx_y);
	h_vx_x[2]->Fill(vx_x);
	h_vPhotonZ[2]->Fill(vPhoton_Z);
	h_vPhotonE[2]->Fill(vPhoton_E);

	h_vx_r[2]->Fill(vx_x, vx_y);
	h_qP_theta2[2]->Fill(qP_traj2, theta_traj2);
	h_theta1_phi1[2]->Fill(theta_traj1, phi_traj1);
	h_theta2_phi2[2]->Fill(theta_traj2, phi_traj2);
	
      }
      else{
	h_MasterTrig[3]->Fill(MasterTrigMask);
	h_theta_traj1[3]->Fill(theta_traj1);
	h_phi_traj1[3]->Fill(phi_traj1);
	h_qP_traj1[3]->Fill(qP_traj1);
	h_theta_traj2[3]->Fill(theta_traj2);
	h_phi_traj2[3]->Fill(phi_traj2);
	h_qP_traj2[3]->Fill(qP_traj2);
	h_theta_trajPIn[3]->Fill(theta_trajPIn);
	h_phi_trajPIn[3]->Fill(phi_trajPIn);
	h_qP_trajPIn[3]->Fill(qP_trajPIn);
	h_vx_z[3]->Fill(vx_z);
	h_vx_y[3]->Fill(vx_y);
	h_vx_x[3]->Fill(vx_x);
	h_vPhotonZ[3]->Fill(vPhoton_Z);
	h_vPhotonE[3]->Fill(vPhoton_E);

	h_vx_r[3]->Fill(vx_x, vx_y);
	h_qP_theta2[3]->Fill(qP_traj2, theta_traj2);
	h_theta1_phi1[3]->Fill(theta_traj1, phi_traj1);
	h_theta2_phi2[3]->Fill(theta_traj2, phi_traj2);
	
      }
    }
  }//event loop


  TCanvas* cTrig = new TCanvas(); cTrig->Divide(2,2);
  for (Int_t i=0; i<nTarg; i++) {
    cTrig->cd(i+1);
    h_MasterTrig[i]->Draw();
    h_MasterTrig[i]->SetLineColor(i+1);
  }

  TCanvas* c1 = new TCanvas(); c1->Divide(3);
  for (Int_t i=0; i<nTarg; i++) {
    c1->cd(1);
    if (i==0) h_theta_traj1[i]->Draw();
    else h_theta_traj1[i]->Draw("sames");

    h_theta_traj1[i]->SetLineColor(i+1);

    c1->cd(2);
    if (i==0) h_phi_traj1[i]->Draw();
    else h_phi_traj1[i]->Draw("sames");

    h_phi_traj1[i]->SetLineColor(i+1);

    c1->cd(3);
    if (i==0) h_qP_traj1[i]->Draw();
    else h_qP_traj1[i]->Draw("sames");

    h_qP_traj1[i]->SetLineColor(i+1);
  }

  TCanvas* c2 = new TCanvas(); c2->Divide(3);
  for (Int_t i=0; i<nTarg; i++) {
    c2->cd(1);
    if (i==0) h_theta_traj2[i]->Draw();
    else h_theta_traj2[i]->Draw("sames");

    h_theta_traj2[i]->SetLineColor(i+1);

    c2->cd(2);
    if (i==0) h_phi_traj2[i]->Draw();
    else h_phi_traj2[i]->Draw("sames");

    h_phi_traj2[i]->SetLineColor(i+1);

    c2->cd(3);
    if (i==0) h_qP_traj2[i]->Draw();
    else h_qP_traj2[i]->Draw("sames");

    h_qP_traj2[i]->SetLineColor(i+1);
  }

  TCanvas* cPIn = new TCanvas(); cPIn->Divide(3);
  for (Int_t i=0; i<nTarg; i++) {
    cPIn->cd(1);
    if (i==0) h_theta_trajPIn[i]->Draw();
    else h_theta_trajPIn[i]->Draw("sames");

    h_theta_trajPIn[i]->SetLineColor(i+1);

    cPIn->cd(2);
    if (i==0) h_phi_trajPIn[i]->Draw();
    else h_phi_trajPIn[i]->Draw("sames");

    h_phi_trajPIn[i]->SetLineColor(i+1);

    cPIn->cd(3);
    if (i==0) h_qP_trajPIn[i]->Draw();
    else h_qP_trajPIn[i]->Draw("sames");

    h_qP_trajPIn[i]->SetLineColor(i+1);
  }

  TCanvas* cvx = new TCanvas(); cvx->Divide(3);
  for (Int_t i=0; i<nTarg; i++) {
    cvx->cd(1);
    if (i==0) h_vx_z[i]->Draw();
    else h_vx_z[i]->Draw("sames");

    h_vx_z[i]->SetLineColor(i+1);

    cvx->cd(2);
    if (i==0) h_vx_y[i]->Draw();
    else h_vx_y[i]->Draw("sames");

    h_vx_y[i]->SetLineColor(i+1);

    cvx->cd(3);
    if (i==0) h_vx_x[i]->Draw();
    else h_vx_x[i]->Draw("sames");

    h_vx_x[i]->SetLineColor(i+1);
  }

  TCanvas* cvPhoton = new TCanvas(); cvPhoton->Divide(2);
    for (Int_t i=0; i<nTarg; i++) {
    cvPhoton->cd(1);
    if (i==0) h_vPhotonZ[i]->Draw();
    else h_vPhotonZ[i]->Draw("sames");

    h_vPhotonZ[i]->SetLineColor(i+1);

    cvPhoton->cd(2);
    if (i==0) h_vPhotonE[i]->Draw();
    else h_vPhotonE[i]->Draw("sames");

    h_vPhotonE[i]->SetLineColor(i+1);
  }

  TCanvas* cRad = new TCanvas(); cRad->Divide(2,2);
  for (Int_t i=0; i<nTarg; i++) {
    cRad->cd(i+1);
    h_vx_r[i]->Draw("colz");
  }

  TCanvas* cqP2 = new TCanvas(); cqP2->Divide(2,2);
  for (Int_t i=0; i<nTarg; i++) {
    cqP2->cd(i+1);
    h_qP_theta2[i]->Draw("colz");
  }

  TCanvas* cthetaphi1 = new TCanvas(); cthetaphi1->Divide(2,2);
  for (Int_t i=0; i<nTarg; i++) {
    cthetaphi1->cd(i+1);
    h_theta1_phi1[i]->Draw("colz");
  }

  TCanvas* cthetaphi2 = new TCanvas(); cthetaphi2->Divide(2,2);
  for (Int_t i=0; i<nTarg; i++) {
    cthetaphi2->cd(i+1);
    h_theta2_phi2[i]->Draw("colz");
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
