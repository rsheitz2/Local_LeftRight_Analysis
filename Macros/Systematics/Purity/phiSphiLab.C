#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Double_t DeterminePhi(Double_t vX, Double_t vY);

Double_t  ShiftPhiSimple (Double_t PhiS_simple);

void phiSphiLab(){

  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/";
  TString fname="slot1WAll_LowM_AMDY.root";
  TFile *f = OpenFile(pathData+fname);
  TTree *tree = (TTree*)f->Get("pT_Weighted");

  Double_t PhiS, PhiS_simple, Spin_0, Polarization, dilutionFactor, Mmumu;
  Double_t vPhoton_X, vPhoton_Y;
  Int_t targetPosition;
  tree->SetBranchAddress("PhiS", &PhiS);
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("Spin_0", &Spin_0);
  tree->SetBranchAddress("Polarization", &Polarization);
  tree->SetBranchAddress("dilutionFactor", &dilutionFactor);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  tree->SetBranchAddress("Mmumu", &Mmumu);
  tree->SetBranchAddress("vPhoton_X", &vPhoton_X);
  tree->SetBranchAddress("vPhoton_Y", &vPhoton_Y);

  TH2D* hphiSphiLab = new TH2D("hphiSphiLab", "hphiSphiLab",
			       100, -3.2, 3.2, 100, -3.2, 3.2);
  SetUp(hphiSphiLab);

  for (Int_t ev=0; ev<tree->GetEntries(); ev++) {
    //for (Int_t ev=0; ev<1000; ev++) { if (ev==0) cout << "debug" << endl;
    tree->GetEntry(ev);

    if (Mmumu > 3.44) continue;
    if (Mmumu < 2.80) continue;

    Double_t phiPhoton = DeterminePhi(vPhoton_X, vPhoton_Y);
    Double_t phiLab = ShiftPhiSimple(PhiS_simple);
    hphiSphiLab->Fill(phiLab, phiPhoton);

    /*if (targetPosition==0){
      if (Spin_0 > 0 ){
	hPhiS_upS_up->Fill(PhiS_simple);

	if (PhiS_simple > 0) left_upS_up++;
	else right_upS_up++;
      }
      else{
	hPhiS_upS_down->Fill(PhiS_simple);

	if (PhiS_simple < 0) left_upS_down++;
	else right_upS_down++;
      }
    }
    else{
      if (Spin_0 > 0 ){
	hPhiS_downS_up->Fill(PhiS_simple);

	if (PhiS_simple > 0) left_downS_up++;
	else right_downS_up++;
      }
      else{
	hPhiS_downS_down->Fill(PhiS_simple);

	if (PhiS_simple < 0) left_downS_down++;
	else right_downS_down++;
      }
    }//*/

  }//event loop

  TCanvas* c1 = new TCanvas();
  hphiSphiLab->Draw("colz");
  DrawLine(hphiSphiLab, -TMath::Pi()/2.0);
  DrawLine(hphiSphiLab, TMath::Pi()/2.0);

  DrawLine(hphiSphiLab, -TMath::Pi()/2.0, "x");
  DrawLine(hphiSphiLab, TMath::Pi()/2.0, "x");
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
