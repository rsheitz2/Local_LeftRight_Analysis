#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Double_t DeterminePhi(Double_t vX, Double_t vY);

Double_t  ShiftPhiSimple (Double_t PhiS_simple);

Double_t MakeAsym(Double_t l_pUp, Double_t r_pUp);

Double_t eAsym(Double_t l_pUp, Double_t r_pUp);

void mcLeftRight(){
  TString pathData="/Volumes/Seagate/DrellYan/Charles_Official/";
  TString whichMC ="LMDY";//"LMDY", "OC", "Psi", "Jpsi"
  TString fname =Form("Charles_W12_%s.root", whichMC.Data());
  TFile *f = OpenFile(pathData+fname);
  TTree *tree = (TTree*)f->Get("pT_Weighted");

  Double_t PhiS, PhiS_simple, Mmumu;
  Double_t vPhoton_X, vPhoton_Y;
  Int_t targetPosition;
  tree->SetBranchAddress("PhiS", &PhiS);
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  tree->SetBranchAddress("Mmumu", &Mmumu);
  tree->SetBranchAddress("vPhoton_X", &vPhoton_X);
  tree->SetBranchAddress("vPhoton_Y", &vPhoton_Y);

  TH1D* hPhiS = new TH1D("hPhiS", "hPhiS", 100, -3.2, 3.2);

  TH1D* hPhiS_upS_up = new TH1D("hPhiS_upS_up", "hPhiS_upS_up", 100, -3.2, 3.2);
  TH1D* hPhiS_downS_up = new TH1D("hPhiS_downS_up", "hPhiS_downS_up",
				  100, -3.2, 3.2);
  SetUp(hPhiS); SetUp(hPhiS_upS_up); SetUp(hPhiS_downS_up);

  Double_t left_upS_up =0., left_downS_up =0.;
  Double_t right_upS_up =0., right_downS_up =0.;

  TString whichCuts ="";
  for (Int_t ev=0; ev<tree->GetEntries(); ev++) {
  /*for (Int_t ev=0; ev<50000; ev++) {
    if (ev==0) {
      cout << "debug" << endl;
      whichCuts +="_debug";
    }//*/
    tree->GetEntry(ev);

    ////////////
    //Cuts
    //Mass cut
    if (Mmumu > 3.44) continue;
    if (Mmumu < 2.80) continue;
    if (ev==0) whichCuts+= "_2.8Mu3.44";//*/
    /*if (Mmumu > 8.5) continue;
    if (Mmumu < 4.3) continue;
    if (ev==0) whichCuts+= "_4.3Mu8.5";//*/

    //forced acceptance cut
    /*Double_t phiPhoton = DeterminePhi(vPhoton_X, vPhoton_Y);
    Double_t phiLab = ShiftPhiSimple(PhiS_simple);
    Bool_t labLeft =false, labRight =false;
    if ( (phiLab > -TMath::Pi()/2) && (phiLab < TMath::Pi()/2) ) labLeft = true;
    else labRight = true;
    Bool_t tfLeft =false, tfRight =false;
    if ( (phiPhoton > -TMath::Pi()/2) && (phiPhoton < TMath::Pi()/2) )
      tfLeft = true;
    else tfRight = true;

    if (labLeft != tfLeft) continue;
    if (labRight != tfRight) continue;
    if (ev==0) whichCuts+= "_forceAccep";//*/

    hPhiS->Fill(PhiS_simple);

    if (targetPosition==0){
      hPhiS_upS_up->Fill(PhiS_simple);

      if (PhiS_simple > 0) left_upS_up++;
      else right_upS_up++;
    }
    else{
      hPhiS_downS_up->Fill(PhiS_simple);

      if (PhiS_simple > 0) left_downS_up++;
      else right_downS_up++;
    }
  }//Event loop
  
  TCanvas* cPhiS = new TCanvas(); cPhiS->Divide(2);
  cPhiS->cd(1); hPhiS_upS_up->Draw();
  cPhiS->cd(2); hPhiS_downS_up->Draw();

  Double_t xvals[] = {0.5}, ex[] ={0.};
  Double_t An_upS = MakeAsym(left_upS_up, right_upS_up);
  Double_t eAn_upS = eAsym(left_upS_up, right_upS_up);
  Double_t An_downS = MakeAsym(left_downS_up, right_downS_up);
  Double_t eAn_downS = eAsym(left_downS_up, right_downS_up);

  TGraphErrors* gUpS = new TGraphErrors(1, xvals, &An_upS, ex, &eAn_upS);
  TGraphErrors* gDownS = new TGraphErrors(1, xvals, &An_downS, ex, &eAn_downS);
  SetUp(gUpS); SetUp(gDownS);

  TCanvas* cAsym = new TCanvas();
  gUpS->Draw("AP"); gUpS->SetMarkerColor(kRed);
  //gUpS->GetYaxis()->SetRangeUser(-0.005, 0.005);
  gUpS->GetYaxis()->SetRangeUser(-0.02, 0.02);
  gDownS->Draw("PSame"); gDownS->SetMarkerColor(kBlue);
  OffSet(gDownS, 0.04);
  DrawLine(gUpS, 0.0);

  //Write output
  TString fOutName = Form("mcLeftRight_%s%s.root", whichMC.Data(),
			  whichCuts.Data());
  TString thisDirPath ="/Users/robertheitz/Documents/Research/DrellYan/\
Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/Purity/Data/\
LeftRight/";
  TFile *fOut = new TFile(thisDirPath+fOutName, "RECREATE");
  hPhiS_upS_up->Write();
  hPhiS_downS_up->Write();
  gUpS->Write("gUpS");
  gDownS->Write("gDownS");
  fOut->Close();
}


Double_t MakeAsym(Double_t l_pUp, Double_t r_pUp){
  Double_t L = l_pUp;

  Double_t R = r_pUp;

  return (L - R)/(L + R);
}


Double_t eAsym(Double_t l_pUp,  Double_t r_pUp){

  return 1.0/TMath::Sqrt(l_pUp+r_pUp);
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
