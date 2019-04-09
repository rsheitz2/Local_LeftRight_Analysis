#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Double_t DeterminePhi(Double_t vX, Double_t vY);

Double_t  ShiftPhiSimple (Double_t PhiS_simple);

Double_t MakeAsym(Double_t l_pUp, Double_t l_pdown,
		  Double_t r_pUp, Double_t r_pdown);

Double_t eAsym(Double_t l_pUp, Double_t l_pdown,
	       Double_t r_pUp, Double_t r_pdown);

void LeftRight(){
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/";
  TString fname="slot1WAll_LowM_AMDY.root";
  //TString fname="slot1W12_LowM_AMDY.root";
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

  TH1D* hPhiS = new TH1D("hPhiS", "hPhiS", 100, -3.2, 3.2);

  TH1D* hPhiS_upS_up = new TH1D("hPhiS_upS_up", "hPhiS_upS_up",
				100, -3.2, 3.2);
  TH1D* hPhiS_upS_down = new TH1D("hPhiS_upS_down", "hPhiS_upS_down",
				  100, -3.2, 3.2);
  TH1D* hPhiS_downS_up = new TH1D("hPhiS_downS_up", "hPhiS_downS_up",
				  100, -3.2, 3.2);
  TH1D* hPhiS_downS_down = new TH1D("hPhiS_downS_down", "hPhiS_downS_down",
				    100, -3.2, 3.2);
  const Int_t nTarg =4;
  TString targNames[nTarg] = {"ups_up", "ups_down", "downs_up", "downs_down"};
  TH1D* hMuMu[nTarg];
  for (Int_t i=0; i<nTarg; i++) {
    hMuMu[i] = new TH1D(Form("hMuMu%s", targNames[i].Data()),
			Form("hMuMu%s", targNames[i].Data()), 100, 2, 4.3);  
  }

  SetUp(hPhiS); SetUp(hPhiS_upS_up); SetUp(hPhiS_upS_down);
  SetUp(hPhiS_downS_up); SetUp(hPhiS_downS_down);

  //Deterined minimum statistics
  TH1D* hpolUpS = new TH1D("hpolUpS", "hpolUpS", 4, -2, 2);
  TH1D* hpolDownS = new TH1D("hpolDownS", "hpolDownS", 4, -2, 2);
  tree->Draw("Spin_0>>hpolUpS", "targetPosition==0&&Mmumu<3.475&&Mmumu>2.80");
  cout << "UpS stats:   " << hpolUpS->GetBinContent(2) << " "
       << hpolUpS->GetBinContent(4) << " " << endl;

  tree->Draw("Spin_0>>hpolDownS", "targetPosition==1&&Mmumu<3.468&&Mmumu>2.73");
  cout << "DownS stats:   " << hpolDownS->GetBinContent(2) << " "
       << hpolDownS->GetBinContent(4) << " " << endl;
  
  Double_t minStat = hpolUpS->GetBinContent(2);
  Double_t upS_downstat = minStat;
  Double_t upS_upstat = hpolUpS->GetBinContent(4);
  Double_t downS_downstat = hpolDownS->GetBinContent(2);
  Double_t downS_upstat = hpolDownS->GetBinContent(4);
  if (hpolUpS->GetBinContent(4) < minStat) minStat =hpolUpS->GetBinContent(4);
  if (hpolDownS->GetBinContent(4) < minStat)
    minStat =hpolDownS->GetBinContent(2);
  if (hpolDownS->GetBinContent(4) < minStat)
    minStat =hpolDownS->GetBinContent(4);
  cout << "min data " << minStat << endl;
  Int_t upSsb_one =0, downSsb_one =0;
  Int_t upSsb_two =0, downSsb_two =0;

  Double_t left_upS_up =0., left_upS_down =0.;
  Double_t left_downS_up =0., left_downS_down =0.;
  Double_t right_upS_up =0., right_upS_down =0.;
  Double_t right_downS_up =0., right_downS_down =0.;

  Double_t alphaP_UpDown = 0.125;

  for (Int_t ev=0; ev<tree->GetEntries(); ev++) {
    //for (Int_t ev=0; ev<1000; ev++) { if (ev==0) cout << "debug" << endl;
    tree->GetEntry(ev);

    /*if (Mmumu > 3.17) continue;
    if (Mmumu < 3.07) continue;//*/
    /*if (Mmumu > 3.22) continue;
    if (Mmumu < 3.02) continue;//*/
    /*if (Mmumu > 3.44) continue;
    if (Mmumu < 2.80) continue;//*/
    /*if (Mmumu > 8.5) continue;
    if (Mmumu < 4.3) continue;//*/

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
    if (labRight != tfRight) continue;//*/

    hPhiS->Fill(PhiS_simple);

    if (targetPosition==0){
      if (Mmumu > 3.475) continue;
      if (Mmumu < 2.80) continue;//*/
      
      if (Spin_0 > 0 ){
	//if (upSsb_one >= minStat) continue;
	if (upSsb_one >= alphaP_UpDown*upS_downstat) continue;
	upSsb_one++;
	
	hPhiS_upS_up->Fill(PhiS_simple);

	if (PhiS_simple > 0) left_upS_up++;
	else right_upS_up++;

	hMuMu[0]->Fill(Mmumu);
      }
      else{
	//if (upSsb_two >= minStat) continue;
	if (upSsb_two*alphaP_UpDown >= upS_upstat) continue;
	upSsb_two++;
	
	hPhiS_upS_down->Fill(PhiS_simple);

	if (PhiS_simple < 0) left_upS_down++;
	else right_upS_down++;

	hMuMu[1]->Fill(Mmumu);
      }
    }
    else{
      if (Mmumu > 3.468) continue;
      if (Mmumu < 2.73) continue;//*/
      
      if (Spin_0 > 0 ){
	//if (downSsb_two > minStat) continue;
	if (downSsb_two > alphaP_UpDown*downS_downstat) continue;
	downSsb_two++;
		  
	hPhiS_downS_up->Fill(PhiS_simple);

	if (PhiS_simple > 0) left_downS_up++;
	else right_downS_up++;

	hMuMu[2]->Fill(Mmumu);
      }
      else{
	//if (downSsb_one > minStat) continue;
	if (downSsb_one*alphaP_UpDown > downS_upstat) continue;
	downSsb_one++;
	
	hPhiS_downS_down->Fill(PhiS_simple);

	if (PhiS_simple < 0) left_downS_down++;
	else right_downS_down++;

	hMuMu[3]->Fill(Mmumu);
      }
    }
  }
  
  TCanvas* cPhiS = new TCanvas(); cPhiS->Divide(2,2);
  cPhiS->cd(1); hPhiS_upS_up->Draw();
  cPhiS->cd(2); hPhiS_upS_down->Draw();
  cPhiS->cd(3); hPhiS_downS_up->Draw();
  cPhiS->cd(4); hPhiS_downS_down->Draw();
  hPhiS_upS_up->Sumw2(); hPhiS_upS_down->Sumw2();
  hPhiS_downS_up->Sumw2(); hPhiS_downS_down->Sumw2();

  TH1D* hPhiS_upS = (TH1D*)hPhiS_upS_up->Clone("hPhiS_upS");
  hPhiS_upS->Divide(hPhiS_upS_down);
  TH1D* hPhiS_downS = (TH1D*)hPhiS_downS_up->Clone("hPhiS_downS");
  hPhiS_downS->Divide(hPhiS_downS_down);
  TCanvas* cRphiS = new TCanvas(); cRphiS->Divide(2);
  cRphiS->cd(1); hPhiS_upS->Draw("e");
  cRphiS->cd(2); hPhiS_downS->Draw("e");

  Double_t xvals[] = {0.5}, ex[] ={0.};
  Double_t An_upS =
    MakeAsym(left_upS_up, left_upS_down, right_upS_up, right_upS_down);
  Double_t eAn_upS =
    eAsym(left_upS_up, left_upS_down, right_upS_up, right_upS_down);
  Double_t An_downS =
    MakeAsym(left_downS_up, left_downS_down, right_downS_up, right_downS_down);
  Double_t eAn_downS =
    eAsym(left_downS_up, left_downS_down, right_downS_up, right_downS_down);

  TGraphErrors* gUpS = new TGraphErrors(1, xvals, &An_upS, ex, &eAn_upS);
  TGraphErrors* gDownS = new TGraphErrors(1, xvals, &An_downS, ex, &eAn_downS);
  SetUp(gUpS); SetUp(gDownS);

  TCanvas* cAsym = new TCanvas();
  gUpS->Draw("AP"); gUpS->SetMarkerColor(kRed);
  //gUpS->GetYaxis()->SetRangeUser(-0.005, 0.005);
  gUpS->GetYaxis()->SetRangeUser(-0.05, 0.05);
  gDownS->Draw("PSame"); gDownS->SetMarkerColor(kBlue);

  TCanvas* cMuMu = new TCanvas(); cMuMu->Divide(2,2);
  for (Int_t i=0; i<nTarg; i++) {
    cMuMu->cd(i+1);
    hMuMu[i]->Draw();
  }

}


Double_t MakeAsym(Double_t l_pUp, Double_t l_pdown,
		  Double_t r_pUp, Double_t r_pdown){

  Double_t L = l_pUp*l_pdown;
  L = TMath::Sqrt(L);

  Double_t R = r_pUp*r_pdown;
  R = TMath::Sqrt(R);

  return (L - R)/(L + R);
}


Double_t eAsym(Double_t l_pUp, Double_t l_pdown,
	       Double_t r_pUp, Double_t r_pdown){

  //Regular geomean asymmetry
  Double_t L = l_pUp*l_pdown;
  L = TMath::Sqrt(L);

  Double_t R = r_pUp*r_pdown;
  R = TMath::Sqrt(R);

  Double_t A = MakeAsym(l_pUp, l_pdown, r_pUp, r_pdown);

  Double_t dL2 = 0.5*L*0.5*L*( 1./l_pUp + 1./l_pdown );
  Double_t dR2 = 0.5*R*0.5*R*( 1./r_pUp + 1./r_pdown );

  Double_t error = ( 1-A )*( 1-A )*dL2 + ( 1+A )*( 1+A )*dR2;
  error = TMath::Sqrt( error );
  error *= 1.0/( L + R );

  if (error < 10e-9) {
    cout << "Error: AN error way too small" << endl;
    exit(EXIT_FAILURE);
  }

  return error;
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
