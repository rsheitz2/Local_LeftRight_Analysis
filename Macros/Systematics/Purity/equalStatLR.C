#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Double_t DeterminePhi(Double_t vX, Double_t vY);

Double_t  ShiftPhiSimple (Double_t PhiS_simple);

Double_t MakeAsym(Double_t l_pUp, Double_t l_pdown,
		  Double_t r_pUp, Double_t r_pdown);

Double_t eAsym(Double_t l_pUp, Double_t l_pdown,
	       Double_t r_pUp, Double_t r_pdown);

void equalStatLR(){
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/";

  const Int_t nPeriod =9;
  TString periods[] =
  {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};//*/
  /*const Int_t nPeriod =7;
  TString periods[] =
    {"W07", "W08", "W10", "W11", "W12", "W13", "W15"};//*/
    
  Double_t PhiS, PhiS_simple, Spin_0, Polarization, dilutionFactor, Mmumu;
  Double_t vPhoton_X, vPhoton_Y, vOpenAngle;
  Double_t vx_z, vx_x, vx_y;
  Int_t targetPosition;

  //Left/right values
  Double_t An_upS[nPeriod], An_downS[nPeriod];
  Double_t eAn_upS[nPeriod], eAn_downS[nPeriod];
  Double_t xvals[nPeriod], ex[nPeriod] ={0.0};
  
  for (Int_t p=0; p<nPeriod; p++) {

    TString fname=Form("slot1%s_LowM_AMDY.root", periods[p].Data() );
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
    tree->SetBranchAddress("vx_z", &vx_z);
    tree->SetBranchAddress("vx_y", &vx_y);
    tree->SetBranchAddress("vx_x", &vx_x);

    TH1D* hpolUpS = new TH1D("hpolUpS", "hpolUpS", 4, -2, 2);
    TH1D* hpolDownS = new TH1D("hpolDownS", "hpolDownS", 4, -2, 2);
    Int_t totalStat =0;
    cout << "\n\n\nPeriod:  " << periods[p] << endl;
    tree->Draw("Spin_0>>hpolUpS", "targetPosition==0&&Mmumu<3.38&&Mmumu>2.87");
    cout << "UpS stats:   " << hpolUpS->GetBinContent(2) << " "
	 << hpolUpS->GetBinContent(4) << " " << endl;
    totalStat +=hpolUpS->GetBinContent(2); totalStat+=hpolUpS->GetBinContent(4);

    tree->Draw("Spin_0>>hpolDownS", "targetPosition==1&&Mmumu<3.38&&Mmumu>2.87");
    cout << "DownS stats:   " << hpolDownS->GetBinContent(2) << " "
	 << hpolDownS->GetBinContent(4) << " " << endl;
    totalStat +=hpolDownS->GetBinContent(2);
    totalStat+=hpolDownS->GetBinContent(4);
    cout << "Final total stats:  " << totalStat << endl;
  
    Double_t minStat = hpolUpS->GetBinContent(2);
    if (hpolUpS->GetBinContent(4) < minStat) minStat =hpolUpS->GetBinContent(4);
    if (hpolDownS->GetBinContent(4) < minStat)
      minStat =hpolDownS->GetBinContent(2);
    if (hpolDownS->GetBinContent(4) < minStat)
      minStat =hpolDownS->GetBinContent(4);
    cout << "min data " << minStat << endl;

    Double_t left_upS_up =0., left_upS_down =0.;
    Double_t left_downS_up =0., left_downS_down =0.;
    Double_t right_upS_up =0., right_upS_down =0.;
    Double_t right_downS_up =0., right_downS_down =0.;

    Int_t upSsb_one =0, downSsb_one =0;
    Int_t upSsb_two =0, downSsb_two =0;
    for (Int_t ev=0; ev<tree->GetEntries(); ev++) {
      //for (Int_t ev=0; ev<1000; ev++) { if (ev==0) cout << "debug" << endl;
      tree->GetEntry(ev);

      ////////////
      //Cuts
      //Mass cut
      if (Mmumu > 3.44) continue;
      if (Mmumu < 2.80) continue;

      //forced acceptance cut
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
      if (vOpenAngle > 0.17) continue;
      if (vOpenAngle < 0.05) continue;//*/

      //Target radius cut
      Double_t vx_rad2 = vx_x*vx_x + vx_y*vx_y;
      if (vx_rad2 > 0.4*0.4) continue;//*/

      if (targetPosition==0){
	if (Spin_0 > 0 ){
	  if (upSsb_one > minStat) continue;
	
	  if (PhiS_simple > 0) left_upS_up++;
	  else right_upS_up++;

	  upSsb_one++;
	}
	else{
	  if (upSsb_two > minStat) continue;

	  if (PhiS_simple < 0) left_upS_down++;
	  else right_upS_down++;

	  upSsb_two++;
	}
      }
      else{
	if (Spin_0 > 0 ){
	  if (downSsb_two > minStat) continue;
	
	  if (PhiS_simple > 0) left_downS_up++;
	  else right_downS_up++;

	  downSsb_two++;
	}
	else{
	  if (downSsb_one > minStat) continue;
	
	  if (PhiS_simple < 0) left_downS_down++;
	  else right_downS_down++;

	  downSsb_one++;
	}
      }
    }//event loop
    cout << upSsb_one << " " << upSsb_two << endl;
    cout << downSsb_one << " " << downSsb_two << endl;
						 
    An_upS[p] =
      MakeAsym(left_upS_up, left_upS_down, right_upS_up, right_upS_down);
    eAn_upS[p] =
      eAsym(left_upS_up, left_upS_down, right_upS_up, right_upS_down);
    An_downS[p] =
      MakeAsym(left_downS_up, left_downS_down, right_downS_up,right_downS_down);
    eAn_downS[p] =
      eAsym(left_downS_up, left_downS_down, right_downS_up, right_downS_down);

    xvals[p] = p+1;
  }//period loop
  
  TGraphErrors* gUpS = new TGraphErrors(nPeriod, xvals, An_upS, ex, eAn_upS);
  TGraphErrors* gDownS = new TGraphErrors(nPeriod, xvals,An_downS,ex,eAn_downS);
  OffSet(gDownS, 0.2);
  SetUp(gUpS); SetUp(gDownS);

  Double_t An_upS_wavg[1], eAn_upS_wavg[1];
  An_upS_wavg[0] = WeightedAvgAndError(gUpS, eAn_upS_wavg);
  Double_t An_downS_wavg[1], eAn_downS_wavg[1];
  An_downS_wavg[0] = WeightedAvgAndError(gDownS, eAn_downS_wavg);

  TGraphErrors* gUpS_wavg =
    new TGraphErrors(1, xvals, An_upS_wavg, ex, eAn_upS_wavg);
  TGraphErrors* gDownS_wavg =
    new TGraphErrors(1, xvals,An_downS_wavg,ex,eAn_downS_wavg);
  OffSet(gDownS_wavg, -0.4); OffSet(gUpS_wavg, -0.2);
  SetUp(gUpS_wavg); SetUp(gDownS_wavg);

  TCanvas* cAsym = new TCanvas();
  gUpS->Draw("AP"); gUpS->SetMarkerColor(kRed);
  gUpS->GetYaxis()->SetRangeUser(-0.013, 0.013);
  //gUpS->GetYaxis()->SetRangeUser(-0.05, 0.05);
  gDownS->Draw("PSame"); gDownS->SetMarkerColor(kBlue);//*/

  gUpS_wavg->Draw("Psame"); gUpS_wavg->SetMarkerStyle(25);
  gUpS_wavg->SetMarkerColor(kRed);
  gDownS_wavg->Draw("Psame"); gDownS_wavg->SetMarkerStyle(25);
  gDownS_wavg->SetMarkerColor(kBlue);

  DrawLine(gUpS, 0.0);

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
