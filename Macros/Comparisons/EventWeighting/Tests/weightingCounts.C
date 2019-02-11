#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

void SetBounds (TString binfile, TString boundVar, Double_t *bounds,
		Double_t *xval, Int_t nBins);

void binVal(Double_t *bounds, Double_t binVar, Double_t *counts,
	    Double_t weight, Int_t nBins);


Double_t  ShiftPhiSimple (Double_t PhiS_simple) {
  //Lab reference frame:
  //    -Pi/2 ->   P/i/2 = left
  //    Pi/2  -> 3*P/i/2 = right
  Double_t phi = TMath::Pi()/2 - PhiS_simple;
  
  return phi;
}//ShiftPhiSimple


void weightingCounts(){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  TString production ="slot1";
  TString period ="W13";
  TString physBinned ="xPi";
  
  Double_t Mmin =4.3;//Not changed for now
  Double_t Mmax =8.5;
  Double_t phiScut =0.0;

  Bool_t debug =false;
  Bool_t toWrite =true;
  //Setup_______________

  //Settings
  cout << "nBins:   " << nBins << endl;
  cout << "Mtype:   " << Mtype << endl;
  cout << "\nData type____" << endl;
  cout << "production:   " << production << endl;
  cout << "period:   " << period << endl;
  cout << "physBinned:   " << physBinned << endl;
  cout << "\nCuts____" << endl;
  cout << Mmin << "  < Mmumu <  " << Mmax << endl;
  cout << "phiScut:  " << phiScut << endl;
  cout << "\ntoWrite:   " << toWrite << endl;
  
  //Get data
  TString pathData="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData/";
  TString fname = Form("%s/%s/%s%s_%s.root", pathData.Data(), Mtype.Data(),
		       production.Data(), period.Data(), Mtype.Data());
  TFile *fIn = OpenFile(fname);
  TTree *tree = (TTree*)fIn->Get("pT_Weighted");

  //Binning setup
  TString binfile =
    Form("%s/%s/BinValues/%sWAll_%s_%ibinsRelease.txt", pathData.Data(),
	 Mtype.Data(), production.Data(), Mtype.Data(), nBins);

  Double_t bounds[nBins+1], xval[nBins];
  if (physBinned=="M") physBinned = "mass";
  SetBounds(binfile, physBinned, bounds, xval, nBins);
  if (physBinned=="mass") physBinned = "M";

  //Loop number setup
  Int_t numEvents = (debug) ? 1000 : tree->GetEntries();
  cout << "\nNumber of events in tree:  " << tree->GetEntries();
  cout << "\nNumber of events to loop through:  " << numEvents << endl;
  cout << "Debug mode is ineffect:   " << debug << endl;
  cout << "" << endl;

  //Tree Variable Setup
  Double_t PhiS_simple, Spin_0;
  Double_t dilutionFactor, Polarization;
  Int_t targetPosition;
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;
  Double_t *binVar;
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("Spin_0", &Spin_0);
  tree->SetBranchAddress("dilutionFactor", &dilutionFactor);
  tree->SetBranchAddress("Polarization", &Polarization);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  tree->SetBranchAddress("x_beam", &x_beam);
  tree->SetBranchAddress("x_target", &x_target);
  tree->SetBranchAddress("x_feynman", &x_feynman);
  tree->SetBranchAddress("q_transverse", &q_transverse);
  tree->SetBranchAddress("Mmumu", &Mmumu);

  if (physBinned=="xN") binVar =&x_target;
  else if (physBinned=="xPi") binVar =&x_beam;
  else if (physBinned=="xF") binVar =&x_feynman;
  else if (physBinned=="pT") binVar =&q_transverse;
  else if (physBinned=="M") binVar =&Mmumu;

  //Varible count setup
  Double_t upS_Pup=0.;//upstream
  Double_t c_upS_Pup=0., c2_upS_Pup=0., c3_upS_Pup=0., c4_upS_Pup=0.;
  Double_t upS_Pup_phys[nBins]={0.};
  Double_t c_upS_Pup_phys[nBins]={0.}, c2_upS_Pup_phys[nBins]={0.};
  Double_t c3_upS_Pup_phys[nBins]={0.}, c4_upS_Pup_phys[nBins]={0.};
  Double_t upS_Pdown=0.;
  Double_t c_upS_Pdown=0., c2_upS_Pdown=0., c3_upS_Pdown=0., c4_upS_Pdown=0.;
  Double_t upS_Pdown_phys[nBins]={0.};
  Double_t c_upS_Pdown_phys[nBins]={0.}, c2_upS_Pdown_phys[nBins]={0.};
  Double_t c3_upS_Pdown_phys[nBins]={0.}, c4_upS_Pdown_phys[nBins]={0.};

  Double_t downS_Pup=0.;//downstream
  Double_t c_downS_Pup=0., c2_downS_Pup=0., c3_downS_Pup=0., c4_downS_Pup=0.;
  Double_t downS_Pup_phys[nBins]={0.};
  Double_t c_downS_Pup_phys[nBins]={0.}, c2_downS_Pup_phys[nBins]={0.};
  Double_t c3_downS_Pup_phys[nBins]={0.}, c4_downS_Pup_phys[nBins]={0.};
  Double_t downS_Pdown=0.;
  Double_t c_downS_Pdown=0., c2_downS_Pdown=0.;
  Double_t c3_downS_Pdown=0., c4_downS_Pdown=0.;
  Double_t downS_Pdown_phys[nBins]={0.};
  Double_t c_downS_Pdown_phys[nBins]={0.}, c2_downS_Pdown_phys[nBins]={0.};
  Double_t c3_downS_Pdown_phys[nBins]={0.}, c4_downS_Pdown_phys[nBins]={0.};

  //Polarizations/dilutions
  Double_t polUp_upS =0.0, polDown_upS =0.0;
  Double_t polUp_downS =0.0, polDown_downS =0.0;
  Double_t pol_upS =0.0, pol_downS =0.0;

  Double_t polUp_upS_phys[nBins] ={0.0}, polDown_upS_phys[nBins] ={0.0};
  Double_t polUp_downS_phys[nBins] ={0.0}, polDown_downS_phys[nBins] ={0.0};
  Double_t pol_upS_phys[nBins] ={0.0}, pol_downS_phys[nBins] ={0.0};

  //Tree Loop
  Int_t n=0;
  for (Int_t ev=0; ev<numEvents; ev++) {
    tree->GetEntry(ev);
    
    //Cuts
    if (Mmumu < Mmin) continue;
    if (Mmumu > Mmax) continue;
    if ((PhiS_simple < phiScut) && (PhiS_simple > -phiScut)) continue;
    if (PhiS_simple < -TMath::Pi() + phiScut) continue;
    if (PhiS_simple > TMath::Pi() - phiScut) continue;

    Double_t phi_photon_lab = ShiftPhiSimple(PhiS_simple);
    Double_t cos = TMath::Cos(phi_photon_lab);

    if (targetPosition==0){//upstream
      Double_t dil = TMath::Abs(dilutionFactor);
      Double_t correct_dil = 0.95*dil;
      Double_t pol = TMath::Abs(Polarization);
      
      if (Spin_0>0){//spin up
	upS_Pup++;
	c_upS_Pup += cos;
	c2_upS_Pup += cos*cos;
	c3_upS_Pup += cos*cos*cos;
	c4_upS_Pup += cos*cos*cos*cos;

	binVal(bounds, *binVar, upS_Pup_phys, 1.0, nBins);
	binVal(bounds, *binVar, c_upS_Pup_phys, cos, nBins);
	binVal(bounds, *binVar, c2_upS_Pup_phys, cos*cos, nBins);
	binVal(bounds, *binVar, c3_upS_Pup_phys, cos*cos*cos, nBins);
	binVal(bounds, *binVar, c4_upS_Pup_phys, cos*cos*cos*cos, nBins);

	polUp_upS += pol*correct_dil;
	binVal(bounds, *binVar, polUp_upS_phys, pol*correct_dil, nBins);
      }
      else {//spin down
	upS_Pdown++;
	c_upS_Pdown += cos;
	c2_upS_Pdown += cos*cos;
	c3_upS_Pdown += cos*cos*cos;
	c4_upS_Pdown += cos*cos*cos*cos;

	binVal(bounds, *binVar, upS_Pdown_phys, 1.0, nBins);
	binVal(bounds, *binVar, c_upS_Pdown_phys, cos, nBins);
	binVal(bounds, *binVar, c2_upS_Pdown_phys, cos*cos, nBins);
	binVal(bounds, *binVar, c3_upS_Pdown_phys, cos*cos*cos, nBins);
	binVal(bounds, *binVar, c4_upS_Pdown_phys, cos*cos*cos*cos, nBins);

	polDown_upS += pol*correct_dil;
	binVal(bounds, *binVar, polDown_upS_phys, pol*correct_dil, nBins);
      }
    }
    else if (targetPosition==1){//downstream
      Double_t dil = TMath::Abs(dilutionFactor);
      Double_t correct_dil = 0.91*dil;
      Double_t pol = TMath::Abs(Polarization);
	
      if (Spin_0>0){//spin up
	downS_Pup++;
	c_downS_Pup += cos;
	c2_downS_Pup += cos*cos;
	c3_downS_Pup += cos*cos*cos;
	c4_downS_Pup += cos*cos*cos*cos;

	binVal(bounds, *binVar, downS_Pup_phys, 1.0, nBins);
	binVal(bounds, *binVar, c_downS_Pup_phys, cos, nBins);
	binVal(bounds, *binVar, c2_downS_Pup_phys, cos*cos, nBins);
	binVal(bounds, *binVar, c3_downS_Pup_phys, cos*cos*cos, nBins);
	binVal(bounds, *binVar, c4_downS_Pup_phys, cos*cos*cos*cos, nBins);

	polUp_downS += pol*correct_dil;
	binVal(bounds, *binVar, polUp_downS_phys, pol*correct_dil, nBins);
      }
      else {//spin down
	downS_Pdown++;
	c_downS_Pdown += cos;
	c2_downS_Pdown += cos*cos;
	c3_downS_Pdown += cos*cos*cos;
	c4_downS_Pdown += cos*cos*cos*cos;

	binVal(bounds, *binVar, downS_Pdown_phys, 1.0, nBins);
	binVal(bounds, *binVar, c_downS_Pdown_phys, cos, nBins);
	binVal(bounds, *binVar, c2_downS_Pdown_phys, cos*cos, nBins);
	binVal(bounds, *binVar, c3_downS_Pdown_phys, cos*cos*cos, nBins);
	binVal(bounds, *binVar, c4_downS_Pdown_phys, cos*cos*cos*cos, nBins);

	polDown_downS += pol*correct_dil;
	binVal(bounds, *binVar, polDown_downS_phys, pol*correct_dil, nBins);
      }
    }
    else{
      cout << "Target does not exist" << endl;
      exit(EXIT_FAILURE);
    }
  }//end tree loop

  //Polarization average
  pol_upS = polUp_upS + polDown_upS;
  pol_upS /= (upS_Pup + upS_Pdown);
  polUp_upS /= upS_Pup; polDown_upS /= upS_Pdown;
  pol_downS = polUp_downS + polDown_downS;
  pol_downS /= (downS_Pup + downS_Pdown);
  polUp_downS /= downS_Pup; polDown_downS /= downS_Pdown;
  for (Int_t i=0; i<nBins; i++) {
    pol_upS_phys[i] = polUp_upS_phys[i] + polDown_upS_phys[i];
    pol_upS_phys[i] /= (upS_Pup_phys[i] + upS_Pdown_phys[i]);
    polUp_upS_phys[i] /= upS_Pup_phys[i];
    polDown_upS_phys[i] /= upS_Pdown_phys[i];

    pol_downS_phys[i] = polUp_downS_phys[i] + polDown_downS_phys[i];
    pol_downS_phys[i] /= (downS_Pup_phys[i] + downS_Pdown_phys[i]);
    polUp_downS_phys[i] /= downS_Pup_phys[i];
    polDown_downS_phys[i] /= downS_Pdown_phys[i];
  }

  //upstream
  TGraph *g_upS_Pup = new TGraph(1, &(xval[0]), &upS_Pup);//pol up
  TGraph *g_c_upS_Pup = new TGraph(1, &(xval[0]), &c_upS_Pup);
  TGraph *g_c2_upS_Pup = new TGraph(1, &(xval[0]), &c2_upS_Pup);
  TGraph *g_c3_upS_Pup = new TGraph(1, &(xval[0]), &c3_upS_Pup);
  TGraph *g_c4_upS_Pup = new TGraph(1, &(xval[0]), &c4_upS_Pup);
  TGraph *g_upS_Pup_phys = new TGraph(nBins, xval, upS_Pup_phys);
  TGraph *g_c_upS_Pup_phys = new TGraph(nBins, xval, c_upS_Pup_phys);
  TGraph *g_c2_upS_Pup_phys = new TGraph(nBins, xval, c2_upS_Pup_phys);
  TGraph *g_c3_upS_Pup_phys = new TGraph(nBins, xval, c3_upS_Pup_phys);
  TGraph *g_c4_upS_Pup_phys = new TGraph(nBins, xval, c4_upS_Pup_phys);

  TGraph *g_upS_Pdown = new TGraph(1, &(xval[0]), &upS_Pdown);//pol down
  TGraph *g_c_upS_Pdown = new TGraph(1, &(xval[0]), &c_upS_Pdown);
  TGraph *g_c2_upS_Pdown = new TGraph(1, &(xval[0]), &c2_upS_Pdown);
  TGraph *g_c3_upS_Pdown = new TGraph(1, &(xval[0]), &c3_upS_Pdown);
  TGraph *g_c4_upS_Pdown = new TGraph(1, &(xval[0]), &c4_upS_Pdown);
  TGraph *g_upS_Pdown_phys = new TGraph(nBins, xval, upS_Pdown_phys);
  TGraph *g_c_upS_Pdown_phys = new TGraph(nBins, xval, c_upS_Pdown_phys);
  TGraph *g_c2_upS_Pdown_phys = new TGraph(nBins, xval, c2_upS_Pdown_phys);
  TGraph *g_c3_upS_Pdown_phys = new TGraph(nBins, xval, c3_upS_Pdown_phys);
  TGraph *g_c4_upS_Pdown_phys = new TGraph(nBins, xval, c4_upS_Pdown_phys);

  TGraph *g_pol_upS = new TGraph(1, &(xval[0]), &pol_upS);
  TGraph *g_polUp_upS = new TGraph(1, &(xval[0]), &polUp_upS);
  TGraph *g_polDown_upS = new TGraph(1, &(xval[0]), &polDown_upS);
  TGraph *g_pol_upS_phys = new TGraph(nBins, xval, pol_upS_phys);
  TGraph *g_polUp_upS_phys = new TGraph(nBins, xval, polUp_upS_phys);
  TGraph *g_polDown_upS_phys = new TGraph(nBins, xval, polDown_upS_phys);

  //downstream
  TGraph *g_downS_Pup = new TGraph(1, &(xval[0]), &downS_Pup);//pol up
  TGraph *g_c_downS_Pup = new TGraph(1, &(xval[0]), &c_downS_Pup);
  TGraph *g_c2_downS_Pup = new TGraph(1, &(xval[0]), &c2_downS_Pup);
  TGraph *g_c3_downS_Pup = new TGraph(1, &(xval[0]), &c3_downS_Pup);
  TGraph *g_c4_downS_Pup = new TGraph(1, &(xval[0]), &c4_downS_Pup);
  TGraph *g_downS_Pup_phys = new TGraph(nBins, xval, downS_Pup_phys);
  TGraph *g_c_downS_Pup_phys = new TGraph(nBins, xval, c_downS_Pup_phys);
  TGraph *g_c2_downS_Pup_phys = new TGraph(nBins, xval, c2_downS_Pup_phys);
  TGraph *g_c3_downS_Pup_phys = new TGraph(nBins, xval, c3_downS_Pup_phys);
  TGraph *g_c4_downS_Pup_phys = new TGraph(nBins, xval, c4_downS_Pup_phys);

  TGraph *g_downS_Pdown = new TGraph(1, &(xval[0]), &downS_Pdown);//pol down
  TGraph *g_c_downS_Pdown = new TGraph(1, &(xval[0]), &c_downS_Pdown);
  TGraph *g_c2_downS_Pdown = new TGraph(1, &(xval[0]), &c2_downS_Pdown);
  TGraph *g_c3_downS_Pdown = new TGraph(1, &(xval[0]), &c3_downS_Pdown);
  TGraph *g_c4_downS_Pdown = new TGraph(1, &(xval[0]), &c4_downS_Pdown);
  TGraph *g_downS_Pdown_phys = new TGraph(nBins, xval, downS_Pdown_phys);
  TGraph *g_c_downS_Pdown_phys = new TGraph(nBins, xval, c_downS_Pdown_phys);
  TGraph *g_c2_downS_Pdown_phys = new TGraph(nBins, xval, c2_downS_Pdown_phys);
  TGraph *g_c3_downS_Pdown_phys = new TGraph(nBins, xval, c3_downS_Pdown_phys);
  TGraph *g_c4_downS_Pdown_phys = new TGraph(nBins, xval, c4_downS_Pdown_phys);

  TGraph *g_pol_downS = new TGraph(1, &(xval[0]), &pol_downS);
  TGraph *g_polUp_downS = new TGraph(1, &(xval[0]), &polUp_downS);
  TGraph *g_polDown_downS = new TGraph(1, &(xval[0]), &polDown_downS);
  TGraph *g_pol_downS_phys = new TGraph(nBins, xval, pol_downS_phys);
  TGraph *g_polUp_downS_phys = new TGraph(nBins, xval, polUp_downS_phys);
  TGraph *g_polDown_downS_phys = new TGraph(nBins, xval, polDown_downS_phys);

  SetUp(g_upS_Pup); SetUp(g_c_upS_Pup); SetUp(g_c2_upS_Pup);SetUp(g_c3_upS_Pup);
  SetUp(g_c4_upS_Pup); SetUp(g_upS_Pup_phys); SetUp(g_c_upS_Pup_phys);
  SetUp(g_c2_upS_Pup_phys); SetUp(g_c3_upS_Pup_phys); SetUp(g_c4_upS_Pup_phys);

  SetUp(g_upS_Pdown); SetUp(g_c_upS_Pdown); SetUp(g_c2_upS_Pdown);
  SetUp(g_c3_upS_Pdown); SetUp(g_c4_upS_Pdown); SetUp(g_upS_Pdown_phys);
  SetUp(g_c_upS_Pdown_phys); SetUp(g_c2_upS_Pdown_phys);
  SetUp(g_c3_upS_Pdown_phys); SetUp(g_c4_upS_Pdown_phys);
  SetUp(g_polUp_upS); SetUp(g_polDown_upS);
  SetUp(g_polUp_upS_phys); SetUp(g_polDown_upS_phys);
  SetUp(g_pol_upS); SetUp(g_pol_upS_phys);

  SetUp(g_downS_Pup); SetUp(g_c_downS_Pup); SetUp(g_c2_downS_Pup);
  SetUp(g_c3_downS_Pup); SetUp(g_c4_downS_Pup); SetUp(g_downS_Pup_phys);
  SetUp(g_c_downS_Pup_phys); SetUp(g_c2_downS_Pup_phys);
  SetUp(g_c3_downS_Pup_phys); SetUp(g_c4_downS_Pup_phys);

  SetUp(g_downS_Pdown); SetUp(g_c_downS_Pdown); SetUp(g_c2_downS_Pdown);
  SetUp(g_c3_downS_Pdown); SetUp(g_c4_downS_Pdown); SetUp(g_downS_Pdown_phys);
  SetUp(g_c_downS_Pdown_phys); SetUp(g_c2_downS_Pdown_phys);
  SetUp(g_c3_downS_Pdown_phys); SetUp(g_c4_downS_Pdown_phys);
  SetUp(g_polUp_downS); SetUp(g_polDown_downS);
  SetUp(g_polUp_downS_phys); SetUp(g_polDown_downS_phys);
  SetUp(g_pol_downS); SetUp(g_pol_downS_phys);

  //Testing
  /*TCanvas* c1 = new TCanvas();
  //g_upS_Pup->Draw("AP");
  g_polUp_upS->Draw("AP");
  g_polDown_upS->Draw("Psame"); g_polDown_upS->SetMarkerColor(kBlue);
  //*/
  
  //Write output
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting/Data/Counts";
  TString fOutput =
    Form("%s/counts_%s%i_%s_%s%s.root", thisDirPath.Data(),
	 physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	 period.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_upS_Pup->Write("upS_Pup");
    g_c_upS_Pup->Write("c_upS_Pup");
    g_c2_upS_Pup->Write("c2_upS_Pup");
    g_c3_upS_Pup->Write("c3_upS_Pup");
    g_c4_upS_Pup->Write("c4_upS_Pup");
    g_upS_Pup_phys->Write("upS_Pup_phys");
    g_c_upS_Pup_phys->Write("c_upS_Pup_phys");
    g_c2_upS_Pup_phys->Write("c2_upS_Pup_phys");
    g_c3_upS_Pup_phys->Write("c3_upS_Pup_phys");
    g_c4_upS_Pup_phys->Write("c4_upS_Pup_phys");

    g_upS_Pdown->Write("upS_Pdown");
    g_c_upS_Pdown->Write("c_upS_Pdown");
    g_c2_upS_Pdown->Write("c2_upS_Pdown");
    g_c3_upS_Pdown->Write("c3_upS_Pdown");
    g_c4_upS_Pdown->Write("c4_upS_Pdown");
    g_upS_Pdown_phys->Write("upS_Pdown_phys");
    g_c_upS_Pdown_phys->Write("c_upS_Pdown_phys");
    g_c2_upS_Pdown_phys->Write("c2_upS_Pdown_phys");
    g_c3_upS_Pdown_phys->Write("c3_upS_Pdown_phys");
    g_c4_upS_Pdown_phys->Write("c4_upS_Pdown_phys");

    g_pol_upS->Write("pol_upS");
    g_polUp_upS->Write("polUp_upS");
    g_polDown_upS->Write("polDown_upS");
    g_pol_upS_phys->Write("pol_upS_phys");
    g_polUp_upS_phys->Write("polUp_upS_phys");
    g_polDown_upS_phys->Write("polDown_upS_phys");

    g_downS_Pup->Write("downS_Pup");
    g_c_downS_Pup->Write("c_downS_Pup");
    g_c2_downS_Pup->Write("c2_downS_Pup");
    g_c3_downS_Pup->Write("c3_downS_Pup");
    g_c4_downS_Pup->Write("c4_downS_Pup");
    g_downS_Pup_phys->Write("downS_Pup_phys");
    g_c_downS_Pup_phys->Write("c_downS_Pup_phys");
    g_c2_downS_Pup_phys->Write("c2_downS_Pup_phys");
    g_c3_downS_Pup_phys->Write("c3_downS_Pup_phys");
    g_c4_downS_Pup_phys->Write("c4_downS_Pup_phys");

    g_downS_Pdown->Write("downS_Pdown");
    g_c_downS_Pdown->Write("c_downS_Pdown");
    g_c2_downS_Pdown->Write("c2_downS_Pdown");
    g_c3_downS_Pdown->Write("c3_downS_Pdown");
    g_c4_downS_Pdown->Write("c4_downS_Pdown");
    g_downS_Pdown_phys->Write("downS_Pdown_phys");
    g_c_downS_Pdown_phys->Write("c_downS_Pdown_phys");
    g_c2_downS_Pdown_phys->Write("c2_downS_Pdown_phys");
    g_c3_downS_Pdown_phys->Write("c3_downS_Pdown_phys");
    g_c4_downS_Pdown_phys->Write("c4_downS_Pdown_phys");

    g_pol_downS->Write("pol_downS");
    g_polUp_downS->Write("polUp_downS");
    g_polDown_downS->Write("polDown_downS");
    g_pol_downS_phys->Write("pol_downS_phys");
    g_polUp_downS_phys->Write("polUp_downS_phys");
    g_polDown_downS_phys->Write("polDown_downS_phys");

    fResults->Close();
  }
}//end weightingCounts.C


void SetBounds (TString binfile, TString boundVar, Double_t *bounds,
		Double_t *xval, Int_t nBins){
  Int_t iter=0;
  if (boundVar=="xN"||boundVar=="xPi"||boundVar=="pT"||
      boundVar=="mass"||boundVar=="rad") { bounds[iter] =(0.0);}
  else if (boundVar=="xF") bounds[iter] =(-1.0);
  else if (boundVar=="vxZ_upstream") bounds[iter] =(-294.5);
  else if (boundVar=="vxZ_downstream") bounds[iter] =(-219.5);
  else {
    cout << "Invalid boundVar: " << boundVar << endl;
    exit(EXIT_FAILURE);
  }

  TString boundaries = boundVar; boundaries += " bin boundaries";
  TString averages = boundVar; averages += " bin averages";
  
  string line;
  bool first = false, found = false;
  ifstream f_bins(binfile);
  if(!f_bins.is_open() ) {
    cout << "\nbinFile did not open from" << endl;
    exit(EXIT_FAILURE);
  }
  
  while (!f_bins.eof()) {
    getline(f_bins,line);
    TString tline (line);

    if (tline == boundaries) {
      found = true;
      first = true;
      iter=1;
      continue;
    }
    else if (tline == averages) {
      first = false;
      iter=0;
      continue;
    }
    else if (!found) continue;
    
    if (first) {
      bounds[iter] =(atof(line.c_str() ) );
      iter++;
    }
    else {
      xval[iter] =(atof(line.c_str() ) );
      iter++;

      if (iter==nBins) break;
    }
  }

  iter =nBins;
  if (boundVar=="xN"||boundVar=="xPi"||boundVar=="xF") bounds[iter] =(1.0);
  else if (boundVar=="pT") bounds[iter] =(5.0);
  else if (boundVar=="mass") bounds[iter] =(12.0);
  else if (boundVar=="rad") bounds[iter] =(2.0);
  else if (boundVar=="vxZ_upstream") bounds[iter] =(-239.3);
  else if (boundVar=="vxZ_downstream") bounds[iter] =(-164.3);
}


void binVal(Double_t *bounds, Double_t binVar, Double_t *counts,
	    Double_t weight, Int_t nBins){
  if (binVar < bounds[0]){
    cout << "binVar is too small..." << endl;
    exit(EXIT_FAILURE);
  }
  
  for (Int_t i=1; i<nBins+1; i++) {
    if (binVar < bounds[i]){
      counts[i-1] += weight;
      return;
    }
  }

  if (binVar > bounds[nBins]){
    cout << "binVar is too big..." << endl;
    exit(EXIT_FAILURE);
  }
}
