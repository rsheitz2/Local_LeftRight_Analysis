#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PhiScut/Include/\
helperFunctions.h"

void realDataCutImpact(TString start=""){
  //Setup__________
  TString Mtype ="HMDY";
  TString period ="WAll";
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
    TString production ="slot1";//""=t3, "slot1"
  Bool_t debug =false;
  
  Bool_t toWrite =false;
  //Setup__________

  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);
  if (start==""){
    cout << "\nCurrent Settings:" << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Debug mode in use:            " << debug << endl;
    cout << "\nTo Write:                   " << toWrite << endl;
    cout << "\n\nUsage: " << endl;
    cout << "root \'resolutionPhiS(1).C\'" << endl;
    exit(EXIT_FAILURE);
  }

  //Get data paths
  TString pathRD ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData";
  TString inputFile = Form("%s/%s/%s%s_%s.root", pathRD.Data(), Mtype.Data(),
			   production.Data(), period.Data(), Mtype.Data() );
  TFile *fRD = OpenFile(inputFile);
  TTree *tRD = (TTree*) fRD->Get("pT_Weighted");

  //Setup Tree variables
  Double_t PhiS_simple, Mmumu, Spin_0;
  tRD->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tRD->SetBranchAddress("Mmumu", &Mmumu);
  tRD->SetBranchAddress("Spin_0", &Spin_0);

  //Distributions
  TH1D* hPhiS = new TH1D("hPhiS", "hPhiS", 100, -TMath::Pi(), TMath::Pi());
  TH1D* hPhiS_polUp = new TH1D("hPhiS_polUp", "hPhiS_polUp",
			       100, -TMath::Pi(), TMath::Pi());
  TH1D* hPhiS_polDown = new TH1D("hPhiS_polDown", "hPhiS_polDown",
				 100, -TMath::Pi(), TMath::Pi());

  //PhiS cut setup
  const Int_t nPhiScut =9;
  Double_t phiScuts[nPhiScut] =
    {0.0, 0.044, 0.088, 0.17, 0.36, 0.53, 0.71, 0.88, 1.07};
  Double_t NotCutData[nPhiScut] = {0.0}, ex[nPhiScut] ={0.0};
  
  //Tree loop
  Int_t tree_entries = (debug) ? 10000 : tRD->GetEntries();
  if (debug) cout << "Debug mode is in effect" << endl;
  else cout << "Debug mode is   NOT   in effect" << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    tRD->GetEntry(ev, 0);

    if (Mmumu < Mmin) continue;//Mmumu cut
    if (Mmumu > Mmax) continue;

    hPhiS->Fill(PhiS_simple);

    if (Spin_0 > 0) {//Spin up
      hPhiS_polUp->Fill(PhiS_simple);
    }
    else if (Spin_0 < 0){//Spin down
      hPhiS_polDown->Fill(PhiS_simple);
    }

    for (Int_t c=0; c<nPhiScut; c++) {
      if ((PhiS_simple > -TMath::Pi() + phiScuts[c])
	  && (PhiS_simple < -phiScuts[c]) ) {
	NotCutData[c]++;
      }
      else if ((PhiS_simple > phiScuts[c])
	       && (PhiS_simple < TMath::Pi() -phiScuts[c]) ){
	NotCutData[c]++;
      }
    }//nPhiSCuts loop
  }//tree loop

  //Draw distributions
  TCanvas* cPhiS = new TCanvas();
  hPhiS->Draw(); Setup(hPhiS);
  hPhiS_polUp->Draw("same"); Setup(hPhiS_polUp);
  hPhiS_polUp->SetLineColor(kGreen);
  hPhiS_polDown->Draw("same"); Setup(hPhiS_polDown);
  hPhiS_polDown->SetLineColor(kRed);

  //PhiS cut Distribution setup
  Double_t resPhiS =0.1776; //Yu rms value 2.0-8.5
  const Int_t nPhiSCuts =6;
  Double_t phiScut[nPhiSCuts] =
    {0., 0.25*resPhiS, 0.5*resPhiS, resPhiS, 2*resPhiS, 3*resPhiS};
  for (Int_t i=0; i<nPhiSCuts; i++) {
    Double_t zero_plus = phiScut[i];
    Double_t zero_minus = -1.0*phiScut[i];

    Double_t pi_plus = -TMath::Pi() + phiScut[i];
    Double_t pi_minus = TMath::Pi() - phiScut[i];
    
    DrawLine(hPhiS, zero_plus, i+1, "v");
    DrawLine(hPhiS, zero_minus, i+1, "v");
    DrawLine(hPhiS, pi_plus, i+1, "v");
    DrawLine(hPhiS, pi_minus, i+1, "v");
  }

  //PhiS not cut
  Double_t eNotCut[nPhiScut];
  Double_t CutPercent[nPhiScut];
  Double_t Pol = 0.13;
  for (Int_t c=0; c<nPhiScut; c++) {
    eNotCut[c] = 1.0/( TMath::Sqrt(NotCutData[c])*Pol );

    CutPercent[c] = (NotCutData[0] - NotCutData[c])/NotCutData[0];
  }
  TGraphErrors* gNotCut =
    new TGraphErrors(nPhiScut, phiScuts, NotCutData, ex, eNotCut);
  Setup(gNotCut);
  TGraphErrors *gCutPercent =
    new TGraphErrors(nPhiScut, phiScuts, CutPercent, ex, eNotCut);
  Setup(gCutPercent);
  
  TCanvas* c2 = new TCanvas(); c2->Divide(2);
  c2->cd(1);
  gNotCut->Draw("AP"); gNotCut->SetTitle("PhiS Cut Remaining Data");
  c2->cd(2);
  gCutPercent->Draw("AP"); gCutPercent->SetTitle("Percent Data Cut");
  gCutPercent->GetXaxis()->SetLimits(-0.2, 1.5);
  
  /*TString outName = Form("Data/ResolutionPhiS/resolution_%s_%0.2f_%0.2f.root",
			 whichMC.Data(), Mmin, Mmax);
  if (toWrite){
    TFile* fOut = new TFile(outName, "RECREATE");
    for (Int_t i=0; i<4; i++) { h1All[i]->Write(); }
    for (Int_t i=0; i<3; i++) { h2All[i]->Write(); }
    hAMDY_Diff->Write();
    hOC_Diff->Write();
    hJPsi_Diff->Write();
    hPsi_Diff->Write();

    gCrossOver->Write("gCrossOver");
    fOut->Close();
    }//*/

  //Final Output
  cout << "\nSettings !!!!" << endl;
  cout << "Debug mode in use:            " << debug << endl;
  //if (toWrite) { cout << outName << "   was written" << endl; }
  //else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
