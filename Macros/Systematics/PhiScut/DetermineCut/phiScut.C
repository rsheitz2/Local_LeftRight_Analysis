#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/Systematics/PhiScut/Include/helperFunctions.h"

void phiScut(TString start=""){
  //Setup_______________
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
  TString whichRD ="slot1"; //"t3", "slot1"

  Double_t resPhiS =0.1776; //Yu rms value 2.0-8.5
  //Double_t resPhiS =0.299; //Yu rms value 4.3-8.5
  
  Bool_t toWrite =true;
  //Setup_______________

  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);
  if (start==""){
    cout << "Macro determines the statistical impact on Real Data ";
    cout << " for a range of PhiS cuts determined from resolutionPhiS.C\n";
    cout << "\nCurrent Settings:" << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Which MC considered:          " << whichRD << endl;
    cout << "PhiS resolution assume:       " << resPhiS << endl;
    cout << "\nTo Write:                   " << toWrite << endl;
    cout << "\n\nUsage: " << endl;
    cout << "root \'phiScut(1).C\'" << endl;
    exit(EXIT_FAILURE);
  }

  //Get data file
  TString nameRD = "/Users/robertheitz/Documents/Research/DrellYan\
/Analysis/TGeant/Presents/DATA/RealData/LowM_AMDY/";
  if (whichRD =="t3") nameRD += "WAll_LowM_AMDY.root";
  else if (whichRD =="slot1") nameRD += "slot1WAll_LowM_AMDY.root";
    
  TFile *fData = TFile::Open(nameRD);
  if (!fData){
    cout << "File do not open" << endl;
    exit(EXIT_FAILURE);
  }

  TTree *t = (TTree*)fData->Get("pT_Weighted");
  Double_t PhiS, Mmumu;
  t->SetBranchAddress("PhiS", &PhiS);
  t->SetBranchAddress("Mmumu", &Mmumu);

  //Distributions
  TH1D* hPhiS = new TH1D("hPhiS", "hPhiS", 100, -TMath::Pi(), TMath::Pi());
  TH1D* hPhiPhoton = new TH1D("hPhiPhoton", "hPhiPhoton", 100,
			      -0.5*TMath::Pi(), 3*0.5*TMath::Pi());

  //PhiS cuts
  const Int_t nPhiScut =12;
  Double_t phiScut[nPhiScut] =
    {0., 0.044, 0.088, 0.17, 0.26, 0.36, 0.53, 0.71, 0.88, 1.07, 1.24, 1.44};
  /*Double_t phiScut[nPhiScut] =
    {0., 0.25*resPhiS, 0.5*resPhiS, resPhiS, 2*resPhiS, 3*resPhiS};//*/
  Double_t dataCut[nPhiScut] ={0.0};
  
  for (Int_t ev=0; ev<t->GetEntries(); ev++) {
  //cout << "\nDebug Mode \n"; for (Int_t ev=0; ev<1000; ev++) { //Debug mode
    t->GetEntry(ev);

    if (Mmumu < Mmin) continue;//Mmumu cut
    if (Mmumu > Mmax) continue;

    Double_t phiPhoton = ShiftPhiSimple(PhiS);
    
    hPhiS->Fill(PhiS);
    hPhiPhoton->Fill(phiPhoton);

    for (Int_t i=0; i<nPhiScut; i++) {
      Bool_t goodRange =true;
      if ( (phiPhoton < TMath::Pi()/2 + phiScut[i]) &&
	   (phiPhoton > TMath::Pi()/2 - phiScut[i])) goodRange =false;
      if ( (phiPhoton < -TMath::Pi()/2 + phiScut[i]) ||
	   (phiPhoton > 3*TMath::Pi()/2 - phiScut[i])) goodRange =false;

      if (!goodRange){ dataCut[i]++; }
    }
  }//End tree loop
  
  TCanvas* cPhi = new TCanvas(); cPhi->Divide(2);
  cPhi->cd(1); hPhiS->Draw();
  cPhi->cd(2); hPhiPhoton->Draw();

  TCanvas *cCut = new TCanvas(); cCut->Divide(2);
  Int_t totalData =hPhiS->GetEntries();
  Double_t totError[nPhiScut], totalPhiScut[nPhiScut];
  for (Int_t i=0; i<nPhiScut; i++) {
    totError[i] = 1.0/TMath::Sqrt(totalData - dataCut[i]);
      
    dataCut[i] *= (100.0/totalData);

    totalPhiScut[i] = 4*phiScut[i];
  }

  TGraph* gCutData = new TGraph(nPhiScut, totalPhiScut, dataCut);
  TGraph* gTotError = new TGraph(nPhiScut, totalPhiScut, totError);
  Setup(gCutData); Setup(gTotError);
  cCut->cd(1); gCutData->Draw("AP"); gCutData->SetTitle("Cut Data Percent");
  cCut->cd(2); gTotError->Draw("AP"); gTotError->SetTitle("Total Error");

  TString outName = Form("Data/PhiScut/phiScut_%s_%0.2f_%0.2f.root",
			 whichRD.Data(), Mmin, Mmax);
  if (toWrite){
    TFile* fOut = new TFile(outName, "RECREATE");
    gTotError->Write("gTotError");
    gCutData->Write("gCutData");
    hPhiS->Write();
    hPhiPhoton->Write();
    fOut->Close();
  }

  //Final Output
  cout << "\nSettings !!!!" << endl;
  cout << "Mass Range is		:  " << MassCut << endl;
  cout << "Which RD considered:          " << whichRD << endl;
  cout << "PhiS resolution assume:       " << resPhiS << endl;
  if (toWrite) { cout << outName << "   was written" << endl; }
  else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
