#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PhiScut/Include/helperFunctions.h"

void sysStatError(TString start=""){
  //Setup__________
  const Double_t Mmin=2.0, Mmax=8.5; //Mass cut
  
  TString whichMC ="Yu"; //"Yu", "Charles"
  TString whichRD ="slot1"; //"t3", "slot1"

  Double_t AN_assumed =0.05;
  
  Bool_t toWrite =false;
  //Setup__________

  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);
  if (start==""){
    cout << "Macro determines total error as a function of PhiS cut" << endl;
    cout << "     Systematic error from resolutionPhiS.C" << endl;
    cout << "     Statistical error from phiScut.C" << endl;
    cout << "\nCurrent Settings:" << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Which MC considered:          " << whichMC << endl;
    cout << "Which RD considered:          " << whichRD << endl;
    cout << "\nTo Write:                   " << toWrite << endl;
    cout << "\n\nUsage: " << endl;
    cout << "root \'sysStatError(1).C\'" << endl;
    exit(EXIT_FAILURE);
  }

  //Get data files
  TString path ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PhiScut/DetermineCut/";
  TString nameSys = Form("Data/ResolutionPhiS/resolution_%s_%0.2f_%0.2f.root",
			 whichMC.Data(), Mmin, Mmax);
  TString nameStat = Form("Data/PhiScut/phiScut_%s_%0.2f_%0.2f.root",
			  whichRD.Data(), Mmin, Mmax);
  TFile *fsys = OpenFile(path+nameSys);
  TFile *fstat = OpenFile(path+nameStat);

  TGraph *gCrossOver = (TGraph*)fsys->Get("gCrossOver");
  TGraph *gstat = (TGraph*)fstat->Get("gTotError");

  Double_t *yCrossOver = gCrossOver->GetY();
  Double_t *ystat= gstat->GetY();
  const Int_t nPhiSCuts =6;
  Double_t ySys[nPhiSCuts], totError[nPhiSCuts];
  for (Int_t i=0; i<nPhiSCuts; i++) {
    ySys[i] = yCrossOver[i]*2.0*AN_assumed;
    
    totError[i] = TMath::Sqrt(ySys[i]*ySys[i] + ystat[i]*ystat[i]);
  }

  Double_t resPhiS =0.303; //rms value
  Double_t phiScut[nPhiSCuts] =
    {0., 0.25*resPhiS, 0.5*resPhiS, resPhiS, 2*resPhiS, 3*resPhiS};
  TGraph *gsys = new TGraph(nPhiSCuts, phiScut, ySys);
  TGraph *gTotError = new TGraph(nPhiSCuts, phiScut, totError);
  Setup(gsys); Setup(gTotError);

  TCanvas* c1 = new TCanvas();
  gstat->SetMarkerColor(kBlue);
  gTotError->SetMarkerColor(kRed);
  gTotError->Draw("AP");
  gTotError->GetYaxis()->SetRangeUser(0.0, AN_assumed/5.0);
  gsys->Draw("Psame");
  gstat->Draw("Psame");

  gTotError->SetTitle("Errors");

  TString outName =Form("Data/SysStatError/sysStatError_%s_%s_%0.2f_%0.2f.root",
			whichMC.Data(), whichRD.Data(), Mmin, Mmax);
  if (toWrite){
    TFile* fOut = new TFile(outName, "RECREATE");
    gsys->Write("sysError");
    gstat->Write("statError");
    gTotError->Write("totalError");

    fOut->Close();
  }

  //Final Output
  cout << "\nSettings !!!!" << endl;
  cout << "Mass Range is		:  " << MassCut << endl;
  cout << "Which MC considered:          " << whichMC << endl;
  cout << "Which RD considered:          " << whichRD << endl;
  if (toWrite) { cout << outName << "   was written" << endl; }
  else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
