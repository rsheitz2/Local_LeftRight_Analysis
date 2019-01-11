#include "Include/helperFunctions.h"

void writeFile_Mmumu(TString start=""){
  //Setup_____
  TString whichProcess ="slot1"; //"Jpsi", "Psi", "OC", "t3", "slot1"
  TString whichMC ="Charles";

  Int_t nBins =200;
  Double_t Mmin =4.3;
  Double_t Mmax =8.5;
	
  Bool_t toWrite =true;
  //Setup_____

  if (start==""){
    cout << "\nMacro combines writes normalized invariant mass distribution ";
    cout << "to file ";
    cout << "Macro only works for Charles MC for now" << endl;
    cout << "Settings" << endl;
    cout << "nBins in ouput hist:   " << nBins << endl;
    cout << "Minimum mass range:   " << Mmin << endl;
    cout << "Maximum mass range:   " << Mmax << endl;
    cout << "Which MC component considered: " << whichProcess << endl;
    cout << "\nOutput is written to a file:   " << toWrite << endl;
    cout << "" << endl;
    exit(EXIT_FAILURE);
  }

  //Basic checks
  if (whichMC != "Charles"){
    cout << "\nOnly works for Charles MC for now....\n" << endl;
    exit(EXIT_FAILURE);
  }

  TString path, nameFile;
  if (whichProcess =="t3" || whichProcess =="slot1"){
    if (whichProcess =="slot1") nameFile ="slot1";
    nameFile +="WAll_LowM_AMDY.root";
    path ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Presents/DATA/RealData/LowM_AMDY/";
  }
  else {
    nameFile =Form("Charles_W12_%s.root", whichProcess.Data());
    path ="/Volumes/Seagate/DrellYan/Charles_Official/";
  }
  TFile *f = TFile::Open(path+nameFile);
  TTree *t = (TTree*) f->Get("pT_Weighted");

  TCanvas *c1 = new TCanvas(); 
  gPad->SetLogy();
  TH1D *h = new TH1D(Form("h_%s", whichProcess.Data()),
		     Form("h_%s", whichProcess.Data()),
		     nBins, Mmin, Mmax);
  t->Draw(Form("Mmumu>>h_%s", whichProcess.Data()),
	  "Mmumu>2.0&&Mmumu<8.5", "0");

  h->Sumw2();
  if (whichProcess =="t3" || whichProcess =="slot1"){
    cout << "Not normalizing real data" << endl;
  }
  else{	h->Scale(1.0/(h->Integral() ) ); }
  h->Draw();
  
  TString outName;
  if (whichProcess =="t3" || whichProcess =="slot1"){
    outName=Form("M_%s_%ibins%0.1f_%0.1f.root", whichProcess.Data(), nBins,
		 Mmin, Mmax);
  }
  else{
    outName=Form("M_%s_%s_%ibins%0.1f_%0.1f.root", whichMC.Data(),
		 whichProcess.Data(), nBins, Mmin, Mmax);
  }
  if (toWrite){
    TString thisPath =
      "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MC_bgDetermination/Data/MassDist/";
    cout << "     Full file path and name:  ";
    cout << thisPath + outName << endl;
    TFile *fOutput = 
      new TFile(thisPath+outName, "RECREATE");
    h->Write();
    fOutput->Close();
  }

  //Final output
  cout << "Settings" << endl;
  cout << "nBins in ouput hist:   " << nBins << endl;
  cout << "Minimum mass range:   " << Mmin << endl;
  cout << "Maximum mass range:   " << Mmax << endl;
  cout << "Which MC component considered: " << whichProcess << endl;
  cout << "\nFile  " << outName << "   was written:   " << toWrite << endl;
  cout << "" << endl;
}
