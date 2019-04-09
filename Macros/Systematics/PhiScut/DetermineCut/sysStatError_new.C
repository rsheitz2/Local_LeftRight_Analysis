#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PhiScut/Include/helperFunctions.h"

void sysStatError(TString start=""){
  //Setup__________
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
  //const Double_t Mmin=2.87, Mmax=3.38; //Mass cut
  
  TString whichMC ="Charles"; //"Yu", "Charles"
  TString whichRD ="slot1"; //"t3", "slot1"
    
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
TGeant/Local_LeftRight_Analysis/Macros/";
  TString pathSys = path+"Systematics/PhiScut/DetermineCut/";
  TString pathStat = path+"FinalAsym/Data/WAvg/";
  TString nameSys =
    Form("Data/ResolutionPhiS/resolution_%s_%0.2f_%0.2f_noDebug.root",
	 whichMC.Data(), Mmin, Mmax);
  TString nameStat =
    Form("wAvg_true_LowM_AMDY_JPsi%0.2f_%0.2f_29_34xN1_150hbin_slot1_phiS0.0.root",
	 Mmin, Mmax);
  cout << "Sys error file: " << nameSys << "    Stat error file: "
       << nameStat << endl;
  TFile *fsys = OpenFile(pathSys+nameSys);
  TFile *fstat = OpenFile(pathStat+nameStat);

  TGraph *gCrossOver = (TGraph*)fsys->Get("gCrossOver");
  TGraph *gTotStatError = (TGraph*)fstat->Get("AN");

  Double_t *yCrossOver = gCrossOver->GetY();
  Double_t *yTotStatError= gTotStatError->GetY();
  const Int_t nPhiScut =1;
  Double_t ySys_PhiS[nPhiScut], totError_PhiS[nPhiScut], yStat_PhiS[nPhiScut];
  Double_t yFinalSys[1], xFinalSys[1] ={1.5};
  Double_t phiScut[nPhiScut] =
    {0.};
  
  for (Int_t i=0; i<nPhiScut; i++) {
    ySys_PhiS[i] =
      yCrossOver[i]*(gTotStatError->GetY()[0]+gTotStatError->GetEY()[0]);
    cout << yCrossOver[i] << " " << gTotStatError->GetY()[0] << " " << gTotStatError->GetEY()[0] << endl;
    
    totError_PhiS[i] =
      TMath::Sqrt(ySys_PhiS[i]*ySys_PhiS[i] + yStat_PhiS[i]*yStat_PhiS[i] );

    if (i==0){
      cout << "\nFinal systematics for phiS cut:  " << phiScut[i] << endl;
      cout << "For file naming, phiS cut should equal:  " << endl;
      cout << " " << endl;
      yFinalSys[0] = ySys_PhiS[i];
    }
  }

  Double_t totalPhiScut[nPhiScut];
  for (Int_t i=0; i<nPhiScut; i++) { totalPhiScut[i] = 4*phiScut[i]; }
  TGraph *gSys_PhiS = new TGraph(nPhiScut, totalPhiScut, ySys_PhiS);
  TGraph *gTotError_PhiS = new TGraph(nPhiScut, totalPhiScut, totError_PhiS);
  TGraph *gStat_PhiS = new TGraph(nPhiScut, totalPhiScut, yStat_PhiS);
  Setup(gSys_PhiS); Setup(gTotError_PhiS); Setup(gStat_PhiS);

  cout << yFinalSys[0] << endl;
  TGraph *gSys = new TGraph(1, xFinalSys, yFinalSys);
  Setup(gSys);

  TCanvas* c1 = new TCanvas(); c1->Divide(2);
  c1->cd(1);
  gStat_PhiS->SetMarkerColor(kBlue);
  gTotError_PhiS->SetMarkerColor(kRed);
  gTotError_PhiS->Draw("AP");
  gSys_PhiS->Draw("Psame");
  gStat_PhiS->Draw("Psame");
  gTotError_PhiS->SetTitle("Errors");

  c1->cd(2);
  gSys->Draw("AP"); gSys->SetTitle("FinalSystematics");

  TString outName =
    Form("Data/SysStatError/sysStatError_%s_%s_%0.2f_%0.2f_phiS0.0.root",
	 whichMC.Data(), whichRD.Data(), Mmin, Mmax );
  if (toWrite){
    TFile* fOut = new TFile(outName, "RECREATE");
    gSys_PhiS->Write("gSys_PhiS");
    gStat_PhiS->Write("gStat_PhiS");
    gTotError_PhiS->Write("gTotError_PhiS");
    gSys->Write("gSys");
    
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
