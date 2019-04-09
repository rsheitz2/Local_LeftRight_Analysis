#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PhiScut/Include/helperFunctions.h"

void resolutionHMDYPhiS(TString start=""){
  //Setup__________
  TString whichMC ="Charles"; //"Yu", "Charles"
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
  Bool_t debug =false;
  
  Bool_t toWrite =true;
  //Setup__________

  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);
  if (start==""){
    cout << "Macro compares " << endl;
    cout << "    reconstruction MC PhiS" << endl;
    cout << "to the same events for";
    cout << "    generated MC PhiS" << endl;
    cout << "to determine left/right cross over" << endl;
    cout << "\nCurrent Settings:" << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Which MC considered:          " << whichMC << endl;
    cout << "Debug mode in use:            " << debug << endl;
    cout << "\nTo Write:                   " << toWrite << endl;
    cout << "\n\nUsage: " << endl;
    cout << "root \'resolutionPhiS(1).C\'" << endl;
    exit(EXIT_FAILURE);
  }

  //Get data paths
  TString dataPathName;
  if (whichMC=="Yu"){
    TString pathMC ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/MC_Data/YuShiangMC/";
    dataPathName = pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root";
  }
  else if (whichMC=="Charles"){
    TString pathMC ="/Volumes/Seagate/DrellYan/Charles_Official/";
    dataPathName = pathMC+"Charles_W12_HMDY.root";
  }
  else {
    cout << "MC specified does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  //Get files and trees
  TFile *f_HMDY = OpenFile(dataPathName);
  TTree *t_HMDY = (TTree*) f_HMDY->Get("pT_Weighted");

  //Setup Tree variables
  Double_t PhiS, Gen_PhiS, Mmumu, q_transverse;
  t_HMDY->SetBranchAddress("PhiS", &PhiS);
  t_HMDY->SetBranchAddress("Gen_PhiS", &Gen_PhiS);
  t_HMDY->SetBranchAddress("Mmumu", &Mmumu);
  t_HMDY->SetBranchAddress("q_transverse", &q_transverse);

  //Distributions
  TH1D* hPhiPhoton = new TH1D("hPhiPhoton", "hPhiPhoton", 100, -TMath::Pi(),
			      TMath::Pi());
  TH1D* hGenPhiPhoton = new TH1D("hGenPhiPhoton", "hGenPhiPhoton", 100,
				 -TMath::Pi(), TMath::Pi());

  TH1D* hDiff = new TH1D("hDiff", "hDiff", 100, -6.28, 6.28);
  TH1D* hHMDY_Diff = new TH1D("hHMDY_Diff", "hHMDY_Diff", 100, -6.28, 6.28);
  TH2D* hPhiDiff = new TH2D("hPhiDiff", "hPhiDiff",
			    100, -TMath::Pi(), TMath::Pi(),
			    100, -6.28, 6.28);
  TH2D* hMmumuDiff = new TH2D("hMmumuDiff", "hMmumuDiff",
			      100, 1.0, 10.0,
			      100, -6.28, 6.28);
  
  TH1D* hLR = new TH1D("hLR", "hLR", 2, 0, 2);
  TH2D* hPhiLR = new TH2D("hPhiLR", "hPhiLR",
			  100, -TMath::Pi(), TMath::Pi(), 2, 0, 2);

  //Setup for LR crossover percentage
  //Double_t resPhiS =0.1776; //Yu rms value 2.0-8.5
  Double_t resPhiS =0.1869; //Charles rms value 4.3-8.5
    
  const Int_t nPhiSCuts =12;
  Double_t phiScut[nPhiSCuts] =
    {0., 0.044, 0.088, 0.17, 0.26, 0.36, 0.53, 0.71, 0.88, 1.07, 1.24, 1.44};
  /*Double_t phiScut[nPhiSCuts] =
    {0., 0.5*resPhiS, resPhiS, 2*resPhiS, 3*resPhiS, 4*resPhiS,
    5*resPhiS, 6*resPhiS};//*/
  Int_t LR_true[nPhiSCuts] ={0}, LR_false[nPhiSCuts] ={0};

  //Treeloops
  Int_t tree_entries = (debug) ? 10000 : t_HMDY->GetEntries();
  if (debug) cout << "\nDebug mode is in effect\n" << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    t_HMDY->GetEntry(ev, 0);

    if (Mmumu < Mmin) continue;//Mmumu cut
    if (Mmumu > Mmax) continue;

    Double_t phiPhoton = ShiftPhiSimple(PhiS);
    Double_t Gen_phiPhoton = ShiftPhiSimple(Gen_PhiS);

    hPhiPhoton->Fill(PhiS);//PhiPhoton
    hGenPhiPhoton->Fill(Gen_PhiS);//Gen_phiPhoton

    //Rec left/right
    Bool_t Left =false, Right =false;
    if ((phiPhoton < TMath::Pi()/2) && (phiPhoton > -TMath::Pi()/2) )
      Left = true;
    else if ((phiPhoton <3*TMath::Pi()/2) && (phiPhoton>TMath::Pi()/2) )
      Right = true;

    //Gen left/right
    Bool_t Gen_Left =false, Gen_Right =false;
    if ((Gen_phiPhoton < TMath::Pi()/2) && (Gen_phiPhoton > -TMath::Pi()/2) )
      Gen_Left = true;
    else if ((Gen_phiPhoton <3*TMath::Pi()/2) && (Gen_phiPhoton>TMath::Pi()/2) )
      Gen_Right = true;

    Int_t LR_consistent;
    if ( (Gen_Left==Left) && (Gen_Right==Right) ) LR_consistent =1;
    else LR_consistent =0;
    hLR->Fill(LR_consistent);//hLR
    hPhiLR->Fill(Gen_PhiS, LR_consistent);//hPhiLR

    Double_t Diff = phiPhoton-Gen_phiPhoton;
    if (Diff > 3.14){ Diff = -2*TMath::Pi() + Diff; }
    else if (Diff < -3.14) { Diff = 2*TMath::Pi() + Diff; }
    hDiff->Fill(Diff);//hDiff
    hHMDY_Diff->Fill(Diff);
    hPhiDiff->Fill(Gen_PhiS, Diff);//hPhiDiff
    hMmumuDiff->Fill(Mmumu, Diff);//hMmumuDiff

    //Crossover percent calculation
    for (Int_t i=0; i<nPhiSCuts; i++) {
      Bool_t goodRange =true;
      if ( (phiPhoton < TMath::Pi()/2 + phiScut[i]) &&
	   (phiPhoton > TMath::Pi()/2 - phiScut[i])) goodRange =false;
      if ( (phiPhoton < -TMath::Pi()/2 + phiScut[i]) ||
	   (phiPhoton > 3*TMath::Pi()/2 - phiScut[i])) goodRange =false;

      if (goodRange){
	if (LR_consistent) LR_true[i]++;
	else LR_false[i]++;
      }
    }
  }//End tree loop
  
  //Draw distributions
  TCanvas* cPhiPhoton = new TCanvas();
  hPhiPhoton->Draw();
  hGenPhiPhoton->Draw("sames");
  Setup(hPhiPhoton); Setup(hGenPhiPhoton); hGenPhiPhoton->SetLineColor(kRed);

  TCanvas* cDiff = new TCanvas(); cDiff->Divide(3);
  cDiff->cd(1); hDiff->Draw();
  Setup(hDiff); hDiff->SetTitle("Rec - Gen");
  cDiff->cd(2); hPhiDiff->Draw("colz");
  cDiff->cd(3); hPhiDiff->ProfileX()->Draw();
  Setup(hMmumuDiff); Setup(hPhiDiff);

  TCanvas *cMdiff = new TCanvas(); cMdiff->Divide(2);
  cMdiff->cd(1); hMmumuDiff->Draw("colz");
  Double_t massVal[6], massRMS[6];
  for (Int_t i=0; i<6; i++) {
    Double_t binWidth = (Mmax-Mmin)/6.0;
    Double_t lowBin =hMmumuDiff->GetXaxis()->FindBin(Mmin+(i-0.5)*binWidth );
    Double_t highBin =hMmumuDiff->GetXaxis()->FindBin(Mmin+(i+0.5)*binWidth );
    
    massRMS[i] = hMmumuDiff->ProjectionY("h1", lowBin, highBin)->GetRMS();
    massVal[i] = Mmin + binWidth*i;
  }
  TGraph* gMdiff = new TGraph(6, massVal, massRMS);
  Setup(gMdiff); gMdiff->SetTitle("Mass vs. PhiS Resolution");
  cMdiff->cd(2); gMdiff->Draw("AP");

  TCanvas* cLR = new TCanvas(); cLR->Divide(2);
  cLR->cd(1); hLR->Draw();
  cLR->cd(2); hPhiLR->Draw("colz");
  hPhiLR->GetXaxis()->SetTitle("PhiS Generated");
  Setup(hLR); Setup(hPhiLR);

  //crossover percent calculation
  Double_t crossOver[nPhiSCuts], totalPhiScut[nPhiSCuts];
  for (Int_t i=0; i<nPhiSCuts; i++) {
    crossOver[i] = 100.0*LR_false[i]/(1.0*LR_true[i] + 1.0*LR_false[i]);
    totalPhiScut[i] = 4*phiScut[i];
  }
  TGraph* gCrossOver = new TGraph(nPhiSCuts, totalPhiScut, crossOver);
  Setup(gCrossOver);
  gCrossOver->GetXaxis()->SetTitle("PhiS Cut");
  gCrossOver->SetTitle("Misidentification_Percent");

  TCanvas* cCrossOver = new TCanvas();
  gCrossOver->Draw("AP");

  TString outName = Form("Data/ResolutionPhiS/resolution_%s_%0.2f_%0.2f",
			 whichMC.Data(), Mmin, Mmax);
  if (debug) outName += ".root";
  else outName += "_noDebug.root";
  if (toWrite){
    TFile* fOut = new TFile(outName, "RECREATE");
    hPhiPhoton->Write();
    hGenPhiPhoton->Write();
    hLR->Write();
    hDiff->Write();

    hPhiLR->Write();
    hPhiDiff->Write();
    hMmumuDiff->Write();
    
    hHMDY_Diff->Write();
        
    gCrossOver->Write("gCrossOver");
    fOut->Close();
 } 

  //Final Output
  cout << "\nSettings !!!!" << endl;
  cout << "Mass Range is		:  " << MassCut << endl;
  cout << "Which MC considered:          " << whichMC << endl;
  cout << "Debug mode in use:            " << debug << endl;
  if (toWrite) { cout << outName << "   was written" << endl; }
  else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
