#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PhiScut/Include/helperFunctions.h"

void TreeLoop(TTree *t, TH1D** h1, TH2D** h2, TH1D *hProcess, Double_t *var[],
	      Double_t Mmin, Double_t Mmax,
	      Int_t *LR_true, Int_t *LR_false, Double_t *phiScut,
	      Int_t nPhiSCuts, Bool_t debug);

void resolutionPhiS(TString start=""){
  //Setup__________
  TString whichMC ="Yu"; //"Yu", "Charles"
  const Double_t Mmin=4.3, Mmax=8.5; //Mass cut
  //const Double_t Mmin=2.87, Mmax=3.38; //Mass cut
  Bool_t debug =false;
  
  Bool_t toWrite =false;
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
  const Int_t nFiles =5;
  TString dataPathNames[nFiles];
  if (whichMC=="Yu"){
    TString pathMC ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/MC_Data/YuShiangMC/";
    dataPathNames[0] = pathMC+"JPsi/Yu_Wall_full_main_JPsi_20bins.root";
    dataPathNames[1] = pathMC+"psi/Yu_Wall_full_main_psi_20bins.root";
    dataPathNames[2] = pathMC+"OC/Yu_Wall_full_main_OC_20bins.root";
    dataPathNames[3] = pathMC+"AMDY/Yu_Wall_full_main_AMDY_20bins.root";
  }
  else if (whichMC=="Charles"){
    TString pathMC ="/Volumes/Seagate/DrellYan/Charles_Official/";
    dataPathNames[0] = pathMC+"Charles_W12_Jpsi.root";
    dataPathNames[1] = pathMC+"Charles_W12_Psi.root";
    dataPathNames[2] = pathMC+"Charles_W12_OC.root";
    dataPathNames[3] = pathMC+"Charles_W12_LMDY.root";
    dataPathNames[4] = pathMC+"Charles_W12_HMDY.root";
  }
  else {
    cout << "MC specified does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  //Get files and trees
  TFile *f_JPsi = TFile::Open(dataPathNames[0]);
  TFile *f_psi = TFile::Open(dataPathNames[1]);
  TFile *f_OC = TFile::Open(dataPathNames[2]);
  TFile *f_AMDY = TFile::Open(dataPathNames[3]);
  TFile *f_HMDY = (whichMC=="Charles") ? TFile::Open(dataPathNames[4]) : NULL;
  TFile *Files[nFiles] = {f_JPsi, f_psi, f_OC, f_AMDY, f_HMDY};
  TString type[nFiles] = {"JPsi", "psi", "OC", "AMDY", "HMDY"};
  for (Int_t i=0; i<nFiles; i++) {
    if ( (i==nFiles-1)&&(whichMC!="Charles") ){ break; }
    
    if (!Files[i]) {
      cout << "File:   " << i << "   does not exist " << endl;
      exit(EXIT_FAILURE); }
  }
  TTree *t_JPsi = (TTree*) f_JPsi->Get("pT_Weighted");
  TTree *t_psi = (TTree*) f_psi->Get("pT_Weighted");
  TTree *t_OC = (TTree*) f_OC->Get("pT_Weighted");
  TTree *t_AMDY = (TTree*) f_AMDY->Get("pT_Weighted");
  TTree *t_HMDY =
    (whichMC=="Charles") ? (TTree*) f_HMDY->Get("pT_Weighted") : NULL;

  //Setup Tree variables
  Double_t PhiS, Gen_PhiS, Mmumu;
  t_JPsi->SetBranchAddress("PhiS", &PhiS);
  t_JPsi->SetBranchAddress("Gen_PhiS", &Gen_PhiS);
  t_JPsi->SetBranchAddress("Mmumu", &Mmumu);
  t_psi->SetBranchAddress("PhiS", &PhiS);
  t_psi->SetBranchAddress("Gen_PhiS", &Gen_PhiS);
  t_psi->SetBranchAddress("Mmumu", &Mmumu);
  t_OC->SetBranchAddress("PhiS", &PhiS);
  t_OC->SetBranchAddress("Gen_PhiS", &Gen_PhiS);
  t_OC->SetBranchAddress("Mmumu", &Mmumu);
  t_AMDY->SetBranchAddress("PhiS", &PhiS);
  t_AMDY->SetBranchAddress("Gen_PhiS", &Gen_PhiS);
  t_AMDY->SetBranchAddress("Mmumu", &Mmumu);
  if (whichMC=="Charles"){
    t_HMDY->SetBranchAddress("PhiS", &PhiS);
    t_HMDY->SetBranchAddress("Gen_PhiS", &Gen_PhiS);
    t_HMDY->SetBranchAddress("Mmumu", &Mmumu);
  }

  //Distributions
  TH1D* hPhiPhoton = new TH1D("hPhiPhoton", "hPhiPhoton", 100, -0.5*TMath::Pi(),
			      3.0*0.5*TMath::Pi());
  TH1D* hGenPhiPhoton = new TH1D("hGenPhiPhoton", "hGenPhiPhoton", 100,
				 -0.5*TMath::Pi(), 3.0*0.5*TMath::Pi());

  TH1D* hDiff = new TH1D("hDiff", "hDiff", 100, -6.28, 6.28);
  TH1D* hAMDY_Diff = new TH1D("hAMDY_Diff", "hAMDY_Diff", 100, -6.28, 6.28);
  TH1D* hHMDY_Diff = new TH1D("hHMDY_Diff", "hHMDY_Diff", 100, -6.28, 6.28);
  TH1D* hJPsi_Diff = new TH1D("hJPsi_Diff", "hJPsi_Diff", 100, -6.28, 6.28);
  TH1D* hPsi_Diff = new TH1D("hPsi_Diff", "hPsi_Diff", 100, -6.28, 6.28);
  TH1D* hOC_Diff = new TH1D("hOC_Diff", "hOC_Diff", 100, -6.28, 6.28);
  TH2D* hPhiDiff = new TH2D("hPhiDiff", "hPhiDiff",
			    100, -0.5*TMath::Pi(), 3.0*0.5*TMath::Pi(),
			    100, -6.28, 6.28);
  TH2D* hMmumuDiff = new TH2D("hMmumuDiff", "hMmumuDiff",
			      100, 1.0, 10.0,
			      100, -6.28, 6.28);
  
  TH1D* hLR = new TH1D("hLR", "hLR", 2, 0, 2);
  TH2D* hPhiLR = new TH2D("hPhiLR", "hPhiLR",
			  100, -0.5*TMath::Pi(), 3.0*0.5*TMath::Pi(), 2, 0, 2);

  //Setup for LR crossover percentage
  //Double_t resPhiS =0.133; //leading gaussian
  //Double_t resPhiS =0.303; //rms value AMDY only

  //Double_t resPhiS =0.2; //Charles rms value 2.0-8.5 debug mode
  Double_t resPhiS =0.1776; //Yu rms value 2.0-8.5
  //Double_t resPhiS =0.299; //Yu rms value 4.3-8.5
  
  const Int_t nPhiSCuts =9;
  Double_t phiScut[nPhiSCuts] =
    {0., 0.25*resPhiS, 0.5*resPhiS, resPhiS, 2*resPhiS, 3*resPhiS, 4*resPhiS,
     5*resPhiS, 6*resPhiS};
  Int_t LR_true[nPhiSCuts] ={0}, LR_false[nPhiSCuts] ={0};

  TH1D *h1All[] =
    {hPhiPhoton, hGenPhiPhoton, hLR, hDiff};
  TH2D *h2All[] =
    {hPhiLR, hPhiDiff, hMmumuDiff};
  Double_t *variables[] = {&PhiS, &Gen_PhiS, &Mmumu};

  //Treeloops
  if (Mmin < 4.3){
    cout << "only JPsi used dut to mass range" << endl;
      TreeLoop(t_JPsi, h1All, h2All, hJPsi_Diff, variables, Mmin, Mmax, LR_true,
	       LR_false, phiScut, nPhiSCuts, debug);
  }
  else{
    TreeLoop(t_AMDY, h1All, h2All, hAMDY_Diff, variables, Mmin, Mmax, LR_true,
	     LR_false, phiScut, nPhiSCuts, debug);
    TreeLoop(t_OC, h1All, h2All, hOC_Diff, variables, Mmin, Mmax, LR_true,
	     LR_false, phiScut, nPhiSCuts, debug);
    TreeLoop(t_psi, h1All, h2All, hPsi_Diff, variables, Mmin, Mmax, LR_true,
	     LR_false, phiScut, nPhiSCuts, debug);
    TreeLoop(t_JPsi, h1All, h2All, hJPsi_Diff, variables, Mmin, Mmax, LR_true,
	     LR_false, phiScut, nPhiSCuts, debug);
    if (whichMC=="Charles"){
      TreeLoop(t_HMDY, h1All, h2All, hHMDY_Diff, variables, Mmin, Mmax, LR_true,
	       LR_false, phiScut, nPhiSCuts, debug);
    }
  }

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
    Double_t lowBin =hMmumuDiff->GetXaxis()->FindBin(2+i-0.5);
    Double_t highBin =hMmumuDiff->GetXaxis()->FindBin(2+i+0.5);
    
    massRMS[i] = hMmumuDiff->ProjectionY("h1", lowBin, highBin)->GetRMS();
    massVal[i] = 2.0 + 1.0*i;
  }
  TGraph* gMdiff = new TGraph(6, massVal, massRMS);
  Setup(gMdiff);
  cMdiff->cd(2); gMdiff->Draw("AP");

  TCanvas* cLR = new TCanvas(); cLR->Divide(2);
  cLR->cd(1); hLR->Draw();
  cLR->cd(2); hPhiLR->Draw("colz");
  hPhiLR->GetXaxis()->SetTitle("PhiS Generated");
  Setup(hLR); Setup(hPhiLR);

  //crossover percent calculation
  Double_t crossOver[nPhiSCuts];
  for (Int_t i=0; i<nPhiSCuts; i++) {
    crossOver[i] = 1.0*LR_false[i]/(1.0*LR_true[i] + 1.0*LR_false[i]);
  }
  TGraph* gCrossOver = new TGraph(nPhiSCuts, phiScut, crossOver);
  Setup(gCrossOver);
  gCrossOver->GetXaxis()->SetTitle("PhiS Cut");

  TCanvas* cCrossOver = new TCanvas();
  gCrossOver->Draw("AP");

  TString outName = Form("Data/ResolutionPhiS/resolution_%s_%0.2f_%0.2f.root",
			 whichMC.Data(), Mmin, Mmax);
  if (toWrite){
    TFile* fOut = new TFile(outName, "RECREATE");
    for (Int_t i=0; i<4; i++) { h1All[i]->Write(); }
    for (Int_t i=0; i<3; i++) { h2All[i]->Write(); }
    hAMDY_Diff->Write();
    hOC_Diff->Write();
    hJPsi_Diff->Write();
    hPsi_Diff->Write();

    if (whichMC=="Charles"){
      hHMDY_Diff->Write();
    }
    
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


void TreeLoop(TTree *t, TH1D** h1, TH2D** h2, TH1D *hProcess, Double_t *var[],
	      Double_t Mmin, Double_t Mmax,
	      Int_t *LR_true, Int_t *LR_false, Double_t *phiScut,
	      Int_t nPhiSCuts, Bool_t debug){
  
  Int_t tree_entries = (debug) ? 10000 : t->GetEntries();
  if (debug) cout << "Debug mode is in effect" << endl;
  else cout << "Debug mode is   NOT   in effect" << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    t->GetEntry(ev, 0);

    if (*var[2] < Mmin) continue;//Mmumu cut
    if (*var[2] > Mmax) continue;

    Double_t phiPhoton = ShiftPhiSimple(*var[0]);
    Double_t Gen_phiPhoton = ShiftPhiSimple(*var[1]);

    h1[0]->Fill(phiPhoton);//PhiPhoton
    h1[1]->Fill(Gen_phiPhoton);//Gen_phiPhoton

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
    if ( Gen_Left==Left && Gen_Right==Right) LR_consistent =1;
    else LR_consistent =0;
    h1[2]->Fill(LR_consistent);//hLR
    h2[0]->Fill(Gen_phiPhoton, LR_consistent);//hPhiLR

    Double_t Diff = phiPhoton-Gen_phiPhoton;
    if (Diff > 3.14){ Diff = -2*TMath::Pi() + Diff; }
    else if (Diff < -3.14) { Diff = 2*TMath::Pi() + Diff; }
    h1[3]->Fill(Diff);//hDiff
    hProcess->Fill(Diff);
    h2[1]->Fill(Gen_phiPhoton, Diff);//hPhiDiff
    h2[2]->Fill(*var[2], Diff);//hMmumuDiff

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
}
