#include "include/helperFunctions.h"

Bool_t BinFill(TH1D *h[], Double_t binVar, Double_t *bounds,
	       Double_t fillVal, Int_t nBins);

Bool_t BinAdd(Double_t *Count, Double_t binVar, Double_t *bounds,
	      Double_t fillVal, Int_t nBins);

Double_t ErrorCal(Int_t n_UpS_up, Int_t n_DownS_up,
		  Int_t n_UpS_down, Int_t n_DownS_down);

void doubleRatio(TString start=""){
  //Setup_______________
  const Int_t nBins =1;//# of physBinned bins
  const Int_t nHbins =8;
  TString period ="WAll";
  TString Mtype ="HMDY";
  TString physBinned ="xPi";//"xN", "xPi", "xF", "pT", "M"
  TString production ="slot1";//"t3", "slot1"
  TString whichTSA ="Trans";//"Siv", "Trans", "Pretz"
  
  Bool_t toWrite =false;
  //Setup_______________

  //More setup
  TString Mrange, binRange;
  if (Mtype =="HMDY"){ Mrange ="4.30_8.50";  binRange ="43_85"; }
  else{ cout << Mtype << "not defined well" << endl; exit(EXIT_FAILURE); }

  //Input file paths/names
  TString pathRD ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData";
  if (production =="t3") production ="";
  TString RDfile =Form("%s/%s/%s%s_%s.root", pathRD.Data(), Mtype.Data(),
		       production.Data(), period.Data(), Mtype.Data() );
  if (production =="") production ="t3";

  TString pathLR ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Data";
  TString LRfile =
    Form("%s/leftRight_byTarget_%s_%s%s_%ibins%s_150hbin_%s_phiS0.0.root",
	 pathLR.Data(), period.Data(), Mtype.Data(), Mrange.Data(), nBins,
	 binRange.Data(), production.Data());

  if (production =="t3") production ="";
  TString binfile =
    Form("/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData/%s/BinValues/%sWAll_%s_%ibins.txt",
	 Mtype.Data(), production.Data(), Mtype.Data(), nBins);
  if (production =="") production ="t3";

  //Settings
  if (start==""){
    cout << "\nCurrent settings:" << endl;
    cout << "Real data path:             " << pathRD << endl;
    cout << "Real data file considered:  " << RDfile << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "Polarizations from file:     " << LRfile << endl;
    cout << "TSA determined:     " << whichTSA << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //Get Input Files 
  TFile *fRD  = OpenFile(RDfile);
  TTree *tRD = (TTree*)fRD->Get("pT_Weighted");

  //Get Polarization
  TFile *fLR = OpenFile(LRfile);
  TGraph* g_Pol =(TGraph*)fLR->Get(Form("%s_Pol", physBinned.Data()));
  TGraph* g_Dil =(TGraph*)fLR->Get(Form("%s_Dil", physBinned.Data()));
    
  Double_t Pol[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);

  //X-axis and bounds
  Double_t *xPhysVal = g_Pol->GetX();
  Double_t bounds[nBins+1];
  if (physBinned =="M") physBinned ="mass";
  GetBounds(physBinned, binfile, bounds);
  if (physBinned =="mass") physBinned ="M";
  
  //Setup tree
  Double_t PhiS_simple, Phi_CS, Theta_CS, Spin_0;
  Int_t targetPosition;
  Double_t x_target, x_beam, x_feynman, q_transverse, Mmumu;
  tRD->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tRD->SetBranchAddress("Phi_CS", &Phi_CS);
  tRD->SetBranchAddress("Theta_CS", &Theta_CS);
  tRD->SetBranchAddress("targetPosition", &targetPosition);
  tRD->SetBranchAddress("Spin_0", &Spin_0);
  tRD->SetBranchAddress("x_target", &x_target);
  tRD->SetBranchAddress("x_beam", &x_beam);
  tRD->SetBranchAddress("x_feynman", &x_feynman);
  tRD->SetBranchAddress("q_transverse", &q_transverse);
  tRD->SetBranchAddress("Mmumu", &Mmumu);

  Double_t *physVal;
  if ( physBinned == "xN" ) physVal = &x_target;
  else if ( physBinned == "xPi" ) physVal = &x_beam;
  else if ( physBinned == "xF" ) physVal = &x_feynman;
  else if ( physBinned == "pT" ) physVal = &q_transverse;
  else if ( physBinned == "M" ) physVal = &Mmumu;

  //Setup histograms
  TH1D *h_UpS_up[nBins], *h_UpS_down[nBins];
  TH1D *h_DownS_up[nBins], *h_DownS_down[nBins];
  TH1D *h_Ratio[nBins];
  for (Int_t i=0; i<nBins; i++) {
    Double_t xmax;
    if ( whichTSA == "Siv" ) xmax = TMath::Pi();
    else if ( whichTSA == "Pretz" ) xmax = 3*TMath::Pi();
    else if ( whichTSA == "Trans" ) xmax = 3*TMath::Pi();
    else{
      cout << "Wrong TSA input" << endl;
      exit(EXIT_FAILURE);
    }
    
    h_UpS_up[i] = new TH1D(Form("UpS_up_%i", i), Form("UpS_up_%i", i),
			   nHbins, -xmax, xmax);
    h_UpS_down[i] = new TH1D(Form("UpS_down_%i", i), Form("UpS_down_%i", i),
			     nHbins, -xmax, xmax);

    h_DownS_up[i] = new TH1D(Form("DownS_up_%i", i), Form("DownS_up_%i", i),
			     nHbins, -xmax, xmax);
    h_DownS_down[i] = new TH1D(Form("DownS_down_%i", i),Form("DownS_down_%i",i),
			       nHbins, -xmax, xmax);

    h_Ratio[i] = new TH1D(Form("Ratio_%i", i),Form("Ratio_%i",i),
			  nHbins, -xmax, xmax);
  }

  //Depolarization factor
  Double_t Depol[nBins] ={0.0}, binCounts[nBins] ={0.0};
  
  //Loop over tree
  Int_t treeEntries = tRD->GetEntries();
  cout << "\nNumber of entries considered:   " << treeEntries << "\n" << endl;
  for (Int_t ev=0; ev<treeEntries; ev++) {
    tRD->GetEntry(ev);

    Double_t depolarization, angle;
    if ( whichTSA == "Siv" ) {
      depolarization = 1.0;
      angle = PhiS_simple;
    }
    else if ( whichTSA == "Trans" ) {
      depolarization =
	TMath::Sin(Theta_CS)/(1+TMath::Cos(Theta_CS)*TMath::Cos(Theta_CS));
      angle = 2*Phi_CS - PhiS_simple;
    }
    else if ( whichTSA == "Pretz" ) {
      depolarization =
	TMath::Sin(Theta_CS)/(1+TMath::Cos(Theta_CS)*TMath::Cos(Theta_CS));
      angle = 2*Phi_CS + PhiS_simple;
    }
    else{
      cout << "wrong TSA input" << endl;
      exit(EXIT_FAILURE);
    }

    //Depolarization average
    BinAdd(Depol, *physVal, bounds, depolarization, nBins);
    BinAdd(binCounts, *physVal, bounds, 1.0, nBins);

    if (targetPosition == 0){//upstream
      if (Spin_0 >0){//spin up
	BinFill(h_UpS_up, *physVal, bounds, angle, nBins);
      }
      else{//spin down
	BinFill(h_UpS_down, *physVal, bounds, angle, nBins);
      }
    }
    else{
      if (Spin_0 >0){//spin up
	BinFill(h_DownS_up, *physVal, bounds, angle, nBins);
      }
      else{//spin down
	BinFill(h_DownS_down, *physVal, bounds, angle, nBins);
      }
    }

  }//tree loop

  //Draw distributions
  for (Int_t bi=0; bi<nBins; bi++) {
    for (Int_t hi=1; hi<nHbins+1; hi++) {
      Int_t n_UpS_up =h_UpS_up[bi]->GetBinContent(hi);
      Int_t n_UpS_down =h_UpS_down[bi]->GetBinContent(hi);
      Int_t n_DownS_up =h_DownS_up[bi]->GetBinContent(hi);
      Int_t n_DownS_down =h_DownS_down[bi]->GetBinContent(hi);

      Double_t Ratio = 1.0*n_UpS_up*n_DownS_up/(n_UpS_down*n_DownS_down);
      h_Ratio[bi]->SetBinContent(hi, Ratio);
      
      Double_t eRatio =ErrorCal(n_UpS_up, n_DownS_up, n_UpS_down, n_DownS_down);
      h_Ratio[bi]->SetBinError(hi, eRatio);
    }
  }
  
  TCanvas* cRatio = new TCanvas(); cRatio->Divide(nBins);
  Double_t Amp[nBins], eAmp[nBins];
  for (Int_t bi=0; bi<nBins; bi++) {
    cRatio->cd(bi+1); h_Ratio[bi]->Draw("e");

    TF1 *fit;
    if ( whichTSA == "Siv" )
      fit = new TF1("fit", "[0]*(1+4*[1]*sin(x))", -TMath::Pi(),TMath::Pi());
    else if ( whichTSA == "Pretz" )
      fit = new TF1("fit", "[0]*(1+4*[1]*sin(x))",-3*TMath::Pi(),3*TMath::Pi());
    else if ( whichTSA == "Trans" )
      fit = new TF1("fit", "[0]*(1+4*[1]*sin(x))",-3*TMath::Pi(),3*TMath::Pi());
    else{
      cout << "Wrong TSA input" << endl;
      exit(EXIT_FAILURE);
    }
    fit->SetParameters(1.0, 0.0);
    
    if ( h_Ratio[bi]->Fit("fit", "RQ") ){
      cout << period << " bin " << bi << " phys binned " << physBinned;
      cout << "  Fit failed" << endl;
      exit(EXIT_FAILURE); 
    }

    //Depolarization average
    Depol[bi] /= binCounts[bi];

    Double_t *pars =fit->GetParameters();
    Amp[bi] = pars[1];
    Amp[bi] /= (Pol[bi]);
    Amp[bi] /= (Depol[bi]);

    const Double_t *e_pars =fit->GetParErrors();
    eAmp[bi] = e_pars[1];
    eAmp[bi] /= (Pol[bi]);
    eAmp[bi] /= (Depol[bi]);
  }

  Double_t ex[nBins] ={0.0};
  TGraphErrors* gAmp = new TGraphErrors(nBins, xPhysVal, Amp, ex, eAmp);
  SetUp(gAmp);

  TCanvas* cAmp = new TCanvas();
  gAmp->Draw("AP"); DrawLine(gAmp, 0.0);
  gAmp->SetTitle("TSA Amplitude");

  //Write output/Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/DoubleRatio/Data/DoubleRatio";
  TString fOutput =
    Form("%s/doubleRatio_%s_%s_%s%i_%ihbins_%s_%s.root", thisDirPath.Data(),
	 period.Data(), Mtype.Data(), physBinned.Data(), nBins, nHbins,
	 production.Data(), whichTSA.Data());
  if (toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    gAmp->Write("Amp");

    for (Int_t bi=0; bi<nBins; bi++) {
      h_UpS_up[bi]->Write(); h_UpS_down[bi]->Write();
      h_DownS_up[bi]->Write(); h_DownS_down[bi]->Write();

      h_Ratio[bi]->Write(Form("Ratio_%i", bi));
    }

  }
  cout << "\nSettings______" << endl;
  cout << "Real data path:             " << pathRD << endl;
  cout << "Real data file considered:  " << RDfile << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Fit in nHbins:              " << nHbins << endl;
  cout << "Mass range considered:      " << Mtype << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
  cout << "Production considered:      " << production << endl;
  cout << "TSA determined:     " << whichTSA << endl;
  cout << "Polarizations from file:     " << LRfile << endl;
  cout << "Bin file used:               " << binfile << endl;
  cout << "\nTo write output file:       " << toWrite << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}

Bool_t BinFill(TH1D *h[], Double_t binVar, Double_t *bounds,
	       Double_t fillVal, Int_t nBins){
  for (Int_t bi=0; bi<nBins; bi++) {
    if (binVar < bounds[bi+1]) {

      h[bi]->Fill(fillVal);
      return true;
    }
  }

  return false;
}//BinFill

Bool_t BinAdd(Double_t *Count, Double_t binVar, Double_t *bounds,
	       Double_t fillVal, Int_t nBins){
  for (Int_t bi=0; bi<nBins; bi++) {
    if (binVar < bounds[bi+1]) {

      Count[bi] += fillVal;
      return true;
    }
  }

  return false;
}//BinAdd

Double_t ErrorCal(Int_t n_UpS_up, Int_t n_DownS_up,
		  Int_t n_UpS_down, Int_t n_DownS_down){
  
  Double_t Ratio = 1.0*n_UpS_up*n_DownS_up/(n_UpS_down*n_DownS_down);
  
  Double_t error = 1.0/n_UpS_up + 1.0/n_DownS_up;
  error += 1.0/n_UpS_down + 1.0/n_DownS_down;
  error = TMath::Sqrt(error);

  return error*Ratio;
}
