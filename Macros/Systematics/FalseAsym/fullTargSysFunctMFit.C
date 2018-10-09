#include "include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/AN_calculation/include/fitFunctions.h"


TF1* FitGetLR(TH1D **h, Int_t bin, Double_t *LR, Double_t *e_LR,
	      Double_t LR_Mmin, Double_t LR_Mmax, TString process,
	      TString whichFit, Double_t Mmin, Double_t Mmax,
	      TCanvas *c1, TH1D *hChi2){

  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_u", 11) == 0)
    c1->cd(2*bin+1);
  else if (strncmp(Form("%s", h[bin]->GetTitle()),"MuMu_left_d", 11)==0)
    c1->cd(2*bin+1);
  else
    c1->cd(2*bin+2);

  gPad->SetLogy();
  TF1 *fitFunc = FitGetLR(h, bin, LR, e_LR,LR_Mmin, LR_Mmax, process,
			  whichFit, Mmin, Mmax);

  //Get reduced Chi2
  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();
  Double_t redChi2 = Chi2/ndf;
  hChi2->Fill(redChi2);

  return fitFunc;
}


void fullTargSysFunctMFit(TString start=""){
  
  //Setup_______________
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;
  TString binRange ="25_43";
  TString physBinned ="xN";//"xN", "xPi", "xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.90;
  Double_t LR_Mmax =3.30;//L/R counts mass range
  Double_t Mmin =2.00;//Fit Mass minimum
  Double_t Mmax =7.50;//Fit Mass maximum
  TString whichFit ="ten";
  
  Bool_t toWrite =false;
  //Setup_______________

  Double_t nominal_Mmin =Mmin, nominal_Mmax =Mmax;
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =Form("systematic_leftRight_%s1.00_8.50_%ibins%s_%ihbin.root",
		       period_Mtype.Data(), nBins, binRange.Data(), hbins);
  
  if (start==""){
    cout<<"Script outputs AN and left/right counts per target and polarization";
    cout << " using functional mass fitting for a given fit" << endl;
    cout << "Outputs are in the formate needed for GeoMean4Targ.C" << endl;
    cout << "Script also outputs information on the fit quality" <<endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'sysFunctMFit.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Real data path:             " << pathRD << endl;
    cout << "Real data file considered:  " << RDfile << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << LR_Mmin << "  -  "
	 << LR_Mmax << endl;
    cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //Get Input Files from Local_leftRight
  TFile *fRD  = TFile::Open(pathRD + RDfile);
  if (!fRD ){
    cout << "RD file does not exist " << endl;
    exit(EXIT_FAILURE);
  }
  
  //Get Input Hist/Determine LR counts for specified process
  TH1D *hRD_upstream_up_L[nBins], *hRD_upstream_up_R[nBins];
  TH1D *hRD_downstream_up_L[nBins], *hRD_downstream_up_R[nBins];
  TH1D *hRD_upstream_down_L[nBins], *hRD_upstream_down_R[nBins];
  TH1D *hRD_downstream_down_L[nBins], *hRD_downstream_down_R[nBins];
  
  //L/R counts for specified process
  Double_t LR_upstream_up_L[nBins], LR_upstream_up_R[nBins];
  Double_t LR_downstream_up_L[nBins], LR_downstream_up_R[nBins];
  Double_t LR_upstream_down_L[nBins], LR_upstream_down_R[nBins];
  Double_t LR_downstream_down_L[nBins], LR_downstream_down_R[nBins];
  
  Double_t e_LR_upstream_up_L[nBins], e_LR_upstream_up_R[nBins];
  Double_t e_LR_downstream_up_L[nBins], e_LR_downstream_up_R[nBins];
  Double_t e_LR_upstream_down_L[nBins], e_LR_upstream_down_R[nBins];
  Double_t e_LR_downstream_down_L[nBins], e_LR_downstream_down_R[nBins];
    
  const Int_t nTargPol =4;
  TCanvas* cFit[nTargPol];
  TString targNames[nTargPol] = {"upstream_up", "downstream_down",
				 "upstream_down", "downstream_up"};
  for (Int_t c=0; c<nTargPol; c++) {
    cFit[c] = new TCanvas(targNames[c]); cFit[c]->Divide(2, nBins); }
  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 30, 0, 3); SetUp(hChi2);
  //Perform Fits and integrate to get L/R
  for (Int_t bi=0; bi<nBins; bi++) {
    hRD_upstream_up_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s%i",
					       physBinned.Data(), bi) );
    hRD_upstream_up_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_up_%s%i",
					       physBinned.Data(), bi) );
    SetUp(hRD_upstream_up_L[bi]); SetUp(hRD_upstream_up_R[bi]);
    
    hRD_upstream_down_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_upstream_down_%s%i",
			   physBinned.Data(), bi) );
    hRD_upstream_down_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_upstream_down_%s%i",
			   physBinned.Data(), bi) );
    SetUp(hRD_upstream_down_L[bi]); SetUp(hRD_upstream_down_R[bi]);
    
    hRD_downstream_up_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_downstream_up_%s%i",
			   physBinned.Data(), bi) );
    hRD_downstream_up_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_downstream_up_%s%i",
			   physBinned.Data(), bi) );
    SetUp(hRD_downstream_up_L[bi]); SetUp(hRD_downstream_up_R[bi]);

    hRD_downstream_down_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_downstream_down_%s%i",
			   physBinned.Data(), bi) );
    hRD_downstream_down_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_downstream_down_%s%i",
			   physBinned.Data(), bi) );
    SetUp(hRD_downstream_down_L[bi]); SetUp(hRD_downstream_down_R[bi]);

    FitGetLR(hRD_upstream_up_L, bi, LR_upstream_up_L, e_LR_upstream_up_L,
	     LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
	     cFit[0], hChi2);
    FitGetLR(hRD_upstream_up_R, bi, LR_upstream_up_R, e_LR_upstream_up_R,
	     LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
	     cFit[0], hChi2);

    FitGetLR(hRD_upstream_down_L, bi, LR_upstream_down_L,
	     e_LR_upstream_down_L, LR_Mmin, LR_Mmax, process, whichFit,
	     Mmin, Mmax, cFit[1], hChi2);
    FitGetLR(hRD_upstream_down_R, bi, LR_upstream_down_R,
	     e_LR_upstream_down_R, LR_Mmin, LR_Mmax, process, whichFit,
	     Mmin, Mmax, cFit[1], hChi2);
    
    FitGetLR(hRD_downstream_up_L, bi, LR_downstream_up_L,
	     e_LR_downstream_up_L, LR_Mmin, LR_Mmax, process, whichFit,
	     Mmin, Mmax, cFit[2], hChi2);
    FitGetLR(hRD_downstream_up_R, bi, LR_downstream_up_R,
	     e_LR_downstream_up_R, LR_Mmin, LR_Mmax, process, whichFit,
	     Mmin, Mmax, cFit[2], hChi2);

    FitGetLR(hRD_downstream_down_L, bi, LR_downstream_down_L,
	     e_LR_downstream_down_L, LR_Mmin, LR_Mmax, process, whichFit,
	     Mmin, Mmax, cFit[3], hChi2);
    FitGetLR(hRD_downstream_down_R, bi, LR_downstream_down_R,
	     e_LR_downstream_down_R, LR_Mmin, LR_Mmax, process, whichFit,
	     Mmin, Mmax, cFit[3], hChi2);
  }//Loop over nBins physics binning
  
  //Draw hChi2
  TCanvas* cChi2 = new TCanvas("Chi2");
  hChi2->Draw();
  
  //L/R count graphs/Drawing
  TGraphErrors *g_count_upstream_up =
    (TGraphErrors*)fRD->Get(Form("%s_left_upSup_upP", physBinned.Data() ));
  Double_t *xvals = g_count_upstream_up->GetX();
  Double_t ex[nBins] = {0.0};
  TGraphErrors* g_Left_upstream_up =
    new TGraphErrors(nBins, xvals, LR_upstream_up_L, ex, e_LR_upstream_up_L);
  TGraphErrors* g_Right_upstream_up =
    new TGraphErrors(nBins, xvals, LR_upstream_up_R, ex, e_LR_upstream_up_R);

  TGraphErrors* g_Left_upstream_down =
    new TGraphErrors(nBins, xvals, LR_upstream_down_L, ex, e_LR_upstream_down_L);
  TGraphErrors* g_Right_upstream_down =
    new TGraphErrors(nBins, xvals, LR_upstream_down_R, ex, e_LR_upstream_down_R);

  TGraphErrors* g_Left_downstream_up =
    new TGraphErrors(nBins, xvals, LR_downstream_up_L, ex, e_LR_downstream_up_L);
  TGraphErrors* g_Right_downstream_up =
    new TGraphErrors(nBins, xvals, LR_downstream_up_R, ex, e_LR_downstream_up_R);

  TGraphErrors* g_Left_downstream_down =
    new TGraphErrors(nBins, xvals, LR_downstream_down_L, ex,
		     e_LR_downstream_down_L);
  TGraphErrors* g_Right_downstream_down =
    new TGraphErrors(nBins, xvals, LR_downstream_down_R, ex,
		     e_LR_downstream_down_R);
  
  SetUp(g_Left_upstream_up); SetUp(g_Right_upstream_up);
  SetUp(g_Left_upstream_down); SetUp(g_Right_upstream_down);
  SetUp(g_Left_downstream_up); SetUp(g_Right_downstream_up);
  SetUp(g_Left_downstream_down); SetUp(g_Right_downstream_down);

  TCanvas* cLR = new TCanvas("LR counts"); cLR->Divide(2);
  cLR->cd(1); g_Left_upstream_up->Draw("AP");//Upstream
  g_Left_upstream_up->SetTitle("L/R upstream");
  g_Right_upstream_up->Draw("Psame");
  g_Left_upstream_down->Draw("Psame");
  g_Right_upstream_down->Draw("Psame");
  g_Right_upstream_up->SetMarkerColor(kRed);
  g_Left_upstream_down->SetMarkerColor(kBlue);
  g_Right_upstream_down->SetMarkerColor(kGreen);

  cLR->cd(2); g_Left_downstream_up->Draw("AP");//Downstream
  g_Left_downstream_up->SetTitle("L/R downstream");
  g_Right_downstream_up->Draw("Psame");
  g_Left_downstream_down->Draw("Psame");
  g_Right_downstream_down->Draw("Psame");
  g_Right_downstream_up->SetMarkerColor(kRed);
  g_Left_downstream_down->SetMarkerColor(kBlue);
  g_Right_downstream_down->SetMarkerColor(kGreen);
    
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
fullTargSysFunctMFit";
  TString fOutput =
    Form("%s/fullTargSysFunctMFit_%s%.2f_%.2f_%s_%s%.2f_%.2f_%s%i_%ihbin.root",
	 thisDirPath.Data(), whichFit.Data(), nominal_Mmin, nominal_Mmax,
	 period_Mtype.Data(), process.Data(), LR_Mmin, LR_Mmax,
	 physBinned.Data(), nBins, hbins);
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    
    g_Left_upstream_up->Write(Form("%s_left_upstream_up", physBinned.Data()));
    g_Right_upstream_up->Write(Form("%s_right_upstream_up", physBinned.Data()));
    g_Left_upstream_down->Write(Form("%s_left_upstream_down",
				     physBinned.Data()));
    g_Right_upstream_down->Write(Form("%s_right_upstream_down",
				      physBinned.Data()));

    g_Left_downstream_up->Write(Form("%s_left_downstream_up",
				     physBinned.Data()));
    g_Right_downstream_up->Write(Form("%s_right_downstream_up",
				      physBinned.Data()));
    g_Left_downstream_down->Write(Form("%s_left_downstream_down",
				       physBinned.Data()));
    g_Right_downstream_down->Write(Form("%s_right_downstream_down",
					physBinned.Data()));
  }

  cout << " " << endl;
  cout << "Settings______" << endl;
  cout << "Real data file considered:              " << RDfile << endl;
  cout <<"Number of histogram bins used in M fit:  " << hbins << endl;
  cout << "Mass Range for fitting is:              " << Mmin
       << " - " << Mmax << endl;
  cout << "physBinned nBins times:                 " << nBins << endl;
  cout << "AN for physical process:                " << process << endl;
  cout << "Physics binning is:                     " << physBinned << endl;
  cout << "LR integral mass range:                 " << LR_Mmin
       << " - " << LR_Mmax << endl;
  cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
  cout << "Which fit considered:       " << whichFit << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
