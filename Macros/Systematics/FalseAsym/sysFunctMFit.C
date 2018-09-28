#include "include/helperFunctions.h"

//2 Gaussian w/ psi' M/W = A*JPsi M/W by target
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_six.h"

//2 Gaussian w/ psi' M/W = A*JPsi M/W by target
//1 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_seven.h"

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_eight.h"

//2 Crystal Ball w/ psi' M/W = A*JPsi M/W by target
//1 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_nine.h"


void DoFit(TH1D *h, TMatrixDSym &cov, Double_t Mmin, Double_t Mmax){
  h->Sumw2(); 
  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ){
    cout << h->GetTitle() << "    Fit failed Once!!" << endl;
    
    status = h->Fit("fitFunc", "RISQ", "", Mmin, Mmax+0.05);
    if (status->Status() ){
      cout << h->GetTitle() << "  Fit failed Twice!!!\n" << endl;
      //exit(EXIT_FAILURE);
    }
  }

  cov = status->GetCovarianceMatrix();
}


Double_t FitGetPars(TH1D **h, Int_t bin,
		    Double_t *LR, Double_t *e_LR,
		    Double_t LR_Mmin, Double_t LR_Mmax, TString process,
		    TString whichFit, Double_t Mmin, Double_t Mmax){
  
  Double_t processPars[8] = {0.0};//{aJPsi,mJPsi,wJPsi,Apsi,Mpsi,Wpsi,aDY,cDY}
  Double_t LR_cov[25] = {0.0};
  Double_t binWidth = h[bin]->GetBinWidth(1);
  
  Bool_t hIsUpS =false;
  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_upstream", 18) == 0)
    hIsUpS =true;
  else if (strncmp(Form("%s", h[bin]->GetTitle() ),"MuMu_right_upstream",19)==0)
    hIsUpS =true;

  //Fit Setups_____
  TF1 *fitFunc = NULL;
  Int_t nPar;
  if (whichFit =="six"){
    fitFunc = SetupFunc_six(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_six(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    TF1 *f_LR =ComponentFuncts_six(processPars, Mmin, Mmax, process);

    IntegrateLR_six(f_LR, processPars, LR_cov, binWidth, LR_Mmin, LR_Mmax,
		    &(LR[bin]), &(e_LR[bin]), Mmin, Mmax, process);
  }
  else if (whichFit =="seven"){
    fitFunc = SetupFunc_seven(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_seven(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);

    TF1 *f_LR =ComponentFuncts_seven(processPars, Mmin, Mmax, process);

    IntegrateLR_seven(f_LR, processPars, LR_cov, binWidth, LR_Mmin, LR_Mmax,
		      &(LR[bin]), &(e_LR[bin]), Mmin, Mmax, process);
  }
  else if (whichFit =="eight"){
    fitFunc = SetupFunc_eight(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);

    ProcessPars_eight(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
    
    Double_t *pars = fitFunc->GetParameters();
    TF1 *f_LR =ComponentFuncts_eight(pars, Mmin, Mmax, process, hIsUpS);
    IntegrateLR_eight(f_LR, processPars, LR_cov, binWidth,
		      (LR_Mmax - LR_Mmin)/2.0, &(LR[bin]), &(e_LR[bin]), hIsUpS,
		      process);
  }
  else if (whichFit =="nine"){
    fitFunc = SetupFunc_nine(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_nine(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);

    Double_t *pars = fitFunc->GetParameters();
    TF1 *f_LR =ComponentFuncts_nine(pars, Mmin, Mmax, process);
    IntegrateLR_nine(f_LR, pars, LR_cov, binWidth, LR_Mmin, LR_Mmax,
		    &(LR[bin]), &(e_LR[bin]), Mmin, Mmax, process);    
  }
  else{
    cout << "Int valid fit type:   " << whichFit << endl;
    exit(EXIT_FAILURE);
  }

  //Get reduced Chi2
  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();
  return Chi2/ndf;
}


void WrapperFitGetPars(TH1D **h, Int_t bin,
		       Double_t *LR, Double_t *e_LR,
		       Double_t LR_Mmin, Double_t LR_Mmax, TString process,
		       TString whichFit, Double_t Mmin, Double_t Mmax,
		       TCanvas *c1, TH1D *hChi2){

  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left", 9) == 0)
    c1->cd(2*bin+1);
  else
    c1->cd(2*bin+2);

  gPad->SetLogy();
  Double_t Chi =FitGetPars(h, bin, LR, e_LR, LR_Mmin, LR_Mmax, process,
			   whichFit, Mmin, Mmax);
  hChi2->Fill(Chi);
}


void Correct_LRerror(Double_t *e_left, Double_t *e_right, Int_t nBins){
  //If errors are 0 set errors = to nonzero value
  for (Int_t bi=0; bi<nBins; bi++) {
    Double_t eL = e_left[bi];
    Double_t eR = e_right[bi];

    if ( (eL==0) && (eR==0) ){
      cout << "Cannot correct left/right errors" << endl;
      exit(EXIT_FAILURE);
    }
    else if (eL==0) e_left[bi] = eR;
    else e_right[bi] = eL;

  }
}


void sysFunctMFit(TString start=""){
  
  //Setup_______________
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;
  TString binRange ="25_43";
  TString physBinned ="pT";//"xN", "xPi", "xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.90;
  Double_t LR_Mmax =3.30;//L/R counts mass range
  Double_t Mmin =1.00;//Fit Mass minimum
  Double_t Mmax =8.50;//Fit Mass maximum
  TString whichFit ="eight";
  
  Bool_t toWrite =false;
  //Setup_______________

  Double_t nominal_Mmin =Mmin, nominal_Mmax =Mmax;
  if (whichFit == "eight"){
    cout << "Mass fitting range preset for fit eight" << endl;
    if (physBinned=="pT") {Mmin = 1.75; Mmax = 8.00;}
    else if (physBinned=="xF") {Mmin = 1.75; Mmax = 7.50;}
    
    cout << "Fit mass range updated to:  " << Mmin << "  -  " << Mmax << endl;
  }
  
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
  TH1D *hRD_upSup_upP_L[nBins], *hRD_upSup_upP_R[nBins];
  TH1D *hRD_upSdown_upP_L[nBins], *hRD_upSdown_upP_R[nBins];
  TH1D *hRD_downSup_upP_L[nBins], *hRD_downSup_upP_R[nBins];
  TH1D *hRD_downSdown_upP_L[nBins], *hRD_downSdown_upP_R[nBins];
  TH1D *hRD_upSup_downP_L[nBins], *hRD_upSup_downP_R[nBins];
  TH1D *hRD_upSdown_downP_L[nBins], *hRD_upSdown_downP_R[nBins];
  TH1D *hRD_downSup_downP_L[nBins], *hRD_downSup_downP_R[nBins];
  TH1D *hRD_downSdown_downP_L[nBins], *hRD_downSdown_downP_R[nBins];
  
  //L/R counts for specified process
  Double_t LR_upSup_upP_L[nBins], LR_upSup_upP_R[nBins];
  Double_t LR_upSdown_upP_L[nBins], LR_upSdown_upP_R[nBins];
  Double_t LR_downSup_upP_L[nBins], LR_downSup_upP_R[nBins];
  Double_t LR_downSdown_upP_L[nBins], LR_downSdown_upP_R[nBins];
  Double_t LR_upSup_downP_L[nBins], LR_upSup_downP_R[nBins];
  Double_t LR_upSdown_downP_L[nBins], LR_upSdown_downP_R[nBins];
  Double_t LR_downSup_downP_L[nBins], LR_downSup_downP_R[nBins];
  Double_t LR_downSdown_downP_L[nBins], LR_downSdown_downP_R[nBins];

  Double_t e_LR_upSup_upP_L[nBins], e_LR_upSup_upP_R[nBins];
  Double_t e_LR_upSdown_upP_L[nBins], e_LR_upSdown_upP_R[nBins];
  Double_t e_LR_downSup_upP_L[nBins], e_LR_downSup_upP_R[nBins];
  Double_t e_LR_downSdown_upP_L[nBins], e_LR_downSdown_upP_R[nBins];
  Double_t e_LR_upSup_downP_L[nBins], e_LR_upSup_downP_R[nBins];
  Double_t e_LR_upSdown_downP_L[nBins], e_LR_upSdown_downP_R[nBins];
  Double_t e_LR_downSup_downP_L[nBins], e_LR_downSup_downP_R[nBins];
  Double_t e_LR_downSdown_downP_L[nBins], e_LR_downSdown_downP_R[nBins];
  
  const Int_t nTargPol =8;
  TCanvas* cFit[nTargPol];
  TString targNames[nTargPol] = {"upSup_upP", "upSdown_upP",       //Sub one
				 "downSup_downP", "downSdown_downP",
				 "upSup_downP", "upSdown_downP",   //Sub two
				 "downSup_upP", "downSdown_upP"};
  for (Int_t c=0; c<nTargPol; c++) {
    cFit[c] = new TCanvas(targNames[c]); cFit[c]->Divide(2, nBins); }
  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 30, 0, 3); SetUp(hChi2);
  //Perform Fits and integrate to get L/R
  for (Int_t bi=0; bi<nBins; bi++) {
    hRD_upSup_upP_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upSup_upP_%s%i",
					       physBinned.Data(), bi) );
    hRD_upSup_upP_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upSup_upP_%s%i",
					       physBinned.Data(), bi) );
    hRD_upSdown_upP_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upSdown_upP_%s%i",
					       physBinned.Data(), bi) );
    hRD_upSdown_upP_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upSdown_upP_%s%i",
					       physBinned.Data(), bi) );
    SetUp(hRD_upSup_upP_L[bi]); SetUp(hRD_upSup_upP_R[bi]);
    SetUp(hRD_upSdown_upP_L[bi]); SetUp(hRD_upSdown_upP_R[bi]); 

    hRD_upSup_downP_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upSup_downP_%s%i",
					       physBinned.Data(), bi) );
    hRD_upSup_downP_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upSup_downP_%s%i",
					       physBinned.Data(), bi) );
    hRD_upSdown_downP_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_upSdown_downP_%s%i",
			   physBinned.Data(), bi) );
    hRD_upSdown_downP_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_upSdown_downP_%s%i",
			   physBinned.Data(), bi) );
    SetUp(hRD_upSup_downP_L[bi]); SetUp(hRD_upSup_downP_R[bi]);
    SetUp(hRD_upSdown_downP_L[bi]); SetUp(hRD_upSdown_downP_R[bi]); 
    
    hRD_downSup_upP_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_downSup_upP_%s%i",
					       physBinned.Data(), bi) );
    hRD_downSup_upP_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_downSup_upP_%s%i",
					       physBinned.Data(), bi) );
    hRD_downSdown_upP_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_downSdown_upP_%s%i",
			   physBinned.Data(), bi) );
    hRD_downSdown_upP_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_downSdown_upP_%s%i",
			   physBinned.Data(), bi) );
    SetUp(hRD_downSup_upP_L[bi]); SetUp(hRD_downSup_upP_R[bi]);
    SetUp(hRD_downSdown_upP_L[bi]); SetUp(hRD_downSdown_upP_R[bi]);

    hRD_downSup_downP_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_downSup_downP_%s%i",
			   physBinned.Data(), bi) );
    hRD_downSup_downP_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_downSup_downP_%s%i",
			   physBinned.Data(), bi) );
    hRD_downSdown_downP_L[bi] =
      (TH1D*)fRD->Get(Form("MuMu_left_downSdown_downP_%s%i",
			   physBinned.Data(), bi) );
    hRD_downSdown_downP_R[bi] =
      (TH1D*)fRD->Get(Form("MuMu_right_downSdown_downP_%s%i",
			   physBinned.Data(), bi) );
    SetUp(hRD_downSup_downP_L[bi]); SetUp(hRD_downSup_downP_R[bi]);
    SetUp(hRD_downSdown_downP_L[bi]); SetUp(hRD_downSdown_downP_R[bi]); 

    WrapperFitGetPars(hRD_upSup_upP_L, bi, LR_upSup_upP_L, e_LR_upSup_upP_L,
		      LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
		      cFit[0], hChi2);
    WrapperFitGetPars(hRD_upSup_upP_R, bi, LR_upSup_upP_R, e_LR_upSup_upP_R,
		      LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
		      cFit[0], hChi2);
    WrapperFitGetPars(hRD_upSdown_upP_L, bi, LR_upSdown_upP_L,
		      e_LR_upSdown_upP_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[1], hChi2);
    WrapperFitGetPars(hRD_upSdown_upP_R, bi, LR_upSdown_upP_R,
		      e_LR_upSdown_upP_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[1], hChi2);

    WrapperFitGetPars(hRD_upSup_downP_L, bi, LR_upSup_downP_L,
		      e_LR_upSup_downP_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[2], hChi2);
    WrapperFitGetPars(hRD_upSup_downP_R, bi, LR_upSup_downP_R,
		      e_LR_upSup_downP_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[2], hChi2);
    WrapperFitGetPars(hRD_upSdown_downP_L, bi, LR_upSdown_downP_L,
		      e_LR_upSdown_downP_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[3], hChi2);
    WrapperFitGetPars(hRD_upSdown_downP_R, bi, LR_upSdown_downP_R,
		      e_LR_upSdown_downP_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[3], hChi2);
    
    WrapperFitGetPars(hRD_downSup_upP_L, bi, LR_downSup_upP_L,
		      e_LR_downSup_upP_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[4], hChi2);
    WrapperFitGetPars(hRD_downSup_upP_R, bi, LR_downSup_upP_R,
		      e_LR_downSup_upP_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[4], hChi2);
    WrapperFitGetPars(hRD_downSdown_upP_L, bi, LR_downSdown_upP_L,
		      e_LR_downSdown_upP_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[5], hChi2);
    WrapperFitGetPars(hRD_downSdown_upP_R, bi, LR_downSdown_upP_R,
		      e_LR_downSdown_upP_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[5], hChi2);

    WrapperFitGetPars(hRD_downSup_downP_L, bi, LR_downSup_downP_L,
		      e_LR_downSup_downP_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[6], hChi2);
    WrapperFitGetPars(hRD_downSup_downP_R, bi, LR_downSup_downP_R,
		      e_LR_downSup_downP_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[6], hChi2);
    WrapperFitGetPars(hRD_downSdown_downP_L, bi, LR_downSdown_downP_L,
		      e_LR_downSdown_downP_L, LR_Mmin, LR_Mmax, process,
		      whichFit, Mmin, Mmax, cFit[7], hChi2);
    WrapperFitGetPars(hRD_downSdown_downP_R, bi, LR_downSdown_downP_R,
		      e_LR_downSdown_downP_R, LR_Mmin, LR_Mmax, process,
		      whichFit, Mmin, Mmax, cFit[7], hChi2);
  }//Loop over nBins physics binning

  //Correct for zero errors
  Correct_LRerror(e_LR_upSup_upP_L, e_LR_upSup_upP_R, nBins);
  Correct_LRerror(e_LR_upSdown_upP_L, e_LR_upSdown_upP_R, nBins);
  Correct_LRerror(e_LR_upSup_downP_L, e_LR_upSup_downP_R, nBins);
  Correct_LRerror(e_LR_upSdown_downP_L, e_LR_upSdown_downP_R, nBins);

  Correct_LRerror(e_LR_downSup_upP_L, e_LR_downSup_upP_R, nBins);
  Correct_LRerror(e_LR_downSdown_upP_L, e_LR_downSdown_upP_R, nBins);
  Correct_LRerror(e_LR_downSup_downP_L, e_LR_downSup_downP_R, nBins);
  Correct_LRerror(e_LR_downSdown_downP_L, e_LR_downSdown_downP_R, nBins);
  
  //Draw hChi2
  TCanvas* cChi2 = new TCanvas("Chi2");
  hChi2->Draw();
  
  //L/R count graphs/Drawing
  TGraphErrors *g_count_upSup_upP =
    (TGraphErrors*)fRD->Get(Form("%s_left_upSup_upP", physBinned.Data() ));
  Double_t *xvals = g_count_upSup_upP->GetX();
  Double_t ex[nBins] = {0.0};
  TGraphErrors* g_Left_upSup_upP =
    new TGraphErrors(nBins, xvals, LR_upSup_upP_L, ex, e_LR_upSup_upP_L);
  TGraphErrors* g_Right_upSup_upP =
    new TGraphErrors(nBins, xvals, LR_upSup_upP_R, ex, e_LR_upSup_upP_R);
  TGraphErrors* g_Left_upSdown_upP =
    new TGraphErrors(nBins, xvals, LR_upSdown_upP_L, ex, e_LR_upSdown_upP_L);
  TGraphErrors* g_Right_upSdown_upP =
    new TGraphErrors(nBins, xvals, LR_upSdown_upP_R, ex, e_LR_upSdown_upP_R);

  TGraphErrors* g_Left_upSup_downP =
    new TGraphErrors(nBins, xvals, LR_upSup_downP_L, ex, e_LR_upSup_downP_L);
  TGraphErrors* g_Right_upSup_downP =
    new TGraphErrors(nBins, xvals, LR_upSup_downP_R, ex, e_LR_upSup_downP_R);
  TGraphErrors* g_Left_upSdown_downP =
    new TGraphErrors(nBins, xvals, LR_upSdown_downP_L, ex,
		     e_LR_upSdown_downP_L);
  TGraphErrors* g_Right_upSdown_downP =
    new TGraphErrors(nBins, xvals, LR_upSdown_downP_R, ex,
		     e_LR_upSdown_downP_R);

  TGraphErrors* g_Left_downSup_upP =
    new TGraphErrors(nBins, xvals, LR_downSup_upP_L, ex, e_LR_downSup_upP_L);
  TGraphErrors* g_Right_downSup_upP =
    new TGraphErrors(nBins, xvals, LR_downSup_upP_R, ex, e_LR_downSup_upP_R);
  TGraphErrors* g_Left_downSdown_upP =
    new TGraphErrors(nBins, xvals, LR_downSdown_upP_L, ex,
		     e_LR_downSdown_upP_L);
  TGraphErrors* g_Right_downSdown_upP =
    new TGraphErrors(nBins, xvals, LR_downSdown_upP_R, ex,
		     e_LR_downSdown_upP_R);

  TGraphErrors* g_Left_downSup_downP =
    new TGraphErrors(nBins, xvals, LR_downSup_downP_L, ex,
		     e_LR_downSup_downP_L);
  TGraphErrors* g_Right_downSup_downP =
    new TGraphErrors(nBins, xvals, LR_downSup_downP_R, ex,
		     e_LR_downSup_downP_R);
  TGraphErrors* g_Left_downSdown_downP =
    new TGraphErrors(nBins, xvals, LR_downSdown_downP_L, ex,
		     e_LR_downSdown_downP_L);
  TGraphErrors* g_Right_downSdown_downP =
    new TGraphErrors(nBins, xvals, LR_downSdown_downP_R, ex,
		     e_LR_downSdown_downP_R);
  
  SetUp(g_Left_upSup_upP); SetUp(g_Right_upSup_upP);
  SetUp(g_Left_upSdown_upP); SetUp(g_Right_upSdown_upP);
  SetUp(g_Left_upSup_downP); SetUp(g_Right_upSup_downP);
  SetUp(g_Left_upSdown_downP); SetUp(g_Right_upSdown_downP);

  SetUp(g_Left_downSup_upP); SetUp(g_Right_downSup_upP);
  SetUp(g_Left_downSdown_upP); SetUp(g_Right_downSdown_upP);
  SetUp(g_Left_downSup_downP); SetUp(g_Right_downSup_downP);
  SetUp(g_Left_downSdown_downP); SetUp(g_Right_downSdown_downP); 

  TCanvas* cLR = new TCanvas("LR counts"); cLR->Divide(2, 2);
  cLR->cd(1); g_Left_upSup_upP->Draw("AP");//UpSup
  g_Left_upSup_upP->SetTitle("L/R upSup");
  g_Right_upSup_upP->Draw("Psame");
  g_Left_upSup_downP->Draw("Psame");
  g_Right_upSup_downP->Draw("Psame");
  g_Right_upSup_upP->SetMarkerColor(kRed);
  g_Left_upSup_downP->SetMarkerColor(kBlue);
  g_Right_upSup_downP->SetMarkerColor(kGreen);

  cLR->cd(2); g_Left_upSdown_upP->Draw("AP");//UpSdown
  g_Left_upSdown_upP->SetTitle("L/R upSdown");
  g_Right_upSdown_upP->Draw("Psame");
  g_Left_upSdown_downP->Draw("Psame");
  g_Right_upSdown_downP->Draw("Psame");
  g_Right_upSdown_upP->SetMarkerColor(kRed);
  g_Left_upSdown_downP->SetMarkerColor(kBlue);
  g_Right_upSdown_downP->SetMarkerColor(kGreen);

  cLR->cd(3); g_Left_downSup_upP->Draw("AP");//DownSup
  g_Left_downSup_upP->SetTitle("L/R downSup");
  g_Right_downSup_upP->Draw("Psame");
  g_Left_downSup_downP->Draw("Psame");
  g_Right_downSup_downP->Draw("Psame");
  g_Right_downSup_upP->SetMarkerColor(kRed);
  g_Left_downSup_downP->SetMarkerColor(kBlue);
  g_Right_downSup_downP->SetMarkerColor(kGreen);

  cLR->cd(4); g_Left_downSdown_upP->Draw("AP");//DownSdown
  g_Left_downSdown_upP->SetTitle("L/R downSdown");
  g_Right_downSdown_upP->Draw("Psame");
  g_Left_downSdown_downP->Draw("Psame");
  g_Right_downSdown_downP->Draw("Psame");
  g_Right_downSdown_upP->SetMarkerColor(kRed);
  g_Left_downSdown_downP->SetMarkerColor(kBlue);
  g_Right_downSdown_downP->SetMarkerColor(kGreen);
    
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/sysFunctMFit";
  TString fOutput =
    Form("%s/sysFunctMFit_%s%.2f_%.2f_%s_%s%.2f_%.2f_%s%i_%ihbin.root",
	 thisDirPath.Data(), whichFit.Data(), nominal_Mmin, nominal_Mmax,
	 period_Mtype.Data(), process.Data(), LR_Mmin, LR_Mmax,
	 physBinned.Data(), nBins, hbins);
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    
    g_Left_upSup_upP->Write(Form("%s_left_upSup_upP", physBinned.Data()));
    g_Right_upSup_upP->Write(Form("%s_right_upSup_upP", physBinned.Data()));
    g_Left_upSdown_upP->Write(Form("%s_left_upSdown_upP", physBinned.Data()));
    g_Right_upSdown_upP->Write(Form("%s_right_upSdown_upP", physBinned.Data()));
    g_Left_upSup_downP->Write(Form("%s_left_upSup_downP", physBinned.Data()));
    g_Right_upSup_downP->Write(Form("%s_right_upSup_downP", physBinned.Data()));
    g_Left_upSdown_downP->Write(Form("%s_left_upSdown_downP", physBinned.Data()));
    g_Right_upSdown_downP->Write(Form("%s_right_upSdown_downP", physBinned.Data()));

    g_Left_downSup_upP->Write(Form("%s_left_downSup_upP", physBinned.Data()));
    g_Right_downSup_upP->Write(Form("%s_right_downSup_upP", physBinned.Data()));
    g_Left_downSdown_upP->Write(Form("%s_left_downSdown_upP", physBinned.Data()));
    g_Right_downSdown_upP->Write(Form("%s_right_downSdown_upP", physBinned.Data()));
    g_Left_downSup_downP->Write(Form("%s_left_downSup_downP", physBinned.Data()));
    g_Right_downSup_downP->Write(Form("%s_right_downSup_downP", physBinned.Data()));
    g_Left_downSdown_downP->Write(Form("%s_left_downSdown_downP", physBinned.Data()));
    g_Right_downSdown_downP->Write(Form("%s_right_downSdown_downP", physBinned.Data()));
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
