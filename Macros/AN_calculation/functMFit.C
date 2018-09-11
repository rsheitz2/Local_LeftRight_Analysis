#include "include/helperFunctions.h"

//2 Gaussian w/ psi' M/W = A*JPsi M/W by target
//2 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_six.h"

//2 Gaussian w/ psi' M/W = A*JPsi M/W by target
//1 Exponential
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/MassFitting/include/Fit_seven.h"


void DoFit(TH1D *h, TMatrixDSym &cov, Double_t Mmin, Double_t Mmax){
  h->Sumw2(); 
  TFitResultPtr status = h->Fit("fitFunc", "RLSQ", "", Mmin, Mmax);
  if (status->Status() ){
    cout << "Fit failed!!" << endl;
    cout << h->GetTitle() << endl;
    exit(EXIT_FAILURE);
  }

  cov = status->GetCovarianceMatrix();
}


Double_t FitGetPars(TH1D **h, Int_t bin, Double_t *LR, Double_t *e_LR,
		    Double_t LR_Mmin, Double_t LR_Mmax, TString process,
		    TString whichFit, Double_t Mmin, Double_t Mmax){
  
  Double_t processPars[8] = {0.0};//{aJPsi,mJPsi,wJPsi,Apsi,Mpsi,Wpsi,aDY,cDY}
  Double_t LR_cov[9];
  
  Bool_t hIsUpS =false;
  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_upstream", 18) == 0)
    hIsUpS =true;
  else if (strncmp(Form("%s", h[bin]->GetTitle() ),"MuMu_right_upstream",19)==0)
    hIsUpS =true;

  //Fit Setup_____
  TF1 *fitFunc = NULL;
  Int_t nPar;
  if (whichFit =="six"){
    fitFunc = SetupFunc_six(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_six(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
  }
  else if (whichFit =="seven"){
    fitFunc = SetupFunc_seven(h[bin], hIsUpS, fitFunc, Mmin, Mmax, &nPar);
    
    TMatrixDSym cov;
    cov.ResizeTo(nPar, nPar);
    DoFit(h[bin], cov, Mmin, Mmax);
    
    ProcessPars_seven(fitFunc, processPars, LR_cov, cov, process,nPar,hIsUpS);
  }
  else{
    cout << "Int valid fit type:   " << whichFit << endl;
    exit(EXIT_FAILURE);
  }

  //Draw physics processes
  Double_t *pars = fitFunc->GetParameters();
  const Double_t *e_pars = fitFunc->GetParErrors();
  TF1 *f_JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    Mmin, Mmax);
  f_JPsi->SetParameters(processPars[0], processPars[1], processPars[2]);
  f_JPsi->SetLineColor(kGreen); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]) )",
	    Mmin, Mmax);
  f_psi->SetParameters(processPars[3], processPars[4], processPars[5]);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");

  Double_t DY_pars[] = {processPars[6], processPars[7], Mmin};
  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  f_DY->SetParameters(DY_pars[0], DY_pars[1], DY_pars[2]);
  f_DY->SetLineColor(kGreen); f_DY->Draw("same");

  //Tmp
  //TF1 *f_Bg = new TF1("f_Bg", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  //f_Bg->SetParameters(pars[4], pars[5], Mmin);
  //f_Bg->SetLineColor(kGreen); f_Bg->Draw("same");
  //Tmp

  //Integrate for L/R counts
  if (process =="JPsi"){
    f_JPsi->SetLineColor(kBlack);
    
    LR[bin] = f_JPsi->Integral(LR_Mmin, LR_Mmax);
    e_LR[bin] = f_JPsi->IntegralError(LR_Mmin, LR_Mmax, processPars, LR_cov);
  }
  else if (process =="psi"){
    f_psi->SetLineColor(kBlack); 
  
    LR[bin] = f_psi->Integral(LR_Mmin, LR_Mmax);
    e_LR[bin] = f_psi->IntegralError(LR_Mmin, LR_Mmax, &processPars[3], LR_cov);
  }
  else if (process =="DY"){
    f_DY->SetLineColor(kBlack);

    LR[bin] = f_DY->Integral(LR_Mmin, LR_Mmax);
    e_LR[bin] = f_DY->IntegralError(LR_Mmin, LR_Mmax, DY_pars, LR_cov);
  }
  
  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();
  return Chi2/ndf;
}


void WrapperFitGetPars(TH1D **h, Int_t bin, Double_t *LR, Double_t *e_LR,
		       Double_t LR_Mmin, Double_t LR_Mmax, TString process,
		       TString whichFit, Double_t Mmin, Double_t Mmax,
		       TCanvas *c1, TH1D *hChi2){

  if (strncmp(Form("%s", h[bin]->GetTitle() ), "MuMu_left_upstream", 11) == 0)
    c1->cd(2*bin+1);
  else if (strncmp(Form("%s", h[bin]->GetTitle()),"MuMu_left_downstream",11)==0)
    c1->cd(2*bin+1);
  else
    c1->cd(2*bin+2);

  gPad->SetLogy();
  Double_t Chi =FitGetPars(h, bin, LR, e_LR, LR_Mmin, LR_Mmax, process,
			   whichFit, Mmin, Mmax);
  hChi2->Fill(Chi);
}


Double_t MakeAsym(Double_t L, Double_t R, Double_t P){
  Double_t A = L - R;
  A /= ( L + R );
  A /= P;

  return A;
}


Double_t MakeAsymError(Double_t L, Double_t R, Double_t e_L, Double_t e_R,
		       Double_t P){
  Double_t dL2 = e_L*e_L;
  Double_t dR2 = e_R*e_R;

  Double_t e = dL2/( L*L )  + dR2/( R*R );
  e = TMath::Sqrt( e );
  e *= 2.0*L*R/( (L+R)*(L+R) );
  e /= P;
  
  return e;
}


void functMFit(Bool_t PolCorr =true, TString start=""){
  //Setup_______________
  TString period_Mtype ="W09_LowM_AMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;//# of bins M is binned (must exist for input files)
      
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =Form("leftRight_byTarget_%s_%ibins_%ihbin.root",
		       period_Mtype.Data(), nBins, hbins);
  TString RDfile_noCorr =Form("leftRight_byTarget_%s_%ibins_noCorr.root",
			      period_Mtype.Data(), nBins);

  TString physBinned ="xN";//"xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.80;
  Double_t LR_Mmax =3.50;//L/R counts mass range
  Double_t Mmin =2.50;//Fit Mass minimum
  Double_t Mmax =8.50;//Fit Mass maximum
  TString whichFit ="seven";
  
  Bool_t toWrite =false;
  //Setup_______________
  
  if (start==""){
    cout<<"Script outputs AN and left/right counts per target and polarization";
    cout << " using functional mass fitting for a given fit" << endl;
    cout << "Outputs are in the formate needed for GeoMean4Targ.C" << endl;
    cout << "Script also outputs information on the fit quality" <<endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'functMFit.C(Bool_t PolCorr =true, 1)\'" << endl;
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
  TFile *fRD_noCorr = TFile::Open(pathRD + RDfile_noCorr);
  if (!fRD || !fRD_noCorr ){
    cout << "RD or RD_noCorr file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Determine polarization factor
  Double_t ex[nBins] = {0.0};
  TGraphErrors* g_asym
    =(TGraphErrors*)fRD->Get(Form("%s_asym",physBinned.Data()));
  TGraphErrors* g_asym_noCorr
    =(TGraphErrors*)fRD_noCorr->Get(Form("%s_asym",physBinned.Data()));
  if (g_asym->GetN() != nBins){
    cout << "nBins not defined well!!!" << endl;
    exit(EXIT_FAILURE);
  }
  Double_t *xvals = g_asym->GetX();
  Double_t *yvals =g_asym->GetY();
  Double_t *yvals_noCorr =g_asym_noCorr->GetY();
  Double_t Pol[nBins];
  if (PolCorr) GetPolarization(yvals_noCorr, yvals, Pol, nBins);
  else {
    for (Int_t bi=0; bi<nBins; bi++) Pol[bi] = 1.0;
  }

  //Get Input Hist/Determine LR counts for specified process
  TH1D *hRD_upS_up_L[nBins], *hRD_upS_up_R[nBins]; //Invariant M dist to be Fit
  TH1D *hRD_upS_down_L[nBins], *hRD_upS_down_R[nBins];
  TH1D *hRD_downS_up_L[nBins], *hRD_downS_up_R[nBins];
  TH1D *hRD_downS_down_L[nBins], *hRD_downS_down_R[nBins];

  //L/R counts for specified process
  Double_t LR_upS_up_L[nBins], LR_upS_up_R[nBins]; 
  Double_t LR_upS_down_L[nBins], LR_upS_down_R[nBins];
  Double_t LR_downS_up_L[nBins], LR_downS_up_R[nBins];
  Double_t LR_downS_down_L[nBins], LR_downS_down_R[nBins];
  
  Double_t e_LR_upS_up_L[nBins], e_LR_upS_up_R[nBins];
  Double_t e_LR_upS_down_L[nBins], e_LR_upS_down_R[nBins];
  Double_t e_LR_downS_up_L[nBins], e_LR_downS_up_R[nBins];
  Double_t e_LR_downS_down_L[nBins], e_LR_downS_down_R[nBins];


  const Int_t nTargPol =4;
  TCanvas* cFit[nTargPol];
  TString targNames[nTargPol] = {"upS_up", "upS_down", "downS_up","downS_down"};
  for (Int_t c=0; c<nTargPol; c++) {
    cFit[c] = new TCanvas(targNames[c]); cFit[c]->Divide(2, nBins); }
  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 30, 0, 3); SetUpHist(hChi2);
  //Perform Fits and integrate to get L/R
  for (Int_t bi=0; bi<nBins; bi++) {
    hRD_upS_up_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s%i", 
					    physBinned.Data(), bi) );
    hRD_upS_up_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_up_%s%i", 
					    physBinned.Data(), bi) );
    hRD_upS_down_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_down_%s%i", 
					      physBinned.Data(), bi) );
    hRD_upS_down_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_down_%s%i", 
					      physBinned.Data(), bi) );
    SetUpHist(hRD_upS_up_L[bi]); SetUpHist(hRD_upS_up_R[bi]); 
    SetUpHist(hRD_upS_down_L[bi]); SetUpHist(hRD_upS_down_R[bi]);
    
    hRD_downS_up_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_up_%s%i", 
					      physBinned.Data(), bi) );
    hRD_downS_up_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_downstream_up_%s%i", 
					      physBinned.Data(), bi) );
    hRD_downS_down_L[bi]= (TH1D*)fRD->Get(Form("MuMu_left_downstream_down_%s%i",
					       physBinned.Data(), bi) );
    hRD_downS_down_R[bi]=(TH1D*)fRD->Get(Form("MuMu_right_downstream_down_%s%i",
					      physBinned.Data(), bi) );
    SetUpHist(hRD_downS_up_L[bi]); SetUpHist(hRD_downS_up_R[bi]);
    SetUpHist(hRD_downS_down_L[bi]); SetUpHist(hRD_downS_down_R[bi]);
    

    WrapperFitGetPars(hRD_upS_up_L, bi, LR_upS_up_L, e_LR_upS_up_L, LR_Mmin,
		      LR_Mmax, process, whichFit, Mmin, Mmax, cFit[0], hChi2);
    WrapperFitGetPars(hRD_upS_up_R, bi, LR_upS_up_R, e_LR_upS_up_R, LR_Mmin,
		      LR_Mmax, process, whichFit, Mmin, Mmax, cFit[0], hChi2);
    WrapperFitGetPars(hRD_upS_down_L, bi, LR_upS_down_L, e_LR_upS_down_L,
		      LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
		      cFit[1], hChi2);
    WrapperFitGetPars(hRD_upS_down_R, bi, LR_upS_down_R, e_LR_upS_down_R,
		      LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
		      cFit[1], hChi2);

    WrapperFitGetPars(hRD_downS_up_L, bi, LR_downS_up_L,e_LR_downS_up_L,LR_Mmin,
		      LR_Mmax, process, whichFit, Mmin, Mmax, cFit[2], hChi2);
    WrapperFitGetPars(hRD_downS_up_R, bi, LR_downS_up_R,e_LR_downS_up_R,LR_Mmin,
		      LR_Mmax, process, whichFit, Mmin, Mmax, cFit[2], hChi2);
    WrapperFitGetPars(hRD_downS_down_L, bi, LR_downS_down_L, e_LR_downS_down_L,
		      LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax,
		      cFit[3], hChi2);
    WrapperFitGetPars(hRD_downS_down_R, bi, LR_downS_down_R, e_LR_downS_down_R,
		      LR_Mmin, LR_Mmax, process, whichFit, Mmin, Mmax, cFit[3],
		      hChi2);
  }//Loop over nBins physics binning
  
  //Draw hChi2
  TCanvas* cChi2 = new TCanvas("Chi2");
  hChi2->Draw();
  
  //L/R count graphs/Drawing
  TGraphErrors* g_Left_upS_up =
    new TGraphErrors(nBins, xvals, LR_upS_up_L, ex, e_LR_upS_up_L);
  TGraphErrors* g_Right_upS_up =
    new TGraphErrors(nBins, xvals, LR_upS_up_R, ex, e_LR_upS_up_R);
  TGraphErrors* g_Left_upS_down =
    new TGraphErrors(nBins, xvals, LR_upS_down_L, ex, e_LR_upS_down_L);
  TGraphErrors* g_Right_upS_down =
    new TGraphErrors(nBins, xvals, LR_upS_down_R, ex, e_LR_upS_down_R);
  TGraphErrors* g_Left_downS_up =
    new TGraphErrors(nBins, xvals, LR_downS_up_L, ex, e_LR_downS_up_L);
  TGraphErrors* g_Right_downS_up =
    new TGraphErrors(nBins, xvals, LR_downS_up_R, ex, e_LR_downS_up_R);
  TGraphErrors* g_Left_downS_down =
    new TGraphErrors(nBins, xvals, LR_downS_down_L, ex, e_LR_downS_down_L);
  TGraphErrors* g_Right_downS_down =
    new TGraphErrors(nBins, xvals, LR_downS_down_R, ex, e_LR_downS_down_R);
  
  SetUpTGraph(g_Left_upS_up); SetUpTGraph(g_Right_upS_up);
  SetUpTGraph(g_Left_upS_down); SetUpTGraph(g_Right_upS_down);
  SetUpTGraph(g_Left_downS_up); SetUpTGraph(g_Right_downS_up);
  SetUpTGraph(g_Left_downS_down); SetUpTGraph(g_Right_downS_down);

  TCanvas* cLR = new TCanvas("LR counts"); cLR->Divide(2, 2);
  cLR->cd(1); g_Left_upS_up->Draw("AP");//Upstream
  g_Left_upS_up->SetTitle("L/R upS_up");
  g_Right_upS_up->Draw("Psame"); g_Right_upS_up->SetMarkerColor(kRed);
  
  cLR->cd(2); g_Left_upS_down->Draw("AP");
  g_Left_upS_down->SetTitle("L/R upS_down");
  g_Right_upS_down->Draw("Psame");
  g_Right_upS_down->SetMarkerColor(kRed);
  
  cLR->cd(3); g_Left_downS_up->Draw("AP");//Downstream
  g_Left_downS_up->SetTitle("L/R downS_up");
  g_Right_downS_up->Draw("Psame");
  g_Right_downS_up->SetMarkerColor(kRed);
  
  cLR->cd(4); g_Left_downS_down->Draw("AP");
  g_Left_downS_down->SetTitle("L/R downS_down");
  g_Right_downS_down->Draw("Psame");
  g_Right_downS_down->SetMarkerColor(kRed);


  //Make AN asymmetery per targ && pol
  Double_t AN_upS_up[nBins], e_AN_upS_up[nBins];
  Double_t AN_upS_down[nBins], e_AN_upS_down[nBins];
  Double_t AN_downS_up[nBins], e_AN_downS_up[nBins];
  Double_t AN_downS_down[nBins], e_AN_downS_down[nBins];
  
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_upS_up[bi] = MakeAsym(LR_upS_up_L[bi], LR_upS_up_R[bi], Pol[bi]);
    e_AN_upS_up[bi]
      = MakeAsymError(LR_upS_up_L[bi], LR_upS_up_R[bi], e_LR_upS_up_L[bi],
		      e_LR_upS_up_R[bi],Pol[bi]);
    
    AN_upS_down[bi] = MakeAsym(LR_upS_down_L[bi], LR_upS_down_R[bi], Pol[bi]);
    e_AN_upS_down[bi]
      =MakeAsymError(LR_upS_down_L[bi], LR_upS_down_R[bi], e_LR_upS_down_L[bi],
		     e_LR_upS_down_R[bi],Pol[bi]);
    
    AN_downS_up[bi] = MakeAsym(LR_downS_up_L[bi], LR_downS_up_R[bi], Pol[bi]);
    e_AN_downS_up[bi] =
      MakeAsymError(LR_downS_up_L[bi], LR_downS_up_R[bi], e_LR_downS_up_L[bi],
		    e_LR_downS_up_R[bi],Pol[bi]);
    
    AN_downS_down[bi]=MakeAsym(LR_downS_down_L[bi],LR_downS_down_R[bi],Pol[bi]);
    e_AN_downS_down[bi]
      =MakeAsymError(LR_downS_down_L[bi],LR_downS_down_R[bi],
		     e_LR_downS_down_L[bi],e_LR_downS_down_R[bi], Pol[bi]);
  }

  //Make/Draw AN graphs per targ && pol
  TGraphErrors *g_AN_upS_up
    = new TGraphErrors(nBins, xvals,AN_upS_up, ex, e_AN_upS_up);
  TGraphErrors *g_AN_upS_down =
    new TGraphErrors(nBins, xvals,AN_upS_down, ex, e_AN_upS_down);
  TGraphErrors *g_AN_downS_up =
    new TGraphErrors(nBins, xvals,AN_downS_up, ex, e_AN_downS_up);
  TGraphErrors *g_AN_downS_down =
    new TGraphErrors(nBins, xvals,AN_downS_down, ex, e_AN_downS_down);
  
  SetUpTGraph(g_AN_upS_up); SetUpTGraph(g_AN_upS_down);
  SetUpTGraph(g_AN_downS_up); SetUpTGraph(g_AN_downS_down);
  
  TCanvas* cAN = new TCanvas("AN by tar && pol"); cAN->Divide(2, 2);
  Double_t yMax;
  if (process =="JPsi")  (PolCorr) ? yMax =0.5 : yMax =0.5;
  else if (process =="psi")  (PolCorr) ? yMax =1.0 : yMax =0.5;
  else if (process =="DY")  (PolCorr) ? yMax =3.0 : yMax =0.5;
  cAN->cd(1);
  g_AN_upS_up->Draw("AP"); g_AN_upS_up->SetTitle("AN_upS_up");
  g_AN_upS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_upS_up, 0.0);

  cAN->cd(2);
  g_AN_upS_down->Draw("AP"); g_AN_upS_up->SetTitle("AN_upS_down");
  g_AN_upS_down->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_upS_down, 0.0);

  cAN->cd(3);
  g_AN_downS_up->Draw("AP"); g_AN_upS_up->SetTitle("AN_downS_up");
  g_AN_downS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_downS_up, 0.0);

  cAN->cd(4);
  g_AN_downS_down->Draw("AP"); g_AN_upS_up->SetTitle("AN_downS_down");
  g_AN_downS_down->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_downS_down, 0.0);

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation";
  TString fOutput
    = Form("%s/Data/functMFit/functMFit_%s%.2f_%.2f_%s_%s%.2f_%.2f_%s%i_%ihbin",
	   thisDirPath.Data(), whichFit.Data(), Mmin, Mmax,
	   period_Mtype.Data(), process.Data(), LR_Mmin, LR_Mmax,
	   physBinned.Data(), nBins, hbins);
  fOutput += (PolCorr) ? "_corr.root" : "_noCorr.root";
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TList *doc = new TList();
    doc->Add((TObject*)(new TObjString(pathRD+"\n"+RDfile+"\n"+RDfile_noCorr)));
    doc->Write("InputData");
    
    g_AN_upS_up->Write("AN_upS_up");
    g_AN_upS_down->Write("AN_upS_down");
    g_AN_downS_up->Write("AN_downS_up");
    g_AN_downS_down->Write("AN_downS_down");

    g_Left_upS_up->Write("Left_upS_up");
    g_Left_upS_down->Write("Left_upS_down");
    g_Left_downS_up->Write("Left_downS_up");
    g_Left_downS_down->Write("Left_downS_down");

    g_Right_upS_up->Write("Right_upS_up");
    g_Right_upS_down->Write("Right_upS_down");
    g_Right_downS_up->Write("Right_downS_up");
    g_Right_downS_down->Write("Right_downS_down");
  }

  cout << " " << endl;
  cout << "Settings______" << endl;
  cout << "Real data file considered:              " << RDfile << endl;
  cout <<"Number of histogram bins used in M fit:  "
       <<hRD_upS_up_L[0]->GetNbinsX()<<endl;
  cout << "Mass Range for fitting is:              " << Mmin
       << " - " << Mmax << endl;
  cout << "physBinned nBins times:                 " << nBins << endl;
  cout << "Polarization was performed:             " << PolCorr << endl;
  cout << "\n";
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
