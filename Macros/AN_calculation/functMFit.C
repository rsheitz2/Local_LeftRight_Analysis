#include "include/helperFunctions.h"
#include "include/fitFunctions.h"

TF1* FitGetLR(TH1D **h, TH1D **hRatio, Int_t bin, Double_t *LR,Double_t *e_LR,
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
  TF1 *fitFunc = FitGetLR(h, hRatio, bin, LR, e_LR,LR_Mmin, LR_Mmax, process,
			  whichFit, Mmin, Mmax);

  //Get reduced Chi2
  Double_t Chi2 =fitFunc->GetChisquare();
  Double_t ndf =fitFunc->GetNDF();
  Double_t redChi2 = Chi2/ndf;
  hChi2->Fill(redChi2);

  return fitFunc;
}


void FitOne(TH1D **h, TH1D **hRatio, Int_t bin, Double_t *LR,Double_t *e_LR,
	    Double_t LR_Mmin, Double_t LR_Mmax, TString process,
	    TString whichFit, Double_t Mmin, Double_t Mmax, TCanvas *c1){
  //Used to draw just one histogram fit (for better visibility)
  c1->cd(1);
  gPad->SetLogy();
  TF1 *fitFunc = FitGetLR(h, hRatio, bin, LR, e_LR, LR_Mmin, LR_Mmax, process,
			  whichFit, Mmin, Mmax);
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


Double_t LocalChi2(TH1D *h, Double_t min, Double_t max){
  Int_t lowerBin = h->FindBin(min);
  Int_t upperBin = h->FindBin(max);

  Double_t Chi2 =0.0;
  for (Int_t bi=lowerBin; bi<upperBin+1; bi++) {
    Double_t val = h->GetBinContent(bi);
    Double_t eVal = h->GetBinError(bi);
    Chi2 += (val-1.0)*(val-1.0)/(eVal*eVal);
  }

  Chi2 = TMath::Sqrt(Chi2);
  Int_t ndf = upperBin - lowerBin;
  Chi2 /= (1.0*ndf+1.0);

  return Chi2;
}


void SetupRatio(TH1D* h, Double_t RatioPercent){
  h->GetYaxis()->SetRangeUser(1-RatioPercent,1+RatioPercent);
  DrawLine(h, 1.0);
}


void functMFit(TString start=""){
  
  //Setup_______________
  /*TString period_Mtype ="W07_HMDY";
    Int_t hbins =150;//# of histogram bins using in mass fitting
    const Int_t nBins =3;
    TString binRange ="43_85";
    TString physBinned ="pT";//"xN", "xPi", "xF", "pT"
    TString process ="DY";//JPsi, psi, DY
    Double_t LR_Mmin =4.30;
    Double_t LR_Mmax =8.50;//L/R counts mass range
    Double_t Mmin =4.30;//Fit Mass minimum
    Double_t Mmax =8.50;//Fit Mass maximum
    TString whichFit ="eight";//*/
  
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;
  TString binRange ="25_43";
  TString physBinned ="pT";//"xN", "xPi", "xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.90;
  Double_t LR_Mmax =3.30;//L/R counts mass range
  Double_t Mmin =2.00;//Fit Mass minimum
  Double_t Mmax =7.50;//Fit Mass maximum
  TString whichFit ="eight";//*/
  
  Bool_t toWrite =false;
  //Setup_______________

  Double_t nominal_Mmin =Mmin, nominal_Mmax =Mmax;
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =Form("leftRight_byTarget_%s1.00_8.50_%ibins%s_%ihbin.root",
		       period_Mtype.Data(), nBins, binRange.Data(), hbins);
  
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
  if (!fRD ){
    cout << "RD file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Determine polarization factor
  Double_t ex[nBins] = {0.0};
  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  TGraph* g_Dil =(TGraph*)fRD->Get(Form("%s_Dil", physBinned.Data()));
  Double_t *xvals = g_Pol->GetX();
  Double_t Pol[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);

  //Get Input Hist/fit ratio/Determine LR counts for specified process
  TH1D *hRD_upS_up_L[nBins], *hRD_upS_up_R[nBins]; //Invariant M dist to be Fit
  TH1D *hRD_upS_down_L[nBins], *hRD_upS_down_R[nBins];
  TH1D *hRD_downS_up_L[nBins], *hRD_downS_up_R[nBins];
  TH1D *hRD_downS_down_L[nBins], *hRD_downS_down_R[nBins];

  TH1D *hRatio_upS_up_L[nBins], *hRatio_upS_up_R[nBins]; //RD/fit ratio
  TH1D *hRatio_upS_down_L[nBins], *hRatio_upS_down_R[nBins];
  TH1D *hRatio_downS_up_L[nBins], *hRatio_downS_up_R[nBins];
  TH1D *hRatio_downS_down_L[nBins], *hRatio_downS_down_R[nBins];

  //Fit parameters to check
  const Int_t nTargPol =4, nSelectPars =5;
  Double_t pars_upS[nBins*nTargPol*nSelectPars];
  Double_t pars_downS[nBins*nTargPol*nSelectPars];
  Double_t e_pars_upS[nBins*nTargPol*nSelectPars];
  Double_t e_pars_downS[nBins*nTargPol*nSelectPars];
  
  //L/R counts for specified process
  Double_t LR_upS_up_L[nBins], LR_upS_up_R[nBins]; 
  Double_t LR_upS_down_L[nBins], LR_upS_down_R[nBins];
  Double_t LR_downS_up_L[nBins], LR_downS_up_R[nBins];
  Double_t LR_downS_down_L[nBins], LR_downS_down_R[nBins];
  
  Double_t e_LR_upS_up_L[nBins], e_LR_upS_up_R[nBins];
  Double_t e_LR_upS_down_L[nBins], e_LR_upS_down_R[nBins];
  Double_t e_LR_downS_up_L[nBins], e_LR_downS_up_R[nBins];
  Double_t e_LR_downS_down_L[nBins], e_LR_downS_down_R[nBins];

  //TargPol setup
  TCanvas* cFit[nTargPol];
  TString targNames[nTargPol] = {"upS_up", "upS_down", "downS_up","downS_down"};
  for (Int_t c=0; c<nTargPol; c++) {
    cFit[c] = new TCanvas(targNames[c]); cFit[c]->Divide(2, nBins); }

  //To draw just one fit/ratio for better visibility
  //TCanvas* cOne = new TCanvas(); cOne->Divide(1, 2);
  
  TH1D* hChi2 = new TH1D("hChi2", "hChi2", 30, 0, 2); SetUpHist(hChi2);
  
  //Perform Fits and integrate to get L/R
  for (Int_t bi=0, iPar_upS=0, iPar_dS=0; bi<nBins; bi++) {
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

    hRatio_upS_up_L[bi] = (TH1D*)hRD_upS_up_L[bi]->Clone();
    hRatio_upS_up_R[bi] = (TH1D*)hRD_upS_up_R[bi]->Clone();
    hRatio_upS_down_L[bi] = (TH1D*)hRD_upS_down_L[bi]->Clone();
    hRatio_upS_down_R[bi] = (TH1D*)hRD_upS_down_R[bi]->Clone();
    SetUpHist(hRatio_upS_up_L[bi]); SetUpHist(hRatio_upS_up_R[bi]); 
    SetUpHist(hRatio_upS_down_L[bi]); SetUpHist(hRatio_upS_down_R[bi]);

    hRatio_downS_up_L[bi] = (TH1D*)hRD_downS_up_L[bi]->Clone();
    hRatio_downS_up_R[bi] = (TH1D*)hRD_downS_up_R[bi]->Clone();
    hRatio_downS_down_L[bi] = (TH1D*)hRD_downS_down_L[bi]->Clone();
    hRatio_downS_down_R[bi] = (TH1D*)hRD_downS_down_R[bi]->Clone();
    SetUpHist(hRatio_downS_up_L[bi]); SetUpHist(hRatio_downS_up_R[bi]); 
    SetUpHist(hRatio_downS_down_L[bi]); SetUpHist(hRatio_downS_down_R[bi]);

    Bool_t hIsUpS =true; 
    TF1 *fitFunc = FitGetLR(hRD_upS_up_L, hRatio_upS_up_L, bi, LR_upS_up_L,
			    e_LR_upS_up_L, LR_Mmin, LR_Mmax, process, whichFit,
			    Mmin, Mmax, cFit[0], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_upS[iPar_upS], &e_pars_upS[iPar_upS],
			 process, whichFit, hIsUpS); iPar_upS +=nSelectPars;
    
    fitFunc = FitGetLR(hRD_upS_up_R, hRatio_upS_up_R, bi, LR_upS_up_R,
		       e_LR_upS_up_R, LR_Mmin, LR_Mmax, process, whichFit, Mmin,
		       Mmax, cFit[0], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_upS[iPar_upS], &e_pars_upS[iPar_upS],
			 process, whichFit, hIsUpS); iPar_upS +=nSelectPars;
    
    fitFunc = FitGetLR(hRD_upS_down_L, hRatio_upS_down_L, bi, LR_upS_down_L,
		       e_LR_upS_down_L, LR_Mmin, LR_Mmax, process,whichFit,Mmin,
		       Mmax, cFit[1], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_upS[iPar_upS], &e_pars_upS[iPar_upS],
			 process, whichFit, hIsUpS); iPar_upS +=nSelectPars;
    
    fitFunc =FitGetLR(hRD_upS_down_R, hRatio_upS_down_R, bi, LR_upS_down_R,
		      e_LR_upS_down_R, LR_Mmin, LR_Mmax, process, whichFit,Mmin,
		      Mmax, cFit[1], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_upS[iPar_upS], &e_pars_upS[iPar_upS],
			 process, whichFit, hIsUpS); iPar_upS +=nSelectPars;//*/
    /*if (bi != 0){ //To draw just one fit/ratio for better visibility
      fitFunc =FitGetLR(hRD_upS_down_R, hRatio_upS_down_R, bi, LR_upS_down_R,
      e_LR_upS_down_R, LR_Mmin, LR_Mmax, process, whichFit,
      Mmin, Mmax, cFit[1], hChi2);
      }
      else {
      FitOne(hRD_upS_down_R, hRatio_upS_down_R, bi, LR_upS_down_R,
      e_LR_upS_down_R, LR_Mmin, LR_Mmax, process, whichFit,
      Mmin, Mmax, cOne);
      cout << "\n\nDraw option currently in use\n\n";
      }//*/

    hIsUpS =false;
    fitFunc =FitGetLR(hRD_downS_up_L, hRatio_downS_up_L, bi, LR_downS_up_L,
		      e_LR_downS_up_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[2], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_downS[iPar_dS], &e_pars_downS[iPar_dS],
			 process, whichFit, hIsUpS); iPar_dS +=nSelectPars;
    
    fitFunc =FitGetLR(hRD_downS_up_R, hRatio_downS_up_R, bi, LR_downS_up_R,
		      e_LR_downS_up_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[2], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_downS[iPar_dS], &e_pars_downS[iPar_dS],
			 process, whichFit, hIsUpS); iPar_dS +=nSelectPars;
    
    fitFunc =FitGetLR(hRD_downS_down_L, hRatio_downS_down_L, bi,LR_downS_down_L,
		      e_LR_downS_down_L, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[3], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_downS[iPar_dS], &e_pars_downS[iPar_dS],
			 process, whichFit, hIsUpS); iPar_dS +=nSelectPars;
    
    fitFunc =FitGetLR(hRD_downS_down_R, hRatio_downS_down_R, bi,LR_downS_down_R,
		      e_LR_downS_down_R, LR_Mmin, LR_Mmax, process, whichFit,
		      Mmin, Mmax, cFit[3], hChi2);
    SelectFitPhysicsPars(fitFunc, &pars_downS[iPar_dS], &e_pars_downS[iPar_dS],
			 process, whichFit, hIsUpS); iPar_dS +=nSelectPars;
  }//Loop over nBins physics binning

  //Draw RD/fit ratios
  TCanvas* cRatio[nTargPol];
  for (Int_t c=0; c<nTargPol; c++) {
    cRatio[c] = new TCanvas("Ratio_"+targNames[c]);cRatio[c]->Divide(2, nBins);}
  Double_t RatioPercent =0.2;
  for (Int_t bi=0; bi<nBins; bi++) {
    cRatio[0]->cd(2*bi + 1); hRatio_upS_up_L[bi]->Draw();
    SetupRatio(hRatio_upS_up_L[bi], RatioPercent);
    cRatio[0]->cd(2*bi + 2); hRatio_upS_up_R[bi]->Draw();
    SetupRatio(hRatio_upS_up_R[bi], RatioPercent);
    cRatio[1]->cd(2*bi + 1); hRatio_upS_down_L[bi]->Draw();
    SetupRatio(hRatio_upS_down_L[bi], RatioPercent);
    cRatio[1]->cd(2*bi + 2); hRatio_upS_down_R[bi]->Draw();
    SetupRatio(hRatio_upS_down_R[bi], RatioPercent);

    cRatio[2]->cd(2*bi + 1); hRatio_downS_up_L[bi]->Draw();
    SetupRatio(hRatio_downS_up_L[bi], RatioPercent);
    cRatio[2]->cd(2*bi + 2); hRatio_downS_up_R[bi]->Draw();
    SetupRatio(hRatio_downS_up_R[bi], RatioPercent);
    cRatio[3]->cd(2*bi + 1); hRatio_downS_down_L[bi]->Draw();
    SetupRatio(hRatio_downS_down_L[bi], RatioPercent);
    cRatio[3]->cd(2*bi + 2); hRatio_downS_down_R[bi]->Draw();
    SetupRatio(hRatio_downS_down_R[bi], RatioPercent);
  }
  
  //Draw hChi2
  TCanvas* cChi2 = new TCanvas("Chi2");
  hChi2->Draw();

  //Make LR Integration region Chi2
  TH1D* hLRchi2 = new TH1D("hLRchi2", "hLRchi2", 20, 0, 1);
  Double_t lrChi2[nBins] = {0.0};
  for (Int_t bi=0; bi<nBins; bi++) {
    Double_t Chi2 =LocalChi2(hRatio_upS_up_L[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    Chi2 =LocalChi2(hRatio_upS_up_R[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    Chi2 =LocalChi2(hRatio_upS_down_L[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    Chi2 =LocalChi2(hRatio_upS_down_R[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;

    Chi2 =LocalChi2(hRatio_downS_up_L[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    Chi2 =LocalChi2(hRatio_downS_up_R[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    Chi2 =LocalChi2(hRatio_downS_down_L[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    Chi2 =LocalChi2(hRatio_downS_down_R[bi], LR_Mmin, LR_Mmax);
    hLRchi2->Fill(Chi2); lrChi2[bi] += Chi2;
    lrChi2[bi] /= 8.0;
  }
  TCanvas* cLRChi2 = new TCanvas(); cLRChi2->Divide(2);
  cLRChi2->cd(1); SetUpHist(hLRchi2); hLRchi2->Draw();
  TGraph *gLRChi2 = new TGraph(nBins, xvals, lrChi2);
  SetUpTGraph(gLRChi2); gLRChi2->SetTitle("LRChi2");
  cLRChi2->cd(2); gLRChi2->Draw("AP");

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

  gStyle->SetOptFit(11);
  TCanvas* cLR = new TCanvas("LR counts"); cLR->Divide(2, 2);
  cLR->cd(1); g_Left_upS_up->Draw("AP");//Upstream
  g_Left_upS_up->SetTitle("L/R upS_up");
  g_Right_upS_up->Draw("Psames"); g_Right_upS_up->SetMarkerColor(kRed);
  //g_Left_upS_up->Fit("pol0", "Q"); g_Right_upS_up->Fit("pol0", "Q"); 
  
  cLR->cd(2); g_Left_upS_down->Draw("AP");
  g_Left_upS_down->SetTitle("L/R upS_down");
  g_Right_upS_down->Draw("Psames");
  g_Right_upS_down->SetMarkerColor(kRed);
  //g_Left_upS_down->Fit("pol0", "Q"); g_Right_upS_down->Fit("pol0", "Q"); 
  
  cLR->cd(3); g_Left_downS_up->Draw("AP");//Downstream
  g_Left_downS_up->SetTitle("L/R downS_up");
  g_Right_downS_up->Draw("Psames");
  g_Right_downS_up->SetMarkerColor(kRed);
  //g_Left_downS_up->Fit("pol0", "Q"); g_Right_downS_up->Fit("pol0", "Q"); 
  
  cLR->cd(4); g_Left_downS_down->Draw("AP");
  g_Left_downS_down->SetTitle("L/R downS_down");
  g_Right_downS_down->Draw("Psames");
  g_Right_downS_down->SetMarkerColor(kRed);
  //g_Left_downS_down->Fit("pol0", "Q"); g_Right_downS_down->Fit("pol0", "Q"); 

  //Setup/Draw select fit parameters
  Double_t ex_pars[nBins*nTargPol*nSelectPars] = {0.0};
  Double_t xvals_pars[nBins*nTargPol*nSelectPars];
  for (Int_t i=0; i<nBins*nTargPol*nSelectPars; i+=nSelectPars) {
    for (Int_t j=0; j<nSelectPars; j++) {
      xvals_pars[i+j] = 1.0 + 1.0*i + j/5.0; } }

  TGraphErrors* g_pars_upS = new
    TGraphErrors(nBins*nTargPol*nSelectPars, xvals_pars, pars_upS, ex_pars,
		 e_pars_upS); SetUpTGraph(g_pars_upS);
  TGraphErrors* g_pars_downS = new
    TGraphErrors(nBins*nTargPol*nSelectPars, xvals_pars, pars_downS, ex_pars,
		 e_pars_downS); SetUpTGraph(g_pars_downS);
  g_pars_downS->SetMarkerColor(kRed);
  TCanvas* cPars = new TCanvas(); cPars->Divide(2);
  
  //Change the following ranges to look at different parameters
  Double_t min_parMrange =3.1, max_parMrange =3.15;//JPsi 
  Double_t min_parWrange =0.15, max_parWrange =0.175;//JPsi
  
  //Double_t min_parMrange =3.5, max_parMrange =3.65;//psi' eight
  //Double_t min_parWrange =0.18, max_parWrange =0.21;//psi' eight
  //Double_t min_parMrange =3.3, max_parMrange =3.5;//psi' ten
  //Double_t min_parWrange =0.22, max_parWrange =0.35;//psi' 
  
  //Double_t min_parMrange =-1.5, max_parMrange =-0.5;//DY slope
  
  cPars->cd(1);
  g_pars_upS->Draw("AP");
  g_pars_upS->GetYaxis()->SetRangeUser(min_parMrange, max_parMrange); 
  g_pars_upS->SetTitle("Fit Mass");
  g_pars_downS->Draw("Psame");
  cPars->cd(2);
  g_pars_downS->Draw("AP");
  g_pars_downS->GetYaxis()->SetRangeUser(min_parWrange, max_parWrange);
  g_pars_downS->SetTitle("Fit Width"); 
  g_pars_upS->Draw("Psame"); 
    
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
  if (process =="JPsi")  yMax =0.45;
  else if (process =="psi")  yMax =1.0;
  else if (process =="DY")  yMax =3.0;
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

  //Draw One fit/Ratio //To draw just one fit/ratio for better visibility
  /*cOne->cd(2);
  hRatio_upS_down_R[0]->Draw();
  hRatio_upS_down_R[0]->GetXaxis()->SetRangeUser(1.0, 8.5);
  DrawLine(hRatio_upS_down_R[0], 1.0);//*/

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation";
  TString fOutput
    = Form("%s/Data/functMFit/functMFit_%s%.2f_%.2f_%s_%s%.2f_%.2f_%s%i_%ihbin",
	   thisDirPath.Data(), whichFit.Data(), nominal_Mmin, nominal_Mmax,
	   period_Mtype.Data(), process.Data(), LR_Mmin, LR_Mmax,
	   physBinned.Data(), nBins, hbins);
  fOutput += "_corr.root";
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TList *doc = new TList();
    doc->Add((TObject*)(new TObjString(pathRD+"\n"+RDfile)));
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
