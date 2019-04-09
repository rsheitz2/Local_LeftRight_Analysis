#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

TH1D *hist[4] = {NULL, NULL, NULL, NULL}; //{JPsi, psi, OC, AMDY}

Double_t Fit_MC(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );

  if (hist[2]->GetBinError(hist[2]->FindBin(x[0])) > 0.5*nOC ) {
    return par[0]*nJPsi +par[1]*npsi +par[3]*nAMDY;
  }
  else{
    return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
  }
}

TH1D* GetMCHist(TFile *f, TString name, Int_t iMC, Int_t hbins);

Int_t BasicChecks(Double_t *pars, const Double_t *epars);

void mcMFit_new(){
  
  //Setup_______________
  TString period_Mtype ="WAll_LowM_AMDY"; 
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;
  TString binRange ="25_43"; //"29_34";
  TString physBinned ="pT";//"xN", "xPi", "xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.5;
  Double_t LR_Mmax =4.3;//L/R counts mass range
  Double_t Mmin =2.00;//Fit Mass minimum
  Double_t Mmax =8.50;//Fit Mass maximum
  
  Bool_t toWrite =false;
  //Setup_______________

  //Basic setup
  gStyle->SetOptFit(11111);
  gStyle->SetOptStat(111111);

  TString pathData = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =
    Form("leftRight_byTarget_%s1.00_8.50_%ibins%s_%ihbin_slot1_phiS0.0.root",
		       period_Mtype.Data(), nBins, binRange.Data(), hbins);
  TString JPsiFile, psiFile, OCFile, AMDYFile, LMDYFile, HMDYFile;
  //Yu MC used
  /*JPsiFile = Form("leftRight_byTarget_Yu_JPsi1.00_8.50_%ibins%s_%ihbin.root",
		  nBins, binRange.Data(), hbins);
  psiFile = Form("leftRight_byTarget_Yu_psi1.00_8.50_%ibins%s_%ihbin.root",
		 nBins, binRange.Data(), hbins);
  OCFile = Form("leftRight_byTarget_Yu_OC1.00_8.50_%ibins%s_%ihbin.root",
		nBins, binRange.Data(), hbins);
  AMDYFile = Form("leftRight_byTarget_Yu_AMDY1.00_8.50_%ibins%s_%ihbin.root",
		  nBins, binRange.Data(), hbins);
  cout << "\n\nYu-Shiang MC used\n\n";//*/

  //Charles MC used
  JPsiFile =
    Form("leftRight_byTarget_Charles_Jpsi1.00_8.50_%ibins%s_%ihbin.root",
	 nBins, binRange.Data(), hbins);
  psiFile =
    Form("leftRight_byTarget_Charles_Psi1.00_8.50_%ibins%s_%ihbin.root",
	 nBins, binRange.Data(), hbins);
  OCFile =
    Form("leftRight_byTarget_Charles_OC1.00_8.50_%ibins%s_%ihbin.root",
	 nBins, binRange.Data(), hbins);
  AMDYFile =
    Form("leftRight_byTarget_Charles_AMDY1.00_8.50_%ibins%s_%ihbin.root",
	 nBins, binRange.Data(), hbins);
  cout << "\n\nCharles MC used\n\n";//*/
    
  /*if (start==""){
    cout<<"Script outputs AN and left/right counts per target and polarization";
    cout << " using functional mass fitting for a given fit" << endl;
    cout << "Outputs are in the formate needed for GeoMean4Targ.C" << endl;
    cout << "Script also outputs information on the fit quality" <<endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'functMFit.C(Bool_t PolCorr =true, 1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Real data path:             " << pathData << endl;
    cout << "Real data file considered:  " << RDfile << endl;
    cout << "Monte Carlo files:  " << endl;
    cout << JPsiFile << "\n" << psiFile << "\n" << OCFile << "\n" << AMDYFile;
    cout << "\nphysBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << LR_Mmin << "  -  "
	 << LR_Mmax << endl;
    cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
    }//*/

  //Get Input Files from Local_leftRight
  TFile *fRD  = OpenFile(pathData + RDfile);
  TFile *fJPsi = OpenFile(pathData + JPsiFile);
  TFile *fpsi = OpenFile(pathData + psiFile);
  TFile *fOC = OpenFile(pathData + OCFile);
  TFile *fAMDY = OpenFile(pathData + AMDYFile);

  //Determine polarization factor
  Double_t ex[nBins] = {0.0};
  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  TGraph* g_Dil =(TGraph*)fRD->Get(Form("%s_Dil", physBinned.Data()));
  Double_t *xvals = g_Pol->GetX();
  Double_t Pol[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);

  //Get Input Hist/fit ratio/Determine LR counts for specified process
  const Int_t nTargPol =4;
  TH1D *hRD[nBins*nTargPol*2], *hjpsi[nBins*nTargPol*2],*hpsi[nBins*nTargPol*2];
  TH1D *hOC[nBins*nTargPol*2], *hAMDY[nBins*nTargPol*2];
  TH1D *hRatio[nBins*nTargPol*2];

  //TargPol setup
  TString targNames[nTargPol] = {"upstream_up", "upstream_down",
				 "downstream_up","downstream_down"};

  //Parameter setup
  Double_t massJPsi =3.13, widthJPsi=0.17;
  Double_t alphaJPsi =1.5, nJPsi =1.5;
  Double_t Bg_slope =-2.3;
  Double_t DY_slope = -0.90, DY_Mmin = 4.5;

  //Drawn fit canvases
  TCanvas* cFit[nTargPol], *cRatio[nTargPol];
  for (Int_t c=0; c<nTargPol; c++) {
    cFit[c] = new TCanvas(targNames[c]);
    cFit[c]->Divide(2, nBins);//cleanup

    cRatio[c] = new TCanvas(Form("ratio_%s", targNames[c].Data() ));
    cRatio[c]->Divide(2, nBins);//cleanup
  }

  //Left/Right
  Double_t leftRight[nBins*nTargPol*2], eLeftRight[nBins*nTargPol*2];

  //JPsi purity
  TH1D* h_purity = new TH1D("h_purity", "h_purity", 30, 0.6,1); SetUp(h_purity);

  //Quality parameters
  Int_t NDF;//cleanup
  Double_t ndfEst = (1.0*hbins-4.0)*(Mmax - Mmin)/(8.5-1.0);
  TH1D *hChi2 = new TH1D("hChi2", "hChi2", 10, 0, 15); SetUp(hChi2);

  //Perform Fits and integrate to get L/R
  for (Int_t targ=0, iter=0; targ<nTargPol; targ++) {
    for (Int_t l=0; l<2; l++) {//*/
      for (Int_t bi=0; bi<nBins; bi++, iter++) {
      
      /*for (Int_t targ=2, iter=0; targ<nTargPol-1; targ++) { cout << "debug" <<endl;
	for (Int_t bi=0; bi<1; bi++) {
	for (Int_t l=1; l<2; l++, iter++) {//*/

	TString histName;
	if (l==0){
	  histName = Form("MuMu_left_%s_%s%i", targNames[targ].Data(),
			  physBinned.Data(), bi);	  
	}
	else{
	  histName = Form("MuMu_right_%s_%s%i", targNames[targ].Data(),
			  physBinned.Data(), bi);
	}
	
	hRD[iter] = (TH1D*)fRD->Get(histName);
	hjpsi[iter] = (TH1D*)GetMCHist(fJPsi, histName, 0, hbins);
	hpsi[iter] = (TH1D*)GetMCHist(fpsi, histName, 1, hbins);
	hOC[iter] = (TH1D*)GetMCHist(fOC, histName, 2, hbins);
	hAMDY[iter] = (TH1D*)GetMCHist(fAMDY, histName, 3, hbins);

	hRD[iter]->Sumw2();
	SetUp(hRD[iter]); SetUp(hjpsi[iter]);
	SetUp(hpsi[iter]); SetUp(hOC[iter]); SetUp(hAMDY[iter]);

	//Get starting amplitudes
	Int_t binJPsi = hjpsi[iter]->GetMaximumBin();
	Double_t A_JPsi = hRD[iter]->GetBinContent(binJPsi);
	A_JPsi /= hjpsi[iter]->GetBinContent(binJPsi);

	Double_t A_psi = hRD[iter]->GetBinContent(binJPsi)/30.0;
	Int_t binpsi = hpsi[iter]->GetMaximumBin();
	A_psi /= hpsi[iter]->GetBinContent(binpsi);

	Int_t binOC = hOC[iter]->GetMaximumBin();
	Double_t A_OC = hRD[iter]->GetBinContent(binOC);
	A_OC /= hOC[iter]->GetBinContent(binOC);

	Int_t binDY = hAMDY[iter]->FindBin(4.7);
	Double_t A_DY = hRD[iter]->GetBinContent(binDY);
	A_DY /= hAMDY[iter]->GetBinContent(binDY);
	
	TF1 *fit = new TF1("Sum", Fit_MC, Mmin, Mmax, 4);
	Double_t factor =100.0;
	fit->SetParameters(A_JPsi, A_psi, A_OC, A_DY);
	fit->SetParLimits(0, 0, A_JPsi*factor);
	fit->SetParLimits(1, 0, A_psi*factor);
	fit->SetParLimits(2, 0, A_OC*factor);
	fit->SetParLimits(3, 0, A_DY*factor);	
	

	cFit[targ]->cd(2*bi + l + 1);
	gPad->SetLogy();
	hRD[iter]->Draw();
	hRD[iter]->Fit("Sum", "WLRQ0");

	//Draw result distributons
	Double_t *pars = fit->GetParameters();
	const Double_t *epars = fit->GetParErrors();
	hjpsi[iter]->Scale(pars[0]);
	hpsi[iter]->Scale(pars[1]);
	hOC[iter]->Scale(pars[2]);
	hAMDY[iter]->Scale(pars[3]);
	
	hjpsi[iter]->Draw("same"); hjpsi[iter]->SetLineColor(kGreen);
	hpsi[iter]->Draw("same"); hpsi[iter]->SetLineColor(kRed);
	hOC[iter]->Draw("same"); hOC[iter]->SetLineColor(kBlack);
	hAMDY[iter]->Draw("same"); hAMDY[iter]->SetLineColor(kYellow);
	TH1D *hTotal = (TH1D*) hjpsi[iter]->Clone("hTotal");
	hTotal->Add(hpsi[iter]);
	hTotal->Add(hOC[iter]);
	hTotal->Add(hAMDY[iter]);
	hTotal->Draw("same"); hTotal->SetLineColor(kRed);

	//Get JPsi Purity
	Int_t lrMinBin = hRD[iter]->FindBin(LR_Mmin);
	Int_t lrMaxBin = hRD[iter]->FindBin(LR_Mmax);
	Double_t lrRangeJPsi = hjpsi[iter]->Integral(lrMinBin, lrMaxBin);
	Double_t lrRangepsi = hpsi[iter]->Integral(lrMinBin, lrMaxBin);
	Double_t lrRangeOC = hOC[iter]->Integral(lrMinBin, lrMaxBin);
	Double_t lrRangeDY = hAMDY[iter]->Integral(lrMinBin, lrMaxBin);
	Double_t lrRangeTotal = lrRangeJPsi + lrRangepsi + lrRangeOC +lrRangeDY;
	h_purity->Fill( lrRangeJPsi/lrRangeTotal );
	
	//Left/Right
	leftRight[iter] = pars[0];
	eLeftRight[iter] = epars[0];

	//Quality checks
	Double_t chi2 = fit->GetChisquare();
	Int_t ndf = fit->GetNDF();

	hChi2->Fill( 1.0*chi2/ndf );

	BasicChecks(pars, epars);

	hRatio[iter] =
	  (TH1D*)hRD[iter]->Clone(Form("ratio_%s", histName.Data()));
	hRatio[iter]->Divide(hTotal);
	hRatio[iter]->GetYaxis()->SetRangeUser(0.5, 1.5);
	cRatio[targ]->cd(2*bi + l + 1);
	hRatio[iter]->Draw();
	DrawLine(hRatio[iter], 1.0);
    }//left/right
    }//targ
  }//Loop nBins

  //Left/right plots
  TGraphErrors* g_Left_upS_up =
    new TGraphErrors(nBins, xvals, leftRight, ex, eLeftRight);
  TGraphErrors* g_Right_upS_up =
    new TGraphErrors(nBins, xvals, &leftRight[nBins], ex, &eLeftRight[nBins]);
  TGraphErrors* g_Left_upS_down =
    new TGraphErrors(nBins, xvals, &leftRight[nBins*2], ex,
		     &eLeftRight[nBins*2]);
  TGraphErrors* g_Right_upS_down =
    new TGraphErrors(nBins, xvals, &leftRight[nBins*3], ex,
		     &eLeftRight[nBins*3]);
  TGraphErrors* g_Left_downS_up =
    new TGraphErrors(nBins, xvals, &leftRight[nBins*4], ex,
		     &eLeftRight[nBins*4]);
  TGraphErrors* g_Right_downS_up =
    new TGraphErrors(nBins, xvals, &leftRight[nBins*5], ex,
		     &eLeftRight[nBins*5]);
  TGraphErrors* g_Left_downS_down =
    new TGraphErrors(nBins, xvals, &leftRight[nBins*6], ex,
		     &eLeftRight[nBins*6]);
  TGraphErrors* g_Right_downS_down =
  new TGraphErrors(nBins, xvals, &leftRight[nBins*7], ex,
		     &eLeftRight[nBins*7]);

  SetUp(g_Left_upS_up); SetUp(g_Right_upS_up);
  SetUp(g_Left_upS_down); SetUp(g_Right_upS_down);
  SetUp(g_Left_downS_up); SetUp(g_Right_downS_up);
  SetUp(g_Left_downS_down); SetUp(g_Right_downS_down);

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

  //Purity
  TCanvas* cPur = new TCanvas();
  h_purity->Draw();
  
  //Quality output
  TCanvas* cChi2 = new TCanvas();
  hChi2->Draw();

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/Purity/Data/Purity";
  TString fOutput
    = Form("%s/mcMFit_%.2f_%.2f_%s_%s%.2f_%.2f_%s%s%i_%ihbin.root",
	   thisDirPath.Data(), Mmin, Mmax,
	   period_Mtype.Data(), process.Data(), LR_Mmin, LR_Mmax,
	   binRange.Data(), physBinned.Data(), nBins, hbins);
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");

    h_purity->Write();
  }

  cout << " " << endl;
  cout << "Settings______" << endl;
  cout << "Real data file considered:              " << RDfile << endl;
  cout <<"Number of histogram bins used in M fit:  " << hbins << endl;
  cout << "Mass Range for fitting is:              " << Mmin
       << " - " << Mmax << endl;
  cout << "physBinned nBins times:                 " << nBins << endl;
  cout << "\n";
  cout << "AN for physical process:                " << process << endl;
  cout << "Physics binning is:                     " << physBinned << endl;
  cout << "LR integral mass range:                 " << LR_Mmin
       << " - " << LR_Mmax << endl;
  cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;//*/
}


TH1D* GetMCHist(TFile *f, TString name, Int_t iMC, Int_t hbins){
  TH1D *h = (TH1D*)f->Get(name);
  if (iMC != 3) //Removes warnings
    h->Sumw2();
  
  h->Scale(1.0/(h->Integral() ) );

  hist[iMC] = (TH1D*)h->Clone(Form("%i%s", iMC, name.Data()));

  //Basic setup check
  if (h->GetNbinsX() != hbins){
    cout << "MC with wrong number of bins" << endl;
    exit(EXIT_FAILURE);
  }

  return h;
}


Int_t BasicChecks(Double_t *pars, const Double_t *epars){
  Int_t failure = 0;
  if ( pars[0] < 0 ) {//Parameter Checks
    cout << pars[0] << " A_JPsi" << endl;
    failure++;
  }
  if ( pars[1] < 0 ){
    cout << pars[1] << " A_psi" << endl;
    failure++;
  }
  if ( pars[2] < 0 ){
    cout << pars[2] << " A_OC" << endl;
    failure++;
  }
  if ( pars[3] < 0 ){
    cout << pars[3] << " A_DY" << endl;
    failure++;
  }

  Double_t factor = 5.0;
  if ( epars[0] > TMath::Sqrt(pars[0])*factor ){//Error checks
    cout << epars[0] << " " << TMath::Sqrt(pars[0]) << " eJPsi" << endl;
    failure++;
  }
  if ( epars[1] > TMath::Sqrt(pars[1])*factor ){
    cout << epars[1] << " " << TMath::Sqrt(pars[1]) << " epsi" << endl;
    failure++;
  }
  if ( epars[2] > TMath::Sqrt(pars[2])*factor ){
    cout << epars[2] << " " << TMath::Sqrt(pars[2]) << " eOC" << endl;
    failure++;
  }
  if ( epars[3] > TMath::Sqrt(pars[3])*factor ){
    cout << epars[3] << " " << TMath::Sqrt(pars[3]) << " eDY" << endl;
    failure++;
  }

  return failure;
}
