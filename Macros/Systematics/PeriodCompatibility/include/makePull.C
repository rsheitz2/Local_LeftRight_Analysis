#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/include/\
helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void PushBackValues(TGraphErrors* gFA, vector<Double_t> &val,
		    vector<Double_t> &e_val, Double_t &wAvg, Double_t &wSigma2,
		    Int_t nBins){
  
  Double_t *yVal=gFA->GetY();
  Double_t *e_yVal=gFA->GetEY();
  for (Int_t bi=0; bi<nBins; bi++) {
    val.push_back(yVal[bi]);
    e_val.push_back(e_yVal[bi]);

    wAvg += yVal[bi]/(e_yVal[bi]*e_yVal[bi]);
    wSigma2 += 1.0/(e_yVal[bi]*e_yVal[bi]);
  }
}


void PushBackValues(TGraphErrors* gFA, Double_t &wAvg, Double_t &wSigma2,
		    Int_t nBins){
  
  Double_t *yVal=gFA->GetY();
  Double_t *e_yVal=gFA->GetEY();
  for (Int_t bi=0; bi<nBins; bi++) {
    wAvg += yVal[bi]/(e_yVal[bi]*e_yVal[bi]);
    wSigma2 += 1.0/(e_yVal[bi]*e_yVal[bi]);
  }
}


void FillPull(TH1D* h, vector<Double_t> &val, vector<Double_t> &e_val,
	      Double_t wAvg, Double_t wSigma2){
  for(vector<Double_t>::iterator it=val.begin(), e_it=e_val.begin();
      it!=val.end(); it++, e_it++){
    Double_t pull = *it - wAvg;
    Double_t denom = (*e_it)*(*e_it) - wSigma2;
    
    if (denom < 0.0) denom *= -1.0;
    denom = TMath::Sqrt(denom);

    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    h->Fill(pull);
  }
}


void FillPull(TH1D** h, vector<Double_t> &val, vector<Double_t> &e_val,
	      Double_t *wAvg, Double_t *wSigma2, Int_t counts){

  Int_t whichPhys=0, count=0;
  for(vector<Double_t>::iterator it=val.begin(), e_it=e_val.begin();
      it!=val.end(); it++, e_it++){
    Double_t pull = *it - wAvg[whichPhys];
    Double_t denom = (*e_it)*(*e_it) - wSigma2[whichPhys];
    
    if (denom < 0.0) denom *= -1.0;
    denom = TMath::Sqrt(denom);

    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    h[whichPhys]->Fill(pull);

    if (count == counts -1){
      count =0;
      whichPhys++;
    }
    else count++;
  }
}


void AddWAvgSigm2(TH1D *h, Double_t &wMean, Double_t &wMeanSig2,
		  Double_t &wSig2, Double_t &wSig2Sig2){
  Double_t mean = h->GetFunction("gaus")->GetParameter(1);
  Double_t sig = h->GetFunction("gaus")->GetParameter(2);

  Double_t eMean = h->GetFunction("gaus")->GetParError(1);
  Double_t eSig = h->GetFunction("gaus")->GetParError(2);

  wMean += mean/(eMean*eMean);
  wMeanSig2 += 1.0/(eMean*eMean);

  wSig2 += sig*sig/(eSig*eSig);
  wSig2Sig2 += 1.0/(eSig*eSig);
}


void FinalWAvgs(Double_t &wMean, Double_t &wMeanSig2,
		Double_t &wSig2, Double_t &wSig2Sig2){
  wMean = wMean/wMeanSig2;
  wSig2 = wSig2/wSig2Sig2;

  wMeanSig2 = 1.0/wMeanSig2;
  wSig2Sig2 = 1.0/wSig2Sig2;

  cout << "\nwMean = " << wMean << " +/- " << TMath::Sqrt(wMeanSig2) << endl;
  cout << "wSig2 = " << wSig2 << " +/- " << TMath::Sqrt(wSig2Sig2) << endl;
}


Double_t CalSysErrorPeriod(TH1D*h){
  Double_t mean = h->GetFunction("gaus")->GetParameter(1);
  Double_t sigma = h->GetFunction("gaus")->GetParameter(2);
  if (mean < 0) mean *= -1.0;
  Double_t sigma2_m1 = sigma*sigma - 1;
  if (sigma2_m1 < 0) sigma2_m1 *= -1.0;
  
  return TMath::Sqrt(sigma2_m1) + mean/2.0;
}


Double_t CalSysErrorPeriod(Double_t mean, Double_t sigma2){
  if (mean < 0) mean *= mean;

  Double_t sigma2_m1 = sigma2 - 1;
  if (sigma2_m1 < 0) sigma2_m1 *= -1.0;

  return TMath::Sqrt(sigma2_m1) + mean/2.0;
}


void SetUpPull(TH1D *h){
  SetUp(h); 
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlue);
  h->SetLineColor(kBlue);
  h->Draw("E X0");
  h->Fit("gaus", "Q");
  h->GetFunction("gaus")->SetLineStyle(7);
  FinalSetup(h);
}


void makePull(TString fileNameStart, TString fileNameMiddle,
	      TString fileNameEnd, TString tgraphName, const Int_t nBins){

  const Int_t nPeriod = 9;
  const Int_t nPhysBinned =5;
  TString period[nPeriod] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "M"};

  //Aesthics
  gStyle->SetOptStat(11); gStyle->SetOptFit(111);
  cout << "Number of periods considered:    " << nPeriod << endl;
  TString xNames[nPhysBinned] =
    {"x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};

  //Get Data by physBinned
  vector<Double_t> pull, e_pull;
  Double_t wAvg_pull =0.0, wSigma2_pull =0.0;
  Double_t wAvgPhys_pull[nPhysBinned] ={0.0};
  Double_t wSigma2Phys_pull[nPhysBinned] ={0.0};
  
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    for (Int_t p=0; p<nPeriod; p++) {
      TFile *f_FA =
	OpenFile(Form("%s%s%s%s%s", fileNameStart.Data(), period[p].Data(),
		      fileNameMiddle.Data(), physBinned[phys].Data(),
		      fileNameEnd.Data()));
      TGraphErrors *g_pull = (TGraphErrors*)f_FA->Get(tgraphName);
      PushBackValues(g_pull, pull, e_pull, wAvg_pull, wSigma2_pull, nBins);
      PushBackValues(g_pull, wAvgPhys_pull[phys], wSigma2Phys_pull[phys],
		     nBins);
    }//period loop
  }//nPhysBinned loop

  wAvg_pull = wAvg_pull/wSigma2_pull;
  wSigma2_pull = 1.0/wSigma2_pull;

  //Kinematic averaging
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    wAvgPhys_pull[phys] = wAvgPhys_pull[phys]/wSigma2Phys_pull[phys];
    wSigma2Phys_pull[phys] = 1.0/wSigma2Phys_pull[phys];
  }

  //Make and fill pull dist
  TH1D* hPulls_pull[nPhysBinned], *hPulls_faSubper[nPhysBinned];;
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    hPulls_pull[phys] =
      new TH1D(Form("hPulls_%s", physBinned[phys].Data()),
	       Form("hPulls_%s", physBinned[phys].Data()),
	       12, -4, 4);
    SetUp(hPulls_pull[phys]);
  }
  
  TH1D* hPull_pull = new TH1D("hPull_pull", "hPull_pull", 20, -4, 4);
  FillPull(hPull_pull, pull, e_pull, wAvg_pull, wSigma2_pull);
  FillPull(hPulls_pull, pull, e_pull, wAvgPhys_pull, wSigma2Phys_pull,
	   nPeriod*nBins);
  
  //Draw Pulls
  TCanvas* cPull = new TCanvas("CorrelatedPulls"); 
  SetUpPull(hPull_pull);

  //Pull by kinematic
  TCanvas* cPulls_pull = new TCanvas("Pullspull");
  cPulls_pull->Divide(nPhysBinned, 1, 0.00001, 0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cPulls_pull->cd(i+1);
    gPad->SetFrameLineWidth(2);
    SetUpPull(hPulls_pull[i]);
    FinalClearTitles(hPulls_pull[i]);
    SetTitleName(hPulls_pull[i], xNames[i]);
  }

  //Systematic errors
  Double_t error_pull = CalSysErrorPeriod(hPull_pull);
  Double_t sysError_pull[nBins];
  Double_t  xvals[nBins];
  for (Int_t i=0; i<nBins; i++) {
    sysError_pull[i] = error_pull;
    
    xvals[i] = 1 + i;
  }

  //Sys error by kinematics
  Double_t sysErrorPhys_pull[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    sysErrorPhys_pull[i] = CalSysErrorPeriod(hPulls_pull[i]);
  }
  
  //Draw systematic error
  TGraph *gSys_pull = new TGraphErrors(nBins, xvals, sysError_pull);
  gSys_pull->GetYaxis()->SetRangeUser(0.0, 1.0);
  SetUp(gSys_pull); 
  
  TCanvas* cSys = new TCanvas("CorrelatedPullSys");
  gSys_pull->Draw("AP");
  gSys_pull->SetTitle("Systematic Error/Statistical Error pull");

  //Drawing sys error by kinematics
  TCanvas* cSysPhys_pull = new TCanvas("pullSys");
  cSysPhys_pull->Divide(nPhysBinned, 1, 0, 0.01);
  Double_t wMean_pull =0.0, wMeanSig_pull =0.0;
  Double_t wSig2_pull =0.0, wSig2Sig_pull =0.0;
  TGraph *gSysPhys_pull[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    gSysPhys_pull[i] = new TGraph(1, xvals, &(sysErrorPhys_pull[i]));
    SetUp(gSysPhys_pull[i]);
    gSysPhys_pull[i]->GetYaxis()->SetRangeUser(0.0, 2.5);

    cSysPhys_pull->cd(i+1);
    gSysPhys_pull[i]->Draw("AP");
    gSysPhys_pull[i]->SetTitle(Form("sysError_%s = %0.2f",
				     physBinned[i].Data(),
				     sysErrorPhys_pull[i]));

    AddWAvgSigm2(hPulls_pull[i], wMean_pull, wMeanSig_pull,
		 wSig2_pull, wSig2Sig_pull);
  }
  FinalWAvgs(wMean_pull, wMeanSig_pull, wSig2_pull, wSig2Sig_pull);

  Double_t wSysErrorPhys_pull = CalSysErrorPeriod(wMean_pull, wSig2_pull);
  TGraph *gWSysPhys_pull = new TGraph(1, xvals, &(wSysErrorPhys_pull));
  SetUp(gWSysPhys_pull);

  TCanvas* wSys = new TCanvas("wAvgSys");
  gWSysPhys_pull->Draw("AP"); gWSysPhys_pull->SetTitle("Weighted Sys Error");
}
