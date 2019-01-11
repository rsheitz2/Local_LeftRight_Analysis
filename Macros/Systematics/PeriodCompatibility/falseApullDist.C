#include "include/helperFunctions.h"


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


Double_t CalSysErrorPeriod(TH1D*h){
  Double_t mean = h->GetFunction("gaus")->GetParameter(1);
  Double_t sigma = h->GetFunction("gaus")->GetParameter(2);
  if (mean < 0) mean *= -1.0;
  Double_t sigma2_m1 = sigma*sigma - 1;
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
}


void falseApullDist(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  /*const Int_t nBins =5;//JPsi
    TString fitMrangeType ="LowM_AMDY";
    Int_t hbins =150;
    TString process ="JPsi";//JPsi, psi
    TString lrMrange ="2.00_5.00";
    TString fitMrange ="2.00_8.50";
    TString binRange ="25_43";
    TString whichFit ="thirteen";
    TString production ="slot1";//"t3", "slot1"
    TString additionalCuts ="phiS0.53";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/TargFlip";
  const Int_t nPeriod = 9;
  const Int_t nPhysBinned =5;
  TString period[nPeriod] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "M"};
    
  if (start==""){
    cout <<"\nThis Script computes the pull distributions for an L/R asymmetry";
    cout <<"\nThe asymmetry is determined ";
    cout << "(# periods)*(# physics kinematics)*(# physics kinematic bins)\n";
    cout << "Input data should be a TGraphErrors for each period and";
    cout << "physical kinematic";
    cout << "\n\nCurrent Setup:" << endl;
    cout << "Data comes from:                 " << pathFA << endl;
    cout << "Number of periods considered:    " << nPeriod << endl;
    cout << "Number of physics kinematic bins " << nBins << endl;
    cout << "Mass type considered:            " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Binned in Mass range:       " << binRange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Fit considered:     " << whichFit << endl;
    cout << "Production considered:   " << production << endl;
    cout << "additional cuts considered:   " << additionalCuts << endl;
    cout << "\nOutput is to be written:     " << toWrite << endl;
    cout << "\nUsage:" << endl;
    cout <<"root \'PullDist(1)\'\n" << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data by physBinned
  TString docName=""; TString physBinnedNames =""; TString fitNames="";
  vector<Double_t> faPol, e_faPol;
  vector<Double_t> fa2Targ_upS, e_fa2Targ_upS, fa2Targ_downS, e_fa2Targ_downS;
  Double_t wAvg_faPol =0.0, wSigma2_faPol =0.0;
  Double_t wAvg_fa2Targ_upS =0.0, wSigma2_fa2Targ_upS =0.0;
  Double_t wAvg_fa2Targ_downS =0.0, wSigma2_fa2Targ_downS =0.0;
  Double_t wAvgPhys_faPol[nPhysBinned] ={0.0};
  Double_t wSigma2Phys_faPol[nPhysBinned] ={0.0};
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    for (Int_t p=0; p<nPeriod; p++) {
      TString fInput;
      if (whichFit == "true" ){
	if (lrMrange != fitMrange){
	  cout << "Error: lrMrange != fitMrange with whichFit==true" << endl;
	  exit(EXIT_FAILURE);
	}
      
	fInput =
	  Form("%s/falseGeoMean4Targ_true_%s_%s_%s%s_%s%s%i_%s_%s.root",
	       pathFA.Data(), period[p].Data(), fitMrangeType.Data(),
	       process.Data(), fitMrange.Data(), binRange.Data(),
	       physBinned[phys].Data(), nBins, production.Data(),
	       additionalCuts.Data());
      }
      else{
	cout << "doesn't work for now" << endl;
	exit(EXIT_FAILURE);
	fInput =
	  Form("%s/GeoMean4Targ/\
GeoMean4Targ_%s%s_%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
	       pathFA.Data(), whichFit.Data(), fitMrange.Data(),
	       period[p].Data(), fitMrangeType.Data(), process.Data(),
	       lrMrange.Data(), physBinned[phys].Data(), nBins, hbins,
	       production.Data(), additionalCuts.Data());
      }
    
      TFile *f_FA = OpenFile(fInput);
      docName += fInput+"\n";
      physBinnedNames += physBinned[phys]+" "; fitNames += whichFit+" ";

      TGraphErrors *g_faPol = (TGraphErrors*)f_FA->Get("falseAN_pol");
      PushBackValues(g_faPol, faPol, e_faPol, wAvg_faPol, wSigma2_faPol, nBins);
      PushBackValues(g_faPol, wAvgPhys_faPol[phys], wSigma2Phys_faPol[phys],
		     nBins);

      TGraphErrors *g_fa2Targ_upS =
	(TGraphErrors*)f_FA->Get("falseAN_2Targ_upS");
      PushBackValues(g_fa2Targ_upS, fa2Targ_upS, e_fa2Targ_upS,
		     wAvg_fa2Targ_upS, wSigma2_fa2Targ_upS, nBins);

      TGraphErrors *g_fa2Targ_downS =
	(TGraphErrors*)f_FA->Get("falseAN_2Targ_downS");
      PushBackValues(g_fa2Targ_downS, fa2Targ_downS, e_fa2Targ_downS,
		     wAvg_fa2Targ_downS, wSigma2_fa2Targ_downS, nBins);
    }//period loop
  }//nPhysBinned loop

  wAvg_faPol = wAvg_faPol/wSigma2_faPol;
  wSigma2_faPol = 1.0/wSigma2_faPol;
  wAvg_fa2Targ_upS = wAvg_fa2Targ_upS/wSigma2_fa2Targ_upS;
  wSigma2_fa2Targ_upS = 1.0/wSigma2_fa2Targ_upS;
  wAvg_fa2Targ_downS = wAvg_fa2Targ_downS/wSigma2_fa2Targ_downS;
  wSigma2_fa2Targ_downS = 1.0/wSigma2_fa2Targ_downS;

  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    wAvgPhys_faPol[phys] = wAvgPhys_faPol[phys]/wSigma2Phys_faPol[phys];
    wSigma2Phys_faPol[phys] = 1.0/wSigma2Phys_faPol[phys];
  }

  //Make and fill pull dist
  gStyle->SetOptStat(11); gStyle->SetOptFit(111);
  TH1D* hPull_faPol = new TH1D("hPull_faPol", "hPull_faPol", 25, -4, 4);
  FillPull(hPull_faPol, faPol, e_faPol, wAvg_faPol, wSigma2_faPol);
  TH1D* hPulls_faPol[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    /*if (phys==nPhysBinned-1){
      hPulls_faPol[phys] = new TH1D("hPull_int", "hPull_int", 5, -4, 4);
    }
    else{
      hPulls_faPol[phys] = new TH1D(Form("hPull_%s", physBinned[phys].Data()),
			      Form("hPull_%s", physBinned[phys].Data()),
			      12, -4, 4);
			      }//*/
    hPulls_faPol[phys] = new TH1D(Form("hPull_%s", physBinned[phys].Data()),
			      Form("hPull_%s", physBinned[phys].Data()),
			      12, -4, 4);
    SetUp(hPulls_faPol[phys]);
  }
  Int_t whichPhys=0, count=0;
  for(vector<Double_t>::iterator it=faPol.begin(), e_it=e_faPol.begin();
      it!=faPol.end(); it++, e_it++){
    Double_t pull = *it - wAvgPhys_faPol[whichPhys]; //Per physics pulls
    Double_t denom = (*e_it)*(*e_it) - wSigma2Phys_faPol[whichPhys];

    if (denom < 0.0) denom *= -1.0;
    denom = TMath::Sqrt(denom);

    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    hPulls_faPol[whichPhys]->Fill(pull);
    
    if (count == nPeriod*nBins -1){
      count =0;
      whichPhys++;
    }
    else count++;
  }
  
  TH1D* hPull_fa2Targ_upS = new TH1D("hPull_fa2Targ_upS", "hPull_fa2Targ_upS",
				     25, -4, 4);
  FillPull(hPull_fa2Targ_upS, fa2Targ_upS, e_fa2Targ_upS, wAvg_fa2Targ_upS,
	   wSigma2_fa2Targ_upS);

  TH1D* hPull_fa2Targ_downS =
    new TH1D("hPull_fa2Targ_downS", "hPull_fa2Targ_downS", 25, -4, 4);
  FillPull(hPull_fa2Targ_downS, fa2Targ_downS, e_fa2Targ_downS,
	   wAvg_fa2Targ_downS, wSigma2_fa2Targ_downS);
  
  //Draw Pulls
  TCanvas* cPull = new TCanvas(); cPull->Divide(3);
  cPull->cd(1); SetUpPull(hPull_faPol);
  cPull->cd(2); SetUpPull(hPull_fa2Targ_upS);
  cPull->cd(3); SetUpPull(hPull_fa2Targ_downS);

  //Pull by kinematic
  TCanvas* cPulls_faPol = new TCanvas();
  cPulls_faPol->Divide(nPhysBinned, 1, 0, 0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cPulls_faPol->cd(i+1);
    gPad->SetFrameLineWidth(2);

    SetUpPull(hPulls_faPol[i]);
  }

  //Systematic errors
  Double_t error_faPol = CalSysErrorPeriod(hPull_faPol);
  Double_t error_fa2Targ_upS = CalSysErrorPeriod(hPull_fa2Targ_upS);
  Double_t error_fa2Targ_downS = CalSysErrorPeriod(hPull_fa2Targ_downS);
  Double_t sysError_faPol[nBins], sysError_fa2Targ_upS[nBins];
  Double_t sysError_fa2Targ_downS[nBins];
  Double_t  xvals[nBins];
  for (Int_t i=0; i<nBins; i++) {
    sysError_faPol[i] = error_faPol;
    sysError_fa2Targ_upS[i] = error_fa2Targ_upS;
    sysError_fa2Targ_downS[i] = error_fa2Targ_downS;
    
    xvals[i] = 1 + i;
  }

  //Sys error by kinematics
  Double_t sysErrorPhys_faPol[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    sysErrorPhys_faPol[i] = CalSysErrorPeriod(hPulls_faPol[i]);
  }
  
  //Draw systematic error
  TGraph *gSys_faPol = new TGraphErrors(nBins, xvals, sysError_faPol);
  SetUp(gSys_faPol);
  gSys_faPol->GetYaxis()->SetRangeUser(0.0, 1.0);
  TGraph *gSys_fa2Targ_upS =
    new TGraphErrors(nBins, xvals, sysError_fa2Targ_upS);
  SetUp(gSys_fa2Targ_upS);
  TGraph *gSys_fa2Targ_downS =
    new TGraphErrors(nBins, xvals, sysError_fa2Targ_downS);
  SetUp(gSys_fa2Targ_downS);
  
  TCanvas* cSys = new TCanvas();
  gSys_faPol->Draw("AP");
  gSys_faPol->SetTitle("Systematic Error/Statistical Error faPol");
  gSys_fa2Targ_upS->Draw("Psame"); gSys_fa2Targ_upS->SetMarkerColor(kRed);
  gSys_fa2Targ_downS->Draw("Psame"); gSys_fa2Targ_downS->SetMarkerColor(kBlue);

  TCanvas* cSysPhys_faPol = new TCanvas();
  cSysPhys_faPol->Divide(nPhysBinned, 1, 0, 0.01);
  TGraph *gSysPhys_faPol[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    gSysPhys_faPol[i] = new TGraphErrors(1, xvals, &(sysErrorPhys_faPol[i]));
    SetUp(gSysPhys_faPol[i]);
    gSysPhys_faPol[i]->GetYaxis()->SetRangeUser(0.0, 2.0);

    cSysPhys_faPol->cd(i+1);
    gSysPhys_faPol[i]->Draw("AP");
    gSysPhys_faPol[i]->SetTitle("");
  }

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Data";
  /*TString fOutput;
  if  (lrMrange == fitMrange){
    fOutput =
      Form("%s/pullDist/pullDist_true_%s_%s%s_%s_%i_%s_%s.root",
	   thisDirPath.Data(), fitMrangeType.Data(), process.Data(),
	   fitMrange.Data(), binRange.Data(), nBins, production.Data(),
	   additionalCuts.Data());
  }
  else {
    cout << "Need to fix this file naming..." << endl;
    exit(EXIT_FAILURE);
    fOutput =
      Form("%s/pullDist/pullDist_%s%s_%s_%s%s%i_%ihbin_%s_%s.root",
	   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
	   fitMrangeType.Data(), process.Data(), lrMrange.Data(), nBins, hbins,
	   production.Data(), additionalCuts.Data());
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TNamed docs ("InputData", docName.Data());
    TNamed pBinNam ("physBinned", physBinnedNames.Data());
    TNamed fitNam ("fitNames", fitNames.Data());
    docs.Write();
    pBinNam.Write();
    fitNam.Write();
    
    hPull_faPol->Write();

    gSys_faPol->Write("gSys_faPol");
  }

  if (start!=1){
    cout << "\nSettings______" << endl;
    cout <<"\nThis Script computes the pull distributions for an L/R asymmetry";
    cout <<"\nThe asymmetry is determined ";
    cout << "(# periods)*(# physics kinematics)*(# physics kinematic bins)\n";
    cout << "Input data should be a TGraphErrors for each period and";
    cout << "physical kinematic";
    cout << "\n\nCurrent Setup:" << endl;
    cout << "Data comes from:                 " << pathFA << endl;
    cout << "Number of periods considered:    " << nPeriod << endl;
    cout << "Number of physics kinematic bins " << nBins << endl;
    cout << "Mass type considered:            " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Binned in Mass range:       " << binRange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Fit considered:     " << whichFit << endl;
    cout << "Production considered:   " << production << endl;
    cout << "additional cuts considered:   " << additionalCuts << endl;
    cout << "\nOutput is to be written:     " << toWrite << endl;
    cout << "\nUsage:" << endl;
    cout <<"root \'PullDist(1)\'\n" << endl;
  }
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;//*/
}
