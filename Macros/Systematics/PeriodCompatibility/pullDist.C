#include "include/helperFunctions.h"

Double_t CalSysErrorPeriod(Double_t mean, Double_t sigma){
  if (mean < 0) mean *= -1.0;
  Double_t sigma2_m1 = sigma*sigma - 1;
  if (sigma2_m1 < 0) sigma2_m1 *= -1.0;
  
  return TMath::Sqrt(sigma2_m1) + mean/2.0;
}


void pullDist(TString start=""){
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

  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data";
  const Int_t nPeriod = 9;
  const Int_t nPhysBinned =6;
  TString period[nPeriod] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "M", "xN"};
  //Last xN for integrated
    
  if (start==""){
    cout <<"\nThis Script computes the pull distributions for an L/R asymmetry";
    cout <<"\nThe asymmetry is determined ";
    cout << "(# periods)*(# physics kinematics)*(# physics kinematic bins)\n";
    cout << "Input data should be a TGraphErrors for each period and";
    cout << "physical kinematic";
    cout << "\n\nCurrent Setup:" << endl;
    cout << "Data comes from:                 " << pathAN << endl;
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
  vector<Double_t> lr, e_lr;
  Double_t wAvg =0.0, wSigma2 =0.0;
  Double_t wAvgPhys[nPhysBinned] ={0.0}, wSigma2Phys[nPhysBinned] ={0.0};
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    for (Int_t p=0; p<nPeriod; p++) {

      Int_t nBinsName =nBins;
      if (phys == nPhysBinned-1) {//used for integrated
	nBinsName =1; }
      TString fInput;
      if (whichFit == "true" ){
	if (lrMrange != fitMrange){
	  cout << "Error: lrMrange != fitMrange with whichFit==true" << endl;
	  exit(EXIT_FAILURE);
	}
      
	fInput =
	  Form("%s/GeoMean4Targ/GeoMean4Targ_true_%s_%s_%s%s_%s%s%i_%s_%s.root",
	       pathAN.Data(), period[p].Data(), fitMrangeType.Data(),
	       process.Data(), fitMrange.Data(), binRange.Data(),
	       physBinned[phys].Data(), nBinsName, production.Data(),
	       additionalCuts.Data());
      }
      else{
	cout << "doesn't work for now" << endl;
	exit(EXIT_FAILURE);
	fInput =
	  Form("%s/GeoMean4Targ/\
GeoMean4Targ_%s%s_%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
	       pathAN.Data(), whichFit.Data(), fitMrange.Data(),
	       period[p].Data(), fitMrangeType.Data(), process.Data(),
	       lrMrange.Data(), physBinned[phys].Data(), nBinsName, hbins,
	       production.Data(), additionalCuts.Data());
      }
    
      TFile *f_AN = TFile::Open(fInput);
      if ( !f_AN){
	cout << "File does not exist:  # " << phys << endl;
	exit(EXIT_FAILURE);
      }
      docName += fInput+"\n";
      physBinnedNames += physBinned[phys]+" "; fitNames += whichFit+" ";

      TGraphErrors *g_AN = (TGraphErrors*)f_AN->Get("AN");
      Double_t *yVal=g_AN->GetY();
      Double_t *e_yVal=g_AN->GetEY();
      
      for (Int_t bi=0; bi<nBinsName; bi++) {
	lr.push_back(yVal[bi]);
	e_lr.push_back(e_yVal[bi]);
	  
	if (phys!=nPhysBinned-1){//Do not include integrated in full Pull
	  wAvg += yVal[bi]/(e_yVal[bi]*e_yVal[bi]);
	  wSigma2 += 1.0/(e_yVal[bi]*e_yVal[bi]);
	}

	wAvgPhys[phys] += yVal[bi]/(e_yVal[bi]*e_yVal[bi]);
	wSigma2Phys[phys] += 1.0/(e_yVal[bi]*e_yVal[bi]);
      }

    }//period loop
  }//nPhysBinned loop

  wAvg = wAvg/wSigma2;
  wSigma2 = 1.0/wSigma2;
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    wAvgPhys[phys] = wAvgPhys[phys]/wSigma2Phys[phys];
    wSigma2Phys[phys] = 1.0/wSigma2Phys[phys];
  }

  //Make and fill pull dist
  //gStyle->SetOptFit(1111);
  gStyle->SetOptFit(111); gStyle->SetOptStat(10);
  TH1D* hPull = new TH1D("hPull", "hPull", 25, -4, 4); SetUp(hPull);
  TH1D* hPulls[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    if (phys==nPhysBinned-1){
      hPulls[phys] = new TH1D("hPull_int", "hPull_int", 5, -4, 4);
    }
    else{
      hPulls[phys] = new TH1D(Form("hPull_%s", physBinned[phys].Data()),
			      Form("hPull_%s", physBinned[phys].Data()),
			      12, -4, 4);
    }
    SetUp(hPulls[phys]);
  }

  Int_t whichPhys=0, count=0;
  for(vector<Double_t>::iterator it=lr.begin(), e_it=e_lr.begin(); it!=lr.end();
      it++, e_it++){
    Double_t pull = *it - wAvg;
    Double_t denom = (*e_it)*(*e_it) - wSigma2;
    
    if (denom < 0.0) denom *= -1.0;
    denom = TMath::Sqrt(denom);

    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    if (whichPhys!=nPhysBinned-1)//Do not include integrated in full Pull
      hPull->Fill(pull);

    //Per physics pulls
    pull = *it - wAvgPhys[whichPhys];
    denom = (*e_it)*(*e_it) - wSigma2Phys[whichPhys];

    if (denom < 0.0) denom *= -1.0;
    denom = TMath::Sqrt(denom);

    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    hPulls[whichPhys]->Fill(pull);
    
    if (count == nPeriod*nBins -1){
      count =0;
      whichPhys++;
    }
    else count++;
  }

  //Draw Pull
  TCanvas* cPull = new TCanvas();
  hPull->SetMarkerStyle(20);
  hPull->SetMarkerColor(kBlue);
  hPull->SetLineColor(kBlue);
  hPull->Draw("E X0");
  hPull->Fit("gaus", "Q");
  hPull->GetFunction("gaus")->SetLineStyle(7);
  
  //Pull by kinematic
  TCanvas* cPulls = new TCanvas(); cPulls->Divide(nPhysBinned, 1, 0, 0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cPulls->cd(i+1);
    gPad->SetFrameLineWidth(2);

    hPulls[i]->SetMarkerStyle(20);
    hPulls[i]->SetMarkerColor(kBlue);
    hPulls[i]->SetLineColor(kBlue);
    hPulls[i]->Draw("E X0");
    hPulls[i]->Fit("gaus", "Q");
    hPulls[i]->GetFunction("gaus")->SetLineStyle(7);
  }
  
  //Systematic Error
  TF1 *f_gaus = hPull->GetFunction("gaus");
  Double_t mean = f_gaus->GetParameter(1);
  Double_t sigma = f_gaus->GetParameter(2);

  Double_t error = CalSysErrorPeriod(mean, sigma);
  Double_t sysError[nBins], xvals[nBins];
  for (Int_t i=0; i<nBins; i++) {
    sysError[i] = error;
    xvals[i] = 1 + i;
  }

  Double_t sysErrorPhys[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    TF1 *fPhys_gaus = hPulls[i]->GetFunction("gaus");
    Double_t meanPhys = fPhys_gaus->GetParameter(1);
    Double_t sigmaPhys = fPhys_gaus->GetParameter(2);

    sysErrorPhys[i] = CalSysErrorPeriod(meanPhys, sigmaPhys);
  }
  
  //Draw systematic error
  TGraph *gSys = new TGraphErrors(nBins, xvals, sysError);
  SetUp(gSys); gSys->GetYaxis()->SetRangeUser(0.0, sysError[0]*1.3);

  TCanvas* cSys = new TCanvas();
  gSys->Draw("AP");
  gSys->SetTitle("Systematic Error/Statistical Error");

  TCanvas* cSysPhys = new TCanvas(); cSysPhys->Divide(nPhysBinned, 1, 0, 0.01);
  TGraph *gSysPhys[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    gSysPhys[i] = new TGraphErrors(1, xvals, &(sysErrorPhys[i]));
    SetUp(gSysPhys[i]);
    gSysPhys[i]->GetYaxis()->SetRangeUser(0.0, 2.0);

    cSysPhys->cd(i+1);
    gSysPhys[i]->Draw("AP");
    gSysPhys[i]->SetTitle("");
  }

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Data";
  TString fOutput;
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
    
    hPull->Write();
    gSys->Write("gSys");
    for (Int_t i=0; i<nPhysBinned; i++) {
      hPulls[i]->Write();

      if (i==nPhysBinned-1){
	gSysPhys[i]->Write("gSys_int");	
      }
      else{
	gSysPhys[i]->Write(Form("gSys_%s", physBinned[i].Data()));
      }
    }
    
    cPulls->Write("cPulls");
    cSysPhys->Write("cSysPhys");
  }

  if (start!=1){
    cout << "\nSettings______" << endl;
    cout <<"\nThis Script computes the pull distributions for an L/R asymmetry";
    cout <<"\nThe asymmetry is determined ";
    cout << "(# periods)*(# physics kinematics)*(# physics kinematic bins)\n";
    cout << "Input data should be a TGraphErrors for each period and";
    cout << "physical kinematic";
    cout << "\n\nCurrent Setup:" << endl;
    cout << "Data comes from:                 " << pathAN << endl;
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
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
