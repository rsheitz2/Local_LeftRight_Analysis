#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void PushBackValues(TGraphErrors* gFA, vector<Double_t> &val,
		    vector<Double_t> &e_val, Double_t &wAvg, Double_t &wSigma2,
		    Int_t nBins);

void PushBackValues(TGraphErrors* gFA, Double_t &wAvg, Double_t &wSigma2,
		    Int_t nBins);

void FillPull(TH1D* h, vector<Double_t> &val, vector<Double_t> &e_val,
	      Double_t wAvg, Double_t wSigma2);

void FillPull(TH1D** h, vector<Double_t> &val, vector<Double_t> &e_val,
	      Double_t *wAvg, Double_t *wSigma2, Int_t counts);

void AddWAvgSigm2(TH1D *h, Double_t &wMean, Double_t &wMeanSig2,
		  Double_t &wSig2, Double_t &wSig2Sig2);

void FinalWAvgs(Double_t &wMean, Double_t &wMeanSig2,
		Double_t &wSig2, Double_t &wSig2Sig2);

Double_t CalSysErrorPeriod(TH1D*h);

Double_t CalSysErrorPeriod(Double_t mean, Double_t sigma2);

void SetUpPull(TH1D *h);

void faPullDist2targWavg(TString start=""){
  //Setup_______________
  /*const Int_t nBins =3;//HMDY
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  const Int_t nBins =3;//JPsi
  TString fitMrangeType ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi
  TString lrMrange ="2.87_3.38";
  TString fitMrange ="2.87_3.38";
  TString binRange ="29_34";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathFA = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/TargFlip";
  const Int_t nPeriod = 9;
  const Int_t nPhysBinned =5;
  TString period[nPeriod] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "M"};

  //Aesthics
  gStyle->SetOptStat(10); gStyle->SetOptFit(11);
  TString xNames[nPhysBinned] =
    {"x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};
    
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
  vector<Double_t> faPol, e_faPol, faSubper, e_faSubper;
  vector<Double_t> fa2Targ_upS, e_fa2Targ_upS, fa2Targ_downS, e_fa2Targ_downS;
  Double_t wAvg_faPol =0.0, wSigma2_faPol =0.0;
  Double_t wAvg_faSubper =0.0, wSigma2_faSubper =0.0;
  Double_t wAvg_fa2Targ_upS =0.0, wSigma2_fa2Targ_upS =0.0;
  Double_t wAvg_fa2Targ_downS =0.0, wSigma2_fa2Targ_downS =0.0;
  Double_t wAvgPhys_faPol[nPhysBinned] ={0.0};
  Double_t wSigma2Phys_faPol[nPhysBinned] ={0.0};
  Double_t wAvgPhys_faSubper[nPhysBinned] ={0.0};
  Double_t wSigma2Phys_faSubper[nPhysBinned] ={0.0};
  Double_t wAvgPhys_fa2Targ_upS[nPhysBinned] ={0.0};
  Double_t wSigma2Phys_fa2Targ_upS[nPhysBinned] ={0.0};
  Double_t wAvgPhys_fa2Targ_downS[nPhysBinned] ={0.0};
  Double_t wSigma2Phys_fa2Targ_downS[nPhysBinned] ={0.0};
  
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

      TGraphErrors *g_faSubper = (TGraphErrors*)f_FA->Get("falseAN_subper");
      PushBackValues(g_faSubper, faSubper, e_faSubper, wAvg_faSubper,
		     wSigma2_faSubper, nBins);
      PushBackValues(g_faSubper, wAvgPhys_faSubper[phys],
		     wSigma2Phys_faSubper[phys], nBins);

      TGraphErrors *g_fa2Targ_upS =
	(TGraphErrors*)f_FA->Get("falseAN_2Targ_upS");
      PushBackValues(g_fa2Targ_upS, fa2Targ_upS, e_fa2Targ_upS,
		     wAvg_fa2Targ_upS, wSigma2_fa2Targ_upS, nBins);
      PushBackValues(g_fa2Targ_upS, wAvgPhys_fa2Targ_upS[phys],
		     wSigma2Phys_fa2Targ_upS[phys], nBins);

      TGraphErrors *g_fa2Targ_downS =
	(TGraphErrors*)f_FA->Get("falseAN_2Targ_downS");
      PushBackValues(g_fa2Targ_downS, fa2Targ_downS, e_fa2Targ_downS,
		     wAvg_fa2Targ_downS, wSigma2_fa2Targ_downS, nBins);
      PushBackValues(g_fa2Targ_downS, wAvgPhys_fa2Targ_downS[phys],
		     wSigma2Phys_fa2Targ_downS[phys], nBins);
    }//period loop
  }//nPhysBinned loop

  wAvg_faPol = wAvg_faPol/wSigma2_faPol;
  wSigma2_faPol = 1.0/wSigma2_faPol;
  wAvg_faSubper = wAvg_faSubper/wSigma2_faSubper;
  wSigma2_faSubper = 1.0/wSigma2_faSubper;
  wAvg_fa2Targ_upS = wAvg_fa2Targ_upS/wSigma2_fa2Targ_upS;
  wSigma2_fa2Targ_upS = 1.0/wSigma2_fa2Targ_upS;
  wAvg_fa2Targ_downS = wAvg_fa2Targ_downS/wSigma2_fa2Targ_downS;
  wSigma2_fa2Targ_downS = 1.0/wSigma2_fa2Targ_downS;

  //Kinematic averaging
  Double_t wAvgPhys_upSdownS[nPhysBinned], wSigma2Phys_upSdownS[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    wAvgPhys_faPol[phys] = wAvgPhys_faPol[phys]/wSigma2Phys_faPol[phys];
    wSigma2Phys_faPol[phys] = 1.0/wSigma2Phys_faPol[phys];

    wAvgPhys_faSubper[phys] =
      wAvgPhys_faSubper[phys]/wSigma2Phys_faSubper[phys];
    wSigma2Phys_faSubper[phys] = 1.0/wSigma2Phys_faSubper[phys];

    wAvgPhys_fa2Targ_upS[phys] =
      wAvgPhys_fa2Targ_upS[phys]/wSigma2Phys_fa2Targ_upS[phys];
    wSigma2Phys_fa2Targ_upS[phys] = 1.0/wSigma2Phys_fa2Targ_upS[phys];

    wAvgPhys_fa2Targ_downS[phys] =
      wAvgPhys_fa2Targ_downS[phys]/wSigma2Phys_fa2Targ_downS[phys];
    wSigma2Phys_fa2Targ_downS[phys] = 1.0/wSigma2Phys_fa2Targ_downS[phys];

    wAvgPhys_upSdownS[phys] =
      WeightedAvg(wAvgPhys_fa2Targ_upS[phys], wAvgPhys_fa2Targ_downS[phys],
		  TMath::Sqrt(wSigma2Phys_fa2Targ_upS[phys]),
		  TMath::Sqrt(wSigma2Phys_fa2Targ_downS[phys]));
    wSigma2Phys_upSdownS[phys] =
      WeightedErr(TMath::Sqrt(wSigma2Phys_fa2Targ_upS[phys]),
		  TMath::Sqrt(wSigma2Phys_fa2Targ_downS[phys]));
    wSigma2Phys_upSdownS[phys] *= wSigma2Phys_upSdownS[phys];
  }

  //Make and fill pull dist
  TH1D* hPulls_faPol[nPhysBinned], *hPulls_faSubper[nPhysBinned];;
  TH1D* hPulls_upS[nPhysBinned], *hPulls_downS[nPhysBinned];
  TH1D* hPulls_upSdownS[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    hPulls_faPol[phys] = new TH1D(Form("hPulls_%s", physBinned[phys].Data()),
			      Form("hPulls_%s", physBinned[phys].Data()),
			      12, -4, 4);
    hPulls_faSubper[phys] = new TH1D(Form("hPulls_%s_subper",
					  physBinned[phys].Data()),
				     Form("hPulls_%s_subper",
					  physBinned[phys].Data()),
				     12, -4, 4);
    hPulls_upS[phys] =
      new TH1D(Form("hPulls_upS_%s", physBinned[phys].Data()),
	       Form("hPulls_upS_%s", physBinned[phys].Data()),
	       6, -4, 4);
    hPulls_downS[phys] =
      new TH1D(Form("hPulls_downS_%s", physBinned[phys].Data()),
	       Form("hPulls_downS_%s", physBinned[phys].Data()),
	       6, -4, 4);
    hPulls_upSdownS[phys] =
      new TH1D(Form("hPulls_upSdownS_%s", physBinned[phys].Data()),
	       Form("hPulls_upSdownS_%s", physBinned[phys].Data()),
	       12, -4, 4);
    
    SetUp(hPulls_faPol[phys]); SetUp(hPulls_faSubper[phys]);
    SetUp(hPulls_upS[phys]); SetUp(hPulls_downS[phys]);
    SetUp(hPulls_upSdownS[phys]);
  }
  
  TH1D* hPull_faPol = new TH1D("hPull_faPol", "hPull_faPol", 20, -4, 4);
  FillPull(hPull_faPol, faPol, e_faPol, wAvg_faPol, wSigma2_faPol);
  FillPull(hPulls_faPol, faPol, e_faPol, wAvgPhys_faPol, wSigma2Phys_faPol,
	   nPeriod*nBins);

  TH1D* hPull_faSubper = new TH1D("hPull_faSubper", "hPull_faSubper",
				  20, -4, 4);
  FillPull(hPull_faSubper, faSubper, e_faSubper, wAvg_faSubper,
	   wSigma2_faSubper);
  FillPull(hPulls_faSubper, faSubper, e_faSubper, wAvgPhys_faSubper,
	   wSigma2Phys_faSubper, nPeriod*nBins);
  
  TH1D* hPull_fa2Targ_upS = new TH1D("hPull_fa2Targ_upS", "hPull_fa2Targ_upS",
				     20, -4, 4);
  FillPull(hPull_fa2Targ_upS, fa2Targ_upS, e_fa2Targ_upS, wAvg_fa2Targ_upS,
	   wSigma2_fa2Targ_upS);
  FillPull(hPulls_upS, fa2Targ_upS, e_fa2Targ_upS, wAvgPhys_fa2Targ_upS,
	   wSigma2Phys_fa2Targ_upS, nPeriod*nBins);

  TH1D* hPull_fa2Targ_downS =
    new TH1D("hPull_fa2Targ_downS", "hPull_fa2Targ_downS", 20, -4, 4);
  FillPull(hPull_fa2Targ_downS, fa2Targ_downS, e_fa2Targ_downS,
	   wAvg_fa2Targ_downS, wSigma2_fa2Targ_downS);
  FillPull(hPulls_downS, fa2Targ_downS, e_fa2Targ_downS,
	   wAvgPhys_fa2Targ_downS, wSigma2Phys_fa2Targ_downS, nPeriod*nBins);

  FillPull(hPulls_upSdownS, fa2Targ_upS, e_fa2Targ_upS,
	   wAvgPhys_upSdownS, wSigma2Phys_upSdownS, nPeriod*nBins);
  FillPull(hPulls_upSdownS, fa2Targ_downS, e_fa2Targ_downS,
	   wAvgPhys_upSdownS, wSigma2Phys_upSdownS, nPeriod*nBins);
  
  //Draw Pulls
  TCanvas* cPull = new TCanvas("CorrelatedPulls"); cPull->Divide(4);
  cPull->cd(1); SetUpPull(hPull_faPol);
  cPull->cd(2); SetUpPull(hPull_faSubper);
  cPull->cd(3); SetUpPull(hPull_fa2Targ_upS);
  cPull->cd(4); SetUpPull(hPull_fa2Targ_downS);

  //Pull by kinematic
  TCanvas* cPulls_faPol = new TCanvas("PullsfaPol");
  cPulls_faPol->Divide(nPhysBinned, 1, 0, 0.01);
  TCanvas* cPulls_faSubper = new TCanvas("PullsfaSubper");
  cPulls_faSubper->Divide(nPhysBinned, 1, 0, 0.01);
  TCanvas* cPulls_upS = new TCanvas("PullsUpS");
  cPulls_upS->Divide(nPhysBinned, 1, 0, 0.01);
  TCanvas* cPulls_downS = new TCanvas("PullsDownS");
  cPulls_downS->Divide(nPhysBinned, 1, 0, 0.01);
  TCanvas* cPulls_upSdownS = new TCanvas("PullsUpSdownS");
  cPulls_upSdownS->Divide(nPhysBinned, 1, 0.001, 0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cPulls_faPol->cd(i+1);
    gPad->SetFrameLineWidth(2);
    SetUpPull(hPulls_faPol[i]);

    cPulls_faSubper->cd(i+1);
    gPad->SetFrameLineWidth(2);
    SetUpPull(hPulls_faSubper[i]);

    cPulls_upS->cd(i+1);
    gPad->SetFrameLineWidth(2);
    SetUpPull(hPulls_upS[i]);

    cPulls_downS->cd(i+1);
    gPad->SetFrameLineWidth(2);
    SetUpPull(hPulls_downS[i]);

    cPulls_upSdownS->cd(i+1);
    gPad->SetFrameLineWidth(2);
    SetUpPull(hPulls_upSdownS[i]);
    FinalClearTitles(hPulls_upSdownS[i]);
    SetTitleName(hPulls_upSdownS[i], xNames[i]);
  }

  //Systematic errors
  Double_t error_faPol = CalSysErrorPeriod(hPull_faPol);
  Double_t error_faSubper = CalSysErrorPeriod(hPull_faSubper);
  Double_t error_fa2Targ_upS = CalSysErrorPeriod(hPull_fa2Targ_upS);
  Double_t error_fa2Targ_downS = CalSysErrorPeriod(hPull_fa2Targ_downS);
  Double_t sysError_faPol[nBins], sysError_faSubper[nBins];
  Double_t sysError_fa2Targ_upS[nBins], sysError_fa2Targ_downS[nBins];
  Double_t  xvals[nBins];
  for (Int_t i=0; i<nBins; i++) {
    sysError_faPol[i] = error_faPol;
    sysError_faSubper[i] = error_faSubper;
    sysError_fa2Targ_upS[i] = error_fa2Targ_upS;
    sysError_fa2Targ_downS[i] = error_fa2Targ_downS;
    
    xvals[i] = 1 + i;
  }

  //Sys error by kinematics
  Double_t sysErrorPhys_faPol[nPhysBinned], sysErrorPhys_faSubper[nPhysBinned];
  Double_t sysErrorPhys_upS[nPhysBinned], sysErrorPhys_downS[nPhysBinned];
  Double_t sysErrorPhys_upSdownS[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    sysErrorPhys_faPol[i] = CalSysErrorPeriod(hPulls_faPol[i]);
    sysErrorPhys_faSubper[i] = CalSysErrorPeriod(hPulls_faSubper[i]);
    sysErrorPhys_upS[i] = CalSysErrorPeriod(hPulls_upS[i]);
    sysErrorPhys_downS[i] = CalSysErrorPeriod(hPulls_downS[i]);
    sysErrorPhys_upSdownS[i] = CalSysErrorPeriod(hPulls_upSdownS[i]);
  }
  
  //Draw systematic error
  TGraph *gSys_faPol = new TGraph(nBins, xvals, sysError_faPol);
  gSys_faPol->GetYaxis()->SetRangeUser(0.0, 1.0);
  TGraph *gSys_faSubper = new TGraph(nBins, xvals, sysError_faSubper);
  gSys_faSubper->GetYaxis()->SetRangeUser(0.0, 1.0);
  TGraph *gSys_fa2Targ_upS =
    new TGraph(nBins, xvals, sysError_fa2Targ_upS);
  TGraph *gSys_fa2Targ_downS =
    new TGraph(nBins, xvals, sysError_fa2Targ_downS);
  SetUp(gSys_faPol); SetUp(gSys_fa2Targ_downS); SetUp(gSys_fa2Targ_upS);
  
  TCanvas* cSys = new TCanvas("CorrelatedPullSys");
  gSys_faPol->Draw("AP");
  gSys_faPol->SetTitle("Systematic Error/Statistical Error faPol");
  gSys_faSubper->Draw("Psame"); gSys_faSubper->SetMarkerColor(kBlue+2);
  gSys_fa2Targ_upS->Draw("Psame"); gSys_fa2Targ_upS->SetMarkerColor(kRed);
  gSys_fa2Targ_downS->Draw("Psame"); gSys_fa2Targ_downS->SetMarkerColor(kBlue);

  //Drawing sys error by kinematics
  TCanvas* cSysPhys_faPol = new TCanvas("faPolSys");
  cSysPhys_faPol->Divide(nPhysBinned, 1, 0, 0.01);
  Double_t wMean_faPol =0.0, wMeanSig_faPol =0.0;
  Double_t wSig2_faPol =0.0, wSig2Sig_faPol =0.0;
  TGraph *gSysPhys_faPol[nPhysBinned];
  Double_t wMean_upS =0.0, wMeanSig_upS =0.0;
  Double_t wSig2_upS =0.0, wSig2Sig_upS =0.0;
  TGraph *gSysPhys_upS[nPhysBinned];
  Double_t wMean_downS =0.0, wMeanSig_downS =0.0;
  Double_t wSig2_downS =0.0, wSig2Sig_downS =0.0;
  TGraph *gSysPhys_downS[nPhysBinned];
  Double_t wMean_upSdownS =0.0, wMeanSig_upSdownS =0.0;
  Double_t wSig2_upSdownS =0.0, wSig2Sig_upSdownS =0.0;
  TGraph *gSysPhys_upSdownS[nPhysBinned];
  for (Int_t i=0; i<nPhysBinned; i++) {
    gSysPhys_faPol[i] = new TGraph(1, xvals, &(sysErrorPhys_faPol[i]));
    SetUp(gSysPhys_faPol[i]);
    gSysPhys_faPol[i]->GetYaxis()->SetRangeUser(0.0, 2.5);

    gSysPhys_upS[i] = new TGraph(1, xvals, &(sysErrorPhys_upS[i]));
    gSysPhys_downS[i] = new TGraph(1, xvals, &(sysErrorPhys_downS[i]));
    SetUp(gSysPhys_upS[i]); SetUp(gSysPhys_downS[i]);

    gSysPhys_upSdownS[i] = new TGraph(1, xvals, &(sysErrorPhys_upSdownS[i]));
    SetUp(gSysPhys_upSdownS[i]);

    cSysPhys_faPol->cd(i+1);
    gSysPhys_faPol[i]->Draw("AP");
    gSysPhys_faPol[i]->SetTitle(Form("sysError_%s = %0.2f",
				     physBinned[i].Data(),
				     sysErrorPhys_faPol[i]));

    gSysPhys_upS[i]->Draw("Psame"); gSysPhys_upS[i]->SetMarkerColor(kRed);
    gSysPhys_downS[i]->Draw("Psame"); gSysPhys_downS[i]->SetMarkerColor(kBlue);
    gSysPhys_upSdownS[i]->Draw("Psame");
    gSysPhys_upSdownS[i]->SetMarkerColor(kGreen);

    AddWAvgSigm2(hPulls_faPol[i], wMean_faPol, wMeanSig_faPol,
		 wSig2_faPol, wSig2Sig_faPol);
    AddWAvgSigm2(hPulls_upS[i], wMean_upS, wMeanSig_upS,
		 wSig2_upS, wSig2Sig_upS);
    AddWAvgSigm2(hPulls_downS[i], wMean_downS, wMeanSig_downS,
		 wSig2_downS, wSig2Sig_downS);
    AddWAvgSigm2(hPulls_upSdownS[i], wMean_upSdownS, wMeanSig_upSdownS,
		 wSig2_upSdownS, wSig2Sig_upSdownS);
  }
  FinalWAvgs(wMean_faPol, wMeanSig_faPol, wSig2_faPol, wSig2Sig_faPol);
  FinalWAvgs(wMean_upS, wMeanSig_upS, wSig2_upS, wSig2Sig_upS);
  FinalWAvgs(wMean_downS, wMeanSig_downS, wSig2_downS, wSig2Sig_downS);
  FinalWAvgs(wMean_upSdownS, wMeanSig_upSdownS, wSig2_upSdownS,
	     wSig2Sig_upSdownS); 

  Double_t wSysErrorPhys_faPol = CalSysErrorPeriod(wMean_faPol, wSig2_faPol);
  TGraph *gWSysPhys_faPol = new TGraph(1, xvals, &(wSysErrorPhys_faPol));
  SetUp(gWSysPhys_faPol);

  Double_t wSysErrorPhys_upS = CalSysErrorPeriod(wMean_upS, wSig2_upS);
  TGraph *gWSysPhys_upS = new TGraph(1, xvals, &(wSysErrorPhys_upS));
  SetUp(gWSysPhys_upS);
  
  Double_t wSysErrorPhys_downS = CalSysErrorPeriod(wMean_downS, wSig2_downS);
  TGraph *gWSysPhys_downS = new TGraph(1, xvals, &(wSysErrorPhys_downS));
  SetUp(gWSysPhys_downS);

  Double_t wSysErrorPhys_upSdownS =
    CalSysErrorPeriod(wMean_upSdownS, wSig2_upSdownS);
  TGraph *gWSysPhys_upSdownS =
    new TGraph(1, xvals, &(wSysErrorPhys_upSdownS));
  SetUp(gWSysPhys_upSdownS);

  TCanvas* wSys = new TCanvas("wAvgSys");
  gWSysPhys_faPol->Draw("AP"); gWSysPhys_faPol->SetTitle("Weighted Sys Error");
  gWSysPhys_upS->Draw("Psame"); gWSysPhys_upS->SetMarkerColor(kRed);
  gWSysPhys_downS->Draw("Psame"); gWSysPhys_downS->SetMarkerColor(kBlue);
  gWSysPhys_upSdownS->Draw("Psame"); gWSysPhys_upSdownS->SetMarkerColor(kGreen);

  Double_t yMaxWavg[nBins] = {0.0};
  if (wSysErrorPhys_faPol > yMaxWavg[0])
    yMaxWavg[0] = wSysErrorPhys_faPol;
  if (wSysErrorPhys_upS > yMaxWavg[0])
    yMaxWavg[0] = wSysErrorPhys_upS;
  if (wSysErrorPhys_downS > yMaxWavg[0])
    yMaxWavg[0] = wSysErrorPhys_downS;
  if (wSysErrorPhys_upSdownS > yMaxWavg[0])
    yMaxWavg[0] = wSysErrorPhys_upSdownS;

  for (Int_t i=1; i<nBins; i++) yMaxWavg[i] = yMaxWavg[0];

  TGraph* gSys_Stat = new TGraph(nBins, xvals, yMaxWavg);
  SetUp(gSys_Stat);
  
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Data";
  TString fOutput;
  if  (whichFit=="true"){
    fOutput =
      Form("%s/faPullDist2targWavg/faPullDist2targWavg_true_%s_%s%s_%s_%i_%s_%s.root",
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
    hPull_faPol->Write();

    gSys_faPol->Write("gSys_faPol");
    gSys_faSubper->Write("gSys_faSubper");
    gSys_fa2Targ_upS->Write("gSys_fa2Targ_upS");
    gSys_fa2Targ_downS->Write("gSys_fa2Targ_downS");

    gWSysPhys_faPol->Write("gWSysPhys_faPol");
    gWSysPhys_upS->Write("gWSysPhys_upS");
    gWSysPhys_downS->Write("gWSysPhys_downS");
    gWSysPhys_upSdownS->Write("gWSysPhys_upSdownS");

    gSys_Stat->Write("gSys_Stat");
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
