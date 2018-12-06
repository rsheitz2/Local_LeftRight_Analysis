#include "include/helperFunctions.h"

Double_t Amp(Double_t NL[][4], Double_t NR[][4],
	     Double_t Pol, Int_t bi){
  //NL[][4] => 4 = nTargPol
  
  Double_t Lup, Rup;
  Lup = NL[bi][upS_up]*NL[bi][upS_down];
  Rup = NR[bi][upS_up]*NR[bi][upS_down];

  Double_t Ldown, Rdown;
  Ldown = NL[bi][downS_up]*NL[bi][downS_down];
  Rdown = NR[bi][downS_up]*NR[bi][downS_down];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);
  
  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;
  
  return A;
}


Double_t e_Amp(Double_t NL[][4], Double_t NR[][4],
	       Double_t e_NL[][4], Double_t e_NR[][4],
	       Double_t Pol, Int_t bi){
  //NL[][4] => 4 = nTargPol
  Int_t nTargPol=4;
  
  Double_t L=1.0, R=1.0;
  Double_t LinvSum2=0.0, RinvSum2=0.0;
  for (Int_t tr=0; tr<nTargPol; tr++) {
    L *= NL[bi][tr]; R *= NR[bi][tr];
    LinvSum2 += e_NL[bi][tr]*e_NL[bi][tr]/(NL[bi][tr]*NL[bi][tr] );
    RinvSum2 += e_NR[bi][tr]*e_NR[bi][tr]/(NR[bi][tr]*NR[bi][tr] );
  }
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  Double_t e = L*R/( 2*(L + R)*(L + R) );
  Double_t error = e*TMath::Sqrt(LinvSum2 + RinvSum2)/Pol;

  return error;
}


void GeoMean4Targ(TString start =""){
  //Setup_______________
  const Int_t nBins =1;//HMDY
  TString period_Mtype ="Charles_HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.53";

  //const Int_t nBins =5;//JPsi
  //TString period_Mtype ="W12_LowM_AMDY";
  //Int_t hbins =150;
  //TString physBinned ="xPi";//xN, xPi, xF, pT, M
  //TString process ="JPsi";//JPsi, psi, DY
  //TString lrMrange ="2.00_5.00";
  //TString fitMrange ="2.00_7.50";
  //TString binRange ="25_43";
  //TString whichFit ="thirteen";
  //TString production ="slot1";//"t3", "slot1"
  //TString additionalCuts ="phiS0.195";

  Bool_t toWrite =false;
  //Setup_______________

    TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/";
  
  if (start==""){//Basic info
    cout << "\nScript calculates AN using the 4 target geometric mean method";
    cout << "\nInput needed is a TGraphErrors of AN per target/pol";
    cout<<" with polarization corrections and without polarization corrections";
    cout << "\nand a TGraphErrors of Left and Right counts per target/pol";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C or trueCount.C ";
    cout <<"-> GeoMean4Targ.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'AN_4TargGeoMean_fromMFit(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Period and Mass type considered:   " << period_Mtype << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    if (whichFit == "true"){
      cout << "L/R mass range:     " << fitMrange << endl;
    }
    else{
      cout << "LR integral mass range:     " << lrMrange << endl;
      cout << "Fit mass range:     " << fitMrange << endl;
    }
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //File name setups && get files
  TString inputFiles, RDfile;
  if (whichFit=="true"){
    if (fitMrange != lrMrange){
      cout << "fit mass range != left/right mass range for true fit" << endl;
      exit(EXIT_FAILURE);
    }
    
    inputFiles =
      Form("trueCount_%s_%s%s_%s%s%i_%s_%s_corr.root", period_Mtype.Data(),
	   process.Data(), fitMrange.Data(), binRange.Data(), physBinned.Data(),
	   nBins, production.Data(), additionalCuts.Data());
    pathAN += "trueCount/";

    RDfile =Form("leftRight_byTarget_%s%s_%ibins%s_%ihbin_%s_%s.root",
		 period_Mtype.Data(), lrMrange.Data(), nBins, binRange.Data(),
		 hbins, production.Data(), additionalCuts.Data());
  }
  else if (whichFit=="MC"){
    cout << "GeoMean currently doesn't work for MC fit" << endl;
    exit(EXIT_FAILURE);
    inputFiles = Form("mcMFit_%s%s_%s_%s%s_%s%i_%ihbin_corr.root",
			 whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
			 process.Data(), lrMrange.Data(), physBinned.Data(),
			 nBins, hbins);
    pathAN += "mcMFit/";

    RDfile =
      Form("leftRight_byTarget_%s%s_%ibins%s_150hbin_%s_%s.root",
	   period_Mtype.Data(), fitMrange.Data(), nBins, binRange.Data(),
	   production.Data(), additionalCuts.Data() );
  }
  else{
    inputFiles = Form("functMFit_%s%s_%s_%s%s_%s%s%i_%ihbin_%s_%s_corr.root",
		      whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
		      process.Data(), lrMrange.Data(), binRange.Data(),
		      physBinned.Data(), nBins, hbins, production.Data(),
		      additionalCuts.Data());
    pathAN += "functMFit/";

    RDfile =
      Form("leftRight_byTarget_%s1.00_8.50_%ibins%s_%ihbin_%s_%s.root",
	   period_Mtype.Data(), nBins, binRange.Data(), hbins,
	   production.Data(), additionalCuts.Data() );
}
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  
  TFile *fAN = TFile::Open(pathAN+inputFiles);
  TFile *fRD  = TFile::Open(pathRD + RDfile);
  if (!fAN || !fRD ){
    cout << "fAN or fRD file does not exist " << endl;
    cout << "fAN = \n"; cout << pathAN+inputFiles << endl;
    cout << "fRD = \n"; cout << pathRD+RDfile << endl;
    exit(EXIT_FAILURE);
  }
  
  //Get data from files
  TString inputAN[4];//{upS_up, upS_down, downS_up, downS_down}
  TString inputLeft[4];//{upS_up, upS_down, downS_up, downS_down}
  TString inputRight[4];//{upS_up, upS_down, downS_up, downS_down}
  
  TString targName[4] = {"upS_up", "upS_down", "downS_up", "downS_down"};  
  for (Int_t i=0; i<4; i++){
    inputAN[i] = Form("AN_%s", targName[i].Data() );
    inputLeft[i] = Form("Left_%s", targName[i].Data() );
    inputRight[i] = Form("Right_%s", targName[i].Data() );
  }

  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  TGraph* g_Dil =(TGraph*)fRD->Get(Form("%s_Dil", physBinned.Data()));

  const Int_t nTargPol =4;
  TGraphErrors *g_Left[nTargPol], *g_Right[nTargPol];
  for (Int_t tr=0; tr<nTargPol; tr++) {
    g_Left[tr] = (TGraphErrors*)fAN->Get(inputLeft[tr]);
    g_Right[tr] = (TGraphErrors*)fAN->Get(inputRight[tr]);
  }
  
  //Get L/R counts
  Double_t LeftCounts[nBins][nTargPol], RightCounts[nBins][nTargPol];
  Double_t e_LeftCounts[nBins][nTargPol], e_RightCounts[nBins][nTargPol];
  
  for (Int_t tr=0; tr<nTargPol; tr++) {
    Double_t *y_Left = g_Left[tr]->GetY();
    Double_t *y_Right = g_Right[tr]->GetY();
    Double_t *ey_Left = g_Left[tr]->GetEY();
    Double_t *ey_Right = g_Right[tr]->GetEY();
    for (Int_t bi=0; bi<nBins; bi++) {
      LeftCounts[bi][tr] = y_Left[bi];
      RightCounts[bi][tr] = y_Right[bi];
      e_LeftCounts[bi][tr] = ey_Left[bi];
      e_RightCounts[bi][tr] = ey_Right[bi];
    }
  }

  //Get Polarization
  Double_t Pol[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);
  
  //Make 4Targ Asym
  Double_t AN_4targ[nBins], e_AN_4targ[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_Pol->GetX();
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_4targ[bi] = Amp(LeftCounts, RightCounts, Pol[bi], bi);
    e_AN_4targ[bi] = e_Amp(LeftCounts, RightCounts,
			   e_LeftCounts, e_RightCounts, Pol[bi], bi);
    if (e_AN_4targ[bi] < 10e-9) {
      cout << "Error: AN error way too small" << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Draw AN
  TGraphErrors *g_AN = new TGraphErrors(nBins, xvals, AN_4targ, ex, e_AN_4targ);
  SetUpTGraph(g_AN);

  Double_t yMax =0.15;
  TCanvas* c1 = new TCanvas();
  g_AN->Draw("AP"); g_AN->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN, 0.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =
      Form("%s/GeoMean4Targ_true_%s_%s%s_%s%s%i_%s_%s.root",
	   thisDirPath.Data(), period_Mtype.Data(), process.Data(),
	   fitMrange.Data(), binRange.Data(), physBinned.Data(), nBins, 
	   production.Data(), additionalCuts.Data());
  }
  else{
    fOutput =
      Form("%s/GeoMean4Targ_%s%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(),binRange.Data(),
	   physBinned.Data(), nBins, hbins, production.Data(),
	   additionalCuts.Data());
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_AN->Write("AN");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathAN << endl;
  cout << "Input P corrected data:        " << inputFiles << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Period and Mass type considered:   " << period_Mtype << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
  cout << "AN physical process:        " << process << endl;
  if (whichFit == "true"){
      cout << "L/R mass range:     " << lrMrange << endl;
    }
    else{
      cout << "LR integral mass range:     " << lrMrange << endl;
      cout << "Fit mass range:     " << fitMrange << endl;
    }
  cout << "Which fit considered:       " << whichFit << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
