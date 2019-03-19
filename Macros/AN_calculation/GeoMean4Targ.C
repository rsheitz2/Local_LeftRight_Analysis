#include "include/helperFunctions.h"

Double_t Amp(Double_t NL[][4], Double_t NR[][4],
	     Double_t Pol, Int_t bi);

Double_t OneTargAmp(Double_t *NL, Double_t *NR, Double_t Pol);

Double_t e_Amp(Double_t NL[][4], Double_t NR[][4],
	       Double_t e_NL[][4], Double_t e_NR[][4],
	       Double_t Pol, Int_t bi);

Double_t e_OneTargAmp(Double_t *NL, Double_t *NR,
		      Double_t *e_NL, Double_t *e_NR,
		      Double_t Pol);

void GeoMean4Targ(TString start =""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";

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
  TGraph* g_Pol_ups =
    (TGraph*)fRD->Get(Form("%s_Pol_upstream", physBinned.Data()));
  TGraph* g_Dil_ups =
    (TGraph*)fRD->Get(Form("%s_Dil_upstream", physBinned.Data()));
  TGraph* g_Pol_downs =
    (TGraph*)fRD->Get(Form("%s_Pol_downstream", physBinned.Data()));
  TGraph* g_Dil_downs =
    (TGraph*)fRD->Get(Form("%s_Dil_downstream", physBinned.Data()));

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
  Double_t Pol[nBins], Pol_ups[nBins], Pol_downs[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);
  GetPolarization(g_Pol_ups, g_Dil_ups, Pol_ups);
  GetPolarization(g_Pol_downs, g_Dil_downs, Pol_downs);
  
  //Make 4Targ Asym
  Double_t AN_4targ[nBins], e_AN_4targ[nBins];
  Double_t AN_4targUncorr[nBins], e_AN_4targUncorr[nBins], PolUncorr[nBins];
  Double_t AN_upS[nBins], e_AN_upS[nBins], AN_downS[nBins], e_AN_downS[nBins];
  Double_t ex[nBins] = {0.0};
  Double_t *xvals = g_Pol->GetX();
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_4targ[bi] = Amp(LeftCounts, RightCounts, Pol[bi], bi);
    e_AN_4targ[bi] = e_Amp(LeftCounts, RightCounts,
			   e_LeftCounts, e_RightCounts, Pol[bi], bi);
    
    PolUncorr[bi] = 1.0;
    AN_4targUncorr[bi] = Amp(LeftCounts, RightCounts, PolUncorr[bi], bi);
    e_AN_4targUncorr[bi] = e_Amp(LeftCounts, RightCounts, e_LeftCounts,
				 e_RightCounts, PolUncorr[bi], bi);

    AN_upS[bi] = OneTargAmp(&(LeftCounts[bi][0]), &(RightCounts[bi][0]),
			    Pol_ups[bi]);
    e_AN_upS[bi] = e_OneTargAmp(&(LeftCounts[bi][0]), &(RightCounts[bi][0]),
				&(e_LeftCounts[bi][0]), &(e_RightCounts[bi][0]),
				Pol_ups[bi]);

    AN_downS[bi] = OneTargAmp(&(LeftCounts[bi][2]), &(RightCounts[bi][2]),
			      Pol_downs[bi]);
    e_AN_downS[bi] =
      e_OneTargAmp(&(LeftCounts[bi][2]), &(RightCounts[bi][2]), 
		   &(e_LeftCounts[bi][2]), &(e_RightCounts[bi][2]),
		   Pol_downs[bi]);
  }

  //Draw AN
  TGraphErrors *g_AN = new TGraphErrors(nBins, xvals, AN_4targ, ex, e_AN_4targ);
  TGraphErrors *g_ANUncorr =
    new TGraphErrors(nBins, xvals, AN_4targUncorr, ex, e_AN_4targUncorr);
  TGraphErrors *g_AN_ups =
    new TGraphErrors(nBins, xvals, AN_upS, ex, e_AN_upS);
  TGraphErrors *g_AN_downs =
    new TGraphErrors(nBins, xvals, AN_downS, ex, e_AN_downS);
  SetUpTGraph(g_AN); SetUpTGraph(g_ANUncorr);
  SetUpTGraph(g_AN_ups); SetUpTGraph(g_AN_downs);

  Double_t yMax =0.5;
  TCanvas* c1 = new TCanvas();
  g_AN->Draw("AP"); g_AN->GetYaxis()->SetRangeUser(-yMax, yMax);
  //g_ANUncorr->Draw("Psame"); g_ANUncorr->SetMarkerColor(kBlue);
  //g_AN_ups->Draw("Psame"); g_AN_ups->SetMarkerColor(kBlue);
  //g_AN_downs->Draw("Psame"); g_AN_downs->SetMarkerColor(kRed);
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
    g_ANUncorr->Write("ANuncorr");
    g_AN_ups->Write("AN_ups");
    g_AN_downs->Write("AN_downs");
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


Double_t OneTargAmp(Double_t *NL, Double_t *NR, Double_t Pol){
  //Regular geomean asymmetry
  
  Double_t L = NL[0]*NL[1];
  Double_t R = NR[0]*NR[1];

  L = TMath::Sqrt(L);
  R = TMath::Sqrt(R);
  
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

  if (error < 10e-9) {
      cout << "Error: AN error way too small" << endl;
      exit(EXIT_FAILURE);
  }

  return error;
}


Double_t e_OneTargAmp(Double_t *NL, Double_t *NR,
		      Double_t *e_NL, Double_t *e_NR,
		      Double_t Pol){
  //Regular geomean asymmetry
  
  Double_t L, R;
  L = NL[0]*NL[1];
  R = NR[0]*NR[1];

  L = TMath::Sqrt(L);
  R = TMath::Sqrt(R);
  
  Double_t A = L - R;
  A /= ( L + R );

  Double_t dL2 = 0.5*L*0.5*L*( e_NL[0]*e_NL[0]/(NL[0]*NL[0]) +
			       e_NL[1]*e_NL[1]/(NL[1]*NL[1]) );
  Double_t dR2 = 0.5*R*0.5*R*( e_NR[0]*e_NR[0]/(NR[0]*NR[0]) +
			       e_NR[1]*e_NR[1]/(NR[1]*NR[1]) );

  Double_t error = ( 1-A )*( 1-A )*dL2 + ( 1+A )*( 1+A )*dR2;
  error = TMath::Sqrt( error );
  error *= 1.0/( L + R );
  error /= Pol;

  if (error < 10e-9) {
    cout << "Error: AN error way too small" << endl;
    exit(EXIT_FAILURE);
  }

  return error;
}
