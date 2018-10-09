#include "include/helperFunctions.h"


Double_t Amp(Double_t L, Double_t R){
  
  Double_t A = L - R;
  A /= ( L + R );

  return A;
}


Double_t e_Amp(Double_t L, Double_t R, Double_t LinvSum, Double_t RinvSum){
  Double_t epsilon = Amp(L, R);
  Double_t dL = 0.25*L*LinvSum;
  Double_t dR = 0.25*R*RinvSum;
  Double_t error =
    (1 - epsilon)*(1-epsilon)*dL*dL + (1 + epsilon)*(1 + epsilon)*dR*dR;
  error = TMath::Sqrt(error);
  error *= 1.0/(L + R);
  
  return error;
}


Double_t AvgPol(Double_t *Pol, Int_t NupSup, Int_t NupSdown,
		Int_t NdownSup, Int_t NdownSdown){
  Double_t avgPol = NupSup*Pol[0] + NupSdown*Pol[1] +
    NdownSup*Pol[2] + NdownSdown*Pol[3];

  avgPol /= (NupSup + NupSdown + NdownSup +NdownSdown);

  return avgPol;
}


Double_t AvgPol(Double_t *Pol, Int_t NupS, Int_t NdownS){
  Double_t avgPol = NupS*(Pol[0]+Pol[1]) + NdownS*(Pol[2]+Pol[3]);
  avgPol *= 0.5;
  avgPol /= (NupS + NdownS);

  return avgPol;
}


Double_t CalLR(Double_t upS_upP, Double_t upS_downP,
	       Double_t downS_upP, Double_t downS_downP){
  Double_t LR = upS_upP*upS_downP*downS_upP*downS_downP;

  return TMath::Power(LR, 0.25);
}


Double_t CalInvSum(Double_t upS_upP, Double_t upS_downP,
		   Double_t downS_upP, Double_t downS_downP){
  Double_t invSum =0.0;
  invSum += 1.0/upS_upP;
  invSum += 1.0/downS_downP;
  invSum += 1.0/upS_downP;
  invSum += 1.0/downS_upP;

  return TMath::Sqrt(invSum);
}


void CalAmp_AmpErr(Double_t *Fasym, Double_t *e_Fasym,
		   Double_t NL[][8], Double_t NR[][8],
		   Double_t Pol[][4], Int_t nBins){
    
  for (Int_t bi=0; bi<nBins; bi++) {
    Long64_t LupS_upP = NL[bi][0] + NL[bi][1];
    Long64_t LdownS_downP = NL[bi][2] + NL[bi][3];
    Long64_t LupS_downP = NL[bi][4] + NL[bi][5];
    Long64_t LdownS_upP = NL[bi][6] + NL[bi][7];
    
    Long64_t RupS_upP = NR[bi][0] + NR[bi][1];
    Long64_t RdownS_downP = NR[bi][2] + NR[bi][3];
    Long64_t RupS_downP = NR[bi][4] + NR[bi][5];
    Long64_t RdownS_upP = NR[bi][6] + NR[bi][7];
    
    Double_t L = CalLR(LupS_upP, LupS_downP, LdownS_upP, LdownS_downP);
    Double_t R = CalLR(RupS_upP, RupS_downP, RdownS_upP, RdownS_downP);

    Double_t LinvSum =
      CalInvSum(LupS_upP, LupS_downP, LdownS_upP, LdownS_downP);
    Double_t RinvSum =
      CalInvSum(RupS_upP, RupS_downP, RdownS_upP, RdownS_downP);
    
    Int_t NupSup = NL[bi][0] + NL[bi][4] + NR[bi][0] + NR[bi][4];
    Int_t NupSdown = NL[bi][1] + NL[bi][5] + NR[bi][1] + NR[bi][5];
    Int_t NdownSup = NL[bi][2] + NL[bi][6] + NR[bi][2] + NR[bi][6];
    Int_t NdownSdown = NL[bi][3] + NL[bi][7] + NR[bi][3] + NR[bi][7];
    Double_t avgPol = AvgPol(Pol[bi], NupSup, NupSdown, NdownSup, NdownSdown);
    
    Fasym[bi] = Amp(L, R)/avgPol;
    e_Fasym[bi] = e_Amp(L, R, LinvSum, RinvSum)/avgPol;
  }
}//CalAmp_AmpErr


void CalAmp_AmpErr4(Double_t *Fasym, Double_t *e_Fasym,
		    Double_t NL[][4], Double_t NR[][4],
		    Double_t Pol[][4], Int_t nBins){
    
  for (Int_t bi=0; bi<nBins; bi++) {
    Double_t L = CalLR(NL[bi][0], NL[bi][2], NL[bi][3], NL[bi][1] );
    Double_t R = CalLR(NR[bi][0], NR[bi][2], NR[bi][3], NR[bi][1] );

    Double_t LinvSum =
      CalInvSum(NL[bi][0], NL[bi][2], NL[bi][3], NL[bi][1] );
    Double_t RinvSum =
      CalInvSum(NR[bi][0], NR[bi][2], NR[bi][3], NR[bi][1] );
    
    Int_t NupS = NL[bi][0] + NR[bi][0] + NL[bi][1] + NR[bi][1];
    Int_t NdownS = NL[bi][2] + NR[bi][2] + NL[bi][3] + NR[bi][3];
    Double_t avgPol = AvgPol(Pol[bi], NupS, NdownS);
    
    Fasym[bi] = Amp(L, R)/avgPol;
    e_Fasym[bi] = e_Amp(L, R, LinvSum, RinvSum)/avgPol;
  }
}//CalAmp_AmpErr


void falseGeoMean4Targ_runNums(TString start =""){
  //Setup_______________
  const Int_t nBins =5;
  TString binRange ="25_43";
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.90_3.30";
  TString fitMrange ="2.00_7.50";
  TString whichFit ="ten";

  Bool_t toWrite =false;
  //Setup_______________

  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/";
  TString pathLR = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/\
fullTargSysFunctMFit/";
  
  if (start==""){//Basic info
    cout << "\nScript calculates false AN asymmetries using the 4 target";
    cout << " geometric mean method" << endl;
    cout << "False asymmetries are determined per subperiod by";
    cout << " splitting each target in half and flipping the spin of the ";
    cout <<"central (inner) two targets and by flipping the upper half of each";
    cout << " target" << endl;
    cout << "\nInput needed is a TGraphErrors of AN per target/pol";
    cout<<" with polarization corrections and without polarization corrections";
    cout << "\nand a TGraphErrors of Left and Right counts per target/pol";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "systematic_leftRight  ->  (for mass fitting only) functMFit.C ";
    cout <<"-> falseGeoMean4Targ_splitTarg.C" << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'falseGeoMean4Targ_splitTarg(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Period and Mass type considered:   " << period_Mtype << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    if (whichFit == "true"){
      if (lrMrange != fitMrange){
	cout << "L/R mass range does not equal fit mass range for true fit\n";
	exit(EXIT_FAILURE);
      }
      cout << "L/R mass range:     " << lrMrange << endl;
    }
    else{
      cout << "LR integral mass range:     " << lrMrange << endl;
      cout << "Fit mass range:     " << fitMrange << endl;
    }
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    cout << "\n\nTo do:" << endl;
    cout << "Setup massing fitting input/output properly" << endl;
    exit(EXIT_FAILURE);
  }

  //File name setup && get file
  TString inputFile, inputLR;
  if (whichFit=="true"){
    inputFile =
      Form("Data/systematic_leftRight_%s%s_%ibins%s_%ihbin_optionER.root",
	   period_Mtype.Data(), lrMrange.Data(), nBins, binRange.Data(), hbins);
  }
  else{
    inputFile =
      Form("Data/systematic_leftRight_%s1.00_8.50_%ibins%s_%ihbin_optionER.root",
	   period_Mtype.Data(), nBins, binRange.Data(), hbins);
    
    inputLR = Form("fullTargSysFunctMFit_%s%s_%s_%s%s_%s%i_%ihbin.root",
		   whichFit.Data(), fitMrange.Data(), period_Mtype.Data(),
		   process.Data(), lrMrange.Data(), physBinned.Data(),
		   nBins, hbins);
  }
  
  TFile *f_sys =TFile::Open(path+inputFile), *f_LR =NULL;; 
  if (whichFit != "true") {
    f_LR = TFile::Open(pathLR+inputLR);
    if (!f_LR){
      cout << "Data file does not exist " << endl;
      exit(EXIT_FAILURE);
    } 
  }
  if (!f_sys){
    cout << "Data file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Get polarization values
  const Int_t nTarg=4;
  TString targName[nTarg] = {"upSup", "upSdown", "downSup", "downSdown"};
  Double_t Pol[nBins][nTarg];
  for (Int_t tr=0; tr<nTarg; tr++) {
    TString Polname =
      Form("%s_avgPolDil_%s",physBinned.Data(),targName[tr].Data());
    TGraph *g_Pol = (TGraph*)f_sys->Get(Polname);

    if (!g_Pol) {
      cout << Polname << endl;
      cout << "g_Pol doesn't exist" << endl;
      exit(EXIT_FAILURE);
    }
    
    Double_t *y_Pol = g_Pol->GetY();
    
    for (Int_t bi=0; bi<nBins; bi++) { Pol[bi][tr] = y_Pol[bi]; }
  }

  //Make asymmetries differently for each fit type
  Double_t FA[nBins], eFA[nBins];
  Double_t xvals[nBins], ex[nBins] ={0.0};
  if (whichFit == "true"){
    //File data names setup
    const Int_t nTargPol =8;    
    TString targPolName[nTargPol] = {"upSup_upP", "upSdown_upP",       //Sub one
				     "downSup_downP", "downSdown_downP",
				     "upSup_downP", "upSdown_downP",   //Sub two
				     "downSup_upP", "downSdown_upP"};
    //Get L/R counts 
    Double_t LeftCounts[nBins][nTargPol], RightCounts[nBins][nTargPol];

    for (Int_t tr=0; tr<nTargPol; tr++) {
      TString inputLeft =
	Form("%s_left_%s", physBinned.Data(), targPolName[tr].Data());
      TString inputRight =
	Form("%s_right_%s", physBinned.Data(), targPolName[tr].Data());

      TGraph *g_Left, *g_Right;
      g_Left = (TGraph*)f_sys->Get(inputLeft);
      g_Right = (TGraph*)f_sys->Get(inputRight);

      if (!g_Left || !g_Right){
	cout << inputLeft << " " << inputRight << endl;
	cout << "Graphs don't exist" << endl;
	exit(EXIT_FAILURE);
      }
            
      Double_t *y_Left = g_Left->GetY();
      Double_t *y_Right = g_Right->GetY();
      Double_t *x_Right = g_Right->GetX();
      for (Int_t bi=0; bi<nBins; bi++) {
	LeftCounts[bi][tr] = y_Left[bi];
	RightCounts[bi][tr] = y_Right[bi];

	if (tr==0) xvals[bi] = x_Right[bi];
      }
    }

    //Make 4Targ False Asymmetrie
    CalAmp_AmpErr(FA, eFA, LeftCounts, RightCounts, Pol, nBins);
  }
  else {//Not true fit
    //File data names setup
    const Int_t nTargPol =4;    
    TString targPolName[nTargPol] = {"upstream_up", "upstream_down",
				     "downstream_up", "downstream_down"};
    
    //Get L/R counts 
    Double_t LeftCounts[nBins][nTargPol], RightCounts[nBins][nTargPol];
    for (Int_t tr=0; tr<nTargPol; tr++) {
      TString inputLeft =
	Form("%s_left_%s", physBinned.Data(), targPolName[tr].Data());
      TString inputRight =
	Form("%s_right_%s", physBinned.Data(), targPolName[tr].Data());

      TGraph *g_Left, *g_Right;
      g_Left = (TGraph*)f_LR->Get(inputLeft);
      g_Right = (TGraph*)f_LR->Get(inputRight);

      if (!g_Left || !g_Right){
	cout << inputLeft << " " << inputRight << endl;
	cout << "Graphs don't exist" << endl;
	exit(EXIT_FAILURE);
      }
            
      Double_t *y_Left = g_Left->GetY();
      Double_t *y_Right = g_Right->GetY();
      Double_t *x_Right = g_Right->GetX();
      for (Int_t bi=0; bi<nBins; bi++) {
	LeftCounts[bi][tr] = y_Left[bi];
	RightCounts[bi][tr] = y_Right[bi];

	if (tr==0) xvals[bi] = x_Right[bi];
      }
    }

    //Make 4Targ False Asymmetrie
    CalAmp_AmpErr4(FA, eFA, LeftCounts, RightCounts, Pol, nBins);
  }
  
  //Draw false AN
  TGraphErrors *g_FA = new TGraphErrors(nBins, xvals, FA, ex, eFA);
  SetUp(g_FA);
  
  TCanvas* c1 = new TCanvas();
  Double_t yMax = 0.3;
  g_FA->Draw("AP"); g_FA->SetMarkerColor(kGreen);
  DrawLine(g_FA, 0.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/RunNums";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/falseGeoMeanRunNums_true_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), period_Mtype.Data(), process.Data(),
		  lrMrange.Data(), physBinned.Data(), nBins);
  }
  else{
    fOutput =Form("%s/falseGeoMeanRunNums_%s%s_%s_%s%s_%s%i_%ihbins.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_FA->Write("FA_RunNums");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << path << endl;
  cout << "Input data:        " << inputFile << endl;
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
