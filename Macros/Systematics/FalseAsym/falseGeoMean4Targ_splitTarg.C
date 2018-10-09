#include "include/helperFunctions.h"


Double_t Amp(Double_t NL_t1, Double_t NR_t1,
	     Double_t NL_t2, Double_t NR_t2,
	     Double_t NL_t3, Double_t NR_t3,
	     Double_t NL_t4, Double_t NR_t4,
	     Double_t Pol){
  
  Double_t L = NL_t1*NL_t2*NL_t3*NL_t4;
  Double_t R = NR_t1*NR_t2*NR_t3*NR_t4;
  
  //L = TMath::Sqrt(L);
  //R = TMath::Sqrt(R);
  
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  Double_t A = L - R;
  A /= ( L + R );

  return A/Pol;
}


Double_t e_Amp(Double_t NL_t1, Double_t NR_t1,
	       Double_t NL_t2, Double_t NR_t2,
	       Double_t NL_t3, Double_t NR_t3,
	       Double_t NL_t4, Double_t NR_t4,
	       Double_t e_NL_t1, Double_t e_NR_t1,
	       Double_t e_NL_t2, Double_t e_NR_t2,
	       Double_t e_NL_t3, Double_t e_NR_t3,
	       Double_t e_NL_t4, Double_t e_NR_t4,
	       Double_t Pol){
  
  Double_t LinvSum2 =
    e_NL_t1*e_NL_t1/(NL_t1*NL_t1) +
    e_NL_t2*e_NL_t2/(NL_t2*NL_t2) +
    e_NL_t3*e_NL_t3/(NL_t3*NL_t3) +
    e_NL_t4*e_NL_t4/(NL_t4*NL_t4);
  LinvSum2 = TMath::Sqrt(LinvSum2);
  Double_t RinvSum2 =
    e_NR_t1*e_NR_t1/(NR_t1*NR_t1) +
    e_NR_t2*e_NR_t2/(NR_t2*NR_t2) +
    e_NR_t3*e_NR_t3/(NR_t3*NR_t3) +
    e_NR_t4*e_NR_t4/(NR_t4*NR_t4);
  RinvSum2 = TMath::Sqrt(RinvSum2);
  
  Double_t L =NL_t1*NL_t2*NL_t3*NL_t4;
  Double_t R = NR_t1*NR_t2*NR_t3*NR_t4;

  //L = TMath::Sqrt(L);
  //R = TMath::Sqrt(R);
  //Double_t dL = 0.5*L*LinvSum2;
  //Double_t dR = 0.5*R*RinvSum2;

  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);
  Double_t dL = 0.25*L*LinvSum2;
  Double_t dR = 0.25*R*RinvSum2;

  Double_t epsilon = (L - R)/(L + R);
  Double_t error =
    (1 - epsilon)*(1-epsilon)*dL*dL + (1 + epsilon)*(1 + epsilon)*dR*dR;
  error = TMath::Sqrt(error);
  error *= 1.0/(L + R);
  
  return error/Pol;
}


Double_t AvgPol(Double_t *Pol, Int_t NupSup, Int_t NupSdown,
		Int_t NdownSup, Int_t NdownSdown){
  Double_t avgPol = NupSup*Pol[0] + NupSdown*Pol[1] +
    NdownSup*Pol[2] + NdownSdown*Pol[3];

  avgPol /= (NupSup + NupSdown + NdownSup +NdownSdown);

  return avgPol;
}


void CalAmp_AmpErr(Double_t *Fasym, Double_t *e_Fasym,
		   Double_t NL[][8], Double_t NR[][8],
		   Double_t e_NL[][8], Double_t e_NR[][8],
		   Double_t Pol[][4], Int_t nBins, TString whichFasym){
    
  for (Int_t bi=0; bi<nBins; bi++) {
    Int_t NupSup = NL[bi][0] + NL[bi][4] + NR[bi][0] + NR[bi][4];
    Int_t NupSdown = NL[bi][1] + NL[bi][5] + NR[bi][1] + NR[bi][5];
    Int_t NdownSup = NL[bi][2] + NL[bi][6] + NR[bi][2] + NR[bi][6];
    Int_t NdownSdown = NL[bi][3] + NL[bi][7] + NR[bi][3] + NR[bi][7];
    
    Double_t avgPol = AvgPol(Pol[bi], NupSup, NupSdown, NdownSup, NdownSdown);
    
    if (whichFasym=="sb1_center"){
      Fasym[bi] = Amp(NL[bi][0], NR[bi][0], //upSup_upP
		      NR[bi][1], NL[bi][1], //upSdown_upP     //flipped
		      NR[bi][2], NL[bi][2], //downSup_downP   //flipped
		      NL[bi][3], NR[bi][3], //downSdown_downP
		      avgPol);

      e_Fasym[bi] = e_Amp(NL[bi][0], NR[bi][0], //upSup_upP
			  NR[bi][1], NL[bi][1], //upSdown_upP     //flipped
			  NR[bi][2], NL[bi][2], //downSup_downP   //flipped
			  NL[bi][3], NR[bi][3], //downSdown_downP
			  e_NL[bi][0], e_NR[bi][0], //upSup_upP
			  e_NR[bi][1], e_NL[bi][1], //upSdown_upP     //flipped
			  e_NR[bi][2], e_NL[bi][2], //downSup_downP   //flipped
			  e_NL[bi][3], e_NR[bi][3], //downSdown_downP
			  avgPol);
    }
    else if (whichFasym=="sb2_center"){
      Fasym[bi] = Amp(NL[bi][4], NR[bi][4], //upSup_downP
		      NR[bi][5], NL[bi][5], //upSdown_downP     //flipped
		      NR[bi][6], NL[bi][6], //downSup_upP       //flipped
		      NL[bi][7], NR[bi][7], //downSdown_upP
		      avgPol);

      e_Fasym[bi] = e_Amp(NL[bi][4], NR[bi][4], //upSup_downP
			  NR[bi][5], NL[bi][5], //upSdown_downP     //flipped
			  NR[bi][6], NL[bi][6], //downSup_upP       //flipped
			  NL[bi][7], NR[bi][7], //downSdown_upP
			  e_NL[bi][4], e_NR[bi][4], //upSup_downP
			  e_NR[bi][5], e_NL[bi][5], //upSdown_downP    //flipped
			  e_NR[bi][6], e_NL[bi][6], //downSup_upP      //flipped
			  e_NL[bi][7], e_NR[bi][7], //downSdown_upP
			  avgPol);
      
    }
    else if (whichFasym=="sb1_Sup"){
      Fasym[bi] = Amp(NR[bi][0], NL[bi][0], //upSup_upP       //flipped
		      NL[bi][1], NR[bi][1], //upSdown_upP   
		      NR[bi][2], NL[bi][2], //downSup_downP   //flipped
		      NL[bi][3], NR[bi][3], //downSdown_downP
		      avgPol);

      e_Fasym[bi] = e_Amp(NR[bi][0], NL[bi][0], //upSup_upP       //flipped
			  NL[bi][1], NR[bi][1], //upSdown_upP   
			  NR[bi][2], NL[bi][2], //downSup_downP   //flipped
			  NL[bi][3], NR[bi][3], //downSdown_downP
			  e_NR[bi][0], e_NL[bi][0], //upSup_upP       //flipped
			  e_NL[bi][1], e_NR[bi][1], //upSdown_upP     
			  e_NR[bi][2], e_NL[bi][2], //downSup_downP   //flipped
			  e_NL[bi][3], e_NR[bi][3], //downSdown_downP
			  avgPol);
      
    }
    else if (whichFasym=="sb2_Sup"){
      Fasym[bi] = Amp(NR[bi][4], NL[bi][4], //upSup_downP       //flipped
		      NL[bi][5], NR[bi][5], //upSdown_downP     
		      NR[bi][6], NL[bi][6], //downSup_upP       //flipped
		      NL[bi][7], NR[bi][7], //downSdown_upP
		      avgPol);

      e_Fasym[bi] = e_Amp(NR[bi][4], NL[bi][4], //upSup_downP       //flipped
			  NL[bi][5], NR[bi][5], //upSdown_downP     
			  NR[bi][6], NL[bi][6], //downSup_upP       //flipped
			  NL[bi][7], NR[bi][7], //downSdown_upP
			  e_NR[bi][4], e_NL[bi][4], //upSup_downP      //flipped
			  e_NL[bi][5], e_NR[bi][5], //upSdown_downP     
			  e_NR[bi][6], e_NL[bi][6], //downSup_upP      //flipped
			  e_NL[bi][7], e_NR[bi][7], //downSdown_upP
			  avgPol);
    }
  }
}//CalAmp_AmpErr


void falseGeoMean4Targ_splitTarg(TString start =""){
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
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/sysFunctMFit/";
  
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
    inputFile = Form("Data/systematic_leftRight_%s%s_%ibins%s_%ihbin.root",
		     period_Mtype.Data(), lrMrange.Data(), nBins,
		     binRange.Data(), hbins);
  }
  else{
    inputFile = Form("Data/systematic_leftRight_%s1.00_8.50_%ibins%s_%ihbin.root",
		     period_Mtype.Data(), nBins, binRange.Data(), hbins);
    
    inputLR = Form("sysFunctMFit_%s%s_%s_%s%s_%s%i_%ihbin.root",
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
    
  //File data names setup
  const Int_t nTargPol =8;  
  TString targPolName[nTargPol] = {"upSup_upP", "upSdown_upP",       //Sub one
				   "downSup_downP", "downSdown_downP",
				   "upSup_downP", "upSdown_downP",   //Sub two
				   "downSup_upP", "downSdown_upP"};

  //Get L/R counts 
  Double_t LeftCounts[nBins][nTargPol], RightCounts[nBins][nTargPol];
  Double_t e_LeftCounts[nBins][nTargPol], e_RightCounts[nBins][nTargPol];
  Double_t xvals[nBins];
  for (Int_t tr=0; tr<nTargPol; tr++) {
    TString inputLeft =
      Form("%s_left_%s",physBinned.Data(),targPolName[tr].Data());
    TString inputRight =
      Form("%s_right_%s",physBinned.Data(),targPolName[tr].Data());

    TGraph *g_Left, *g_Right;
    if (whichFit == "true"){
      g_Left = (TGraph*)f_sys->Get(inputLeft);
      g_Right = (TGraph*)f_sys->Get(inputRight);
    }
    else{
      g_Left = (TGraph*)f_LR->Get(inputLeft);
      g_Right = (TGraph*)f_LR->Get(inputRight);
    }

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
      e_LeftCounts[bi][tr] = TMath::Sqrt( y_Left[bi] );
      e_RightCounts[bi][tr] = TMath::Sqrt( y_Right[bi] );

      if (tr==0) xvals[bi] = x_Right[bi];
    }
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
    
    for (Int_t bi=0; bi<nBins; bi++) {
      Pol[bi][tr] = y_Pol[bi];
    }
  }

  //Make 4Targ False Asymmetries
  //sb1 = subperiod One: upstream Pup, downstream Pdown
  //sb2 = subperiod Two: upstream Pdown, downstream Pup
  Double_t FA_sb1_center[nBins], eFA_sb1_center[nBins];
  Double_t FA_sb2_center[nBins], eFA_sb2_center[nBins];
  Double_t FA_sb1_Sup[nBins], eFA_sb1_Sup[nBins];
  Double_t FA_sb2_Sup[nBins], eFA_sb2_Sup[nBins];
  Double_t ex[nBins] = {0.0};

  CalAmp_AmpErr(FA_sb1_center, eFA_sb1_center, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, Pol, nBins, "sb1_center");
  CalAmp_AmpErr(FA_sb2_center, eFA_sb2_center, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, Pol, nBins, "sb2_center");

  CalAmp_AmpErr(FA_sb1_Sup, eFA_sb1_Sup, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, Pol, nBins, "sb1_Sup");
  CalAmp_AmpErr(FA_sb2_Sup, eFA_sb2_Sup, LeftCounts, RightCounts,
		e_LeftCounts, e_RightCounts, Pol, nBins, "sb2_Sup");
    
  //Draw false AN
  TGraphErrors *g_sb1_center =
    new TGraphErrors(nBins, xvals, FA_sb1_center, ex, eFA_sb1_center);
  TGraphErrors *g_sb2_center =
    new TGraphErrors(nBins, xvals, FA_sb2_center, ex, eFA_sb2_center);

  TGraphErrors *g_sb1_Sup =
    new TGraphErrors(nBins, xvals, FA_sb1_Sup, ex, eFA_sb1_Sup);
  TGraphErrors *g_sb2_Sup =
    new TGraphErrors(nBins, xvals, FA_sb2_Sup, ex, eFA_sb2_Sup);
  
  TCanvas* c1 = new TCanvas();
  Double_t yMax = 0.3;
  //Double_t offset =0.01;
  Double_t offset =0.0;
  SetUp(g_sb1_center); SetUp(g_sb2_center, offset);
  SetUp(g_sb1_Sup, offset*2); SetUp(g_sb2_Sup, offset*3);
  g_sb1_center->Draw("AP"); g_sb1_center->SetMarkerColor(36);
  g_sb2_center->Draw("Psame"); g_sb2_center->SetMarkerColor(kRed);
  g_sb1_Sup->Draw("Psame"); g_sb1_Sup->SetMarkerColor(kBlue);
  g_sb2_Sup->Draw("Psame"); g_sb2_Sup->SetMarkerColor(38);
  g_sb1_center->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_sb1_center, 0.0);
  
  //Write output/final settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/SplitTarg";
  TString fOutput;
  if (whichFit=="true"){
    fOutput =Form("%s/falseGeoMeanSplitTarg_true_%s_%s%s_%s%i.root",
		  thisDirPath.Data(), period_Mtype.Data(), process.Data(),
		  lrMrange.Data(), physBinned.Data(), nBins);
  }
  else{
    fOutput =Form("%s/falseGeoMeanSplitTarg_%s%s_%s_%s%s_%s%i_%ihbins.root",
		  thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		  period_Mtype.Data(), process.Data(), lrMrange.Data(),
		  physBinned.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_sb1_center->Write("sb1_center");
    g_sb2_center->Write("sb2_center");
    g_sb1_Sup->Write("sb1_Sup");
    g_sb2_Sup->Write("sb2_Sup");
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
