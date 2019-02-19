#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void upSdownSCompare(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//# of physBinned bins
  TString Mtype ="HMDY";
  TString physBinned ="xN";//"xN", "xPi", "xF", "pT", "M"
  TString production ="slot1";//"t3", "slot1"
  TString minimizer ="MIGRAD";//MIGRAD, HESSE, MINOS
  
  Bool_t toWrite =false;
  //Setup_______________
  
  //Basic Setup
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting/Data";
  const Int_t nPer =9;
  TString period[nPer] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};

  const Int_t nGraph = 14;
  TString upSNames[nGraph] =
    {"g_A_upS_phys", "g_L_upS_Pup_phys", "g_L_upS_Pdown_phys",
     "g_a1_upS_phys", "g_a2_upS_phys", "g_a3_upS_phys",
     "g_A_upS_phys_expected",
     "g_A_upS_int", "g_L_upS_Pup_int", "g_L_upS_Pdown_int",
     "g_a1_upS_int", "g_a2_upS_int", "g_a3_upS_int",
     "g_A_upS_int_expected"};
  TString downSNames[nGraph] =
    {"g_A_downS_phys", "g_L_downS_Pup_phys", "g_L_downS_Pdown_phys",
     "g_a1_downS_phys", "g_a2_downS_phys", "g_a3_downS_phys",
     "g_A_downS_phys_expected",
     "g_A_downS_int", "g_L_downS_Pup_int", "g_L_downS_Pdown_int",
     "g_a1_downS_int", "g_a2_downS_int", "g_a3_downS_int",
     "g_A_downS_int_expected"};

  Double_t c1Pup_c2Pdown[nPer], c1Pdown_c2Pup[nPer];
  Double_t e_c1Pup_c2Pdown[nPer], e_c1Pdown_c2Pup[nPer];
  Double_t ex[nPer] = {0.0}, xvals[nPer];
  
  //Get Data and make ratio
  for (Int_t p=0; p<nPer; p++) {
    xvals[p] = p+0.5;
    
    TString n_Wper =
      Form("%s/Weight/weight_%s%i_%s_%s%s_%s.root", thisDirPath.Data(),
	   physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	   period[p].Data(), minimizer.Data());
    TFile *f_Wper = OpenFile(n_Wper);

    TGraphErrors *gPer_upS[nGraph], *gPer_downS[nGraph];
    for (Int_t g=0; g<nGraph; g++) {
      gPer_upS[g] = (TGraphErrors*) f_Wper->Get(upSNames[g]);
      gPer_downS[g] = (TGraphErrors*) f_Wper->Get(downSNames[g]);
    }//graph loop

    //Integrated
    Double_t *y_c1Pup=gPer_upS[8]->GetY();
    Double_t *y_c1Pdown=gPer_upS[9]->GetY();
    Double_t *y_c2Pup=gPer_downS[8]->GetY();
    Double_t *y_c2Pdown=gPer_downS[9]->GetY();

    Double_t *ey_c1Pup=gPer_upS[8]->GetEY();
    Double_t *ey_c1Pdown=gPer_upS[9]->GetEY();
    Double_t *ey_c2Pup=gPer_downS[8]->GetEY();
    Double_t *ey_c2Pdown=gPer_downS[9]->GetEY();

    c1Pup_c2Pdown[p] = y_c1Pup[0]/y_c2Pdown[0];
    c1Pdown_c2Pup[p] = y_c1Pdown[0]/y_c2Pup[0];

    cout << y_c1Pup[0] << " " << y_c2Pdown[0] << endl;
    cout << y_c1Pdown[0] << " " << y_c2Pup[0] << endl;

    e_c1Pup_c2Pdown[p] = RatioError(y_c1Pup[0], y_c2Pdown[0],
				    ey_c1Pup[0], ey_c2Pdown[0]);

    e_c1Pdown_c2Pup[p] = RatioError(y_c1Pdown[0], y_c2Pup[0],
				    ey_c1Pdown[0], ey_c2Pup[0]);

    cout << "\n\n"<< endl;
    cout << c1Pup_c2Pdown[p] << "  +/-  " << e_c1Pup_c2Pdown[p] << endl;
    cout << c1Pdown_c2Pup[p] << "  +/-  " << e_c1Pdown_c2Pup[p] << endl;
  }//p period loop
  
  TCanvas* c1 = new TCanvas(); 
  TGraphErrors* g_c1Pup_c2Pdown =
    new TGraphErrors(nPer, xvals, c1Pup_c2Pdown, ex, e_c1Pup_c2Pdown);
  TGraphErrors* g_c1Pdown_c2Pup =
    new TGraphErrors(nPer, xvals, c1Pdown_c2Pup, ex, e_c1Pdown_c2Pup);
  SetUp(g_c1Pup_c2Pdown); SetUp(g_c1Pdown_c2Pup);

  g_c1Pup_c2Pdown->Draw("AP");
  g_c1Pup_c2Pdown->SetMarkerColor(kRed);
  g_c1Pdown_c2Pup->Draw("Psame");
  g_c1Pdown_c2Pup->SetMarkerColor(kBlue);
  g_c1Pdown_c2Pup->SetMarkerStyle(24);
  OffSet(g_c1Pdown_c2Pup);

  /*//Write output/final settings
  TString fOutput =
    Form("%s/WAvg/wAvg_%s%i_%s_%s_%s.root", thisDirPath.Data(),
	 physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	 minimizer.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");

    for (Int_t g=0; g<nGraph; g++) {
      g_upS[g]->Write(upSNames[g]);
      g_downS[g]->Write(downSNames[g]);
    }
  }

  cout << "\nSettings______" << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Mass range considered:      " << Mtype << endl;
  cout << "Binned in which DY physics: " << physBinned << endl;
  cout << "Production considered:      " << production << endl;
  cout << "Minimizer algorithmn:      " << minimizer << endl;
  cout << "\nTo write output file:       " << toWrite << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;//*/
}
