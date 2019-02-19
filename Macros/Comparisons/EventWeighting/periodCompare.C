#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void MakePull(TH1D *h, vector<Double_t> &val, vector<Double_t> &e_val);

void periodCompare(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//# of physBinned bins
  TString Mtype ="HMDY";
  TString physBinned ="xN";//"xN", "xPi", "xF", "pT", "M"
  TString production ="slot1";//"t3", "slot1"
  TString minimizer ="MINOS";//MIGRAD, HESSE, MINOS
  
  Bool_t toWrite =false;
  //Setup_______________
  
  //Basic Setup
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting/Data";
  const Int_t nPer =9;
  TString period[nPer] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};
  Int_t icolor[nPer] = {kBlue+2, kRed+2, kGreen+2, kMagenta+2, kCyan+2,
			    kBlue, kRed, kGreen, kMagenta};
  Int_t imarker[nPer] = {20, 21, 22, 23, 24, 25, 26, 27, 28};
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetNColumns(3);

  const Int_t nGraph = 14;
  TString upSNames[nGraph] =
    {"g_A_upS_phys", "g_L_upS_Pup_phys", "g_L_upS_Pup_phys",
     "g_a1_upS_phys", "g_a2_upS_phys", "g_a3_upS_phys",
     "g_A_upS_phys_expected",
     "g_A_upS_int", "g_L_upS_Pup_int", "g_L_upS_Pup_int",
     "g_a1_upS_int", "g_a2_upS_int", "g_a3_upS_int",
     "g_A_upS_int_expected"};
  TString downSNames[nGraph] =
    {"g_A_downS_phys", "g_L_downS_Pup_phys", "g_L_downS_Pup_phys",
     "g_a1_downS_phys", "g_a2_downS_phys", "g_a3_downS_phys",
     "g_A_downS_phys_expected",
     "g_A_downS_int", "g_L_downS_Pup_int", "g_L_downS_Pup_int",
     "g_a1_downS_int", "g_a2_downS_int", "g_a3_downS_int",
     "g_A_downS_int_expected"};
  vector<Double_t> A_upS, eA_upS;
  vector<Double_t> A_downS, eA_downS;
  vector<Double_t> A_all, eA_all;
  
  //Get Data and add for wAvg
  TCanvas* cPeriod = new TCanvas(); cPeriod->Divide(2);
  for (Int_t p=0; p<nPer; p++) {
    TString n_Wper =
      Form("%s/Weight/weight_%s%i_%s_%s%s_%s.root", thisDirPath.Data(),
	   physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	   period[p].Data(), minimizer.Data());
    TFile *f_Wper = OpenFile(n_Wper);

    TGraphErrors *gPer_upS[nGraph], *gPer_downS[nGraph];
    for (Int_t g=0; g<nGraph; g++) {
      gPer_upS[g] = (TGraphErrors*) f_Wper->Get(upSNames[g]);
      gPer_downS[g] = (TGraphErrors*) f_Wper->Get(downSNames[g]);

      FinalSetup(gPer_upS[g]);
      FinalSetup(gPer_downS[g]);
      OffSetPercent(gPer_upS[g], p*0.02);
      OffSetPercent(gPer_downS[g], p*0.02);
      gPer_upS[g]->SetMarkerColor(icolor[p]);
      gPer_upS[g]->SetMarkerStyle(imarker[p]);

      gPer_downS[g]->SetMarkerColor(icolor[p]);
      gPer_downS[g]->SetMarkerStyle(imarker[p]);
    }//graph loop

    Double_t *yVal_upS=gPer_upS[0]->GetY();
    Double_t *e_yVal_upS=gPer_upS[0]->GetEY();
    Double_t *yVal_downS=gPer_downS[0]->GetY();
    Double_t *e_yVal_downS=gPer_downS[0]->GetEY();
    for (Int_t bi=0; bi<nBins; bi++) {
      A_upS.push_back(yVal_upS[bi]);
      eA_upS.push_back(e_yVal_upS[bi]);

      A_downS.push_back(yVal_downS[bi]);
      eA_downS.push_back(e_yVal_downS[bi]);

      A_all.push_back(yVal_upS[bi]);
      eA_all.push_back(e_yVal_upS[bi]);
      A_all.push_back(yVal_downS[bi]);
      eA_all.push_back(e_yVal_downS[bi]);
    }

    legend->AddEntry(gPer_upS[p], Form("%s", period[p].Data() ), "p");
    if (p==0) {
      cPeriod->cd(1);
      gPer_upS[7]->Draw("AP");
      gPer_upS[7]->GetYaxis()->SetRangeUser(-1.0, 1.0);
      gPer_upS[7]->GetXaxis()->SetLimits(0, 0.5);
      DrawLine(gPer_upS[7], 0.0);

      cPeriod->cd(2);
      gPer_downS[7]->Draw("AP");
      gPer_downS[7]->GetYaxis()->SetRangeUser(-1.0, 1.0);
      gPer_downS[7]->GetXaxis()->SetLimits(0, 0.5);
      DrawLine(gPer_downS[7], 0.0);
    }
    else{
      cPeriod->cd(1);
      gPer_upS[7]->Draw("Psame");

      cPeriod->cd(2);
      gPer_downS[7]->Draw("Psame");
    }
    
  }//p period loop
  
  cPeriod->cd(1);
  legend->Draw("same");

  //Make pull
  TH1D* hPull_upS = new TH1D("hPull_upS", "hPull_upS", 12, -4, 4);
  MakePull(hPull_upS, A_upS, eA_upS);
  TH1D* hPull_downS = new TH1D("hPull_downS", "hPull_downS", 12, -4, 4);
  MakePull(hPull_downS, A_downS, eA_downS);
  TH1D* hPull_all = new TH1D("hPull_all", "hPull_all", 12, -4, 4);
  MakePull(hPull_all, A_all, eA_all);
  
  TCanvas* cPulls = new TCanvas(); cPulls->Divide(3);
  cPulls->cd(1);
  hPull_upS->Draw("E1");

  cPulls->cd(2);
  hPull_downS->Draw("E1");

  cPulls->cd(3);
  hPull_all->Draw("E1");

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


void MakePull(TH1D *h, vector<Double_t> &val, vector<Double_t> &e_val){
  Double_t wSigma2;
  Double_t wAvg = WeightedAvgAndError(val, e_val, wSigma2);
  wSigma2 *= wSigma2;
  
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
