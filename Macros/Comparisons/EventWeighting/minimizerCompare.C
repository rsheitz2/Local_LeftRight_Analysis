#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void FinalSetupLocal(TGraphErrors *g, TString xName);

void FinalLocalIntegrated(TGraphErrors *g);

void minimizerCompare(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//# of physBinned bins
  TString Mtype ="HMDY";
  TString production ="slot1";//"t3", "slot1"
  
  Bool_t toWrite =false;
  //Setup_______________
  
  //Basic Setup
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting/Data";
  TString geoMeanPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/WAvg";
  const Int_t nMin =2, nPhys=5;
  TString minimizer[nMin] ={"MINOS", "MIGRAD"};//MIGRAD, HESSE, MINOS
  TString physBinned[nPhys] =
    {"xN", "xPi", "xF", "pT", "M"};
  TString xNames[nPhys] =
    {"x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};

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
  TCanvas* cA_upS = new TCanvas();
  cA_upS->SetLeftMargin(0.2); cA_upS->Divide(nPhys+1, 1, 0.0, 0.01);
  TCanvas* cA_downS = new TCanvas();
  cA_downS->SetLeftMargin(0.2); cA_downS->Divide(nPhys+1, 1, 0.0, 0.01);
  
  //Get Data and Draw by physics binning
  for (Int_t p=0; p<nPhys; p++) {
    for (Int_t m=0; m<nMin; m++) {

      //Get Event weighting
      TString n_evW =
	Form("%s/WAvg/wAvg_%s%i_%s_%s_%s.root", thisDirPath.Data(),
	     physBinned[p].Data(), nBins, Mtype.Data(), production.Data(),
	     minimizer[m].Data());
      TFile *f_evW = OpenFile(n_evW);

      TGraphErrors *gUpS_evW[nGraph], *gDownS_evW[nGraph];
      for (Int_t g=0; g<nGraph; g++) {
	gUpS_evW[g] = (TGraphErrors*) f_evW->Get(upSNames[g]);
	gDownS_evW[g] = (TGraphErrors*) f_evW->Get(downSNames[g]);

	FinalSetupLocal(gUpS_evW[g], xNames[p]);
	FinalSetupLocal(gDownS_evW[g], xNames[p]);
	if (m != 0){
	  gUpS_evW[g]->SetMarkerStyle(26); OffSet(gUpS_evW[g]);
	  gUpS_evW[g]->SetMarkerColor(kBlack);
	  gDownS_evW[g]->SetMarkerStyle(26);  OffSet(gDownS_evW[g]);
	  gDownS_evW[g]->SetMarkerColor(kBlack);
	}
      }//graph loop
      
      if (p==0) {//Integrate
	cA_upS->cd(p+1);//UpS  
	if (m==0) {
	  gUpS_evW[7]->Draw("AP");
	  gUpS_evW[7]->GetYaxis()->SetRangeUser(-0.6, 0.6);
	  DrawLine(gUpS_evW[7], 0.0);
	  FinalLocalIntegrated(gUpS_evW[7]);
	}
	else gUpS_evW[7]->Draw("Psame");

	cA_downS->cd(p+1);//DownS
	if (m==0) {
	  gDownS_evW[7]->Draw("AP");
	  gDownS_evW[7]->GetYaxis()->SetRangeUser(-0.6, 0.6);
	  DrawLine(gDownS_evW[7], 0.0);
	  FinalLocalIntegrated(gDownS_evW[7]);
	}
	else gDownS_evW[7]->Draw("Psame");
      }//p==0
      cA_upS->cd(p+2);//UpS
      if (m==0){
	gUpS_evW[0]->Draw("AP");
	gUpS_evW[0]->GetYaxis()->SetRangeUser(-0.6, 0.6);
	DrawLine(gUpS_evW[0], 0.0);
      }
      else gUpS_evW[0]->Draw("Psame");

      cA_downS->cd(p+2);//DownS
      if (m==0){
	gDownS_evW[0]->Draw("AP");
	gDownS_evW[0]->GetYaxis()->SetRangeUser(-0.6, 0.6);
	DrawLine(gDownS_evW[0], 0.0);
      }
      else gDownS_evW[0]->Draw("Psame");

    }//minimizer loop
  }//p physics loop

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


void FinalSetupLocal(TGraphErrors *g, TString xName){
  FinalSetup(g);
  FinalClearTitles(g);

  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.8);
  g->SetFillStyle(0);
  g->SetLineWidth(2);
  g->SetTitle("");

  SetTitleName(g, xName, "x");
}


void FinalLocalIntegrated(TGraphErrors *g){
  g->GetXaxis()->SetLimits(0.0, 0.35);
  g->GetXaxis()->SetLabelSize(0.0);
  g->GetXaxis()->SetTickSize(0.0);
}
