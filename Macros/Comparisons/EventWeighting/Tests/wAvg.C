#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

void wAvg(TString start=""){
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

  const Int_t nGraph = 14;
  Double_t upS[nGraph][nBins] ={0.0};
  Double_t e_upS[nGraph][nBins] ={0.0};
  Double_t downS[nGraph][nBins] ={0.0};
  Double_t e_downS[nGraph][nBins] ={0.0};
  Double_t ex[nBins] ={0.0};
  Double_t *xvals;

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

  TGraphErrors *g_upS[nGraph], *g_downS[nGraph];  
  
  //Get Data and add for wAvg
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
      
      Double_t *yval_upS = gPer_upS[g]->GetY();
      Double_t *e_yval_upS = gPer_upS[g]->GetEY();
      Double_t *yval_downS = gPer_downS[g]->GetY();
      Double_t *e_yval_downS = gPer_downS[g]->GetEY();
      for (Int_t b=0; b<gPer_upS[g]->GetN(); b++) {
	upS[g][b] += yval_upS[b]/(e_yval_upS[b]*e_yval_upS[b]);
	e_upS[g][b] += 1.0/(e_yval_upS[b]*e_yval_upS[b]);

	downS[g][b] += yval_downS[b]/(e_yval_downS[b]*e_yval_downS[b]);
	e_downS[g][b] += 1.0/(e_yval_downS[b]*e_yval_downS[b]);
      }//bin loop
    }//graph loop
    
    if (p==0) xvals = gPer_upS[0]->GetX();
    else if (p==nPer-1){//Final wAvg calculation
      for (Int_t g=0; g<nGraph; g++) {
	for (Int_t b=0; b<gPer_upS[g]->GetN(); b++) {
	  upS[g][b] /= e_upS[g][b];
	  e_upS[g][b] = TMath::Sqrt(1.0/e_upS[g][b]);

	  downS[g][b] /= e_downS[g][b];
	  e_downS[g][b] = TMath::Sqrt(1.0/e_downS[g][b]);      
	}//bin loop

	//Make all graphs
	g_upS[g] =
	  new TGraphErrors(gPer_upS[g]->GetN(), xvals, upS[g], ex, e_upS[g]);
	g_downS[g] = new TGraphErrors(gPer_downS[g]->GetN(), xvals, downS[g],
				      ex, e_downS[g]);
	SetUp(g_upS[g]); SetUp(g_downS[g]);
	g_upS[g]->SetTitle(upSNames[g]);
	g_downS[g]->SetTitle(downSNames[g]);
      }//graph loop
    }
    
  }//period loop

  //Draw data
  TCanvas* cupS = new TCanvas(); cupS->Divide(3, 2);
  TCanvas* cdownS = new TCanvas(); cdownS->Divide(3, 2);
  for (Int_t i=0; i<6; i++) {
    cupS->cd(i+1);
    g_upS[i]->Draw("AP");

    cdownS->cd(i+1);
    g_downS[i]->Draw("AP");
  }
  
  cupS->cd(1);
  OffSet(g_upS[6]); g_upS[6]->SetMarkerColor(kRed);
  g_upS[6]->Draw("Psame");

  cdownS->cd(1);
  OffSet(g_downS[6]); g_downS[6]->SetMarkerColor(kRed);
  g_downS[6]->Draw("Psame");

  //Write output/final settings
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
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
