#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"

Int_t CheckErrors(Double_t eY){
  if (eY < 0.15){
    return 0;
  }

  return 1;
}


void wAvg(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//# of physBinned bins
  TString Mtype ="HMDY";
  TString physBinned ="xN";//"xN", "xPi", "xF", "pT", "M"
  TString production ="slot1";//"t3", "slot1"
  
  Bool_t toWrite =false;
  //Setup_______________
  
  //Basic Setup
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting/Data";
  const Int_t nPer =9;
  TString period[nPer] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};
  
  Double_t A_upS_phys[nBins] ={0.0}, e_A_upS_phys[nBins] ={0.0};
  Double_t A_downS_phys[nBins] ={0.0}, e_A_downS_phys[nBins] ={0.0};
  Double_t ex[nBins] ={0.0};
  Double_t *xvals;

  //Get Data and add for wAvg
  for (Int_t p=0; p<nPer; p++) {
    TString n_Wper =
      Form("%s/Weight/weight_%s%i_%s_%s%s.root", thisDirPath.Data(),
	   physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	   period[p].Data());
    TFile *f_Wper = OpenFile(n_Wper);
    TGraphErrors *gPer_A_upS = (TGraphErrors*) f_Wper->Get("g_A_upS_phys");
    TGraphErrors *gPer_A_downS = (TGraphErrors*) f_Wper->Get("g_A_downS_phys");

    Double_t *yval_upS = gPer_A_upS->GetY();
    Double_t *e_yval_upS = gPer_A_upS->GetEY();
    Double_t *yval_downS = gPer_A_downS->GetY();
    Double_t *e_yval_downS = gPer_A_downS->GetEY();
    for (Int_t i=0; i<nBins; i++) {
      if ( CheckErrors(e_yval_upS[i]) ){
	A_upS_phys[i] += yval_upS[i]/(e_yval_upS[i]*e_yval_upS[i]);
	e_A_upS_phys[i] += 1.0/(e_yval_upS[i]*e_yval_upS[i]);
      }
      else{
	cout << "Upstream period:  " << period[p] << ", bin:  " << i
	     << "   was skipped" << endl;
      }

      if ( CheckErrors(e_yval_downS[i]) ) {
	A_downS_phys[i] += yval_downS[i]/(e_yval_downS[i]*e_yval_downS[i]);
	e_A_downS_phys[i] += 1.0/(e_yval_downS[i]*e_yval_downS[i]);
      }
      else{
	cout << "Downstream period:  " << period[p] << ", bin:  " << i
	     << "   was skipped" << endl;
      }
    }

    if (p==0) xvals = gPer_A_upS->GetX();
  }

  //Final wAvg calculation
  for (Int_t i=0; i<nBins; i++) {
    A_upS_phys[i] /= e_A_upS_phys[i];
    e_A_upS_phys[i] = TMath::Sqrt(1.0/e_A_upS_phys[i]);

    A_downS_phys[i] /= e_A_downS_phys[i];
    e_A_downS_phys[i] = TMath::Sqrt(1.0/e_A_downS_phys[i]);
  }

  TGraphErrors* g_A_upS_phys =
    new TGraphErrors(nBins, xvals, A_upS_phys, ex, e_A_upS_phys);
  TGraphErrors* g_A_downS_phys =
    new TGraphErrors(nBins, xvals, A_downS_phys, ex, e_A_downS_phys);
  SetUp(g_A_upS_phys); SetUp(g_A_downS_phys); 

  //Draw data
  TCanvas* c1 = new TCanvas();
  g_A_upS_phys->Draw("AP");
  g_A_downS_phys->Draw("Psame"); g_A_downS_phys->SetMarkerColor(kBlue);
  if (Mtype=="HMDY"){ g_A_upS_phys->GetYaxis()->SetRangeUser(-0.5, 0.5); }
  else { g_A_upS_phys->GetYaxis()->SetRangeUser(-0.08, 0.08); }
  g_A_upS_phys->SetTitle("Final Amplitude from wAveraging");
  g_A_upS_phys->GetXaxis()->SetTitle(physBinned);
  DrawLine(g_A_upS_phys, 0.0);

  //Write output/final settings
  TString fOutput =
    Form("%s/WAvg/wAvg_%s%i_%s_%s.root", thisDirPath.Data(),
	 physBinned.Data(), nBins, Mtype.Data(), production.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_A_upS_phys->Write("g_A_upS_phys");
    g_A_downS_phys->Write("g_A_downS_phys");
  }

  cout << "\nSettings______" << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Mass range considered:      " << Mtype << endl;
  cout << "Binned in which DY physics: " << physBinned << endl;
  cout << "Production considered:      " << production << endl;
  cout << "\nTo write output file:       " << toWrite << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
