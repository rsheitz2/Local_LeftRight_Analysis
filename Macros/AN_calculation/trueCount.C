#include "include/helperFunctions.h"

Double_t MakeAsym(Double_t L, Double_t R, Double_t P);

Double_t MakeAsymError(Double_t L, Double_t R, Double_t e_L, Double_t e_R,
		       Double_t P);

void GetLRerror(Double_t *leftCounts, Double_t *rightCounts,
		Double_t *e_leftCounts, Double_t *e_rightCounts, Int_t nBins);

void trueCount(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//# of physBinned bins
  TString period_Mtype ="WAll_HMDY";
  TString binRange ="43_85";
  Double_t Mmin =4.30;//LR Mass minimum
  Double_t Mmax =8.50;//LR Mass maximum
  TString physBinned ="M";//"xF", "pT"
  TString process ="DY";//JPsi, psi, DY
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";
    
  Bool_t toWrite =false;
  //Setup_______________
  
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =
    Form("leftRight_byTarget_%s%.2f_%.2f_%ibins%s_150hbin_%s_%s.root",
	 period_Mtype.Data(), Mmin, Mmax, nBins, binRange.Data(),
	 production.Data(), additionalCuts.Data() );
  
  if (start==""){
    cout<<"Script outputs AN and left/right counts per target and polarization";
    cout << " using functional mass fitting for a given fit" << endl;
    cout << "Outputs are in the formate needed for GeoMean4Targ.C" << endl;
    cout << "Script also outputs information on the fit quality" <<endl;
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  trueCount.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'trueCount.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Real data path:             " << pathRD << endl;
    cout << "Real data file considered:  " << RDfile << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //Get Input Files from Local_leftRight
  TFile *fRD  = TFile::Open(pathRD + RDfile);
  if (!fRD ){
    cout << "RD file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Determine polarization factor
  Double_t ex[nBins] = {0.0};
  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  TGraph* g_Dil =(TGraph*)fRD->Get(Form("%s_Dil", physBinned.Data()));
  Double_t *xvals = g_Pol->GetX();
  Double_t Pol[nBins];
  GetPolarization(g_Pol, g_Dil, Pol);
  
  //Get g_left/right by target and pol
  TGraph *gRD_left_upS_up = (TGraph*)
    fRD->Get(Form("%s_left_upstream_up", physBinned.Data()));
  TGraph *gRD_left_upS_down = (TGraph*)
    fRD->Get(Form("%s_left_upstream_down", physBinned.Data()));
  TGraph *gRD_left_downS_up = (TGraph*)
    fRD->Get(Form("%s_left_downstream_up", physBinned.Data()));
  TGraph *gRD_left_downS_down = (TGraph*)
    fRD->Get(Form("%s_left_downstream_down", physBinned.Data()));

  TGraph *gRD_right_upS_up = (TGraph*)
    fRD->Get(Form("%s_right_upstream_up", physBinned.Data()));
  TGraph *gRD_right_upS_down = (TGraph*)
    fRD->Get(Form("%s_right_upstream_down", physBinned.Data()));
  TGraph *gRD_right_downS_up = (TGraph*)
    fRD->Get(Form("%s_right_downstream_up", physBinned.Data()));
  TGraph *gRD_right_downS_down = (TGraph*)
    fRD->Get(Form("%s_right_downstream_down", physBinned.Data()));

  //L/R counts
  Double_t *LR_upS_up_L = gRD_left_upS_up->GetY();
  Double_t *LR_upS_up_R = gRD_right_upS_up->GetY();
  Double_t *LR_upS_down_L = gRD_left_upS_down->GetY();
  Double_t *LR_upS_down_R = gRD_right_upS_down->GetY();

  Double_t *LR_downS_up_L = gRD_left_downS_up->GetY();
  Double_t *LR_downS_up_R = gRD_right_downS_up->GetY();
  Double_t *LR_downS_down_L = gRD_left_downS_down->GetY();
  Double_t *LR_downS_down_R = gRD_right_downS_down->GetY();
  
  Double_t e_LR_upS_up_L[nBins], e_LR_upS_up_R[nBins];
  Double_t e_LR_upS_down_L[nBins], e_LR_upS_down_R[nBins];
  Double_t e_LR_downS_up_L[nBins], e_LR_downS_up_R[nBins];
  Double_t e_LR_downS_down_L[nBins], e_LR_downS_down_R[nBins];

  //Determine L/R errors
  GetLRerror(LR_upS_up_L, LR_upS_up_R, e_LR_upS_up_L, e_LR_upS_up_R, nBins);
  GetLRerror(LR_upS_down_L, LR_upS_down_R, e_LR_upS_down_L, e_LR_upS_down_R,
	      nBins);
  GetLRerror(LR_downS_up_L, LR_downS_up_R, e_LR_downS_up_L, e_LR_downS_up_R,
	      nBins);
  GetLRerror(LR_downS_down_L, LR_downS_down_R, e_LR_downS_down_L,
	      e_LR_downS_down_R, nBins);
  
  //L/R count graphs/Drawing
  TGraphErrors* g_Left_upS_up =
    new TGraphErrors(nBins, xvals, LR_upS_up_L, ex, e_LR_upS_up_L);
  TGraphErrors* g_Right_upS_up =
    new TGraphErrors(nBins, xvals, LR_upS_up_R, ex, e_LR_upS_up_R);
  TGraphErrors* g_Left_upS_down =
    new TGraphErrors(nBins, xvals, LR_upS_down_L, ex, e_LR_upS_down_L);
  TGraphErrors* g_Right_upS_down =
    new TGraphErrors(nBins, xvals, LR_upS_down_R, ex, e_LR_upS_down_R);
  TGraphErrors* g_Left_downS_up =
    new TGraphErrors(nBins, xvals, LR_downS_up_L, ex, e_LR_downS_up_L);
  TGraphErrors* g_Right_downS_up =
    new TGraphErrors(nBins, xvals, LR_downS_up_R, ex, e_LR_downS_up_R);
  TGraphErrors* g_Left_downS_down =
    new TGraphErrors(nBins, xvals, LR_downS_down_L, ex, e_LR_downS_down_L);
  TGraphErrors* g_Right_downS_down =
    new TGraphErrors(nBins, xvals, LR_downS_down_R, ex, e_LR_downS_down_R);
  
  SetUpTGraph(g_Left_upS_up); SetUpTGraph(g_Right_upS_up);
  SetUpTGraph(g_Left_upS_down); SetUpTGraph(g_Right_upS_down);
  SetUpTGraph(g_Left_downS_up); SetUpTGraph(g_Right_downS_up);
  SetUpTGraph(g_Left_downS_down); SetUpTGraph(g_Right_downS_down);

  TCanvas* cLR = new TCanvas("LR counts"); cLR->Divide(2, 2);
  cLR->cd(1); g_Left_upS_up->Draw("AP");//Upstream
  g_Left_upS_up->SetTitle("L/R upS_up");
  g_Right_upS_up->Draw("Psame"); g_Right_upS_up->SetMarkerColor(kRed);
  
  cLR->cd(2); g_Left_upS_down->Draw("AP");
  g_Left_upS_down->SetTitle("L/R upS_down");
  g_Right_upS_down->Draw("Psame");
  g_Right_upS_down->SetMarkerColor(kRed);
  
  cLR->cd(3); g_Left_downS_up->Draw("AP");//Downstream
  g_Left_downS_up->SetTitle("L/R downS_up");
  g_Right_downS_up->Draw("Psame");
  g_Right_downS_up->SetMarkerColor(kRed);
  
  cLR->cd(4); g_Left_downS_down->Draw("AP");
  g_Left_downS_down->SetTitle("L/R downS_down");
  g_Right_downS_down->Draw("Psame");
  g_Right_downS_down->SetMarkerColor(kRed);

  //Make AN asymmetery per targ && pol
  Double_t AN_upS_up[nBins], e_AN_upS_up[nBins];
  Double_t AN_upS_down[nBins], e_AN_upS_down[nBins];
  Double_t AN_downS_up[nBins], e_AN_downS_up[nBins];
  Double_t AN_downS_down[nBins], e_AN_downS_down[nBins];
  
  for (Int_t bi=0; bi<nBins; bi++) {
    AN_upS_up[bi] = MakeAsym(LR_upS_up_L[bi], LR_upS_up_R[bi], Pol[bi]);
    e_AN_upS_up[bi]
      = MakeAsymError(LR_upS_up_L[bi], LR_upS_up_R[bi], e_LR_upS_up_L[bi],
		      e_LR_upS_up_R[bi],Pol[bi]);
    
    AN_upS_down[bi] = MakeAsym(LR_upS_down_L[bi], LR_upS_down_R[bi], Pol[bi]);
    e_AN_upS_down[bi]
      =MakeAsymError(LR_upS_down_L[bi], LR_upS_down_R[bi], e_LR_upS_down_L[bi],
		     e_LR_upS_down_R[bi],Pol[bi]);
    
    AN_downS_up[bi] = MakeAsym(LR_downS_up_L[bi], LR_downS_up_R[bi], Pol[bi]);
    e_AN_downS_up[bi] =
      MakeAsymError(LR_downS_up_L[bi], LR_downS_up_R[bi], e_LR_downS_up_L[bi],
		    e_LR_downS_up_R[bi],Pol[bi]);
    
    AN_downS_down[bi]=MakeAsym(LR_downS_down_L[bi],LR_downS_down_R[bi],Pol[bi]);
    e_AN_downS_down[bi]
      =MakeAsymError(LR_downS_down_L[bi],LR_downS_down_R[bi],
		     e_LR_downS_down_L[bi],e_LR_downS_down_R[bi], Pol[bi]);
  }

  //Make/Draw AN graphs per targ && pol
  TGraphErrors *g_AN_upS_up
    = new TGraphErrors(nBins, xvals,AN_upS_up, ex, e_AN_upS_up);
  TGraphErrors *g_AN_upS_down =
    new TGraphErrors(nBins, xvals,AN_upS_down, ex, e_AN_upS_down);
  TGraphErrors *g_AN_downS_up =
    new TGraphErrors(nBins, xvals,AN_downS_up, ex, e_AN_downS_up);
  TGraphErrors *g_AN_downS_down =
    new TGraphErrors(nBins, xvals,AN_downS_down, ex, e_AN_downS_down);
  
  SetUpTGraph(g_AN_upS_up); SetUpTGraph(g_AN_upS_down);
  SetUpTGraph(g_AN_downS_up); SetUpTGraph(g_AN_downS_down);
  
  TCanvas* cAN = new TCanvas("AN by tar && pol"); cAN->Divide(2, 2);
  Double_t yMax;
  if (process =="JPsi")  yMax =0.5;
  else if (process =="psi")  yMax =0.5;
  else if (process =="DY")  yMax =0.5;
  cAN->cd(1);
  g_AN_upS_up->Draw("AP"); g_AN_upS_up->SetTitle("AN_upS_up");
  g_AN_upS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_upS_up, 0.0);

  cAN->cd(2);
  g_AN_upS_down->Draw("AP"); g_AN_upS_up->SetTitle("AN_upS_down");
  g_AN_upS_down->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_upS_down, 0.0);

  cAN->cd(3);
  g_AN_downS_up->Draw("AP"); g_AN_upS_up->SetTitle("AN_downS_up");
  g_AN_downS_up->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_downS_up, 0.0);

  cAN->cd(4);
  g_AN_downS_down->Draw("AP"); g_AN_upS_up->SetTitle("AN_downS_down");
  g_AN_downS_down->GetYaxis()->SetRangeUser(-yMax, yMax);
  DrawLine(g_AN_downS_down, 0.0);

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation";
  TString fOutput
    = Form("%s/Data/trueCount/trueCount_%s_%s%.2f_%.2f_%s%s%i",
	   thisDirPath.Data(), period_Mtype.Data(), process.Data(), Mmin, Mmax,
	   binRange.Data(), physBinned.Data(), nBins);
  fOutput += Form("_%s_%s_corr.root", production.Data(), additionalCuts.Data());
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TList *doc = new TList();
    doc->Add((TObject*)(new TObjString(pathRD+"\n"+RDfile)));
    doc->Write("InputData");
    
    g_AN_upS_up->Write("AN_upS_up");
    g_AN_upS_down->Write("AN_upS_down");
    g_AN_downS_up->Write("AN_downS_up");
    g_AN_downS_down->Write("AN_downS_down");

    g_Left_upS_up->Write("Left_upS_up");
    g_Left_upS_down->Write("Left_upS_down");
    g_Left_downS_up->Write("Left_downS_up");
    g_Left_downS_down->Write("Left_downS_down");

    g_Right_upS_up->Write("Right_upS_up");
    g_Right_upS_down->Write("Right_upS_down");
    g_Right_downS_up->Write("Right_downS_up");
    g_Right_downS_down->Write("Right_downS_down");
  }

  cout << " " << endl;
  cout << "Settings______" << endl;
  cout << "Real data file considered:              " << RDfile << endl;
  cout << "Mass Range for LR asymmetry is:         " << Mmin
       << " - " << Mmax << endl;
  cout << "physBinned nBins times:                 " << nBins << endl;
  cout << "\n";
  cout << "AN for physical process:                " << process << endl;
  cout << "Physics binning is:                     " << physBinned << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}


Double_t MakeAsym(Double_t L, Double_t R, Double_t P){
  Double_t A = L - R;
  A /= ( L + R );
  A /= P;

  return A;
}


Double_t MakeAsymError(Double_t L, Double_t R, Double_t e_L, Double_t e_R,
		       Double_t P){
  Double_t dL2 = e_L*e_L;
  Double_t dR2 = e_R*e_R;

  Double_t e = dL2/( L*L )  + dR2/( R*R );
  e = TMath::Sqrt( e );
  e *= 2.0*L*R/( (L+R)*(L+R) );
  e /= P;
  
  return e;
}


void GetLRerror(Double_t *leftCounts, Double_t *rightCounts,
		Double_t *e_leftCounts, Double_t *e_rightCounts, Int_t nBins){

  for (Int_t i=0; i<nBins; i++) {
    e_leftCounts[i] = TMath::Sqrt(leftCounts[i]);
    e_rightCounts[i] = TMath::Sqrt(rightCounts[i]);
  }
}
