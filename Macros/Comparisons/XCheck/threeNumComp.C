#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void threeNumComp(){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString physBinned ="xN";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/
  
  TString period = "W07";
  TString corrected ="uncorr"; //"uncorr"=NOT corrected or ""=corrected
  //Setup_______________

  //Basic Setup
  TString localPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/XCheck";
  TString pathAN="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";

  TString fname =
    Form("%s/GeoMean4Targ_%s_%s_%s_%s%s_%s%s%i_%s_%s.root", pathAN.Data(),
	 whichFit.Data(), period.Data(), Mtype.Data(), process.Data(),
	 lrMrange.Data(), binRange.Data(), physBinned.Data(), nBins,
	 production.Data(), additionalCuts.Data());
  TFile *f_me = OpenFile(fname);
  TGraphErrors *g_me =
    (TGraphErrors*) f_me->Get(Form("AN%s", corrected.Data()));

  Double_t *xvals = g_me->GetX();
  Double_t ex[nBins] = {0.0};
    
  Double_t An_xc[] = {-0.00864868, -0.00510964, -0.0273411};
  Double_t An_xc_flip[] = {0.00864868, 0.00510964, 0.0273411};
  Double_t eAn_xc[] = {0.0287779, 0.0277795, 0.0288889};

  TGraphErrors* g_xc = new TGraphErrors(nBins, xvals, An_xc, ex, eAn_xc);
  SetUp(g_xc); g_xc->SetMarkerColor(kRed);

  TGraphErrors* g_xc_flip =
    new TGraphErrors(nBins, xvals, An_xc_flip, ex, eAn_xc);
  SetUp(g_xc_flip); g_xc_flip->SetMarkerColor(kBlue);

  TCanvas* c1 = new TCanvas();
  g_me->Draw("AP"); g_me->GetYaxis()->SetRangeUser(-0.075, 0.075);
  g_xc->Draw("Psame");
  OffSetPercent(g_xc, 0.04);
  g_xc_flip->Draw("Psame");
  OffSetPercent(g_xc_flip);

  DrawLine(g_me, 0.0);
}
