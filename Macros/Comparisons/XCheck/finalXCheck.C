#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void finalXCheck(){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  
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

  //Micheal uncorrected values
  Double_t An_qT[] = {0.0410412, 0.0175735, -0.024854};
  Double_t An_xPi[] = {-0.0090485, 0.0282803, 0.0199568};
  Double_t An_xF[] = {0.00304508, 0.0221673, 0.0134429}; 
  Double_t An_xN[] = {0.00864868, 0.00510964,  0.0273411};
  
  Double_t eAn_qT[] = {0.0270867, 0.0292951, 0.0292208};
  Double_t eAn_xPi[] = {0.028301,  0.0284889,  0.0286681};
  Double_t eAn_xF[] = {0.0286906, 0.0286612,  0.0280643};
  Double_t eAn_xN[] = {0.0287778, 0.0277795, 0.0288882};
  
  const Int_t nPhysBinned =4;
  TString physBinned[] ={"xN", "xPi", "xF", "pT"};
  TString xNames[] = {"x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)"};  

  TCanvas* cAsym = new TCanvas(); cAsym->Divide(nPhysBinned, 1, 0, 0);

  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    //Get my data
    TString fname =
      Form("%s/GeoMean4Targ_%s_%s_%s_%s%s_%s%s%i_%s_%s.root", pathAN.Data(),
	   whichFit.Data(), period.Data(), Mtype.Data(), process.Data(),
	   lrMrange.Data(), binRange.Data(), physBinned[phys].Data(), nBins,
	   production.Data(), additionalCuts.Data());
    TFile *f_me = OpenFile(fname);
    TGraphErrors *g_me =
      (TGraphErrors*) f_me->Get(Form("AN%s", corrected.Data()));

    //Draw my values
    cAsym->cd(phys + 1); gPad->SetFrameLineWidth(2);
    g_me->Draw("AP"); g_me->GetYaxis()->SetRangeUser(-0.075, 0.075);
    FinalSetup(g_me); g_me->SetMarkerSize(1.3);
    DrawLine(g_me, 0.0); g_me->SetTitle(""); g_me->SetMarkerStyle(20);
    SetTitleName(g_me, xNames[phys], "x");

    //Draw micheals values
    Double_t *xvals = g_me->GetX();
    Double_t ex[nBins] = {0.0};
    cout << xNames[phys] << endl;
    for (Int_t i=0; i<nBins; i++) {
      cout << xvals[i] << endl;
    }
    cout << " " << endl;


    Double_t *An_xc, *eAn_xc; 
    if (physBinned[phys] == "pT"){
      An_xc = An_qT; eAn_xc = eAn_qT;
    }
    else if (physBinned[phys] == "xN"){
      An_xc = An_xN; eAn_xc = eAn_xN;
    }
    else if (physBinned[phys] == "xPi"){
      An_xc = An_xPi; eAn_xc = eAn_xPi;
    }
    else if (physBinned[phys] == "xF"){
      An_xc = An_xF; eAn_xc = eAn_xF;
    }
    else{
      cout << "Problems" << endl;
      exit(EXIT_FAILURE);
    }

    TGraphErrors* g_xc = new TGraphErrors(nBins, xvals, An_xc, ex, eAn_xc);
    FinalSetup(g_xc); g_xc->SetMarkerColor(kRed); g_xc->SetMarkerStyle(22);
    g_xc->SetMarkerSize(1.3); OffSetPercent(g_xc, 0.05);
    g_xc->Draw("Psame");
  }//phys loop

}
