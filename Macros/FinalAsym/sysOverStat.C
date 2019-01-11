#include "include/helperFunctions.h"

void CalSysOverStat(TGraph *gSys, TGraphErrors *gStat, Double_t *sys_stat,
		    Int_t nBins){
  Double_t *ySys = gSys->GetY();
  Double_t *eyStat = gStat->GetEY();

  for (Int_t i=0; i<nBins; i++) {
    if (gSys->GetN() < nBins) sys_stat[i] = ySys[0]/eyStat[i];
    else sys_stat[i] = ySys[i]/eyStat[i];
  }
}


void sysOverStat(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  TString period_Mtype ="WAll_HMDY";//cleanup
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  /*const Int_t nBins =5;//JPsi
    TString Mtype ="LowM_AMDY";
    TString period_Mtype ="WAll_LowM_AMDY";//cleanup
    Int_t hbins =150;
    TString process ="JPsi";//JPsi, psi, DY
    TString lrMrange ="2.00_5.00";
    TString fitMrange ="2.00_8.50";
    TString binRange ="25_43";
    TString whichFit ="thirteen";
    TString production ="slot1";
    TString additionalCuts ="phiS0.53";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString pathAN = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/WAvg";
  TString pathSys = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics";

  const Int_t nPhysBinned =5;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "xN"};
  //Last xN used for integrated values
  
  if (start==""){
    cout <<"Script draws Final AN asymmetry with systematics binned in physics";
    cout << "\nSystematics are input as a TGraph for each physics bin";
    cout << "\n\n\nUsage:" << endl;
    cout << "root \'binned4TargGeoMean.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Asym data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " <<  whichFit << endl;
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  //Setups
  TGraph *gSysStat_alphaAcc[nPhysBinned];
  TGraph *gSysStat_lrCrossOver[nPhysBinned];

  //Get Data file/Get graphs and plot
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName, alphaAccName, lrCrossOverName;

    Int_t nBinsName =nBins;
    if (phys == nPhysBinned-1) {//used for integrated
      nBinsName =1; }
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }

      AsymName =
	Form("%s/wAvg_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), Mtype.Data(), process.Data(),
	     lrMrange.Data(), binRange.Data(), physBinned[phys].Data(),
	     nBinsName, hbins, production.Data(), additionalCuts.Data());

      TString specficPath =
	"FalseAsym/Data/acceptanceFourTargRatio/SystematicError";
      alphaAccName =
	Form("%s/%s/accSys4TargRatio_true_%s_%s%s_%s%s%i_%s_%s.root",
	     pathSys.Data(), specficPath.Data(), period_Mtype.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),
	     physBinned[phys].Data(), nBinsName, production.Data(),
	     additionalCuts.Data());

      specficPath =
	"PhiScut/DetermineCut/Data/SysStatError";
      lrCrossOverName =
	Form("%s/%s/sysStatError_Charles_%s_%s_%s.root",
	     pathSys.Data(), specficPath.Data(), production.Data(),
	     fitMrange.Data(), additionalCuts.Data());
      if (phys ==0){
	cout << "\nSpecfic lrCrossOver Settings:" << endl;
	cout << "Charles MC" << endl;
	cout << " " << endl;
      }
    }
    else {
      cout << "Only works for true fit for now..." << endl;
      exit(EXIT_FAILURE);
      AsymName =
	Form("%s/wAvg_%s%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(), Mtype.Data(),
	     process.Data(), lrMrange.Data(), physBinned[phys].Data(),
	     nBinsName, hbins, production.Data(), additionalCuts.Data());

      TString specficPath =
	"FalseAsym/Data/acceptanceFourTargRatio/SystematicError";
    }

    //True Asymmetry
    TFile *f_AN = OpenFile(AsymName);
    TGraphErrors *g_AN =(TGraphErrors*)f_AN->Get("AN");
    Double_t *xvals = g_AN->GetX();
    
    //Systematic error graphs
    ///////////////
    //alphaAcc
    TFile *f_alphaAcc = OpenFile(alphaAccName);
    TGraph *g_alphaAcc = (TGraph*)f_alphaAcc->Get("gSys");
    Double_t alphaAcc[nBins];
    if (phys == nPhysBinned-1){
      CalSysOverStat(g_alphaAcc, g_AN, alphaAcc, 1);//Integrated
      gSysStat_alphaAcc[phys] = new TGraph(1, xvals, alphaAcc);
    }
    else{
      CalSysOverStat(g_alphaAcc, g_AN, alphaAcc, nBins);
      gSysStat_alphaAcc[phys] = new TGraph(nBins, xvals, alphaAcc);
    }
    SetUp(gSysStat_alphaAcc[phys]);

    //lrCrossOver
    TFile *f_lrCrossOver = OpenFile(lrCrossOverName);
    TGraph *g_lrCrossOver = (TGraph*)f_lrCrossOver->Get("gSys");
    Double_t lrCrossOver[nBins];
    if (phys == nPhysBinned-1){
      CalSysOverStat(g_lrCrossOver, g_AN, lrCrossOver, 1);//Integrated
      gSysStat_lrCrossOver[phys] = new TGraph(1, xvals, lrCrossOver);
    }
    else{
      CalSysOverStat(g_lrCrossOver, g_AN, lrCrossOver, nBins);
      gSysStat_lrCrossOver[phys] = new TGraph(nBins, xvals, lrCrossOver);
    }
    SetUp(gSysStat_lrCrossOver[phys]);
  }//phys binned loop
  
  //Draw Sys over stat
  TCanvas* cAlphaAcc = new TCanvas();
  cAlphaAcc->Divide(nPhysBinned,1,0,0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cAlphaAcc->cd(i+1);
    gSysStat_alphaAcc[i]->Draw("AP");
    gSysStat_alphaAcc[i]->Fit("pol0", "0Q");
    Double_t avg =gSysStat_alphaAcc[i]->GetFunction("pol0")->GetParameter(0);
    gSysStat_alphaAcc[i]->SetTitle(Form("AlphaAcceptance %0.2f", avg));
    gSysStat_alphaAcc[i]->GetYaxis()->SetRangeUser(0, 0.6);
  }

  TCanvas* clrCrossOver = new TCanvas();
  clrCrossOver->Divide(nPhysBinned,1,0,0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    clrCrossOver->cd(i+1);
    gSysStat_lrCrossOver[i]->Draw("AP");
    gSysStat_lrCrossOver[i]->Fit("pol0", "0Q");
    Double_t avg =gSysStat_lrCrossOver[i]->GetFunction("pol0")->GetParameter(0);
    gSysStat_lrCrossOver[i]->SetTitle(Form("lrCrossOvereptance %0.2f", avg));
    gSysStat_lrCrossOver[i]->GetYaxis()->SetRangeUser(0, 0.1);
  }
  
}
