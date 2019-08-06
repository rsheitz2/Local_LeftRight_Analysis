#include "include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void FinalSetupLocal(TGraphErrors *g);

void DrawLegend(TGraphErrors *g);

void FinalSetupLocal(TGraph *g);

void DrawLegend(TGraph *g, TString whichName="SysStat");

void CalTotalSysErrors(Double_t *sysErr, vector<TGraph*> &vec,
		       Int_t nBins, TGraphErrors *gAN);

void CalTotalSysErrors(Double_t *sysErr, vector<TGraph*> &vec,
		       Int_t nBins);

void CalTotalSysOverStat(Double_t *eSys, TGraphErrors* gStat,
			 Double_t *sysOverStat);

void SetUpSys(TGraphAsymmErrors *g);

void sysWphysBinnedData(TString start=""){
  //Setup_______________
  /*const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  const Int_t nBins =4;//JPsi
  TString Mtype ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.87_3.38";
  TString fitMrange ="2.87_3.38";
  TString binRange ="29_34";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  Bool_t b_AnCorrected =true;
  //Setup_______________

  //Included Systematics
  Bool_t b_alphaAcc, b_lrCrossOver, b_faPull, b_additional, b_period;
  Bool_t b_impurity;
  if (Mtype=="HMDY"){
    b_alphaAcc =true; b_lrCrossOver =true; b_faPull =true;
    b_additional =true; b_period =false; b_impurity =false;
  }
  else if (Mtype=="LowM_AMDY"){
    b_alphaAcc =true; b_lrCrossOver =true; b_faPull =true;
      b_additional =true; b_period =true; b_impurity =true;//*/

      /*b_alphaAcc =false; b_lrCrossOver =false; b_faPull =false;
    b_additional =false; b_period =false; b_impurity =false;//*/
  }
  cout << "Systematic effects included" << endl;
  cout << "alphaAcc     " << b_alphaAcc << endl;
  cout << "lrCrossOver     " << b_lrCrossOver << endl;
  cout << "faPull     " << b_faPull << endl;
  cout << "additional     " << b_additional << endl;
  cout << "period     " << b_period << endl;
  cout << "impurity     " << b_impurity << endl;

  TString pathLR ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis";
  TString pathAN = pathLR+"/Macros/FinalAsym/Data/WAvg";
  TString pathSys = pathLR+"/Macros/Systematics";
  const Int_t nPhysBinned =6;
  TString physBinned[nPhysBinned] ={"xN", "xN", "xPi", "xF", "pT", "M"};
  TString xNames[nPhysBinned] =
    {"", "x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};
  //first xN used for integrated values
  
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
  
  //Aesthetics setup
  TCanvas* cAsym = new TCanvas();
  cAsym->SetLeftMargin(0.2); cAsym->Divide(nPhysBinned, 1, 0, 0);
  Double_t yMax =(process=="DY") ? 0.25 : 0.044;
  Double_t ysys =(process=="DY") ? -0.2 : -0.04;
  if (b_AnCorrected) { yMax *= TMath::Pi()/2.0; ysys *= TMath::Pi()/2.0;}
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphErrors *g_AN[nPhysBinned];
  TGraph *g_sys_stat[nPhysBinned];
  TGraphAsymmErrors *g_sys[nPhysBinned];
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString AsymName, alphaAccName, lrCrossOverName, faPullName, periodName;
    TString polDilName;
    Int_t nBinsName =nBins;
    if (phys == 0) {//used for integrated
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
      
      TString specficPath;
      specficPath ="FalseAsym/Data/WAvg";
      alphaAccName =
	Form("%s/%s/wAvg_true_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathSys.Data(), specficPath.Data(), Mtype.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),
	     physBinned[phys].Data(), nBinsName, hbins, production.Data(),
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

      specficPath =
	"PeriodCompatibility/Data/faPullDist2targWavg/";
      faPullName = (Mtype=="HMDY") 
	? Form("%s/%s/faPullDist2targWavg_true_%s_%s%s_%s_%i_%s_%s.root",
	     pathSys.Data(), specficPath.Data(), Mtype.Data(), process.Data(),
	     fitMrange.Data(), binRange.Data(), nBins, production.Data(),
	     additionalCuts.Data())
	: Form("%s/%s/falseApullDist_true_%s_%s%s_%s_%i_%s_%s.root",
	     pathSys.Data(), specficPath.Data(), Mtype.Data(), process.Data(),
	     fitMrange.Data(), binRange.Data(), nBins, production.Data(),
	     additionalCuts.Data());

      specficPath =
	"PeriodCompatibility/Data/pullDist/";
      periodName =
	Form("%s/%s/pullDist_true_%s_%s%s_%s_%i_%s_%s.root", pathSys.Data(),
	     specficPath.Data(), Mtype.Data(), process.Data(),
	     fitMrange.Data(), binRange.Data(), nBins, production.Data(),
	     additionalCuts.Data());

      specficPath =
	"Data";
      polDilName =
	Form("%s/%s/leftRight_byTarget_WAll_%s%s_1bins%s_150hbin_%s_%s.root",
	     pathLR.Data(), specficPath.Data(), Mtype.Data(), fitMrange.Data(),
	     binRange.Data(), production.Data(), additionalCuts.Data());
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
      alphaAccName = //cleanup
	Form("%s/%s/accSys4TargRatio_ten2.00_7.50_%s_%s2.90_3.30_%s%i_%ihbins\
.root",
	     pathSys.Data(), specficPath.Data(), Mtype.Data(),
	     process.Data(),
	     physBinned[phys].Data(), nBinsName, hbins);
    }

    //True Asymmetry
    TFile *f_AN = OpenFile(AsymName);

    cAsym->cd(phys+1);
    gPad->SetFrameLineWidth(2);
    g_AN[phys] =(TGraphErrors*)f_AN->Get("AN");

    if (b_AnCorrected){
      if (phys==0) cout << "\n!!!AN correction pi/2!!!\n" << endl;
      for (Int_t bi=0; bi<nBinsName; bi++) {
	g_AN[phys]->GetY()[bi] *= TMath::Pi()/2.0;
	g_AN[phys]->GetEY()[bi] *= TMath::Pi()/2.0;
      }
    }
    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    
    if (nBins==4){
      if (physBinned[phys] == "xN")
	g_AN[phys]->GetXaxis()->SetLimits(0.02, 0.17);
      else if (physBinned[phys] == "xPi")
	g_AN[phys]->GetXaxis()->SetLimits(0.11, 0.55);
      else if (physBinned[phys] == "xF")
	g_AN[phys]->GetXaxis()->SetLimits(0.0, 0.5);
      else if (physBinned[phys] == "pT")
	g_AN[phys]->GetXaxis()->SetLimits(0.3, 2.2);
      else if (physBinned[phys] == "M")
	g_AN[phys]->GetXaxis()->SetLimits(2.87, 3.4);
    }
    
    FinalSetupLocal(g_AN[phys]); DrawLegend(g_AN[phys]);
    SetTitleName(g_AN[phys], xNames[phys], "x");
    
    //Integrated
    if (phys == 0){
      g_AN[phys]->GetXaxis()->SetLimits(0.11, 0.245);
      SetTitleName(g_AN[phys], "A_{lr}", "y");
      g_AN[phys]->GetXaxis()->SetLabelSize(0.0);
      g_AN[phys]->GetXaxis()->SetTickSize(0.0);

      cout << "\nIntegrated Asym:  " << g_AN[phys]->GetY()[0] << endl;
      cout << "Integrated error:  " << g_AN[phys]->GetEY()[0] << endl;
      cout << "Integrated x value:  " << physBinned[phys] << " = "
	   << g_AN[phys]->GetX()[0] << "\n"<<endl;
    }
    DrawLine(g_AN[phys], 0.0);
    Double_t *xvals = g_AN[phys]->GetX();

    //Systematic error graphs
    ///////////////
    //Calculate and Draw systematic error bars
    //alphaAcc
    TFile *f_alphaAcc = (b_alphaAcc) ? OpenFile(alphaAccName) : NULL;
    TGraph *g_alphaAcc = (b_alphaAcc) ? (TGraph*)f_alphaAcc->Get("gSys") : NULL;

    //lrCrossOver
    TFile *f_lrCrossOver = (b_lrCrossOver) ? OpenFile(lrCrossOverName) : NULL;
    TGraph *g_lrCrossOver =
      (b_lrCrossOver) ? (TGraph*)f_lrCrossOver->Get("gSys") : NULL;

    //faPull
    TFile *f_faPull = (b_faPull) ? OpenFile(faPullName) : NULL; 
    TGraph *g_faPull = (b_faPull) ? (TGraph*)f_faPull->Get("gSys_Stat") : NULL;

    //Dilution && polarization errors
    Double_t dilError =0.05, polError =0.05;
    TFile *f_polDil = (b_additional) ? OpenFile(polDilName) : NULL;
    TGraph *g_Pol = (b_additional) ? (TGraph*)f_polDil->Get("xN_Pol") : NULL;
    TGraph *g_Dil = (b_additional) ? (TGraph*)f_polDil->Get("xN_Dil") : NULL;
    Double_t Pol_1Bin = (b_additional) ? g_Pol->GetY()[0] : 0.0;
    Double_t Dil_1Bin = (b_additional) ? g_Dil->GetY()[0] : 0.0;

    //Impurity
    Double_t impurity =0.0;
    if (b_impurity) impurity = 0.09529;

    //Changable settings (dilution, polarization, impurity)
    if (phys==0){
      cout << "dilution error = " << dilError << endl;
      cout << "polarization error = " << polError << endl;
      cout << "impurity error = " << impurity*100 << " %" << endl;
    }
    Double_t additionalErrors[nBinsName];
    Double_t impurityErrors[nBinsName];
    for (Int_t bi=0; bi<nBinsName; bi++) {
      additionalErrors[bi] = dilError*dilError + polError*polError;

      additionalErrors[bi] = TMath::Sqrt(additionalErrors[bi]);
      impurityErrors[bi] = impurity;
    }
    TGraph *g_Additional =
      (b_additional) ? new TGraph(nBinsName, xvals, additionalErrors) : NULL;
    TGraph *g_Impurity =
      (b_impurity) ? new TGraph(nBinsName, xvals, impurityErrors) : NULL;

    //period
    TFile *f_period = (b_period) ? OpenFile(periodName) : NULL; 
    TGraph *g_period = (b_period) ? (TGraph*)f_period->Get("gSys_Stat") : NULL;
    
    //Correct for sys/stat values
    for (Int_t bi=0; bi<nBinsName; bi++) {
      Double_t deltaAN = g_AN[phys]->GetEY()[bi];
      Double_t abs_AN = g_AN[phys]->GetY()[bi];
      if (abs_AN < 0)  abs_AN *= -1.0;
      if (b_faPull)
	g_faPull->GetY()[bi] *= deltaAN;

      if (b_additional)
	g_Additional->GetY()[bi] *= abs_AN;

      if (b_period)
	g_period->GetY()[bi] *= deltaAN;

      if (b_impurity)
	g_Impurity->GetY()[bi] *= deltaAN;
    }
    
    vector<TGraph*> vec_sys;
    vec_sys.push_back(g_alphaAcc);
    vec_sys.push_back(g_lrCrossOver);
    vec_sys.push_back(g_faPull);
    vec_sys.push_back(g_Additional);
    vec_sys.push_back(g_period);
    vec_sys.push_back(g_Impurity);

    if (phys == 0){//Integrated
      Double_t yvals[1] = {ysys}, eyb[1] = {0.0};
      Double_t exR[1] = {0.2}, exL[1] ={0.2};

      Double_t sysErr[1] = {0.0};
      CalTotalSysErrors(sysErr, vec_sys, 1);
    
      g_sys[phys] = new TGraphAsymmErrors(1, xvals, yvals, exL, exR,
					  eyb, sysErr);
      SetUpSys(g_sys[phys]);
    
      Double_t sysOverStat[1];
      CalTotalSysOverStat(sysErr, g_AN[phys], sysOverStat);
      g_sys_stat[phys] = new TGraph(1, xvals, sysOverStat);

      Double_t tot_int_error = sysErr[0]*sysErr[0]
	+ g_AN[phys]->GetEY()[0]*g_AN[phys]->GetEY()[0];
      cout << "\nIntegrated sys error:  " << sysErr[0] << endl;
      cout << "Integrated stat + sys error:  " << TMath::Sqrt(tot_int_error)
	   << "\n"<<endl;
    }
    else{//Not integrated
      Double_t yvals[nBins], eyb[nBins] = {0.0};
      Double_t exR[nBins], exL[nBins];
      for (Int_t i=0; i<nBins; i++) {
	yvals[i] = ysys;
	if (i==0) {exL[i] = 0.0; }
	else{ exL[i] = (xvals[i] - xvals[i-1])/2.0; }
      
	if (i==nBins-1){exR[i] = 0.0;}
	else{ exR[i] = (xvals[i+1] - xvals[i])/2.0; }
      }

      Double_t sysErr[nBins] = {0.0};
      CalTotalSysErrors(sysErr, vec_sys, nBins);
    
      g_sys[phys] = new TGraphAsymmErrors(nBins, xvals, yvals, exL, exR,
					  eyb, sysErr);
      SetUpSys(g_sys[phys]);
    
      Double_t sysOverStat[nBins];
      CalTotalSysOverStat(sysErr, g_AN[phys], sysOverStat);
      g_sys_stat[phys] = new TGraph(nBins, xvals, sysOverStat);
      SetUp(g_sys_stat[phys]);
    }
  }//phys binned loop

  TCanvas* cSysStat = new TCanvas(); 
  cSysStat->SetLeftMargin(0.2); cSysStat->Divide(nPhysBinned, 1, 0, 0);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cSysStat->cd(i+1);
    g_sys_stat[i]->GetYaxis()->SetRangeUser(0, 0.7);
    g_sys_stat[i]->Fit("pol0", "0Q");
    Double_t avg =g_sys_stat[i]->GetFunction("pol0")->GetParameter(0);
    g_sys_stat[i]->SetTitle(Form("%0.3f", avg));
    g_sys_stat[i]->Draw("AP");
    
    FinalSetupLocal(g_sys_stat[i]); DrawLegend(g_sys_stat[i]);
    SetTitleName(g_sys_stat[i], xNames[i], "x");
    if (i==0){//intgrated
      g_sys_stat[i]->GetXaxis()->SetLimits(0.11, 0.245);
      g_sys_stat[i]->GetXaxis()->SetLabelSize(0.0);
      g_sys_stat[i]->GetXaxis()->SetTickSize(0.0);
    }
  }

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/PhysBinnedData";
  TString fOutput;
  if (whichFit=="true"){
    fOutput = Form("%s/physBinnedData_%s_%s_%s%s_%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), Mtype.Data(),
		   process.Data(), lrMrange.Data(), nBins,
		   hbins, production.Data(), additionalCuts.Data());
  }
  else{
    fOutput = Form("%s/physBinnedData_%s%s_%s_%s%s_%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		   Mtype.Data(), process.Data(), lrMrange.Data(),
		   nBins, hbins, production.Data(), 
		   additionalCuts.Data());
  }

  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    for (Int_t i=0; i<nPhysBinned; i++) {
      g_AN[i]->Write(Form("AN_%s", physBinned[i].Data() ));
    }

    cAsym->Write();
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Asymmetry data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " <<  whichFit << endl;
  }
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
}


void FinalSetupLocal(TGraphErrors *g){
  SetUp(g);
  FinalSetup(g);

  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(20);
  g->SetTitle("");
}


void DrawLegend(TGraphErrors *g){
  g->SetName("AN_phys");
  Double_t sigma;
  Double_t avg = WeightedAvgAndError(g, &sigma);

  TLegend *leg = new TLegend(0.1,0.9,0.7,0.99);
  leg->AddEntry("AN_phys", Form("#bar{A} = %0.2f #pm %0.2f", avg, sigma),
		"p");
  leg->SetBorderSize(0); leg->SetTextFont(132); leg->SetTextSize(0.08);
  leg->Draw("same");
}


void FinalSetupLocal(TGraph *g){
  SetUp(g);
  FinalSetup(g);

  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(20);
  g->SetTitle("");
}


void DrawLegend(TGraph *g, TString whichName="SysStat"){
  g->SetName("gr");
  Double_t avg = Avg(g);

  TLegend *leg = new TLegend(0.1,0.9,0.7,0.99);
  leg->SetBorderSize(0); leg->SetTextFont(132); leg->SetTextSize(0.08);
  
  if(whichName=="SysStat"){
    leg->AddEntry("gr",
		  Form("#bar{#sigma_{systemaic}/#sigma_{statistical}} = %0.2f",
		       avg), "p");
  }
  
  leg->Draw("same");
}


void CalTotalSysErrors(Double_t *sysErr, vector<TGraph*> &vec,
		       Int_t nBins, TGraphErrors *gAN){

  cout << "Function Needs to be fixed" << endl;
  exit(EXIT_FAILURE);
  
  Double_t *e_yStat = gAN->GetEY();
  for (vector<TGraph*>::iterator it=vec.begin(); it!=vec.end(); it++){
    Double_t *yvals = (*it)->GetY();
    for (Int_t bi=0; bi<nBins; bi++) {
      if ((*it)->GetN() < nBins) sysErr[bi] += yvals[0]*yvals[0];
      else sysErr[bi] += yvals[bi]*yvals[bi];
    }
  }

  for (Int_t bi=0; bi<nBins; bi++) {
    sysErr[bi] += 0.28*0.28*e_yStat[bi]*e_yStat[bi];
    sysErr[bi] = TMath::Sqrt(sysErr[bi]);
  }
}


void CalTotalSysErrors(Double_t *sysErr, vector<TGraph*> &vec,
		       Int_t nBins){
  for (vector<TGraph*>::iterator it=vec.begin(); it!=vec.end(); it++){
    if (!(*it)) continue;
    Double_t *yvals = (*it)->GetY();
    for (Int_t bi=0; bi<nBins; bi++) {
      if ((*it)->GetN() < nBins) sysErr[bi] += yvals[0]*yvals[0];
      else sysErr[bi] += yvals[bi]*yvals[bi];
    }
  }

  for (Int_t bi=0; bi<nBins; bi++) {
    sysErr[bi] = TMath::Sqrt(sysErr[bi]);
  }
}


void CalTotalSysOverStat(Double_t *eSys, TGraphErrors* gStat,
			 Double_t *sysOverStat){
  Double_t *eStat = gStat->GetEY();
  for (Int_t i=0; i<gStat->GetN(); i++) {
    sysOverStat[i] = eSys[i]/eStat[i];
  }
}


void SetUpSys(TGraphAsymmErrors *g){
  SetUp(g);
  g->Draw("2same");
  g->SetFillColor(kRed);
  g->SetFillStyle(3002);
}
