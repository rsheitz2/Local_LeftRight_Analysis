#include "include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

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
  leg->AddEntry("AN_phys", Form("#bar{A}_{N} = %0.2f #pm %0.2f", avg, sigma),
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


void sysWphysBinnedData(TString start=""){
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

  /*const Int_t nBins =5;//JPsi
    TString Mtype ="LowM_AMDY";
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
  const Int_t nPhysBinned =6;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT", "M", "xN"};
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
  
  //Aesthetics setup
  TCanvas* cAsym = new TCanvas(); cAsym->Divide(nPhysBinned, 1, 0, 0.01);
  Double_t yMax =(process=="DY") ? 0.25 : 0.1;
  Double_t ysys =(process=="DY") ? -0.15 : -0.075;
  
  //Get Data file/Get graphs and plot
  TString physBinnedNames ="", fitNames="";
  TGraphErrors *g_AN[nPhysBinned];
  TGraph *g_sys_stat[nPhysBinned];
  TGraphAsymmErrors *g_sys[nPhysBinned];
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
	"FalseAsym/Data/WAvg";
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
    g_AN[phys]->Draw("AP"); g_AN[phys]->GetYaxis()->SetRangeUser(-yMax, yMax);
    FinalSetupLocal(g_AN[phys]); DrawLegend(g_AN[phys]);
    DrawLine(g_AN[phys], 0.0); 

    //Integrated
    if (phys == nPhysBinned-1){ g_AN[phys]->GetXaxis()->SetLimits(0.11, 0.245);}

    //Systematic error graphs
    ///////////////
    //Calculate and Draw systematic error bars
    //alphaAcc
    TFile *f_alphaAcc = OpenFile(alphaAccName);
    TGraph *g_alphaAcc = (TGraph*)f_alphaAcc->Get("gSys");

    //lrCrossOver
    TFile *f_lrCrossOver = OpenFile(lrCrossOverName);
    TGraph *g_lrCrossOver = (TGraph*)f_lrCrossOver->Get("gSys");
    
    vector<TGraph*> vec_sys;
    vec_sys.push_back(g_alphaAcc);
    vec_sys.push_back(g_lrCrossOver);

    Double_t *xvals = g_AN[phys]->GetX();
    if (phys == nPhysBinned-1){//Integrated
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

  TCanvas* cSysStat = new TCanvas(); cSysStat->Divide(nPhysBinned, 1, 0, 0.01);
  for (Int_t i=0; i<nPhysBinned; i++) {
    cSysStat->cd(i+1);
    g_sys_stat[i]->GetYaxis()->SetRangeUser(0, 0.7);
    g_sys_stat[i]->Fit("pol0", "0Q");
    Double_t avg =g_sys_stat[i]->GetFunction("pol0")->GetParameter(0);
    g_sys_stat[i]->SetTitle(Form("%0.3f", avg));
    g_sys_stat[i]->Draw("AP");
    FinalSetupLocal(g_sys_stat[i]); DrawLegend(g_sys_stat[i]);
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

