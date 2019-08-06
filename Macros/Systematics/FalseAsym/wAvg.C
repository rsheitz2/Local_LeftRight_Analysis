#include "include/helperFunctions.h"

void AddWAvg(TGraphErrors *g, Double_t *yWavg, Double_t *eyWavg){

  if (!g){
    cout << "TGraph does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t *yval = g->GetY();
  Double_t *e_yval = g->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    yWavg[i] += yval[i]/(e_yval[i]*e_yval[i]);
    eyWavg[i] += 1.0/(e_yval[i]*e_yval[i]);
  }
}


void AddErrorWAvg(TGraph *g, Double_t *yWavg){

  if (!g){
    cout << "TGraph does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t *yval = g->GetY();
  for (Int_t i=0; i<g->GetN(); i++) {
    yWavg[i] += 1.0/(yval[i]*yval[i]);
  }
}


void TableFormat(TString p, TString physBinned, TGraphErrors *g);

void wAvg(){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  /*const Int_t nBins =4;//JPsi
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString physBinned ="xN";//xN, xPi, xF, pT, M
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  //Setup_______________

  //cleanup
  cout << "\n\n!!!!!\nDo something about period canvases....\n!!!!!\n" << endl;
  
  //Basic Setup
  TString pathFA="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data";
  const Int_t nPer =9;
  TString period[nPer] ={"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14",
			 "W15"};
  Int_t icolor[nPer] = {kBlue+2, kRed+2, kGreen+2, kMagenta+2, kCyan+2,
			kBlue, kRed, kGreen, kMagenta};

  //TGraph axis's
  Double_t yFAsubper[nBins] ={0.0}, e_yFAsubper[nBins] ={0.0};
  Double_t yFApol[nBins] ={0.0}, e_yFApol[nBins] ={0.0};
  Double_t yAcc[nBins] ={0.0}, e_yAcc[nBins] ={0.0};
  Double_t ySys[nBins] ={0.0};
  Double_t ex[nBins] ={0.0};
  Double_t *xvals;

  //Period specific
  TCanvas* cPerFApol = new TCanvas();
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9); legend->SetNColumns(3);

  TCanvas* cPer = new TCanvas(); cPer->Divide(3, 3, 0, 0);

  //Get Data and add for wAvg
  for (Int_t p=0; p<nPer; p++) {
    TString perFAsubper, perAcc, perSysErr;
    if (whichFit=="true"){
      perFAsubper =
	Form("%s/TargFlip/falseGeoMean4Targ_%s_%s_%s_%s%s_%s%s%i_%s_%s.root",
	     pathFA.Data(), whichFit.Data(), period[p].Data(), fitMrangeType.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),physBinned.Data(),
	     nBins, production.Data(), additionalCuts.Data());
      
      perAcc =
	Form("%s/acceptanceFourTargRatio/\
acceptanceFourTargRatio_%s_%s_%s_%s%s_%s%s%i_%s_%s.root",
	     pathFA.Data(), whichFit.Data(), period[p].Data(), fitMrangeType.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),physBinned.Data(),
	     nBins, production.Data(), additionalCuts.Data());

      perSysErr =
	Form("%s/acceptanceFourTargRatio/SystematicError/\
accSys4TargRatio_%s_%s_%s_%s%s_%s%s%i_%s_%s.root",
	     pathFA.Data(), whichFit.Data(), period[p].Data(), fitMrangeType.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),physBinned.Data(),
	     nBins, production.Data(), additionalCuts.Data());
    }
    else{
      cout << "Only works for true fit for now" << endl;
    exit(EXIT_FAILURE);
    
      perFAsubper =
	Form("%s/falseGeoMean4Targ_%s%s_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathFA.Data(), whichFit.Data(), fitMrange.Data(), period[p].Data(),
	     fitMrangeType.Data(), process.Data(), lrMrange.Data(), binRange.Data(),
	     physBinned.Data(), nBins, hbins, production.Data(),
	     additionalCuts.Data());
    }
    TFile *f_FA = OpenFile(perFAsubper); TFile *f_acc = OpenFile(perAcc);
    TFile *f_sys = OpenFile(perSysErr);
    TGraphErrors *g_FAsubper = (TGraphErrors*) f_FA->Get("falseAN_subper");
    TGraphErrors *g_FApol = (TGraphErrors*) f_FA->Get("falseAN_pol");
    TGraphErrors *g_upS = (TGraphErrors*) f_FA->Get("falseAN_2Targ_upS");
    TGraphErrors *g_downS = (TGraphErrors*) f_FA->Get("falseAN_2Targ_downS");
    TGraphErrors *g_Acc = (TGraphErrors*) f_acc->Get("alpha_subper");
    TGraphErrors *g_sys = (TGraphErrors*) f_sys->Get("gSys");
    
    AddWAvg(g_FAsubper, yFAsubper, e_yFAsubper);
    AddWAvg(g_FApol, yFApol, e_yFApol);
    AddWAvg(g_Acc, yAcc, e_yAcc);
    AddErrorWAvg(g_sys, ySys);
    
    if (p==0) xvals = g_FAsubper->GetX();

    //period only
    if (p==0){
      cPerFApol->cd();
      g_FApol->Draw("AP");
      DrawLine(g_FApol, 0.0);
    }
    else {
      cPerFApol->cd();
      g_FApol->Draw("Psame");
    }
    OffSet(g_FApol, p*0.001);
    g_FApol->SetMarkerStyle(p+20); g_FApol->SetMarkerColor(icolor[p]);
    legend->AddEntry(g_FApol, Form("%s", period[p].Data()), "p");

    if (p==nPer-1) {
      cPerFApol->cd();
      legend->Draw("same");
    }

    cPer->cd(p+1);//changed for thesis plot
    //g_FApol->Draw("AP"); DrawLine(g_FApol, 0.0);
    g_upS->Draw("AP"); DrawLine(g_upS, 0.0);
    g_upS->GetYaxis()->SetRangeUser(-0.9, 0.9);
    g_upS->GetYaxis()->SetTitle("A_{N,false}");
    g_upS->SetTitle(Form("%s", period[p].Data()));
    OffSet(g_downS, 0.005);
    g_downS->Draw("Psame"); 
  }

  //Final wAvg calculation
  Double_t ySys_stat[nBins] ={0.0};
  for (Int_t i=0; i<nBins; i++) {
    yFAsubper[i] /= e_yFAsubper[i];
    e_yFAsubper[i] = TMath::Sqrt(1.0/e_yFAsubper[i]);

    yFApol[i] /= e_yFApol[i];
    e_yFApol[i] = TMath::Sqrt(1.0/e_yFApol[i]);

    yAcc[i] /= e_yAcc[i];
    e_yAcc[i] = TMath::Sqrt(1.0/e_yAcc[i]);

    ySys[i] = 1.0/TMath::Sqrt(ySys[i]);

    ySys_stat[i] = ySys[i]/(e_yFAsubper[i]);
  }

  TGraphErrors* g_WAvg_FAsubper =
    new TGraphErrors(nBins, xvals, yFAsubper, ex, e_yFAsubper);

  TGraphErrors* g_WAvg_FApol =
    new TGraphErrors(nBins, xvals, yFApol, ex, e_yFApol);

  TGraphErrors* g_WAvg_Acc =
    new TGraphErrors(nBins, xvals, yAcc, ex, e_yAcc);

  TGraph* g_WAvg_Sys = new TGraph(nBins, xvals, ySys);
  TGraph* g_WAvg_Sys_Stat = new TGraph(nBins, xvals, ySys_stat);
  
  SetUp(g_WAvg_FAsubper); SetUp(g_WAvg_Acc); SetUp(g_WAvg_FApol);
  SetUp(g_WAvg_Sys); SetUp(g_WAvg_Sys_Stat);

  TableFormat("False Asymmetry", physBinned, g_WAvg_FAsubper);
  TableFormat("Acceptance Ratio", physBinned, g_WAvg_Acc);

  //Draw data
  TCanvas* cFA = new TCanvas();
  g_WAvg_FAsubper->Draw("AP");
  g_WAvg_FAsubper->SetTitle("falseAN_subper/pol(blue)");
  if (fitMrangeType=="HMDY"){ g_WAvg_FAsubper->GetYaxis()->SetRangeUser(-0.25, 0.25); }
  else { g_WAvg_FAsubper->GetYaxis()->SetRangeUser(-0.08, 0.08); }
  DrawLine(g_WAvg_FAsubper, 0.0);
  g_WAvg_FApol->Draw("Psame"); g_WAvg_FApol->SetMarkerColor(kBlue);

  TCanvas* cAcc = new TCanvas();
  g_WAvg_Acc->Draw("AP"); g_WAvg_Acc->SetTitle("acc_subper");
  if (fitMrangeType=="HMDY"){ g_WAvg_Acc->GetYaxis()->SetRangeUser(1-0.1, 1+0.1); }
  else { g_WAvg_Acc->GetYaxis()->SetRangeUser(1-0.08, 1+0.08); }
  DrawLine(g_WAvg_Acc, 1.0);

  TCanvas* cSys = new TCanvas();
  g_WAvg_Sys->Draw("AP"); g_WAvg_Sys->SetTitle("gSys");

  TCanvas* cSysStat = new TCanvas();
  g_WAvg_Sys_Stat->Draw("AP"); g_WAvg_Sys_Stat->SetTitle("gSys_Stat");
  g_WAvg_Sys_Stat->GetYaxis()->SetRangeUser(0, 1);
  
  //Write output/final settings
  TString fOutput;
  if (whichFit=="true"){
    fOutput = Form("%s/WAvg/wAvg_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
		   pathFA.Data(), whichFit.Data(), fitMrangeType.Data(),
		   process.Data(), lrMrange.Data(), binRange.Data(),
		   physBinned.Data(), nBins, hbins, production.Data(),
		   additionalCuts.Data());
  }
  else{
    fOutput = Form("%s/WAvg/wAvg_%s%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
		   pathFA.Data(), whichFit.Data(), fitMrange.Data(),
		   fitMrangeType.Data(), process.Data(), lrMrange.Data(),
		   physBinned.Data(), binRange.Data(), nBins, hbins,
		   production.Data(), additionalCuts.Data());
  }
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    g_WAvg_FAsubper->Write("falseAN_subper");
    g_WAvg_FApol->Write("falseAN_pol");
    g_WAvg_Acc->Write("acc_subper");
    g_WAvg_Sys->Write("gSys");
    g_WAvg_Sys_Stat->Write("gSys_Stat");
  }

  cout << " " << endl;
  cout << "Settings________" << endl;
  cout << "Data coming from:            " << pathFA << endl;
  cout << "physBinned nBins times:     " << nBins << endl;
  cout << "Mass type considered:   " << fitMrangeType << endl;
  cout << "Binned in which DY physics:  " << physBinned << endl;
  cout << "AN physical process:        " << process << endl;
  if (whichFit == "true"){
      cout << "L/R mass range:     " << lrMrange << endl;
    }
    else{
      cout << "LR integral mass range:     " << lrMrange << endl;
      cout << "Fit mass range:     " << fitMrange << endl;
    }
  cout << "Which fit considered:       " << whichFit << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;
}


void TableFormat(TString p, TString physBinned, TGraphErrors *g){
  Double_t *yval = g->GetY();
  Double_t *e_yval = g->GetEY();
  for (Int_t i=0; i<g->GetN(); i++) {
    cout << p << ", " << physBinned << ", "
	 << yval[i] << ", " << e_yval[i] << endl;
  }
}
