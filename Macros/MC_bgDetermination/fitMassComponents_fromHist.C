#include "Include/helperFunctions.h"

const Int_t nFiles = 5;
TH1D *hist[nFiles]; //{JPsi, psi, OC, AMDY, Real}

Double_t FitMCs(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );

  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}


void fitMassComponents_fromHist(TString start=""){
  //Setup__________
  TString whichMC ="Charles"; //"Yu", "Charles"
  const Double_t Mmin=2.0, Mmax=8.5; //Mass cut
  TString whichRD ="slot1"; //"t3", "slot1"
  Int_t binDivisor =1;

  Bool_t toWrite =false;
  //Setup__________

  TString MassCut = Form("Mmumu>%0.2f&&Mmumu<%0.2f", Mmin, Mmax);

  if (start=="") {
    cout << "Script determines the amount of each MC component in real data";
    cout << "\nThis is done fitting MC distributions as a function of Mass to";
    cout << " to get the contribution from each physical process" << endl;
    cout << "\nCurrent Settings:" << endl;
    cout << "Mass ranged considered:       " << MassCut << endl;
    cout << "Which MC considered:          " << whichMC << endl;
    cout << "Which reald data considered:  " << whichRD << endl;
    cout << "binDivisor:                   " << binDivisor << endl;
    cout << "\nTo Write:                   " << toWrite << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'fitMassComponents(1)\'" << endl;
    exit(EXIT_FAILURE);
  }

  //Aesthetic setup
  gStyle->SetOptFit(1111);
  Int_t icolor[nFiles] = {6,3,9,4,1};//{JPsi, psi, OC, AMDY, RealData}

  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/MC_bgDetermination/";
  TString pathMC;
  TString dataPathNames[nFiles];
  if (whichMC=="Charles"){
    pathMC +=path+"Data/MassDist/";
    dataPathNames[0] = pathMC+"M_Charles_Jpsi.root";
    dataPathNames[1] = pathMC+"M_Charles_Psi.root";
    dataPathNames[2] = pathMC+"M_Charles_OC.root";
    dataPathNames[3] = pathMC+"M_Charles_AMDY.root";
    dataPathNames[4] = pathMC+Form("M_%s.root", whichRD.Data());
  }
  else {
    cout << "MC specified does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  //Open data files && get MC mass histograms
  TFile *fJPsi = TFile::Open(dataPathNames[0]);
  TFile *fpsi = TFile::Open(dataPathNames[1]);
  TFile *fOC = TFile::Open(dataPathNames[2]);
  TFile *fAMDY = TFile::Open(dataPathNames[3]);
  TFile *fReal = TFile::Open(dataPathNames[4]);
  TFile *Files[nFiles] = {fJPsi, fpsi, fOC, fAMDY, fReal};
  TString type[nFiles] = {"Jpsi", "Psi", "OC", "AMDY", whichRD};
  for (Int_t i=0; i<nFiles; i++) {
    if (!Files[i]) {
      cout << "File:   " << i << "   does not exist " << endl;
      exit(EXIT_FAILURE); }

    hist[i] = (TH1D*)Files[i]->Get(Form("h_%s", type[i].Data()));
    hist[i]->Rebin(binDivisor);
    Setup(hist[i]);
  }

  //Clone MC histograms to draw on added up graph
	TCanvas *c1 = new TCanvas();
  hist[nFiles-1]->Draw();
	hist[nFiles-1]->GetXaxis()->SetRangeUser(Mmin, Mmax);
  TH1D *histResult[nFiles];//{JPsi, psi, OC, AMDY, Sum}
  Int_t nBins=hist[0]->GetXaxis()->GetNbins();
  histResult[nFiles-1] = new TH1D("hResult_Sum", "hResult_Sum",nBins,2.0,8.5);
  for (Int_t i=0; i<nFiles-1; i++) {
    hist[i]->Scale(1/(hist[i]->Integral()) );
    hist[i]->SetLineColor(icolor[i]);

    histResult[i] = (TH1D*)hist[i]->Clone();
  }

  //Do fit
  TF1 *fitFunc = new TF1("fitFunc", FitMCs, Mmin, Mmax, nFiles-1);
  hist[nFiles-1]->Fit("fitFunc", "RWL", "", Mmin, Mmax);
  Double_t chi =fitFunc->GetChisquare();
  Int_t ndf =fitFunc->GetNDF();
  cout << "\n\nChi2:  " << chi << "  / ndf:  " << ndf << endl;
  cout << "reduced chi2:  " << chi/(1.0*ndf) << endl;

  c1->cd(); gPad->SetLogy();
  for (Int_t i=0; i<nFiles-1; i++) {
    histResult[i]->Scale(fitFunc->GetParameter(i) );
    histResult[i]->Draw("same");

    histResult[nFiles-1]->Add(histResult[i]);
  }
  
  //Make parameters tgraph
  Double_t *pars =fitFunc->GetParameters();
  const Double_t *e_pars =fitFunc->GetParErrors();
  Double_t xvals[nFiles-1], ex[nFiles-1] = {0.0};
  for (Int_t i=0; i<nFiles-1; i++) { xvals[i] = i+1; }
  TGraphErrors* gPar = new TGraphErrors(nFiles-1, xvals, pars, ex, e_pars);
  Setup(gPar);
  TCanvas* cPars = new TCanvas();
  gPar->Draw("AP"); gPar->SetTitle("Components Coefficents");

  //Make Ratio plot
  TH1D *hRatio =(TH1D*)histResult[nFiles-1]->Clone();
  hRatio->Divide(hist[nFiles-1]);
  TCanvas* cRatio = new TCanvas();
  hRatio->Draw();

  //Write output file
  TString outName =Form("fitComponents_%s_%s_%0.2f_%0.2f.root",
			whichMC.Data(), whichRD.Data(), Mmin, Mmax);
  if (toWrite){
    TFile* fOut = new TFile("Data/FitMassComponents/"+outName, "RECREATE");
    gPar->Write("gPar");
    for (Int_t i=0; i<nFiles-1; i++) {
      hist[i]->Write();
    }

    fOut->Close();
  }

  //Final Output
  cout << "\nSettings !!!!" << endl;
  cout << "Mass Range is		:  " << MassCut << endl;
  cout << "Which MC considered:          " << whichMC << endl;
  cout << "Which reald data considered:  " << whichRD << endl;
  cout << "binDivisor:                   " << binDivisor << endl;
  if (toWrite) { cout << outName << "   was written" << endl; }
  else { cout << outName << "    was NOT written" << endl; }
  cout << " " << endl;
}
