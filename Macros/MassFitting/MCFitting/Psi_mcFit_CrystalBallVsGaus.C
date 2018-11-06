#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/MassFitting/include/helperFunctions.h"


Double_t CrystalBall(Double_t *x, Double_t *par){
  //f(x; alpha, n, xBar, sigma, Amp)
  //   Amp = par[0], xBar = par[1], sigma = par[2]
  //   alpha = par[3], n = par[4]
  Double_t A =TMath::Power(par[4]/par[3],par[4])*TMath::Exp(-par[3]*par[3]/2.0);
  Double_t B =par[4]/par[3] - par[3];
  Double_t C =(par[4]/par[3])*(1.0/(par[4]-1.0))*TMath::Exp(-par[3]*par[3]/2.0);
  Double_t D =TMath::Sqrt( TMath::Pi()/2.0 )*(1+TMath::Erf( par[3]/TMath::Sqrt(2) ) );

  Double_t Norm = 1.0/(par[2]*(C+D) );

  Double_t arg = (x[0] - par[1])/par[2];
  if (arg > -par[3] ) return par[0]*Norm*TMath::Exp(-0.5*arg*arg);
  else return par[0]*Norm*A*TMath::Power((B - arg), -par[4]);
}


void ParaCrystalBall(TH1D *h, TF1 *fitFunc, Double_t mass){
  Double_t Amp = h->GetBinContent(h->FindBin(mass) );
  if (Amp == 0) Amp = 200.0;
  
  fitFunc->SetParameter(0, Amp);
  fitFunc->SetParameter(1, mass);
  fitFunc->SetParameter(2, 0.16);//Width
  fitFunc->SetParameter(3, 1.7);
  //fitFunc->SetParameter(4, 6.5);
  fitFunc->SetParameter(4, 3.5);

  Double_t factor =10.0;
  fitFunc->SetParLimits(0, Amp/factor, Amp*factor/2.0);
  fitFunc->SetParLimits(1, mass-0.4, mass+0.4);
  //fitFunc->SetParLimits(2, 0, 0.25);//Width
  fitFunc->SetParLimits(2, 0, 0.4);//Width
  fitFunc->SetParLimits(3, 0.5, 3.0);
  fitFunc->SetParLimits(4, 1.0, 10.0);
}


Double_t DoFit(TH1D *h, TF1* f, TString fitName, Double_t Mmin, Double_t Mmax){
  //Fit with unified options with either
  //        input function or predefined function fitName

  SetUp(h);
  h->SetLineWidth(3);
  TString fitOptions ="RSQ+";
  //TString fitOptions ="RLSQ+";
  //TString fitOptions ="RWLSQ+";
  TFitResultPtr status;
  if (f) //Fit with input function
    status = h->Fit(f, fitOptions, "", Mmin, Mmax);
  else if (fitName != "NULL")
    status = h->Fit(fitName, fitOptions, "", Mmin, Mmax);
  else {
    cout << "Invalid fit input" << endl;
    exit(EXIT_FAILURE);
  }

  if (status->Status()) {//Check for fit fails
    cout << "\n\n" << h->GetTitle() << "  Fit failed:  \n" << endl;
    //exit(EXIT_FAILURE);
  }

  Int_t ndf = status->Ndf();
  Double_t Chi2 = status->Chi2();

  return Chi2/ndf;
}


void DoFit(TH1D *h, TF1* f, TH1D* hChi2, Double_t Mmin, Double_t Mmax){
  Double_t redChi2 = DoFit(h, f, "NULL", Mmin, Mmax);
  hChi2->Fill(redChi2);
}


TF1* DoFit(TH1D *h, TString fitName, TH1D* hChi2, Double_t Mmin, Double_t Mmax){
  Double_t redChi2 = DoFit(h, NULL, fitName, Mmin, Mmax);
  hChi2->Fill(redChi2);

  TF1 *fgaus = (TF1*) h->GetFunction(fitName);
  fgaus->SetLineColor(kBlue); fgaus->Draw("same");

  return fgaus;
}


void FillFromFit(TF1* f, TH1D* h, Double_t Mmin, Double_t Mmax){
  Int_t minBin = h->FindBin(Mmin);
  Int_t maxBin = h->FindBin(Mmax);
  
  for (Int_t bi=1; bi<h->GetNbinsX()+1; bi++) {
    if (bi < minBin){ h->SetBinContent(bi, 0.0); }
    else if (bi > maxBin){ h->SetBinContent(bi, 0.0); }
    else{
      Double_t xval =h->GetBinCenter(bi);
      Double_t fitVal = f->Eval(xval);
      h->SetBinContent(bi, fitVal);
    }
  }
  h->Sumw2();
}


void Psi_mcFit_CrystalBallVsGaus(TString start=""){
  //Setup_____
  const Int_t nBins =5;
  TString physBinned ="pT";
  Double_t Mmin =3.0;
  Double_t Mmax =4.0;
  TString whichPsi="Psi"; //"JPsi", "Psi";
  //Setup_____

  if (start==""){
    cout << "Macro compares Monte-Carlo data with a Crystal Ball fit";
    cout << " and with Guassian a fit" << endl;
    cout << "\nCharles MC from Blue Waters is being used" << endl;
    exit(EXIT_FAILURE);
  }

  //Decide starting mass
  Double_t M_Psi = (whichPsi=="JPsi") ? 3.12 : 3.60; //"3.12"=JPsi, "3.6"=Psi

  //Get Data file
  TString lrPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/";
  TString namePsi =
    Form("Data/leftRight_byTarget_Charles_%s1.00_8.50_5bins25_43_150hbin.root",
	 whichPsi.Data());
  TFile *fPsi = TFile::Open(lrPath+namePsi);

  if (!fPsi){
    cout << "File does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  //Canvas Setups
  TCanvas* cUpS = new TCanvas(); cUpS->Divide(2, 3);
  TCanvas* cDS = new TCanvas(); cDS->Divide(2, 3);
  TCanvas* cRUpS_cry = new TCanvas(); cRUpS_cry->Divide(2, 3);
  TCanvas* cRUpS_gaus = new TCanvas(); cRUpS_gaus->Divide(2, 3);
  TCanvas* cRDS_cry = new TCanvas(); cRDS_cry->Divide(2, 3);
  TCanvas* cRDS_gaus = new TCanvas(); cRDS_gaus->Divide(2, 3);
  
  //Distribution (h), Fit (hF), and Ratio (r) histograms
  TH1D *h_UpS[nBins], *h_DS[nBins];
  TH1D *hF_UpS_cry[nBins], *hF_DS_cry[nBins];
  TH1D *hF_UpS_gaus[nBins], *hF_DS_gaus[nBins];
  TRatioPlot *r_UpS_cry[nBins], *r_DS_cry[nBins];
  TRatioPlot *r_UpS_gaus[nBins], *r_DS_gaus[nBins];

  //reduced chi2 distributions
  TH1D *hChi2_gaus = new TH1D("chi2_gaus", "chi2_gaus", 100, 100, 400);
  TH1D *hChi2_cry = new TH1D("chi2_cry", "chi2_cry", 40, 0, 150);
  SetUp(hChi2_gaus); SetUp(hChi2_cry);
  hChi2_gaus->SetLineWidth(3); hChi2_cry->SetLineWidth(3);

  //Integral histograms
  TH1D* hInt_cry = new TH1D("hInt_cry", "hInt_cry", 20, 0.96, 1.01);
  TH1D* hInt_gaus = new TH1D("hInt_gaus", "hInt_gaus", 20, 0.96, 1.01);
  SetUp(hInt_cry); SetUp(hInt_gaus);
  hInt_cry->SetLineWidth(3); hInt_gaus->SetLineWidth(3);

  for (Int_t i=0; i<nBins; i++) {
    //////////
    //UpS
    //distribution
    cUpS->cd(i+1); gPad->SetLogy();
    h_UpS[i] =
      (TH1D*)fPsi->Get(Form("MuMu_left_upstream_up_%s%i", physBinned.Data(),i));
    h_UpS[i]->Sumw2(); 
    h_UpS[i]->Draw();

    //fit crystal ball
    TF1 *f_UpS_cry = new TF1("crystal", CrystalBall, Mmin, Mmax, 5);
    ParaCrystalBall(h_UpS[i], f_UpS_cry, M_Psi);
    DoFit(h_UpS[i], f_UpS_cry, hChi2_cry, Mmin, Mmax);
    Double_t *paras = f_UpS_cry->GetParameters();
    cout << "UpS crystal ball Parameter alpha:  "<< paras[3]
	 <<"    Paramter n:  " <<paras[4] << endl;

    //fit gaus 
    TF1* fgaus_UpS =DoFit(h_UpS[i], "gaus", hChi2_gaus, Mmin, Mmax);
    cUpS->Update();

    //Basic histogram needed parameters
    Int_t nhbin = h_UpS[i]->GetNbinsX();
    Double_t xmin = h_UpS[i]->GetXaxis()->GetXmin();
    Double_t xmax = h_UpS[i]->GetXaxis()->GetXmax();

    //crystal ball ratios
    cRUpS_cry->cd(i+1); gPad->SetLogy();
    hF_UpS_cry[i] = new TH1D(Form("hF_UpS_cry_%i", i), Form("hF_UpS_cry_%i", i),
			   nhbin, xmin, xmax);
    FillFromFit(f_UpS_cry, hF_UpS_cry[i], Mmin, Mmax);
    r_UpS_cry[i] = new TRatioPlot(h_UpS[i], hF_UpS_cry[i]);
    SetUp(r_UpS_cry[i]); cRUpS_cry->Update();

    //gaus ratio
    cRUpS_gaus->cd(i+1); gPad->SetLogy();
    hF_UpS_gaus[i] = new TH1D(Form("hR_UpS_gaus_%i", i),
				Form("hR_UpS_gaus_%i", i),
				nhbin, xmin, xmax);
    FillFromFit(fgaus_UpS, hF_UpS_gaus[i], Mmin, Mmax);
    r_UpS_gaus[i] = new TRatioPlot(h_UpS[i], hF_UpS_gaus[i]);
    SetUp(r_UpS_gaus[i]); cRUpS_gaus->Update();

    //Intgrals
    Int_t minBin = h_UpS[i]->FindBin(Mmin);
    Int_t maxBin = h_UpS[i]->FindBin(Mmax);
    Double_t er_UpS, er_UpS_cry, er_UpS_gaus;//cleanup
    Double_t Int_UpS =h_UpS[i]->IntegralAndError(minBin, maxBin, er_UpS);
    Double_t Int_UpS_cry =
      hF_UpS_cry[i]->IntegralAndError(minBin, maxBin, er_UpS_cry);
    Double_t Int_UpS_gaus =
      hF_UpS_gaus[i]->IntegralAndError(minBin, maxBin, er_UpS_gaus);
    hInt_cry->Fill(Int_UpS_cry/Int_UpS);
    hInt_gaus->Fill(Int_UpS_gaus/Int_UpS);

    //////////
    //DS
    //distribution
    cDS->cd(i+1); gPad->SetLogy();
    h_DS[i] =
      (TH1D*)fPsi->Get(Form("MuMu_left_upstream_up_%s%i", physBinned.Data(),i));
    h_DS[i]->Sumw2();
    h_DS[i]->Draw();

    //fit crystal ball
    TF1 *f_DS_cry = new TF1("crystal", CrystalBall, Mmin, Mmax, 5);
    ParaCrystalBall(h_DS[i], f_DS_cry, M_Psi);
    DoFit(h_DS[i], f_DS_cry, hChi2_cry, Mmin, Mmax);

    //fit gaus
    TF1* fgaus_DS =DoFit(h_UpS[i], "gaus", hChi2_gaus, Mmin, Mmax);
    cDS->Update();

    //crystal ball ratios
    cRDS_cry->cd(i+1); gPad->SetLogy();
    hF_DS_cry[i] = new TH1D(Form("hF_DS_cry_%i", i), Form("hF_DS_cry_%i", i),
			   nhbin, xmin, xmax);
    FillFromFit(f_DS_cry, hF_DS_cry[i], Mmin, Mmax);
    r_DS_cry[i] = new TRatioPlot(h_DS[i], hF_DS_cry[i]);
    SetUp(r_DS_cry[i]); cRDS_cry->Update();

    //gaus ratio
    cRDS_gaus->cd(i+1); gPad->SetLogy();
    hF_DS_gaus[i] = new TH1D(Form("hR_DS_gaus_%i", i),
				Form("hR_DS_gaus_%i", i),
				nhbin, xmin, xmax);
    FillFromFit(fgaus_DS, hF_DS_gaus[i], Mmin, Mmax);
    r_DS_gaus[i] = new TRatioPlot(h_DS[i], hF_DS_gaus[i]);
    SetUp(r_DS_gaus[i]); cRDS_gaus->Update();

    //Intgrals
    minBin = h_DS[i]->FindBin(Mmin);
    maxBin = h_DS[i]->FindBin(Mmax);
    Double_t er_DS, er_DS_cry, er_DS_gaus;//cleanup
    Double_t Int_DS =h_DS[i]->IntegralAndError(minBin, maxBin, er_DS);
    Double_t Int_DS_cry =
      hF_DS_cry[i]->IntegralAndError(minBin, maxBin, er_DS_cry);
    Double_t Int_DS_gaus =
      hF_DS_gaus[i]->IntegralAndError(minBin, maxBin, er_DS_gaus);
    hInt_cry->Fill(Int_DS_cry/Int_DS);
    hInt_gaus->Fill(Int_DS_gaus/Int_DS);
  }//Bins loop

  //One ratio Side by side comparison
  TCanvas* c1 = new TCanvas(); c1->Divide(2);
  c1->cd(1); gPad->SetLogy();
  hF_UpS_gaus[0]->SetFillColor(kRed); r_UpS_gaus[0]->SetH2DrawOpt("bar");
  r_UpS_gaus[0]->Draw();
  c1->cd(2); gPad->SetLogy();
  hF_UpS_cry[0]->SetFillColor(kBlue); r_UpS_cry[0]->SetH2DrawOpt("bar");
  r_UpS_cry[0]->Draw();

  //Chi2 Distributions
  TCanvas* cChi2 = new TCanvas(); cChi2->Divide(2);
  cChi2->cd(1);
  hChi2_cry->Draw();
  cChi2->cd(2);
  hChi2_gaus->Draw(); hChi2_gaus->SetLineColor(kRed);

  //Integrals
  TCanvas* cInt = new TCanvas();
  hInt_cry->Draw(); 
  hInt_gaus->Draw("sames"); hInt_gaus->SetLineColor(kRed);


  //Ouput settings
  cout << "\nSettings are:" << endl;
  cout << "Physics binning in:    " << physBinned << endl;
  cout << "Integration Min:       " << Mmin << endl;
  cout << "Integration Max:       " << Mmax << endl;
  cout << "\nWhich resonance fit:   " << whichPsi << endl;
  cout << "\n" << endl;
}
