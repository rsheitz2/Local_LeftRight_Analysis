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


void FillFromFit(TF1* f, TH1D* h, Double_t Mmin, Double_t Mmax){
  Double_t minBin = h->FindBin(Mmin);
  Double_t maxBin = h->FindBin(Mmax);
  
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


void histVfit_integralCompare(TString start=""){
  //Setup_____
  const Int_t nBins =5;
  TString physBinned ="pT";
  Double_t Mmin =3.0;
  Double_t Mmax =4.0;
  TString whichPsi="Psi"; //"JPsi", "Psi";
  TString fitOptions ="RSQ+";
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
  
  //Distribution (h), Fit (hF), and Ratio (r) histograms
  TH1D *h_UpS[nBins];
  TH1D *hF_UpS_cry[nBins];
  TH1D *hF_UpS_gaus[nBins];
  TRatioPlot *r_UpS_cry[nBins];
  TRatioPlot *r_UpS_gaus[nBins];

  //reduced chi2 distributions
  TH1D *hChi2_gaus = new TH1D("chi2_gaus", "chi2_gaus", 100, 0, 500);
  TH1D *hChi2_cry = new TH1D("chi2_cry", "chi2_cry", 100, 0, 500);
  SetUp(hChi2_gaus); SetUp(hChi2_cry);
  hChi2_gaus->SetLineWidth(3); hChi2_cry->SetLineWidth(3);

  //Integral histograms
  TH1D* hInt_cry[nBins];

  //Intgral tgraph setup
  Double_t fInt_cry[nBins], ey_fInt_cry[nBins];
  Double_t histInt_cry[nBins], ey_histInt_cry[nBins];
  Double_t fInt_xval[nBins], ex[nBins] ={0.0};
  Double_t functInt[nBins], histInt[nBins];
  Double_t efunctInt[nBins], ehistInt[nBins];
  
  for (Int_t i=0; i<nBins; i++) {
    //////////
    //UpS
    //distribution
    cUpS->cd(i+1); gPad->SetLogy();
    h_UpS[i] =
      (TH1D*)fPsi->Get(Form("MuMu_left_upstream_up_%s%i", physBinned.Data(),i));
    h_UpS[i]->Sumw2(); 
    h_UpS[i]->Draw();

    //Basic histogram needed parameters
    Double_t width = h_UpS[i]->GetXaxis()->GetBinWidth(1);
    Int_t minBin =h_UpS[i]->FindBin(Mmin), maxBin =h_UpS[i]->FindBin(Mmax);
    Int_t nhbin = h_UpS[i]->GetNbinsX();
    Double_t xmin = h_UpS[i]->GetXaxis()->GetXmin();
    Double_t xmax = h_UpS[i]->GetXaxis()->GetXmax();

    //fit crystal ball
    TF1 *f_UpS_cry = new TF1("crystal", CrystalBall, Mmin, Mmax, 5);
    ParaCrystalBall(h_UpS[i], f_UpS_cry, M_Psi);
    h_UpS[i]->Fit(f_UpS_cry, fitOptions, "", Mmin, Mmax);
    hInt_cry[i] = new TH1D(Form("hF_UpS_cry_%i", i), Form("hF_UpS_cry_%i", i),
			   nhbin, xmin, xmax);
    FillFromFit(f_UpS_cry, hInt_cry[i], Mmin, Mmax);
    
    //Intgrals
    Double_t er_UpS, er_hist_cry;
    Double_t Int_UpS =h_UpS[i]->IntegralAndError(minBin, maxBin, er_UpS);
    Double_t functInt_cry = f_UpS_cry->Integral(Mmin, Mmax)/width;
    fInt_cry[i] = functInt_cry/(Int_UpS);
        
    TVirtualFitter *fitter = TVirtualFitter::GetFitter();
    Double_t *covMatrix = fitter->GetCovarianceMatrix();
    Double_t e_functInt_cry =f_UpS_cry->IntegralError(Mmin, Mmax)/width;
    ey_fInt_cry[i] = RatioError(functInt_cry, Int_UpS, e_functInt_cry, er_UpS);

    Double_t hist_cry =
      hInt_cry[i]->IntegralAndError(minBin, maxBin, er_hist_cry);
    histInt_cry[i] = hist_cry/Int_UpS;
    ey_histInt_cry[i] = RatioError(hist_cry, Int_UpS, er_hist_cry, er_UpS);

    functInt[i] = functInt_cry; efunctInt[i] = e_functInt_cry;
    histInt[i] = hist_cry; ehistInt[i] = er_hist_cry;

    fInt_xval[i] = i+1;
  }//Bins loop

  
  //Integrals
  TCanvas* cInt = new TCanvas();
  TGraphErrors *gfunctInt = new TGraphErrors(nBins, fInt_xval, fInt_cry,
					ex, ey_fInt_cry);
  
  SetUp(gfunctInt); 
  gfunctInt->Draw("AP"); gfunctInt->SetLineColor(kBlue);
  gfunctInt->GetYaxis()->SetRangeUser(0.95, 1.0);
  gfunctInt->SetTitle("Ratio: Int/True integral");

  TGraphErrors *ghistInt = new TGraphErrors(nBins, fInt_xval, histInt_cry,
					    ex, ey_histInt_cry);
  SetUp(ghistInt);
  ghistInt->SetMarkerColor(kBlack); ghistInt->Draw("Psame"); 


  TCanvas *cFull = new TCanvas();
  TGraphErrors *gfunct = new TGraphErrors(nBins, fInt_xval, functInt,
					  ex, efunctInt);
  SetUp(gfunct);
  TGraphErrors *ghist = new TGraphErrors(nBins, fInt_xval, histInt,
					  ex, ehistInt);
  SetUp(ghist);
  gfunct->Draw("AP"); ghist->Draw("Psame"); ghist->SetMarkerColor(kBlack);
  gfunct->SetTitle("Integral Counts");

  //Ouput settings
  cout << "\nSettings are:" << endl;
  cout << "Physics binning in:    " << physBinned << endl;
  cout << "Integration Min:       " << Mmin << endl;
  cout << "Integration Max:       " << Mmax << endl;
  cout << "fit options used:      " << fitOptions << endl;
  cout << "\nWhich resonance fit:   " << whichPsi << endl;
  cout << "\n" << endl;
}
