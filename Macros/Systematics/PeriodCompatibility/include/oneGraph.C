void oneGraph(){
  //Setup
  const Int_t nBins=4;
  Double_t y_sys_stat = 1.19893;

  //Setup
  
  Double_t yvals[nBins], xvals[nBins];
  for (Int_t bi=0; bi<nBins; bi++) {
    yvals[bi] = y_sys_stat;
    xvals[bi] = bi+0.5;
  }

  TGraph *g = new TGraph(nBins, xvals, yvals);
  TCanvas* c1 = new TCanvas();
  g->SetMarkerStyle(21);
  g->Draw("AP");


  /*//Write file
  TString fitMrangeType ="LowM_AMDY";
  TString process ="JPsi";//JPsi, psi
  TString lrMrange ="2.87_3.38";
  TString fitMrange ="2.87_3.38";
  TString binRange ="29_34"; //"25_43";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";
  TString thisDirPath=
    "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/Data/\
faPullDist2targWavg";
  TString fname =Form("%s/falseApullDist_%s_%s_%s%s_%s_%i_%s_%s.root",
		      thisDirPath.Data(), whichFit.Data(), fitMrangeType.Data(),
		      process.Data(), fitMrange.Data(), binRange.Data(),
		      nBins, production.Data(), additionalCuts.Data());
  cout << fname << " was written" << endl;
  TFile *fOut = new TFile(fname, "RECREATE");
  g->Write("gSys_Stat");

  fOut->Close();//*/
}
