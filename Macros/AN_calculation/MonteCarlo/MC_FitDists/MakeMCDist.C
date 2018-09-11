void MakeMCDist(){
  //Setup_______________
  Bool_t toWrite=true;
  
  //TString MCtype ="JPsi"; Int_t icolor =6;
  //TString MCtype ="psi"; Int_t icolor =3;
  //TString MCtype ="OC"; Int_t icolor =4;
  TString MCtype ="AMDY"; Int_t icolor =2;
  
  const Int_t hbins=200;
  const Double_t Mmin=2.5, Mmax=8.5;
  //const Double_t Mmin=4.3, Mmax=8.5;
  //const Double_t Mmin=3.0, Mmax=3.26;

  TString MassCut = Form("Mmumu>%f&&Mmumu<%f", Mmin, Mmax);
  //Setup_______________  

  
  TString pathMC = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/MC_Data/YuShiangMC";
  TFile *fMC = TFile::Open(Form("%s/%s/Yu_Wall_full_main_%s_20bins.root",
				pathMC.Data(), MCtype.Data(), MCtype.Data()) );

  if ( !fMC ){
    cout << "File does not exist " << endl;
    exit(EXIT_FAILURE);
  }
  
  TTree *tMC = (TTree*)fMC->Get("pT_Weighted");
  TH1D *hist = new TH1D(Form("h_%s", MCtype.Data()),Form("h_%s", MCtype.Data()),
			hbins, Mmin, Mmax);

  TCanvas* c1 = new TCanvas();
  tMC->Draw(Form("Mmumu>>h_%s", MCtype.Data() ), MassCut, "0");
  hist->Sumw2();
  hist->Scale(1.0/(hist->Integral()) );
  hist->SetLineColor(icolor);
  hist->Draw();

  TString fOutput = Form("MC_%s_%.2f_%.2f.root", MCtype.Data(), Mmin, Mmax);
  if(toWrite){
    TFile *fResults =
      new TFile(fOutput, "RECREATE");
    hist->Write();
  }

  cout << " " << endl;
  cout << "Settings !!!!" << endl;
  cout << "Mass Range is: " << Mmin << "  -  " << Mmax << endl;
  cout << "Fit data for:  " << MCtype << endl;
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File:  " << fOutput << "   was  NOT  written" << endl;

}
