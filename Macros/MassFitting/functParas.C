#include "include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/AN_calculation/include/fitFunctions.h"

void functParas(TString start=""){
  
  //Setup_______________
  /*TString period_Mtype ="W07_HMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =3;
  TString binRange ="43_85";
  TString physBinned ="pT";//"xN", "xPi", "xF", "pT"
  TString process ="DY";//JPsi, psi, DY
  Double_t LR_Mmin =4.30;
  Double_t LR_Mmax =8.50;//L/R counts mass range
  Double_t Mmin =4.30;//Fit Mass minimum
  Double_t Mmax =8.50;//Fit Mass maximum
  TString whichFit ="eight";//*/
  
  TString period_Mtype ="WAll_LowM_AMDY";
  Int_t hbins =150;//# of histogram bins using in mass fitting
  const Int_t nBins =5;
  TString binRange ="25_43";
  TString physBinned ="xN";//"xN", "xPi", "xF", "pT"
  TString process ="JPsi";//JPsi, psi, DY
  Double_t LR_Mmin =2.90;
  Double_t LR_Mmax =3.30;//L/R counts mass range
  Double_t Mmin =2.00;//Fit Mass minimum
  Double_t Mmax =6.50;//Fit Mass maximum
  TString whichFit ="eight";//*/
  
  Bool_t toWrite =false;
  //Setup_______________

  Double_t nominal_Mmin =Mmin, nominal_Mmax =Mmax;
  TString pathRD = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_leftRight_Analysis/Data/";
  TString RDfile =Form("leftRight_byTarget_%s1.00_8.50_%ibins%s_%ihbin.root",
		       period_Mtype.Data(), nBins, binRange.Data(), hbins);
  if (start==""){
    cout<<"Script outputs the parameters of the specified fit" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'functParas.C(Bool_t PolCorr =true, 1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Real data path:             " << pathRD << endl;
    cout << "Real data file considered:  " << RDfile << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Binned in which DY physics:  " << physBinned << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << LR_Mmin << "  -  "
	 << LR_Mmax << endl;
    cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
    cout << "Which fit considered:       " << whichFit << endl;
    cout << "\nTo write output file:       " << toWrite << endl;
    exit(EXIT_FAILURE);
  }

  //Get Input Files from Local_leftRight
  TFile *fRD  = TFile::Open(pathRD + RDfile);
  if (!fRD ){
    cout << "RD file does not exist " << endl;
    exit(EXIT_FAILURE);
  }

  //Determine x-binning
  Double_t ex[nBins] = {0.0};
  TGraph* g_Pol =(TGraph*)fRD->Get(Form("%s_Pol", physBinned.Data()));
  Double_t *xvals = g_Pol->GetX();

  //Get Input Hist/fit ratio/Determine LR counts for specified process
  TH1D *hRD_upS_up_L[nBins], *hRD_upS_up_R[nBins]; //Invariant M dist to be Fit
  TH1D *hRD_upS_down_L[nBins], *hRD_upS_down_R[nBins];
  TH1D *hRD_downS_up_L[nBins], *hRD_downS_up_R[nBins];
  TH1D *hRD_downS_down_L[nBins], *hRD_downS_down_R[nBins];

  //Output fit parameters 
  const Int_t nSelectPars =10;
  vector <vector<Double_t> > pars_upS_up(nSelectPars);
  vector <vector<Double_t> > pars_upS_down(nSelectPars);
  vector <vector<Double_t> > pars_downS_up(nSelectPars);
  vector <vector<Double_t> > pars_downS_down(nSelectPars);
  vector <vector<Double_t> > e_pars_upS_up(nSelectPars);
  vector <vector<Double_t> > e_pars_upS_down(nSelectPars);
  vector <vector<Double_t> > e_pars_downS_up(nSelectPars);
  vector <vector<Double_t> > e_pars_downS_down(nSelectPars);
  
  //L/R counts for specified process //needed for consistency
  Double_t LR_upS_up_L[nBins], LR_upS_up_R[nBins]; 
  Double_t LR_upS_down_L[nBins], LR_upS_down_R[nBins];
  Double_t LR_downS_up_L[nBins], LR_downS_up_R[nBins];
  Double_t LR_downS_down_L[nBins], LR_downS_down_R[nBins];
  
  Double_t e_LR_upS_up_L[nBins], e_LR_upS_up_R[nBins];
  Double_t e_LR_upS_down_L[nBins], e_LR_upS_down_R[nBins];
  Double_t e_LR_downS_up_L[nBins], e_LR_downS_up_R[nBins];
  Double_t e_LR_downS_down_L[nBins], e_LR_downS_down_R[nBins];
  
  //TargPol setup
  const Int_t nTargPol =4;
  TString targNames[nTargPol] = {"upS_up", "upS_down", "downS_up","downS_down"};

  //Perform Fits and integrate to get L/R
  for (Int_t bi=0, iPar_upS=0, iPar_dS=0; bi<nBins; bi++) {
    hRD_upS_up_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_up_%s%i", 
					    physBinned.Data(), bi) );
    hRD_upS_up_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_up_%s%i", 
					    physBinned.Data(), bi) );
    hRD_upS_down_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_upstream_down_%s%i", 
					      physBinned.Data(), bi) );
    hRD_upS_down_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_upstream_down_%s%i", 
					      physBinned.Data(), bi) );
    SetUp(hRD_upS_up_L[bi]); SetUp(hRD_upS_up_R[bi]); 
    SetUp(hRD_upS_down_L[bi]); SetUp(hRD_upS_down_R[bi]);
    
    hRD_downS_up_L[bi] = (TH1D*)fRD->Get(Form("MuMu_left_downstream_up_%s%i", 
					      physBinned.Data(), bi) );
    hRD_downS_up_R[bi] = (TH1D*)fRD->Get(Form("MuMu_right_downstream_up_%s%i", 
					      physBinned.Data(), bi) );
    hRD_downS_down_L[bi]= (TH1D*)fRD->Get(Form("MuMu_left_downstream_down_%s%i",
					       physBinned.Data(), bi) );
    hRD_downS_down_R[bi]=(TH1D*)fRD->Get(Form("MuMu_right_downstream_down_%s%i",
					      physBinned.Data(), bi) );
    SetUp(hRD_downS_up_L[bi]); SetUp(hRD_downS_up_R[bi]);
    SetUp(hRD_downS_down_L[bi]); SetUp(hRD_downS_down_R[bi]);

    Bool_t hIsUpS =true; 
    TF1 *fitFunc = FitGetLR(hRD_upS_up_L, bi, LR_upS_up_L, e_LR_upS_up_L,
			    LR_Mmin, LR_Mmax, process, whichFit,
			    Mmin, Mmax);
    ProcessPars_eight(fitFunc, pars_upS_up, e_pars_upS_up, hIsUpS);

    hIsUpS =false; 
    fitFunc = FitGetLR(hRD_downS_up_L, bi, LR_downS_up_L, e_LR_downS_up_L,
			    LR_Mmin, LR_Mmax, process, whichFit,
			    Mmin, Mmax);
    ProcessPars_eight(fitFunc, pars_downS_up, e_pars_downS_up, hIsUpS);
  }//Loop over nBins physics binning


  TCanvas *cPars_upS_up = new TCanvas("pars_upS_up", "pars_upS_up");
  cPars_upS_up->Divide(2, nSelectPars/2);
  TCanvas *cPars_downS_up = new TCanvas("pars_downS_up", "pars_downS_up");
  cPars_downS_up->Divide(2, nSelectPars/2);
  for (Int_t i=0; i<nSelectPars; i++) {
    //UpS
    TGraphErrors *gr_upS_up
      = new TGraphErrors(nBins, xvals, &(pars_upS_up[i][0]),
			 ex, &(e_pars_upS_up[i][0]));
    SetUp(gr_upS_up);

    cPars_upS_up->cd(i+1);
    gr_upS_up->Draw("AP");
    gr_upS_up->Fit("pol0", "WQ");
    TF1 *f_upS_up = gr_upS_up->GetFunction("pol0");
    cout << f_upS_up->GetParameter(0) << " +/- " <<
      f_upS_up->GetParError(0) << "          ";

    //DownS
    TGraphErrors *gr_downS_up
      = new TGraphErrors(nBins, xvals, &(pars_downS_up[i][0]),
			 ex, &(e_pars_downS_up[i][0]));
    SetUp(gr_downS_up);

    cPars_downS_up->cd(i+1);
    gr_downS_up->Draw("AP");
    gr_downS_up->Fit("pol0", "WQ");
    TF1 *f_downS_up = gr_downS_up->GetFunction("pol0");
    cout << f_downS_up->GetParameter(0) << " +/- " <<
      f_downS_up->GetParError(0) << endl;
  }
    
  /*//Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/AN_calculation";
  TString fOutput
    = Form("%s/Data/functMFit/functMFit_%s%.2f_%.2f_%s_%s%.2f_%.2f_%s%i_%ihbin",
	   thisDirPath.Data(), whichFit.Data(), nominal_Mmin, nominal_Mmax,
	   period_Mtype.Data(), process.Data(), LR_Mmin, LR_Mmax,
	   physBinned.Data(), nBins, hbins);
  fOutput += "_corr.root";
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TList *doc = new TList();
    doc->Add((TObject*)(new TObjString(pathRD+"\n"+RDfile+"\n")));
    doc->Write("InputData");
    
    g_AN_upS_up->Write("AN_upS_up");
    g_AN_upS_down->Write("AN_upS_down");
    g_AN_downS_up->Write("AN_downS_up");
    g_AN_downS_down->Write("AN_downS_down");

    g_Left_upS_up->Write("Left_upS_up");
    g_Left_upS_down->Write("Left_upS_down");
    g_Left_downS_up->Write("Left_downS_up");
    g_Left_downS_down->Write("Left_downS_down");

    g_Right_upS_up->Write("Right_upS_up");
    g_Right_upS_down->Write("Right_upS_down");
    g_Right_downS_up->Write("Right_downS_up");
    g_Right_downS_down->Write("Right_downS_down");
  }

  cout << " " << endl;
  cout << "Settings______" << endl;
  cout << "Real data file considered:              " << RDfile << endl;
  cout <<"Number of histogram bins used in M fit:  "
       <<hRD_upS_up_L[0]->GetNbinsX()<<endl;
  cout << "Mass Range for fitting is:              " << Mmin
       << " - " << Mmax << endl;
  cout << "physBinned nBins times:                 " << nBins << endl;
  cout << "\n";
  cout << "AN for physical process:                " << process << endl;
  cout << "Physics binning is:                     " << physBinned << endl;
  cout << "LR integral mass range:                 " << LR_Mmin
       << " - " << LR_Mmax << endl;
  cout << "Fit mass range:     " << Mmin << "  -  " << Mmax << endl;
  cout << "Which fit considered:       " << whichFit << endl;
  cout << " " << endl;
  if (toWrite){
    cout << "File:  " << fOutput << "   was written" << endl;
  }
  else cout << "File: " << fOutput << " was NOT written" << endl;//*/
}
