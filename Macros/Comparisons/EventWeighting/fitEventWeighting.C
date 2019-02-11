#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "TMinuit.h"

Double_t sums[6], cov[36], Pol[2]={1.0, 1.0};

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
  //calculate chisquare
  Double_t Lup =par[0], Ldown =par[1];
  Double_t a1 =par[2], a2 =par[3], a3 =par[4];
  Double_t A =par[5];
  Double_t model[6] =
    {Lup*( 1+0.5*a1*Pol[0]*A ),
     0.5*Lup*( Pol[0]*A*( 1+0.5*a2 ) + a1 ),
     0.5*Lup*( 1 + 0.5*a2 + 0.25*(3*a1+a3)*Pol[0]*A ),
     Ldown*( 1-0.5*a1*Pol[1]*A ),
     0.5*Ldown*( a1 - Pol[1]*A*( 1+0.5*a2 ) ),
     0.5*Ldown*( 1 + 0.5*a2 - 0.25*(3*a1+a3)*Pol[1]*A ) };
  TVectorD yobs; yobs.Use(6, sums);
  TVectorD ymodel; ymodel.Use(6, model);
  TMatrixDSym C(6); C.SetMatrixArray(cov);
  C.Invert();
  
  TVectorD diff = yobs - ymodel;
  
  TVectorD mult = C*diff;
  f = diff*mult;
}

void SetUpSumsAndCov(TGraph **g_Pup, TGraph **g_Pdown, Int_t bin=0);

void InitializePars(TMinuit *gMinuit, Double_t *vstart, Double_t *step,
		    Int_t ierflg);

void CorrectForP(Double_t &A_uncorr, Double_t &eA_uncorr,
		 Double_t &A, Double_t &eA, TGraph *g_pol, Int_t bin=0);

void CheckFit(Int_t icstat, TString targ="none", Int_t bin=0){
  if (icstat != 3){
    cout << "Some fit problems target:  " << targ << "  bin:  " << bin << endl;
    exit(EXIT_FAILURE);
  }
}

void ExpectedErrorA(Double_t &eA_expected, TGraph* g_Pup, TGraph *g_Pdown,
		    TGraph *g_pol, Int_t bin=0);

void CheckPhysicality(Double_t L_Pup, Double_t L_Pdown,
		      Double_t a1, Double_t a2, Double_t a3, Double_t A);
  
//______________________________________________________________________________
void fitEventWeighting()
{
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  TString production ="slot1";
  TString period ="W08";
  TString physBinned ="xN";
  TString minimizer ="MINOS";//MIGRAD, HESSE, MINOS

  Bool_t toWrite =false;
  //Setup_______________

  //Minuit setup
  TMinuit *gMinuit = new TMinuit(6);  //initialize TMinuit 
  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1.0;//number of std in error bar (for chi2 fit)
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  arglist[0] = 10000;//500, 10000 Number of calls
  arglist[1] = 1.;

  // Set starting values and step sizes for parameters
  static Double_t vstart[6] = {1000, 1000 , 0.01, 0.01, 0.01, 0.05};
  static Double_t step[6] = {10, 10, 0.1, 0.1, 0.1, 0.001};

  //Get data
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Comparisons/EventWeighting/Data";
  TString fname =
    Form("%s/Counts/counts_%s%i_%s_%s%s.root", thisDirPath.Data(),
	 physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	 period.Data());
  TFile *fIn = OpenFile(fname);

  const Int_t nWeights =5;
  TString whichWeight[nWeights] ={"", "c_", "c2_", "c3_", "c4_"};
  TGraph *g_upS_Pup[nWeights], *g_upS_Pdown[nWeights];
  TGraph *g_upS_Pup_phys[nWeights], *g_upS_Pdown_phys[nWeights];
  TGraph *g_downS_Pup[nWeights], *g_downS_Pdown[nWeights];
  TGraph *g_downS_Pup_phys[nWeights], *g_downS_Pdown_phys[nWeights];
  for (Int_t i=0; i<nWeights; i++) {
    g_upS_Pup[i] =
      (TGraph*)fIn->Get(Form("%supS_Pup", whichWeight[i].Data()));
    g_upS_Pdown[i] =
      (TGraph*)fIn->Get(Form("%supS_Pdown", whichWeight[i].Data()));
    g_upS_Pup_phys[i] =
      (TGraph*)fIn->Get(Form("%supS_Pup_phys", whichWeight[i].Data()));
    g_upS_Pdown_phys[i] =
      (TGraph*)fIn->Get(Form("%supS_Pdown_phys", whichWeight[i].Data()));

    g_downS_Pup[i] =
      (TGraph*)fIn->Get(Form("%sdownS_Pup", whichWeight[i].Data()));
    g_downS_Pdown[i] =
      (TGraph*)fIn->Get(Form("%sdownS_Pdown", whichWeight[i].Data()));
    g_downS_Pup_phys[i] =
      (TGraph*)fIn->Get(Form("%sdownS_Pup_phys", whichWeight[i].Data()));
    g_downS_Pdown_phys[i] =
      (TGraph*)fIn->Get(Form("%sdownS_Pdown_phys", whichWeight[i].Data()));
  }
  TGraph *g_pol_upS = (TGraph*)fIn->Get("pol_upS");
  TGraph *g_pol_downS = (TGraph*)fIn->Get("pol_downS");
  TGraph *g_pol_upS_phys = (TGraph*)fIn->Get("pol_upS_phys");
  TGraph *g_pol_downS_phys = (TGraph*)fIn->Get("pol_downS_phys");
  Double_t *xval = g_pol_upS_phys->GetX();
  Double_t ex[nBins] ={0.0};

  //Upstream
  ///////////////
  //Integrated
  cout << "\n\nUpstream Integrated" << endl;
  SetUpSumsAndCov(g_upS_Pup, g_upS_Pdown);
  InitializePars(gMinuit, vstart, step, ierflg);
  gMinuit->mnexcm(minimizer, arglist ,2,ierflg);
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  CheckFit(icstat, "upS");
  
  Double_t A_upS_int, eA_upS_int;
  gMinuit->GetParameter(5, A_upS_int, eA_upS_int);
  Double_t A_upS_int_corr, eA_upS_int_corr;
  CorrectForP(A_upS_int, eA_upS_int, A_upS_int_corr, eA_upS_int_corr,
	      g_pol_upS);

  Double_t L_upS_Pup_int, L_upS_Pdown_int;
  Double_t a1_upS_int, a2_upS_int, a3_upS_int;
  Double_t eL_upS_Pup_int, eL_upS_Pdown_int;
  Double_t ea1_upS_int, ea2_upS_int, ea3_upS_int;
  Double_t eA_upS_int_expected;
  gMinuit->GetParameter(0, L_upS_Pup_int, eL_upS_Pup_int);
  gMinuit->GetParameter(1, L_upS_Pdown_int, eL_upS_Pdown_int);
  gMinuit->GetParameter(2, a1_upS_int, ea1_upS_int);
  gMinuit->GetParameter(3, a2_upS_int, ea2_upS_int);
  gMinuit->GetParameter(4, a3_upS_int, ea3_upS_int);
  
  ExpectedErrorA(eA_upS_int_expected, g_upS_Pup[0], g_upS_Pdown[0], g_pol_upS);
  CheckPhysicality(L_upS_Pup_int, L_upS_Pdown_int,
		   a1_upS_int, a2_upS_int, a3_upS_int,
		   A_upS_int);
  
  //Phys Binned
  Double_t A_upS_phys[nBins], eA_upS_phys[nBins];
  Double_t A_upS_phys_corr[nBins], eA_upS_phys_corr[nBins];
  Double_t L_upS_Pup_phys[nBins], L_upS_Pdown_phys[nBins];
  Double_t a1_upS_phys[nBins], a2_upS_phys[nBins], a3_upS_phys[nBins];
  Double_t eL_upS_Pup_phys[nBins], eL_upS_Pdown_phys[nBins];
  Double_t ea1_upS_phys[nBins], ea2_upS_phys[nBins], ea3_upS_phys[nBins];
  Double_t eA_upS_phys_expected[nBins];
  for (Int_t i=0; i<nBins; i++) {
    cout << "\n\nUpstream Bin:  " << i << endl;
    SetUpSumsAndCov(g_upS_Pup_phys, g_upS_Pdown_phys, i);
    InitializePars(gMinuit, vstart, step, ierflg);
    gMinuit->mnexcm(minimizer, arglist ,2,ierflg);
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    CheckFit(icstat, "upS", i);
  
    gMinuit->GetParameter(5, A_upS_phys[i], eA_upS_phys[i]);
    gMinuit->GetParameter(0, L_upS_Pup_phys[i], eL_upS_Pup_phys[i]);
    gMinuit->GetParameter(1, L_upS_Pdown_phys[i], eL_upS_Pdown_phys[i]);
    gMinuit->GetParameter(2, a1_upS_phys[i], ea1_upS_phys[i]);
    gMinuit->GetParameter(3, a2_upS_phys[i], ea2_upS_phys[i]);
    gMinuit->GetParameter(4, a3_upS_phys[i], ea3_upS_phys[i]);

    CheckPhysicality(L_upS_Pup_phys[i], L_upS_Pdown_phys[i],
		     a1_upS_phys[i], a2_upS_phys[i], a3_upS_phys[i],
		     A_upS_phys[i]);
    
    CorrectForP(A_upS_phys[i], eA_upS_phys[i],
		A_upS_phys_corr[i], eA_upS_phys_corr[i], g_pol_upS_phys, i);

    ExpectedErrorA(eA_upS_phys_expected[i],
		   g_upS_Pup_phys[0], g_upS_Pdown_phys[0], g_pol_upS_phys, i);
  }

  //Graphs
  TGraphErrors *g_A_upS_int =
    new TGraphErrors(1, &(xval[0]), &A_upS_int_corr, ex, &eA_upS_int_corr);
  TGraphErrors *g_L_upS_Pup_int =
    new TGraphErrors(1, &(xval[0]), &L_upS_Pup_int, ex, &eL_upS_Pup_int);
  TGraphErrors *g_L_upS_Pdown_int =
    new TGraphErrors(1, &(xval[0]), &L_upS_Pdown_int, ex, &eL_upS_Pdown_int);
  TGraphErrors *g_a1_upS_int =
    new TGraphErrors(1, &(xval[0]), &a1_upS_int, ex, &ea1_upS_int);
  TGraphErrors *g_a2_upS_int =
    new TGraphErrors(1, &(xval[0]), &a2_upS_int, ex, &ea2_upS_int);
  TGraphErrors *g_a3_upS_int =
    new TGraphErrors(1, &(xval[0]), &a3_upS_int, ex, &ea3_upS_int);
  TGraphErrors *g_A_upS_int_expected =
    new TGraphErrors(1, &(xval[0]), &A_upS_int_corr, ex, &eA_upS_int_expected);
  
  TGraphErrors *g_A_upS_phys =
    new TGraphErrors(nBins, xval, A_upS_phys_corr, ex, eA_upS_phys_corr);
  TGraphErrors *g_L_upS_Pup_phys =
    new TGraphErrors(nBins, xval, L_upS_Pup_phys, ex, eL_upS_Pup_phys);
  TGraphErrors *g_L_upS_Pdown_phys =
    new TGraphErrors(nBins, xval, L_upS_Pdown_phys, ex, eL_upS_Pdown_phys);
  TGraphErrors *g_a1_upS_phys =
    new TGraphErrors(nBins, xval, a1_upS_phys, ex, ea1_upS_phys);
  TGraphErrors *g_a2_upS_phys =
    new TGraphErrors(nBins, xval, a2_upS_phys, ex, ea2_upS_phys);
  TGraphErrors *g_a3_upS_phys =
    new TGraphErrors(nBins, xval, a3_upS_phys, ex, ea3_upS_phys);
  TGraphErrors *g_A_upS_phys_expected =
    new TGraphErrors(nBins, xval, A_upS_phys_corr, ex, eA_upS_phys_expected);
  
  SetUp(g_A_upS_int); SetUp(g_L_upS_Pup_int); SetUp(g_L_upS_Pdown_int);
  SetUp(g_a1_upS_int); SetUp(g_a2_upS_int); SetUp(g_a3_upS_int);
  
  SetUp(g_A_upS_phys); SetUp(g_L_upS_Pup_phys); SetUp(g_L_upS_Pdown_phys);
  SetUp(g_a1_upS_phys); SetUp(g_a2_upS_phys); SetUp(g_a3_upS_phys);
  SetUp(g_A_upS_int_expected); SetUp(g_A_upS_phys_expected);

  TString parTitles[] = {"A_{N}", "L Pup", "L Pdown", "a1", "a2", "a3"};
  TGraphErrors *g_upS_phys[] =
    {g_A_upS_phys, g_L_upS_Pup_phys, g_L_upS_Pdown_phys,
     g_a1_upS_phys, g_a2_upS_phys, g_a3_upS_phys};
  TGraphErrors *g_upS_int[] =
    {g_A_upS_int, g_L_upS_Pup_int, g_L_upS_Pdown_int,
     g_a1_upS_int, g_a2_upS_int, g_a3_upS_int};

  TCanvas* cupS_phys = new TCanvas(); cupS_phys->Divide(3, 2);
  TCanvas* cupS_int = new TCanvas(); cupS_int->Divide(3, 2);
  for (Int_t i=0; i<6; i++) {
    cupS_phys->cd(i+1);
    g_upS_phys[i]->Draw("AP");
    g_upS_phys[i]->SetTitle(Form("upstream %s", parTitles[i].Data()));

    cupS_int->cd(i+1);
    g_upS_int[i]->Draw("AP");
    g_upS_int[i]->SetTitle(Form("upstream %s", parTitles[i].Data()));
  }

  cupS_phys->cd(1);
  OffSet(g_A_upS_phys_expected);
  g_A_upS_phys_expected->Draw("Psame");
  g_A_upS_phys_expected->SetMarkerColor(kRed);

  cupS_int->cd(1);
  OffSet(g_A_upS_int_expected);
  g_A_upS_int_expected->Draw("Psame");
  g_A_upS_int_expected->SetMarkerColor(kRed);

  //Downstream
  ///////////////
  //Integrated
  cout << "\n\nDownstream Integrated" << endl;
  SetUpSumsAndCov(g_downS_Pup, g_downS_Pdown);
  InitializePars(gMinuit, vstart, step, ierflg);
  gMinuit->mnexcm(minimizer, arglist ,2,ierflg);
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  CheckFit(icstat, "downS");
  
  Double_t A_downS_int, eA_downS_int;
  gMinuit->GetParameter(5, A_downS_int, eA_downS_int);
  Double_t A_downS_int_corr, eA_downS_int_corr;
  CorrectForP(A_downS_int, eA_downS_int, A_downS_int_corr, eA_downS_int_corr,
	      g_pol_downS);

  Double_t L_downS_Pup_int, L_downS_Pdown_int;
  Double_t a1_downS_int, a2_downS_int, a3_downS_int;
  Double_t eL_downS_Pup_int, eL_downS_Pdown_int;
  Double_t ea1_downS_int, ea2_downS_int, ea3_downS_int;
  Double_t eA_downS_int_expected;
  gMinuit->GetParameter(0, L_downS_Pup_int, eL_downS_Pup_int);
  gMinuit->GetParameter(1, L_downS_Pdown_int, eL_downS_Pdown_int);
  gMinuit->GetParameter(2, a1_downS_int, ea1_downS_int);
  gMinuit->GetParameter(3, a2_downS_int, ea2_downS_int);
  gMinuit->GetParameter(4, a3_downS_int, ea3_downS_int);

  ExpectedErrorA(eA_downS_int_expected, g_downS_Pup[0], g_downS_Pdown[0],
		 g_pol_downS);
  CheckPhysicality(L_downS_Pup_int, L_downS_Pdown_int,
		   a1_downS_int, a2_downS_int, a3_downS_int,
		   A_downS_int);

  //Phys Binned
  Double_t A_downS_phys[nBins], eA_downS_phys[nBins];
  Double_t A_downS_phys_corr[nBins], eA_downS_phys_corr[nBins];
  Double_t L_downS_Pup_phys[nBins], L_downS_Pdown_phys[nBins];
  Double_t a1_downS_phys[nBins], a2_downS_phys[nBins], a3_downS_phys[nBins];
  Double_t eL_downS_Pup_phys[nBins], eL_downS_Pdown_phys[nBins];
  Double_t ea1_downS_phys[nBins], ea2_downS_phys[nBins], ea3_downS_phys[nBins];
  Double_t eA_downS_phys_expected[nBins];
  for (Int_t i=0; i<nBins; i++) {
    SetUpSumsAndCov(g_downS_Pup_phys, g_downS_Pdown_phys, i);
    InitializePars(gMinuit, vstart, step, ierflg);
    gMinuit->mnexcm(minimizer, arglist ,2,ierflg);
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    CheckFit(icstat, "downS", i);
  
    gMinuit->GetParameter(5, A_downS_phys[i], eA_downS_phys[i]);
    gMinuit->GetParameter(0, L_downS_Pup_phys[i], eL_downS_Pup_phys[i]);
    gMinuit->GetParameter(1, L_downS_Pdown_phys[i], eL_downS_Pdown_phys[i]);
    gMinuit->GetParameter(2, a1_downS_phys[i], ea1_downS_phys[i]);
    gMinuit->GetParameter(3, a2_downS_phys[i], ea2_downS_phys[i]);
    gMinuit->GetParameter(4, a3_downS_phys[i], ea3_downS_phys[i]);

    CheckPhysicality(L_downS_Pup_phys[i], L_downS_Pdown_phys[i],
		     a1_downS_phys[i], a2_downS_phys[i], a3_downS_phys[i],
		     A_downS_phys[i]);
     
    CorrectForP(A_downS_phys[i], eA_downS_phys[i],
		A_downS_phys_corr[i], eA_downS_phys_corr[i],g_pol_downS_phys,i);

    ExpectedErrorA(eA_downS_phys_expected[i],
		   g_downS_Pup_phys[0], g_downS_Pdown_phys[0],
		   g_pol_downS_phys, i);
  }

  //Graphs
  TGraphErrors *g_A_downS_int =
    new TGraphErrors(1, &(xval[0]), &A_downS_int_corr, ex, &eA_downS_int_corr);
  TGraphErrors *g_L_downS_Pup_int =
    new TGraphErrors(1, &(xval[0]), &L_downS_Pup_int, ex, &eL_downS_Pup_int);
  TGraphErrors *g_L_downS_Pdown_int =
    new TGraphErrors(1, &(xval[0]), &L_downS_Pdown_int, ex,&eL_downS_Pdown_int);
  TGraphErrors *g_a1_downS_int =
    new TGraphErrors(1, &(xval[0]), &a1_downS_int, ex, &ea1_downS_int);
  TGraphErrors *g_a2_downS_int =
    new TGraphErrors(1, &(xval[0]), &a2_downS_int, ex, &ea2_downS_int);
  TGraphErrors *g_a3_downS_int =
    new TGraphErrors(1, &(xval[0]), &a3_downS_int, ex, &ea3_downS_int);
  TGraphErrors *g_A_downS_int_expected =
    new TGraphErrors(1, &(xval[0]),&A_downS_int_corr,ex,&eA_downS_int_expected);
  
  TGraphErrors *g_A_downS_phys =
    new TGraphErrors(nBins, xval, A_downS_phys_corr, ex, eA_downS_phys_corr);
  TGraphErrors *g_L_downS_Pup_phys =
    new TGraphErrors(nBins, xval, L_downS_Pup_phys, ex, eL_downS_Pup_phys);
  TGraphErrors *g_L_downS_Pdown_phys =
    new TGraphErrors(nBins, xval, L_downS_Pdown_phys, ex, eL_downS_Pdown_phys);
  TGraphErrors *g_a1_downS_phys =
    new TGraphErrors(nBins, xval, a1_downS_phys, ex, ea1_downS_phys);
  TGraphErrors *g_a2_downS_phys =
    new TGraphErrors(nBins, xval, a2_downS_phys, ex, ea2_downS_phys);
  TGraphErrors *g_a3_downS_phys =
    new TGraphErrors(nBins, xval, a3_downS_phys, ex, ea3_downS_phys);
  TGraphErrors *g_A_downS_phys_expected =
    new TGraphErrors(nBins, xval, A_downS_phys_corr,ex,eA_downS_phys_expected);

  SetUp(g_A_downS_int); SetUp(g_L_downS_Pup_int); SetUp(g_L_downS_Pdown_int);
  SetUp(g_a1_downS_int); SetUp(g_a2_downS_int); SetUp(g_a3_downS_int);
  
  SetUp(g_A_downS_phys); SetUp(g_L_downS_Pup_phys); SetUp(g_L_downS_Pdown_phys);
  SetUp(g_a1_downS_phys); SetUp(g_a2_downS_phys); SetUp(g_a3_downS_phys);
  SetUp(g_A_downS_int_expected); SetUp(g_A_downS_phys_expected);

  TGraphErrors *g_downS_phys[] =
    {g_A_downS_phys, g_L_downS_Pup_phys, g_L_downS_Pdown_phys,
     g_a1_downS_phys, g_a2_downS_phys, g_a3_downS_phys};
  TGraphErrors *g_downS_int[] =
    {g_A_downS_int, g_L_downS_Pup_int, g_L_downS_Pdown_int,
     g_a1_downS_int, g_a2_downS_int, g_a3_downS_int};

  TCanvas* cdownS_phys = new TCanvas(); cdownS_phys->Divide(3, 2);
  TCanvas* cdownS_int = new TCanvas(); cdownS_int->Divide(3, 2);
  for (Int_t i=0; i<6; i++) {
    cdownS_phys->cd(i+1);
    g_downS_phys[i]->Draw("AP");
    g_downS_phys[i]->SetTitle(Form("downstream %s", parTitles[i].Data()));

    cdownS_int->cd(i+1);
    g_downS_int[i]->Draw("AP");
    g_downS_int[i]->SetTitle(Form("downstream %s", parTitles[i].Data()));
  }

  cdownS_phys->cd(1);
  OffSet(g_A_downS_phys_expected);
  g_A_downS_phys_expected->Draw("Psame");
  g_A_downS_phys_expected->SetMarkerColor(kRed);

  cdownS_int->cd(1);
  OffSet(g_A_downS_int_expected);
  g_A_downS_int_expected->Draw("Psame");
  g_A_downS_int_expected->SetMarkerColor(kRed);
  
  //Write output
  TString fOutput =
    Form("%s/Weight/weight_%s%i_%s_%s%s_%s.root", thisDirPath.Data(),
	 physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	 period.Data(), minimizer.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    
    g_A_upS_int->Write("g_A_upS_int");
    g_L_upS_Pup_int->Write("g_L_upS_Pup_int");
    g_L_upS_Pdown_int->Write("g_L_upS_Pdown_int");
    g_a1_upS_int->Write("g_a1_upS_int");
    g_a2_upS_int->Write("g_a2_upS_int");
    g_a3_upS_int->Write("g_a3_upS_int");
    g_A_upS_int_expected->Write("g_A_upS_int_expected");
    
    g_A_upS_phys->Write("g_A_upS_phys");
    g_L_upS_Pup_phys->Write("g_L_upS_Pup_phys");
    g_L_upS_Pdown_phys->Write("g_L_upS_Pdown_phys");
    g_a1_upS_phys->Write("g_a1_upS_phys");
    g_a2_upS_phys->Write("g_a2_upS_phys");
    g_a3_upS_phys->Write("g_a3_upS_phys");
    g_A_upS_phys_expected->Write("g_A_upS_phys_expected");

    g_A_downS_int->Write("g_A_downS_int");
    g_L_downS_Pup_int->Write("g_L_downS_Pup_int");
    g_L_downS_Pdown_int->Write("g_L_downS_Pdown_int");
    g_a1_downS_int->Write("g_a1_downS_int");
    g_a2_downS_int->Write("g_a2_downS_int");
    g_a3_downS_int->Write("g_a3_downS_int");
    g_A_downS_int_expected->Write("g_A_downS_int_expected");
    
    g_A_downS_phys->Write("g_A_downS_phys");
    g_L_downS_Pup_phys->Write("g_L_downS_Pup_phys");
    g_L_downS_Pdown_phys->Write("g_L_downS_Pdown_phys");
    g_a1_downS_phys->Write("g_a1_downS_phys");
    g_a2_downS_phys->Write("g_a2_downS_phys");
    g_a3_downS_phys->Write("g_a3_downS_phys");
    g_A_downS_phys_expected->Write("g_A_downS_phys_expected");

    fResults->Close();
  }
}


void SetUpSumsAndCov(TGraph **g_Pup, TGraph **g_Pdown, Int_t bin=0){

  // The sums values
  for (Int_t i=0; i<3; i++) {
    sums[i] = g_Pup[i]->GetY()[bin];
    sums[i+3] = g_Pdown[i]->GetY()[bin];
  }

  // the covariance
  for (Int_t i=0; i<36; i++) {
    cov[i] = 0.0;
  }
  
  cov[0] = g_Pup[0]->GetY()[bin];
  cov[1] = g_Pup[1]->GetY()[bin];
  cov[2] = g_Pup[2]->GetY()[bin];
  cov[8] = g_Pup[3]->GetY()[bin];
  cov[14] = g_Pup[4]->GetY()[bin];

  cov[21] = g_Pdown[0]->GetY()[bin];
  cov[22] = g_Pdown[1]->GetY()[bin];
  cov[23] = g_Pdown[2]->GetY()[bin];
  cov[29] = g_Pdown[3]->GetY()[bin];
  cov[35] = g_Pdown[4]->GetY()[bin];

  //make cov symmetric
  cov[6]= cov[1]; cov[7] = cov[2]; cov[12]= cov[2]; cov[13]=cov[8];
  cov[27]= cov[22]; cov[28] = cov[23]; cov[33]= cov[23]; cov[34]=cov[29];
}


void InitializePars(TMinuit *gMinuit, Double_t *vstart, Double_t *step,
		    Int_t ierflg){

  gMinuit->mnparm(0, "a1", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "a2", vstart[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "a3", vstart[2], step[2], 0,0,ierflg);
  gMinuit->mnparm(3, "a4", vstart[3], step[3], 0,0,ierflg);
  gMinuit->mnparm(4, "a5", vstart[4], step[4], -1,1,ierflg);
  gMinuit->mnparm(5, "a6", vstart[5], step[5], 0,0,ierflg);
}


void CorrectForP(Double_t &A_uncorr, Double_t &eA_uncorr,
		 Double_t &A, Double_t &eA, TGraph *g_pol, Int_t bin=0){

  Double_t P = g_pol->GetY()[bin];
  A = A_uncorr/P;
  eA = eA_uncorr/P;
}


void ExpectedErrorA(Double_t &eA_expected, TGraph* g_Pup, TGraph *g_Pdown,
		    TGraph *g_pol, Int_t bin=0){
  
  eA_expected = g_Pup->GetY()[bin] + g_Pdown->GetY()[bin];
  eA_expected = 1.0/TMath::Sqrt( eA_expected );
  eA_expected = eA_expected/(g_pol->GetY()[bin]);
}


void CheckPhysicality(Double_t L_Pup, Double_t L_Pdown,
		      Double_t a1, Double_t a2, Double_t a3, Double_t A){
  Int_t error = 0;
  if (L_Pup < 0) {
    cout << "L_Pup unphysical" << endl;
    error =1;
  }
  if (L_Pdown < 0) {
    cout << "L_Pdown unphysical" << endl;
    error =1;
  }

  if ((a1 > 1) || (a1 < -1)) {
    cout << "a1 unphysical" << endl;
    error =1;
  }
  if ((a2 > 1) || (a2 < -1)) {
    cout << "a2 unphysical" << endl;
    error =1;
  }
  if ((a3 > 1) || (a3 < -1)) {
    cout << "a3 unphysical" << endl;
    error =1;
  }

  if ((A > 1) || (A < -1)) {
    cout << "A unphysical" << endl;
    error =1;
  }

  if (error == 1){
    exit(EXIT_FAILURE);
  }
}
