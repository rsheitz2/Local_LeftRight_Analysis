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

void CorrectForP(Double_t A_uncorr, Double_t eA_uncorr,
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

Int_t DrawContour(TCanvas *c, TString whichBin="int");
  
//______________________________________________________________________________
void byTargFitEventWeighting()
{
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  TString production ="slot1";
  TString period ="WAll";
  TString physBinned ="xN";
  TString whichTarg ="downS";//"upS", "downS"
  TString minimizer ="MIGRAD";//MIGRAD, HESSE, MINOS

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
  //arglist[1] = 6.; //MINOS does errors for this par

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
  TGraph *g_Pup[nWeights], *g_Pdown[nWeights];
  TGraph *g_Pup_phys[nWeights], *g_Pdown_phys[nWeights];
  TGraph *g_downS_Pup[nWeights], *g_downS_Pdown[nWeights];
  TGraph *g_downS_Pup_phys[nWeights], *g_downS_Pdown_phys[nWeights];
  for (Int_t i=0; i<nWeights; i++) {
    g_Pup[i] = (TGraph*)fIn->Get(Form("%s%s_Pup",
				      whichWeight[i].Data(), whichTarg.Data()));
    g_Pdown[i] = (TGraph*)fIn->Get(Form("%s%s_Pdown", whichWeight[i].Data(),
					whichTarg.Data()));
    g_Pup_phys[i] =
      (TGraph*)fIn->Get(Form("%s%s_Pup_phys",
			     whichWeight[i].Data(), whichTarg.Data()));
    g_Pdown_phys[i] =
      (TGraph*)fIn->Get(Form("%s%s_Pdown_phys", whichWeight[i].Data(),
			     whichTarg.Data()));
  }
  TGraph *g_pol = (TGraph*)fIn->Get(Form("pol_%s", whichTarg.Data()));
  TGraph *g_pol_phys = (TGraph*)fIn->Get(Form("pol_%s_phys", whichTarg.Data()));
  Double_t *xval = g_pol_phys->GetX();
  Double_t ex[nBins] ={0.0};

  ///////////////
  //Integrated
  SetUpSumsAndCov(g_Pup, g_Pdown);
  InitializePars(gMinuit, vstart, step, ierflg);
  gMinuit->mnexcm(minimizer, arglist ,2,ierflg);
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  CheckFit(icstat, whichTarg);
  
  Double_t A_int, eA_int, eA_int_up, eA_int_down;
  gMinuit->GetParameter(5, A_int, eA_int);

  Double_t A_int_corr, eA_int_corr, eA_int_up_corr, eA_int_down_corr;
  CorrectForP(A_int, eA_int, A_int_corr, eA_int_corr, g_pol);

  Double_t eparab, gcc;
  if (minimizer=="MINOS"){
    gMinuit->mnerrs(5, eA_int_up, eA_int_down, eparab, gcc);
    eA_int_down *= -1;
    
    CorrectForP(A_int, eA_int_up, A_int_corr, eA_int_up_corr, g_pol);
    CorrectForP(A_int, eA_int_down, A_int_corr, eA_int_down_corr, g_pol);
  }
  
  Double_t L_Pup_int, L_Pdown_int;
  Double_t a1_int, a2_int, a3_int;
  Double_t eL_Pup_int, eL_Pdown_int;
  Double_t ea1_int, ea2_int, ea3_int;
  Double_t eA_int_expected;
  gMinuit->GetParameter(0, L_Pup_int, eL_Pup_int);
  gMinuit->GetParameter(1, L_Pdown_int, eL_Pdown_int);
  gMinuit->GetParameter(2, a1_int, ea1_int);
  gMinuit->GetParameter(3, a2_int, ea2_int);
  gMinuit->GetParameter(4, a3_int, ea3_int);
    
  ExpectedErrorA(eA_int_expected, g_Pup[0], g_Pdown[0], g_pol);
  CheckPhysicality(L_Pup_int, L_Pdown_int,
		   a1_int, a2_int, a3_int,
		   A_int);
    
  //Phys Binned
  Double_t A_phys[nBins], eA_phys[nBins];
  Double_t eA_phys_up[nBins],eA_phys_down[nBins];
  Double_t A_phys_corr[nBins], eA_phys_corr[nBins];
  Double_t eA_phys_up_corr[nBins],eA_phys_down_corr[nBins];
  Double_t L_Pup_phys[nBins], L_Pdown_phys[nBins];
  Double_t a1_phys[nBins], a2_phys[nBins], a3_phys[nBins];
  Double_t eL_Pup_phys[nBins], eL_Pdown_phys[nBins];
  Double_t ea1_phys[nBins], ea2_phys[nBins], ea3_phys[nBins];
  Double_t eA_phys_expected[nBins];
    
  for (Int_t i=0; i<nBins; i++) {
    cout << "\n\nBin:  " << i << endl;
    
    SetUpSumsAndCov(g_Pup_phys, g_Pdown_phys, i);
    InitializePars(gMinuit, vstart, step, ierflg);
    gMinuit->mnexcm(minimizer, arglist ,2,ierflg);
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    CheckFit(icstat, whichTarg, i);
  
    gMinuit->GetParameter(5, A_phys[i], eA_phys[i]);
    gMinuit->GetParameter(0, L_Pup_phys[i], eL_Pup_phys[i]);
    gMinuit->GetParameter(1, L_Pdown_phys[i], eL_Pdown_phys[i]);
    gMinuit->GetParameter(2, a1_phys[i], ea1_phys[i]);
    gMinuit->GetParameter(3, a2_phys[i], ea2_phys[i]);
    gMinuit->GetParameter(4, a3_phys[i], ea3_phys[i]);

    CheckPhysicality(L_Pup_phys[i], L_Pdown_phys[i],
		     a1_phys[i], a2_phys[i], a3_phys[i],
		     A_phys[i]);
    
    CorrectForP(A_phys[i], eA_phys[i],
		A_phys_corr[i], eA_phys_corr[i], g_pol_phys, i);

    if (minimizer=="MINOS"){
      gMinuit->mnerrs(5, eA_phys_up[i], eA_phys_down[i], eparab, gcc);
      eA_phys_down[i] *= -1;
      
      CorrectForP(A_phys[i], eA_phys_up[i], A_phys_corr[i], eA_phys_up_corr[i],
		  g_pol_phys, i);
      CorrectForP(A_phys[i], eA_phys_down[i], A_phys_corr[i],
		  eA_phys_down_corr[i], g_pol_phys, i);
    }

    ExpectedErrorA(eA_phys_expected[i],
		   g_Pup_phys[0], g_Pdown_phys[0], g_pol_phys, i);

    
  }

  //Draw contours
  TCanvas* cCont_int = new TCanvas(); 
  DrawContour(cCont_int);
  TCanvas *cCont_phy[nBins];
  for (Int_t i=0; i<nBins; i++) {
    cCont_phy[i] = new TCanvas();
    DrawContour(cCont_phy[i], Form("%i", i));  
  }

  //Graphs
  TGraphErrors *g_A_int =
    new TGraphErrors(1, &(xval[0]), &A_int_corr, ex, &eA_int_corr);
  TGraphErrors *g_L_Pup_int =
    new TGraphErrors(1, &(xval[0]), &L_Pup_int, ex, &eL_Pup_int);
  TGraphErrors *g_L_Pdown_int =
    new TGraphErrors(1, &(xval[0]), &L_Pdown_int, ex, &eL_Pdown_int);
  TGraphErrors *g_a1_int =
    new TGraphErrors(1, &(xval[0]), &a1_int, ex, &ea1_int);
  TGraphErrors *g_a2_int =
    new TGraphErrors(1, &(xval[0]), &a2_int, ex, &ea2_int);
  TGraphErrors *g_a3_int =
    new TGraphErrors(1, &(xval[0]), &a3_int, ex, &ea3_int);
  TGraphErrors *g_A_int_expected =
    new TGraphErrors(1, &(xval[0]), &A_int_corr, ex, &eA_int_expected);
  
  TGraphErrors *g_A_phys =
    new TGraphErrors(nBins, xval, A_phys_corr, ex, eA_phys_corr);
  TGraphErrors *g_L_Pup_phys =
    new TGraphErrors(nBins, xval, L_Pup_phys, ex, eL_Pup_phys);
  TGraphErrors *g_L_Pdown_phys =
    new TGraphErrors(nBins, xval, L_Pdown_phys, ex, eL_Pdown_phys);
  TGraphErrors *g_a1_phys =
    new TGraphErrors(nBins, xval, a1_phys, ex, ea1_phys);
  TGraphErrors *g_a2_phys =
    new TGraphErrors(nBins, xval, a2_phys, ex, ea2_phys);
  TGraphErrors *g_a3_phys =
    new TGraphErrors(nBins, xval, a3_phys, ex, ea3_phys);
  TGraphErrors *g_A_phys_expected =
    new TGraphErrors(nBins, xval, A_phys_corr, ex, eA_phys_expected);
  
  SetUp(g_A_int); SetUp(g_L_Pup_int); SetUp(g_L_Pdown_int);
  SetUp(g_a1_int); SetUp(g_a2_int); SetUp(g_a3_int);
  SetUp(g_A_int_expected);
  
  SetUp(g_A_phys); SetUp(g_L_Pup_phys); SetUp(g_L_Pdown_phys);
  SetUp(g_a1_phys); SetUp(g_a2_phys); SetUp(g_a3_phys);
  SetUp(g_A_phys_expected);

  TGraphAsymmErrors *g_A_int_minos =
    new TGraphAsymmErrors(1, &(xval[0]), &A_int_corr, ex, ex,
			  &eA_int_down_corr, &eA_int_up_corr);

  TGraphAsymmErrors *g_A_phys_minos =
    new TGraphAsymmErrors(nBins, xval, A_phys_corr, ex, ex,
			  eA_phys_down_corr, eA_phys_up_corr);
  SetUp(g_A_int_minos); SetUp(g_A_phys_minos);

  TString parTitles[] = {"A_{N}", "L Pup", "L Pdown", "a1", "a2", "a3"};
  TGraphErrors *g_phys[] =
    {g_A_phys, g_L_Pup_phys, g_L_Pdown_phys,
     g_a1_phys, g_a2_phys, g_a3_phys};
  TGraphErrors *g_int[] =
    {g_A_int, g_L_Pup_int, g_L_Pdown_int,
     g_a1_int, g_a2_int, g_a3_int};

  TCanvas* cphys = new TCanvas(); cphys->Divide(3, 2);
  TCanvas* cint = new TCanvas(); cint->Divide(3, 2);
  for (Int_t i=0; i<6; i++) {
    cphys->cd(i+1);
    g_phys[i]->Draw("AP");
    g_phys[i]->SetTitle(Form("%stream %s", whichTarg.Data(),
			     parTitles[i].Data()));

    cint->cd(i+1);
    g_int[i]->Draw("AP");
    g_int[i]->SetTitle(Form("%stream %s", whichTarg.Data(),
			    parTitles[i].Data()));
  }

  cphys->cd(1);
  OffSetPercent(g_A_phys_expected, 0.05);
  g_A_phys_expected->Draw("Psame");
  g_A_phys_expected->SetMarkerColor(kRed);

  if (minimizer=="MINOS"){
    OffSet(g_A_phys_minos);
    g_A_phys_minos->Draw("Psame");
    g_A_phys_minos->SetMarkerColor(kBlue);
  }

  cint->cd(1);
  OffSetPercent(g_A_int_expected, 0.05); 
  g_A_int_expected->Draw("Psame");
  g_A_int_expected->SetMarkerColor(kRed);

  if (minimizer == "MINOS"){
    OffSet(g_A_int_minos);
    g_A_int_minos->Draw("Psame");
    g_A_int_minos->SetMarkerColor(kBlue);
  }
  
  //Write output
  TString fOutput =
    Form("%s/Weight/weight_%s%i_%s_%s%s_%s_%s.root", thisDirPath.Data(),
	 physBinned.Data(), nBins, Mtype.Data(), production.Data(),
	 period.Data(), minimizer.Data(), whichTarg.Data());
  
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    
    g_A_int->Write("g_A_int");
    g_L_Pup_int->Write("g_L_Pup_int");
    g_L_Pdown_int->Write("g_L_Pdown_int");
    g_a1_int->Write("g_a1_int");
    g_a2_int->Write("g_a2_int");
    g_a3_int->Write("g_a3_int");
    g_A_int_expected->Write("g_A_int_expected");
    
    g_A_phys->Write("g_A_phys");
    g_L_Pup_phys->Write("g_L_Pup_phys");
    g_L_Pdown_phys->Write("g_L_Pdown_phys");
    g_a1_phys->Write("g_a1_phys");
    g_a2_phys->Write("g_a2_phys");
    g_a3_phys->Write("g_a3_phys");
    g_A_phys_expected->Write("g_A_phys_expected");

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
  //gMinuit->mnparm(4, "a5", vstart[4], step[4], -1,1,ierflg); cout << "par4 limits" << endl;//cleanup
  gMinuit->mnparm(4, "a5", vstart[4], step[4], 0,0,ierflg);  //cleanup
  gMinuit->mnparm(5, "a6", vstart[5], step[5], 0,0,ierflg);
}


void CorrectForP(Double_t A_uncorr, Double_t eA_uncorr,
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
  /*if ((a3 > 1) || (a3 < -1)) {
    cout << "a3 unphysical" << endl;
    error =1;
    }//*///cleanup
  cout << "par4 limits not checked" << endl;//cleanup

  if ((A > 1) || (A < -1)) {
    cout << "A unphysical" << endl;
    error =1;
  }

  if (error == 1){
    exit(EXIT_FAILURE);
  }
}


Int_t DrawContour(TCanvas *c, TString whichBin="int"){
  c->Divide(3, 2);
  
  Int_t status = 0;
  for (Int_t i=0; i<5; i++) {
    (i<2) ? c->cd(i+1) : c->cd(i+2);
    cout << "\nParameter correlation:  " << i << endl;
    
    //gMinuit->SetErrorDef(4);
    gMinuit->SetErrorDef(0.04);
    //gMinuit->SetErrorDef(1);
    TGraph *g_cont_2sig = (TGraph*)gMinuit->Contour(80, 5, i);
    cout << "\nContour status 2sigma:  " << gMinuit->GetStatus() << endl;
    status += gMinuit->GetStatus();
    if (!status){
      SetUp(g_cont_2sig);
      g_cont_2sig->Draw("alf");
      g_cont_2sig->SetFillColor(38);
      g_cont_2sig->SetTitle(Form("%s:  A vs. par_%i", whichBin.Data(), i));
    }
    else{ return status; }

    //gMinuit->SetErrorDef(1);
    //gMinuit->SetErrorDef(0.25);
    gMinuit->SetErrorDef(.01);
    TGraph *g_cont = (TGraph*)gMinuit->Contour(80, 5, i);
    cout << "\nContour status 1sigma:  " << gMinuit->GetStatus() << endl;
    status += gMinuit->GetStatus();
    if (!status){
      SetUp(g_cont);
      g_cont->Draw("lf");
      g_cont->SetFillColor(42);
    }
    else{ return status; }
  }

  return status;
}
