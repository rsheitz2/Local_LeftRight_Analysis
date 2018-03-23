#include "common.h"
#include "functions.h"
#include "setup.h"

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "" << endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t wflag=0, Qflag=0, fflag=0;
  Int_t c;
  TString fname = "", outFile = "";
  
  while ((c = getopt (argc, argv, "wf:Q:")) != -1) {
    switch (c) {
    case 'w':
      wflag = 1;
      break;
    case 'Q':
      Qflag = 1;
      outFile += optarg;
      break;
    case 'f':
      fflag = 1;
      fname += optarg;
      cout << fname << endl;
      break;
    case '?':
      if (optopt == 'u')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Unknown option character `\\x%x'.\n",
		 optopt);
      return 1;
    default:
      abort ();
    }
  }
  if (!fflag) {
    cout << "Please enter an input file" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    TFile *f1 = TFile::Open(fname);
    if(!f1){
      cout << fname << " does not exist" <<endl;
      exit(EXIT_FAILURE);
    }
    f1->Close();
  }
  // }}}
  
  //Opening data files/setting up tree
  ///////////////
  // {{{
  TFile *fdata = TFile::Open(fname);
  TTree *tree = (TTree*)fdata->Get("pT_Weighted");
  
  ////Binned
  //////////
  TVectorD tv_xN_bounds( *( (TVectorD*)fdata->Get("tv_xN_bounds") ) );
  TVectorD tv_xPi_bounds( *( (TVectorD*)fdata->Get("tv_xPi_bounds") ) );
  TVectorD tv_xF_bounds( *( (TVectorD*)fdata->Get("tv_xF_bounds") ) );
  TVectorD tv_pT_bounds( *( (TVectorD*)fdata->Get("tv_pT_bounds") ) );
  TVectorD tv_M_bounds( *( (TVectorD*)fdata->Get("tv_M_bounds") ) );
  Int_t nBounds = tv_xN_bounds.GetNoElements();
  Double_t xN_bounds[nBounds], xPi_bounds[nBounds], xF_bounds[nBounds];
  Double_t pT_bounds[nBounds], M_bounds[nBounds];
  //Double_t pT_bounds[nBounds];
  for (Int_t i=0; i<nBounds; i++) {
    xN_bounds[i] = tv_xN_bounds[i];
    xPi_bounds[i] = tv_xPi_bounds[i];
    xF_bounds[i] = tv_xF_bounds[i];
    pT_bounds[i] = tv_pT_bounds[i];
    M_bounds[i] = tv_M_bounds[i];
  }
  //Double_t M_bounds[] = {0.0, 0.06, 0.08, 0.10, 0.11, 0.12, 0.15, 0.35};
  
  //Vertex specific
  Double_t vx_z, Spin;
  Double_t vx_zVar, vx_xVar, vx_yVar;
  Int_t targetPosition;
  //Drell-Yan Angles
  Double_t PhiS_simple, Theta_CS, Gen_PhiS_simple, vOpenAngle;
  //Virtual Photon
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
  Double_t gen_vPhoton_X, gen_vPhoton_Y, gen_vPhoton_Z, gen_vPhoton_E;
  //Event
  Int_t trigMask;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;

  Int_t Errors = 0;
  //Vertex specific
  Errors += tree->SetBranchAddress("vx_z", &vx_z);
  Errors += tree->SetBranchAddress("Spin", &Spin); //UpStream=1, DownStream=-1, neiter=-10
  Errors += tree->SetBranchAddress("vx_zVar", &vx_zVar);
  Errors += tree->SetBranchAddress("vx_xVar", &vx_xVar);
  Errors += tree->SetBranchAddress("vx_yVar", &vx_yVar);
  Errors += tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  Errors += tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  Errors += tree->SetBranchAddress("Theta_CS", &Theta_CS);
  Errors += tree->SetBranchAddress("Gen_PhiS_simple", &Gen_PhiS_simple);
  Errors += tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Virtual Photon
  Errors += tree->SetBranchAddress("vPhoton_X", &vPhoton_X);
  Errors += tree->SetBranchAddress("vPhoton_Y", &vPhoton_Y);
  Errors += tree->SetBranchAddress("vPhoton_Z", &vPhoton_Z);
  Errors += tree->SetBranchAddress("vPhoton_E", &vPhoton_E);
  Errors += tree->SetBranchAddress("gen_vPhoton_X", &gen_vPhoton_X);
  Errors += tree->SetBranchAddress("gen_vPhoton_Y", &gen_vPhoton_Y);
  Errors += tree->SetBranchAddress("gen_vPhoton_Z", &gen_vPhoton_Z);
  Errors += tree->SetBranchAddress("gen_vPhoton_E", &gen_vPhoton_E);
  //Event
  Errors += tree->SetBranchAddress("trigMask", &trigMask);
  //DY-variables
  Errors += tree->SetBranchAddress("x_beam", &x_beam);
  Errors += tree->SetBranchAddress("x_target", &x_target);
  Errors += tree->SetBranchAddress("x_feynman", &x_feynman);
  Errors += tree->SetBranchAddress("q_transverse", &q_transverse);
  Errors += tree->SetBranchAddress("Mmumu", &Mmumu);

  if(Errors) {
   cout << " " << endl;
   cout << "Errors in Setting tree branch addresses" << endl;
   cout << "Probably some variables are not defined in your file" << endl;
   cout << " " << endl;
   exit(EXIT_FAILURE);
  }
  // }}}
  
  //Kinematic Counting Setup
  ///////////////
  // {{{
  Int_t nBins = nBounds-1;
  unsigned long long xN_Left_UpStream[nBins], xN_Right_UpStream[nBins];
  unsigned long long xN_Left_DownStream[nBins], xN_Right_DownStream[nBins];
  unsigned long long xPi_Left_UpStream[nBins], xPi_Right_UpStream[nBins];
  unsigned long long xPi_Left_DownStream[nBins], xPi_Right_DownStream[nBins];
  unsigned long long xF_Left_UpStream[nBins], xF_Right_UpStream[nBins];
  unsigned long long xF_Left_DownStream[nBins], xF_Right_DownStream[nBins];
  unsigned long long pT_Left_UpStream[nBins], pT_Right_UpStream[nBins];
  unsigned long long pT_Left_DownStream[nBins], pT_Right_DownStream[nBins];
  unsigned long long M_Left_UpStream[nBins], M_Right_UpStream[nBins];
  unsigned long long M_Left_DownStream[nBins], M_Right_DownStream[nBins];
  
  for (Int_t i=0; i<nBins; i++) {
    xN_Left_UpStream[i]=0; xN_Right_UpStream[i]=0;
    xN_Left_DownStream[i]=0; xN_Right_DownStream[i]=0;
    xPi_Left_UpStream[i]=0; xPi_Right_UpStream[i]=0;
    xPi_Left_DownStream[i]=0; xPi_Right_DownStream[i]=0;
    xF_Left_UpStream[i]=0; xF_Right_UpStream[i]=0;
    xF_Left_DownStream[i]=0; xF_Right_DownStream[i]=0;
    pT_Left_UpStream[i]=0; pT_Right_UpStream[i]=0;
    pT_Left_DownStream[i]=0; pT_Right_DownStream[i]=0;
    M_Left_UpStream[i]=0; M_Right_UpStream[i]=0;
    M_Left_DownStream[i]=0; M_Right_DownStream[i]=0;
  }

  const Int_t num_1Dh = 1;
  TH1D *h1[num_1Dh];
  //h1[0] = new TH1D("h_vPhoton", "h_vPhoton", 100, -TMath::Pi(), TMath::Pi() );
  h1[0] = new TH1D("h_vPhoton", "h_vPhoton", 100, -TMath::Pi(), 2*TMath::Pi() );
  // }}}
  
  //Tree loop
  Bool_t first = true;
  Int_t tree_entries = tree->GetEntries();//Tree Loop
  cout << "Number of entries in tree: " << tree_entries << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    //cout << "Debug mode" << endl; for (Int_t ev=0; ev<1000; ev++) {
    tree->GetEntry(ev, 0);
    
    if (first || ev==tree_entries-1){//first or last
      cout << " " << endl;
      cout << "Warning!!!!!!!!!!!!!!!" << endl;
      cout << "True spectrometer left/right (no spin influcence)" << endl;
      cout << "This is always the case for MC data" << endl;
      cout << " " << endl;
      cout << "!!!!!!!!!!!!!!!" << endl;
      //cout << "Gen PhiS Simple currently in use" << endl;
      //cout << "Gen photon tlorentz vector currently in use" << endl;
      //cout << "vPhoton tlorentz vector currently in use" << endl;
      cout << " " << endl;
      //cout << "Trig Mask settings:" << endl;
      //cout << "Only Last-Outer" << endl;
      //cout << "Only Last-Last" << endl;
      //cout << " " << endl;
      //cout << "vOpenAngle < 0.15" << endl;
      //cout << " " << endl;
      //cout << "M bins are binned in vOpenAngle" << endl;
      //cout << " " << endl;
      //cout << "Left = Top, Right = Bottom" << endl;
      //cout << " " << endl;
      
      first = false;
    }

    //Trig Mask
    //524 == Last-Outer
    //768 == Last-Last
    //780 == Last-Last && Last-Outer
    //if (trigMask == 768 || trigMask == 780) continue; //Only Last-Outer
    //if (trigMask == 524 || trigMask == 780) continue; //Only Last-Last

    //if (vOpenAngle <0.15) continue;

    //Choose Left/Right
    Double_t phi_photon_lab = ShiftPhiSimple(PhiS_simple);
    //Double_t phi_photon_lab = ShiftPhiSimple(Gen_PhiS_simple);
    
    TLorentzVector vPhoton(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);
    TLorentzVector gen_vPhoton(gen_vPhoton_X, gen_vPhoton_Y, gen_vPhoton_Z,
			       gen_vPhoton_E);
    //Double_t phi_photon_lab = gen_vPhoton.Phi();
    //Double_t phi_photon_lab = vPhoton.Phi();
    //if (phi_photon_lab < -TMath::Pi()/2) phi_photon_lab += 2*TMath::Pi();
    h1[0]->Fill(phi_photon_lab);
    
    Bool_t Left=false, Right=false;
    if (phi_photon_lab < TMath::Pi()/2 && phi_photon_lab > -TMath::Pi()/2> 0){
      //if (phi_photon_lab < TMath::Pi() && phi_photon_lab > 0.0){//Top 
      Left = true;
    }
    else if (phi_photon_lab < 3*TMath::Pi()/2 && phi_photon_lab>TMath::Pi()/2){
      //else if (phi_photon_lab > TMath::Pi() || phi_photon_lab < 0.0){//Bottom
      Right = true;
    }
    else {
      cout << "No Left or Right choosen" << endl;
      cout << "phi_photon_lab = " << phi_photon_lab << endl;
      cout << " " << endl;
    }

    
    //Bin data
    if (targetPosition == 0) {//UpStream target
      
      if (Left){//Left
	BinDataCounts(xN_Left_UpStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Left_UpStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Left_UpStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Left_UpStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Left_UpStream, nBins, Mmumu, M_bounds);
	//BinDataCounts(M_Left_UpStream, nBins, vOpenAngle, M_bounds);

      }//End Left
      else if (Right){//Right
	BinDataCounts(xN_Right_UpStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Right_UpStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Right_UpStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Right_UpStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Right_UpStream, nBins, Mmumu, M_bounds);
	//BinDataCounts(M_Right_UpStream, nBins, vOpenAngle, M_bounds);
      }//End Right

    }//End UpStream target
    else if (targetPosition == 1) {//DownStream target
      
      if (Left){//Left
	BinDataCounts(xN_Left_DownStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Left_DownStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Left_DownStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Left_DownStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Left_DownStream, nBins, Mmumu, M_bounds);
	//BinDataCounts(M_Left_DownStream, nBins, vOpenAngle, M_bounds);
      }//End Left
      else if (Right){//Right
	BinDataCounts(xN_Right_DownStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Right_DownStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Right_DownStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Right_DownStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Right_DownStream, nBins, Mmumu, M_bounds);
	//BinDataCounts(M_Right_DownStream, nBins, vOpenAngle, M_bounds);
      }//End Right

    }//End DownStream target

  }//End tree loop

  //Draw Basic Dist
  ///////////////
  // {{{
  TCanvas* cDist = new TCanvas();
  h1[0]->Draw();
  // }}}

  //Asymmetries
  ///////////////
  // {{{
  Double_t xN_Asym_UpStream[nBins], xN_Asym_DownStream[nBins], xN_Asym[nBins];
  Double_t e_xN_Asym_UpStream[nBins], e_xN_Asym_DownStream[nBins];
  Double_t e_xN_Asym[nBins];
  BinnedLeftRight(xN_Left_UpStream, xN_Right_UpStream, xN_Left_DownStream,
		  xN_Right_DownStream, xN_Asym_UpStream, xN_Asym_DownStream,
		  xN_Asym, e_xN_Asym_UpStream, e_xN_Asym_DownStream, e_xN_Asym,
		  nBins);
  
  Double_t xPi_Asym_UpStream[nBins],xPi_Asym_DownStream[nBins], xPi_Asym[nBins];
  Double_t e_xPi_Asym_UpStream[nBins], e_xPi_Asym_DownStream[nBins];
  Double_t e_xPi_Asym[nBins];
  BinnedLeftRight(xPi_Left_UpStream, xPi_Right_UpStream, xPi_Left_DownStream,
		  xPi_Right_DownStream, xPi_Asym_UpStream, xPi_Asym_DownStream,
		  xPi_Asym, e_xPi_Asym_UpStream, e_xPi_Asym_DownStream,
		  e_xPi_Asym,
		  nBins);
  
  Double_t xF_Asym_UpStream[nBins], xF_Asym_DownStream[nBins], xF_Asym[nBins];
  Double_t e_xF_Asym_UpStream[nBins], e_xF_Asym_DownStream[nBins];
  Double_t e_xF_Asym[nBins];
  BinnedLeftRight(xF_Left_UpStream, xF_Right_UpStream, xF_Left_DownStream,
		  xF_Right_DownStream, xF_Asym_UpStream, xF_Asym_DownStream,
		  xF_Asym, e_xF_Asym_UpStream, e_xF_Asym_DownStream, e_xF_Asym,
		  nBins);

  Double_t pT_Asym_UpStream[nBins], pT_Asym_DownStream[nBins], pT_Asym[nBins];
  Double_t e_pT_Asym_UpStream[nBins], e_pT_Asym_DownStream[nBins];
  Double_t e_pT_Asym[nBins];
  BinnedLeftRight(pT_Left_UpStream, pT_Right_UpStream, pT_Left_DownStream,
  		  pT_Right_DownStream, pT_Asym_UpStream, pT_Asym_DownStream,
  		  pT_Asym, e_pT_Asym_UpStream, e_pT_Asym_DownStream, e_pT_Asym,
  		  nBins);

  Double_t M_Asym_UpStream[nBins], M_Asym_DownStream[nBins], M_Asym[nBins];
  Double_t e_M_Asym_UpStream[nBins], e_M_Asym_DownStream[nBins];
  Double_t e_M_Asym[nBins];
  BinnedLeftRight(M_Left_UpStream, M_Right_UpStream, M_Left_DownStream,
		  M_Right_DownStream, M_Asym_UpStream, M_Asym_DownStream,
		  M_Asym, e_M_Asym_UpStream, e_M_Asym_DownStream, e_M_Asym,
		  nBins);
  // }}}
  
  //TGraphs
  ///////////////
  // {{{
  TVectorD tv_xN_xval( *( (TVectorD*)fdata->Get("tv_xN_xval") ) );
  TVectorD tv_xPi_xval( *( (TVectorD*)fdata->Get("tv_xPi_xval") ) );
  TVectorD tv_xF_xval( *( (TVectorD*)fdata->Get("tv_xF_xval") ) );
  TVectorD tv_pT_xval( *( (TVectorD*)fdata->Get("tv_pT_xval") ) );
  TVectorD tv_M_xval( *( (TVectorD*)fdata->Get("tv_M_xval") ) );
  Double_t xval_xN[nBins], xval_xPi[nBins], xval_xF[nBins];
  Double_t xval_pT[nBins], xval_M[nBins];
  //Double_t xval_pT[nBins];
  Double_t ex[nBins];
  for (Int_t i=0; i<nBins; i++) {
    xval_xN[i] = tv_xN_xval[i]; xval_xPi[i] = tv_xPi_xval[i];
    xval_xF[i]=tv_xF_xval[i]; xval_pT[i]=tv_pT_xval[i], xval_M[i]=tv_M_xval[i];
    ex[i] = 0.0;
  }

  //Double_t xval_M[] = {0.03, 0.07, 0.1, 0.115, 0.13, 0.155, 0.2};
  
  TGraphErrors* gr_xN_LR_UpStream = new TGraphErrors(nBins, xval_xN,
						     xN_Asym_UpStream, ex,
						     e_xN_Asym_UpStream);
  TGraphErrors* gr_xN_LR_DownStream = new TGraphErrors(nBins, xval_xN,
						       xN_Asym_DownStream, ex,
						       e_xN_Asym_DownStream);
  TGraphErrors* gr_xN_LR = new TGraphErrors(nBins, xval_xN,
					    xN_Asym, ex,
					    e_xN_Asym);
  
  TGraphErrors* gr_xPi_LR_UpStream = new TGraphErrors(nBins, xval_xPi,
						      xPi_Asym_UpStream, ex,
						      e_xPi_Asym_UpStream);
  TGraphErrors* gr_xPi_LR_DownStream = new TGraphErrors(nBins, xval_xPi,
							xPi_Asym_DownStream, ex,
							e_xPi_Asym_DownStream);
  TGraphErrors* gr_xPi_LR = new TGraphErrors(nBins, xval_xPi,
					     xPi_Asym, ex,
					     e_xPi_Asym);
  
  TGraphErrors* gr_xF_LR_UpStream = new TGraphErrors(nBins, xval_xF,
						     xF_Asym_UpStream, ex,
						     e_xF_Asym_UpStream);
  TGraphErrors* gr_xF_LR_DownStream = new TGraphErrors(nBins, xval_xF,
						       xF_Asym_DownStream, ex,
						       e_xF_Asym_DownStream);
  TGraphErrors* gr_xF_LR = new TGraphErrors(nBins, xval_xF,
					    xF_Asym, ex,
					    e_xF_Asym);
  
  TGraphErrors* gr_pT_LR_UpStream = new TGraphErrors(nBins, xval_pT,
  						     pT_Asym_UpStream, ex,
  						     e_pT_Asym_UpStream);
  TGraphErrors* gr_pT_LR_DownStream = new TGraphErrors(nBins, xval_pT,
  						       pT_Asym_DownStream, ex,
  						       e_pT_Asym_DownStream);
  TGraphErrors* gr_pT_LR = new TGraphErrors(nBins, xval_pT,
  					    pT_Asym, ex,
  					    e_pT_Asym);
  
  TGraphErrors* gr_M_LR_UpStream = new TGraphErrors(nBins, xval_M,
						    M_Asym_UpStream, ex,
						    e_M_Asym_UpStream);
  TGraphErrors* gr_M_LR_DownStream = new TGraphErrors(nBins, xval_M,
						      M_Asym_DownStream, ex,
						      e_M_Asym_DownStream);
  TGraphErrors* gr_M_LR = new TGraphErrors(nBins, xval_M,
					   M_Asym, ex,
					   e_M_Asym);

  const Int_t nGr = 3;
  TGraphErrors* gr_xN[nGr] = {gr_xN_LR_UpStream, gr_xN_LR_DownStream, gr_xN_LR};
  TGraphErrors* gr_xPi[nGr] = {gr_xPi_LR_UpStream, gr_xPi_LR_DownStream,
			       gr_xPi_LR};
  TGraphErrors* gr_xF[nGr] = {gr_xF_LR_UpStream, gr_xF_LR_DownStream, gr_xF_LR};
  TGraphErrors* gr_pT[nGr] = {gr_pT_LR_UpStream, gr_pT_LR_DownStream, gr_pT_LR};
  TGraphErrors* gr_M[nGr] = {gr_M_LR_UpStream, gr_M_LR_DownStream, gr_M_LR};
  // }}}

  TString Titles[nGr] = {"UpStream", "DownStream", "Combined Targets"};
  TString xTitles[] = {"xN", "xPi", "xF", "pT", "M"};
  
  ////////////////
  //Draw and pretty up graphs
  ////////////////
  //Binned physics values
  TCanvas* cPhys = new TCanvas();
  //cPhys->Divide(2);
  cPhys->Divide(2, 5);
  Double_t ySym = 10.0;
  //Double_t ymin = -1.0*ySym, ymax = ySym;
  //Double_t ymin = -1.0, ymax = 1.0;
  for (Int_t i=0; i<3; i++) {
    SetupTGraph(gr_xN[i]	, Titles[i], xTitles[0], ySym);
    SetupTGraph(gr_xPi[i]	, Titles[i], xTitles[1], ySym);
    SetupTGraph(gr_xF[i]	, Titles[i], xTitles[2], ySym);
    SetupTGraph(gr_pT[i]	, Titles[i], xTitles[3], ySym);
    SetupTGraph(gr_M[i]		, Titles[i], xTitles[4], ySym);
  }

  TLine *l_xN[nGr];
  for (Int_t i=0, ipad=1; i<2; i++, ipad++) {
    cPhys->cd(ipad);
    gr_xN[i]->Draw("AP");
    
    Double_t xmin = gr_xN[i]->GetXaxis()->GetXmin();
    Double_t xmax = gr_xN[i]->GetXaxis()->GetXmax();
    l_xN[i] = new TLine(xmin,0.0,xmax,0.0);
    SetupTLine(l_xN[i]);
    l_xN[i]->Draw("same");  
  }

  TLine *l_xPi[nGr];
  for (Int_t i=0, ipad=3; i<2; i++, ipad++) {
    cPhys->cd(ipad);
    gr_xPi[i]->Draw("AP");
    
    Double_t xmin = gr_xPi[i]->GetXaxis()->GetXmin();
    Double_t xmax = gr_xPi[i]->GetXaxis()->GetXmax();
    l_xPi[i] = new TLine(xmin,0.0,xmax,0.0);
    SetupTLine(l_xPi[i]);
    l_xPi[i]->Draw("same");  
  }

  TLine *l_xF[nGr];
  for (Int_t i=0, ipad=5; i<2; i++, ipad++) {
    cPhys->cd(ipad);
    gr_xF[i]->Draw("AP");
    
    Double_t xmin = gr_xF[i]->GetXaxis()->GetXmin();
    Double_t xmax = gr_xF[i]->GetXaxis()->GetXmax();
    l_xF[i] = new TLine(xmin,0.0,xmax,0.0);
    SetupTLine(l_xF[i]);
    l_xF[i]->Draw("same");  
  }

  TLine *l_pT[nGr];
  for (Int_t i=0, ipad=7; i<2; i++, ipad++) {
    cPhys->cd(ipad);
    gr_pT[i]->Draw("AP");
    
    Double_t xmin = gr_pT[i]->GetXaxis()->GetXmin();
    Double_t xmax = gr_pT[i]->GetXaxis()->GetXmax();
    l_pT[i] = new TLine(xmin,0.0,xmax,0.0);
    SetupTLine(l_pT[i]);
    l_pT[i]->Draw("same");  
    }

  TLine *l_M[nGr];
  for (Int_t i=0, ipad=9; i<2; i++, ipad++) {
  //for (Int_t i=0, ipad=1; i<2; i++, ipad++) {
    cPhys->cd(ipad);
    gr_M[i]->Draw("AP");
    
    Double_t xmin = gr_M[i]->GetXaxis()->GetXmin();
    Double_t xmax = gr_M[i]->GetXaxis()->GetXmax();
    l_M[i] = new TLine(xmin,0.0,xmax,0.0);
    SetupTLine(l_M[i]);
    l_M[i]->Draw("same");  
  }
    
  //Write output
  ///////////////
  // {{{
  if (Qflag || wflag){
    TFile* fout = (Qflag ? new TFile(outFile, "RECREATE")
		   : new TFile("Output.root", "RECREATE") );

    cout << " " << endl;
    if (outFile != "")  cout << "File: " << outFile << " was written" << endl;
    else cout << "File: Output.root was written" << endl;
    cout << " " << endl;
    
    
    gr_xN_LR_UpStream->Write("gr_xN_LR_UpStream");
    gr_xN_LR_DownStream->Write("gr_xN_LR_DownStream");
    gr_xN_LR->Write("gr_xN_LR");

    gr_xPi_LR_UpStream->Write("gr_xPi_LR_UpStream");
    gr_xPi_LR_DownStream->Write("gr_xPi_LR_DownStream");
    gr_xPi_LR->Write("gr_xPi_LR");

    gr_xF_LR_UpStream->Write("gr_xF_LR_UpStream");
    gr_xF_LR_DownStream->Write("gr_xF_LR_DownStream");
    gr_xF_LR->Write("gr_xF_LR");

    gr_pT_LR_UpStream->Write("gr_pT_LR_UpStream");
    gr_pT_LR_DownStream->Write("gr_pT_LR_DownStream");
    gr_pT_LR->Write("gr_pT_LR");

    gr_M_LR_UpStream->Write("gr_M_LR_UpStream");
    gr_M_LR_DownStream->Write("gr_M_LR_DownStream");
    gr_M_LR->Write("gr_M_LR");

    for (Int_t i=0; i<num_1Dh; i++) h1[i]->Write();
      
  
    fout->Close();
  }
  // }}}
    
  theApp.Run();//Needed to make root graphics work on C++
}
