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
    cout << "Option:  -B binVar	(Variable to bin in)" << endl;
    cout << "     (Current options: HM, xN, xPi, xF, pT, vxZ_upstream, ";
    cout << "vxZ_downstream)" << endl;
    cout << "     (Default option is HM)" << endl;
    cout << "Option   -B binVar (Variables that require a bin count)" << endl;
    cout << "     (Current options: rapidity, phi_vPh, theta_vPh, openAngle"
	 << endl;
    cout << "Option   -N bounds,   bounds=nBins+1"
	 << " (For varibles that require a bin count)"
	 << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t wflag=0, Qflag=0, fflag=0, Bflag=0, Nflag=0;
  Int_t c;
  TString fname = "", outFile = "", binVar="";
  Int_t NVar=0;
  
  while ((c = getopt (argc, argv, "wf:Q:B:N:")) != -1) {
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
    case 'B':
      Bflag = 1;
      binVar += optarg;
      break;
    case 'N':
      Nflag = 1;
      NVar = atoi(optarg);
      break;
    case '?':
      if (optopt == 'u')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'B')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'N')
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

  if (!Bflag){
    cout << " " << endl;
    cout << "Default binning variable = Invariant mass" << endl;
  }
  if (Nflag && !NVar){
    cout << " " << endl;
    cout << "Please enter and integer with -N option" << endl;
    exit(EXIT_FAILURE);
  }
  else if (Nflag && NVar < 2){
    cout << " " << endl;
    cout << "Please enter a number >= 2 for -N option" << endl;
    exit(EXIT_FAILURE);
  }
  if (!Nflag && (binVar!="HM"&&binVar!="xN"&&binVar!="xPi"&&binVar!="xF"&&
		 binVar!="pT"&&binVar!="vxZ_upstream"&&
		 binVar!="vxZ_downstream" ) ){
    cout << " " << endl;
    cout << "Please enter an -N option with -B" << binVar << endl;
    exit(EXIT_FAILURE);
  }
  else if(Nflag && binVar!="rapidity"&&binVar!="phi_vPh"&&binVar!="theta_vPh"&&
	  binVar!="openAngle"){
    cout << " " << endl;
    cout << "Please do not enter an -N option with -B" << binVar << endl;
    exit(EXIT_FAILURE);
  }
  // }}}
  
  //Opening data files
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
  TVectorD tv_rad_bounds( *( (TVectorD*)fdata->Get("tv_rad_bounds") ) );
  TVectorD tv_vxZ_upstream_bounds(*((TVectorD*)fdata->Get("tv_vxZ_upstream_bounds")));
  TVectorD tv_vxZ_downstream_bounds(*((TVectorD*)fdata->Get("tv_vxZ_downstream_bounds")));
  Int_t nBounds = (!Nflag) ? tv_xN_bounds.GetNoElements() : NVar;
  Double_t xN_bounds[nBounds], xPi_bounds[nBounds], xF_bounds[nBounds];
  Double_t pT_bounds[nBounds], M_bounds[nBounds];
  Double_t rad_bounds[nBounds];
  Double_t vxZ_upstream_bounds[nBounds], vxZ_downstream_bounds[nBounds];
  Double_t bounds[nBounds];
  for (Int_t i=0; i<nBounds; i++) {
    xN_bounds[i] = tv_xN_bounds[i];
    xPi_bounds[i] = tv_xPi_bounds[i];
    xF_bounds[i] = tv_xF_bounds[i];
    pT_bounds[i] = tv_pT_bounds[i];
    M_bounds[i] = tv_M_bounds[i];
    rad_bounds[i]=tv_rad_bounds[i];
    vxZ_upstream_bounds[i]=tv_vxZ_upstream_bounds[i];
    vxZ_downstream_bounds[i]=tv_vxZ_downstream_bounds[i];

    //Generic bounds
    if (Nflag) continue;//bounds made after tree loop 1
    if (!Bflag || binVar=="HM") bounds[i] = M_bounds[i];
    else if (binVar=="xN") bounds[i] = xN_bounds[i];
    else if (binVar=="xPi") bounds[i] = xPi_bounds[i];
    else if (binVar=="xF") bounds[i] = xF_bounds[i];
    else if (binVar=="pT") bounds[i] = pT_bounds[i];
    else if (binVar=="rad") bounds[i] = rad_bounds[i];
    else if (binVar=="vxZ_upstream") bounds[i] = vxZ_upstream_bounds[i];
    else if (binVar=="vxZ_downstream") bounds[i] = vxZ_downstream_bounds[i];
    else {cout << "Wrong binning variable option" << endl;
	exit(EXIT_FAILURE);}
  }

  //Dilution and polarization corrections
  TVectorD Dil_int( *( (TVectorD*)fdata->Get("Dil_int") ) );
  TVectorD Pol_int( *( (TVectorD*)fdata->Get("Pol_int") ) );

  ////Binned Amplitudes
  //Dilution and Polarization corrections
  TVectorD Dil_xN( *( (TVectorD*)fdata->Get("Dil_xN") ) );
  TVectorD Dil_xPi( *( (TVectorD*)fdata->Get("Dil_xPi") ) );
  TVectorD Dil_xF( *( (TVectorD*)fdata->Get("Dil_xF") ) );
  TVectorD Dil_pT( *( (TVectorD*)fdata->Get("Dil_pT") ) );
  TVectorD Dil_M( *( (TVectorD*)fdata->Get("Dil_M") ) );

  TVectorD Dil_xN_UpStream( *( (TVectorD*)fdata->Get("Dil_xN_UpStream") ) );
  TVectorD Dil_xPi_UpStream( *( (TVectorD*)fdata->Get("Dil_xPi_UpStream") ) );
  TVectorD Dil_xF_UpStream( *( (TVectorD*)fdata->Get("Dil_xF_UpStream") ) );
  TVectorD Dil_pT_UpStream( *( (TVectorD*)fdata->Get("Dil_pT_UpStream") ) );
  TVectorD Dil_M_UpStream( *( (TVectorD*)fdata->Get("Dil_M_UpStream") ) );

  TVectorD Dil_xN_UpStream_Up( *( (TVectorD*)fdata->Get("Dil_xN_UpStream_Up")));
  TVectorD Dil_xPi_UpStream_Up(*((TVectorD*)fdata->Get("Dil_xPi_UpStream_Up")));
  TVectorD Dil_xF_UpStream_Up( *( (TVectorD*)fdata->Get("Dil_xF_UpStream_Up")));
  TVectorD Dil_pT_UpStream_Up( *( (TVectorD*)fdata->Get("Dil_pT_UpStream_Up")));
  TVectorD Dil_M_UpStream_Up( *( (TVectorD*)fdata->Get("Dil_M_UpStream_Up") ));

  TVectorD Dil_xN_UpStream_Down( *( (TVectorD*)fdata->Get("Dil_xN_UpStream_Down")));
  TVectorD Dil_xPi_UpStream_Down(*((TVectorD*)fdata->Get("Dil_xPi_UpStream_Down")));
  TVectorD Dil_xF_UpStream_Down( *( (TVectorD*)fdata->Get("Dil_xF_UpStream_Down")));
  TVectorD Dil_pT_UpStream_Down( *( (TVectorD*)fdata->Get("Dil_pT_UpStream_Down")));
  TVectorD Dil_M_UpStream_Down( *( (TVectorD*)fdata->Get("Dil_M_UpStream_Down") ));

  TVectorD Dil_xN_DownStream( *( (TVectorD*)fdata->Get("Dil_xN_DownStream") ) );
  TVectorD Dil_xPi_DownStream( *( (TVectorD*)fdata->Get("Dil_xPi_DownStream")));
  TVectorD Dil_xF_DownStream( *( (TVectorD*)fdata->Get("Dil_xF_DownStream") ) );
  TVectorD Dil_pT_DownStream( *( (TVectorD*)fdata->Get("Dil_pT_DownStream") ) );
  TVectorD Dil_M_DownStream( *( (TVectorD*)fdata->Get("Dil_M_DownStream") ) );

  TVectorD Dil_xN_DownStream_Up( *( (TVectorD*)fdata->Get("Dil_xN_DownStream_Up")));
  TVectorD Dil_xPi_DownStream_Up(*((TVectorD*)fdata->Get("Dil_xPi_DownStream_Up")));
  TVectorD Dil_xF_DownStream_Up( *( (TVectorD*)fdata->Get("Dil_xF_DownStream_Up")));
  TVectorD Dil_pT_DownStream_Up( *( (TVectorD*)fdata->Get("Dil_pT_DownStream_Up")));
  TVectorD Dil_M_DownStream_Up( *( (TVectorD*)fdata->Get("Dil_M_DownStream_Up") ));

  TVectorD Dil_xN_DownStream_Down( *( (TVectorD*)fdata->Get("Dil_xN_DownStream_Down")));
  TVectorD Dil_xPi_DownStream_Down(*((TVectorD*)fdata->Get("Dil_xPi_DownStream_Down")));
  TVectorD Dil_xF_DownStream_Down( *( (TVectorD*)fdata->Get("Dil_xF_DownStream_Down")));
  TVectorD Dil_pT_DownStream_Down( *( (TVectorD*)fdata->Get("Dil_pT_DownStream_Down")));
  TVectorD Dil_M_DownStream_Down( *( (TVectorD*)fdata->Get("Dil_M_DownStream_Down") ));
  
  TVectorD Pol_xN( *( (TVectorD*)fdata->Get("Pol_xN") ) );
  TVectorD Pol_xPi( *( (TVectorD*)fdata->Get("Pol_xPi") ) );
  TVectorD Pol_xF( *( (TVectorD*)fdata->Get("Pol_xF") ) );
  TVectorD Pol_pT( *( (TVectorD*)fdata->Get("Pol_pT") ) );
  TVectorD Pol_M( *( (TVectorD*)fdata->Get("Pol_M") ) );

  TVectorD Pol_xN_UpStream( *( (TVectorD*)fdata->Get("Pol_xN_UpStream") ) );
  TVectorD Pol_xPi_UpStream( *( (TVectorD*)fdata->Get("Pol_xPi_UpStream") ) );
  TVectorD Pol_xF_UpStream( *( (TVectorD*)fdata->Get("Pol_xF_UpStream") ) );
  TVectorD Pol_pT_UpStream( *( (TVectorD*)fdata->Get("Pol_pT_UpStream") ) );
  TVectorD Pol_M_UpStream( *( (TVectorD*)fdata->Get("Pol_M_UpStream") ) );

  TVectorD Pol_xN_UpStream_Up( *( (TVectorD*)fdata->Get("Pol_xN_UpStream_Up") ) );
  TVectorD Pol_xPi_UpStream_Up( *( (TVectorD*)fdata->Get("Pol_xPi_UpStream_Up") ) );
  TVectorD Pol_xF_UpStream_Up( *( (TVectorD*)fdata->Get("Pol_xF_UpStream_Up") ) );
  TVectorD Pol_pT_UpStream_Up( *( (TVectorD*)fdata->Get("Pol_pT_UpStream_Up") ) );
  TVectorD Pol_M_UpStream_Up( *( (TVectorD*)fdata->Get("Pol_M_UpStream_Up") ) );

  TVectorD Pol_xN_UpStream_Down( *( (TVectorD*)fdata->Get("Pol_xN_UpStream_Down") ) );
  TVectorD Pol_xPi_UpStream_Down( *( (TVectorD*)fdata->Get("Pol_xPi_UpStream_Down") ) );
  TVectorD Pol_xF_UpStream_Down( *( (TVectorD*)fdata->Get("Pol_xF_UpStream_Down") ) );
  TVectorD Pol_pT_UpStream_Down( *( (TVectorD*)fdata->Get("Pol_pT_UpStream_Down") ) );
  TVectorD Pol_M_UpStream_Down( *( (TVectorD*)fdata->Get("Pol_M_UpStream_Down") ) );

  TVectorD Pol_xN_DownStream( *( (TVectorD*)fdata->Get("Pol_xN_DownStream") ) );
  TVectorD Pol_xPi_DownStream( *( (TVectorD*)fdata->Get("Pol_xPi_DownStream")));
  TVectorD Pol_xF_DownStream( *( (TVectorD*)fdata->Get("Pol_xF_DownStream") ) );
  TVectorD Pol_pT_DownStream( *( (TVectorD*)fdata->Get("Pol_pT_DownStream") ) );
  TVectorD Pol_M_DownStream( *( (TVectorD*)fdata->Get("Pol_M_DownStream") ) );

  TVectorD Pol_xN_DownStream_Up( *( (TVectorD*)fdata->Get("Pol_xN_DownStream_Up") ) );
  TVectorD Pol_xPi_DownStream_Up( *( (TVectorD*)fdata->Get("Pol_xPi_DownStream_Up") ) );
  TVectorD Pol_xF_DownStream_Up( *( (TVectorD*)fdata->Get("Pol_xF_DownStream_Up") ) );
  TVectorD Pol_pT_DownStream_Up( *( (TVectorD*)fdata->Get("Pol_pT_DownStream_Up") ) );
  TVectorD Pol_M_DownStream_Up( *( (TVectorD*)fdata->Get("Pol_M_DownStream_Up") ) );

  TVectorD Pol_xN_DownStream_Down( *( (TVectorD*)fdata->Get("Pol_xN_DownStream_Down") ) );
  TVectorD Pol_xPi_DownStream_Down( *( (TVectorD*)fdata->Get("Pol_xPi_DownStream_Down") ) );
  TVectorD Pol_xF_DownStream_Down( *( (TVectorD*)fdata->Get("Pol_xF_DownStream_Down") ) );
  TVectorD Pol_pT_DownStream_Down( *( (TVectorD*)fdata->Get("Pol_pT_DownStream_Down") ) );
  TVectorD Pol_M_DownStream_Down( *( (TVectorD*)fdata->Get("Pol_M_DownStream_Down") ) );
  // }}}

  //Getting tree
  // {{{
  Int_t errors = 0;
  //Vertex specific
  Double_t vx_z;
  Int_t targetPosition;
  //Drell-Yan Angles
  Double_t PhiS_simple, Theta_CS, rapidity, vOpenAngle;
  //Muons
  Double_t phi_traj1, phi_traj2;
  Double_t theta_traj1, theta_traj2;
  Double_t qP_traj1, qP_traj2;
  //Virtual Photon
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
  //Event
  Int_t trigMask, MasterTrigMask;
  //Target values
  Double_t Spin_0, Spin_1, Spin_2, Spin_3, Spin_4, Spin_5, Spin_6;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;
  //Generic binning
  Double_t *boundValue;
  
  //Vertex specific
  errors += tree->SetBranchAddress("vx_z", &vx_z);
  errors += tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  errors += tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  errors += tree->SetBranchAddress("Theta_CS", &Theta_CS);
  errors += tree->SetBranchAddress("rapidity", &rapidity);
  errors += tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Virtual Photon
  errors += tree->SetBranchAddress("vPhoton_X", &vPhoton_X);
  errors += tree->SetBranchAddress("vPhoton_Y", &vPhoton_Y);
  errors += tree->SetBranchAddress("vPhoton_Z", &vPhoton_Z);
  errors += tree->SetBranchAddress("vPhoton_E", &vPhoton_E);
  //Muons
  errors += tree->SetBranchAddress("phi_traj1", &phi_traj1);
  errors += tree->SetBranchAddress("phi_traj2", &phi_traj2);
  errors += tree->SetBranchAddress("theta_traj1", &theta_traj1);
  errors += tree->SetBranchAddress("theta_traj2", &theta_traj2);
  errors += tree->SetBranchAddress("qP_traj1", &qP_traj1);
  errors += tree->SetBranchAddress("qP_traj2", &qP_traj2);
  //Event
  errors += tree->SetBranchAddress("trigMask", &trigMask);
  errors += tree->SetBranchAddress("MasterTrigMask", &MasterTrigMask);
  //Target values
  errors += tree->SetBranchAddress("Spin_0", &Spin_0);
  errors += tree->SetBranchAddress("Spin_1", &Spin_1);
  errors += tree->SetBranchAddress("Spin_2", &Spin_2);
  errors += tree->SetBranchAddress("Spin_3", &Spin_3);
  errors += tree->SetBranchAddress("Spin_4", &Spin_4);
  errors += tree->SetBranchAddress("Spin_5", &Spin_5);
  errors += tree->SetBranchAddress("Spin_6", &Spin_6);
  //DY-variables
  errors += tree->SetBranchAddress("x_beam", &x_beam);
  errors += tree->SetBranchAddress("x_target", &x_target);
  errors += tree->SetBranchAddress("x_feynman", &x_feynman);
  errors += tree->SetBranchAddress("q_transverse", &q_transverse);
  errors += tree->SetBranchAddress("Mmumu", &Mmumu);
  //Generic bound value
  if (!Bflag || binVar=="HM") boundValue = &Mmumu;
  else if (binVar=="xN") boundValue = &x_target;
  else if (binVar=="xPi") boundValue = &x_beam;
  else if (binVar=="xF") boundValue = &x_feynman;
  else if (binVar=="pT") boundValue = &q_transverse;
  //else if (binVar=="rad") boundValue = &q_transverse; //done in loop
  else if (binVar=="vxZ_upstream"||binVar=="vxZ_downstream") boundValue = &vx_z;
  else if (binVar=="rapidity") boundValue = &rapidity;
  else if (binVar=="openAngle") boundValue = &vOpenAngle;

  if (errors){
    cout << " " << endl;
    cout << "Errors opening trees variables" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  // }}}

  //Sort bin boundries if -N option is given
  ///////////////
  // {{{
  Int_t tree_entries = tree->GetEntries();//Tree Loop 1
  if(Nflag) {
    vector<Double_t> sort_val;
    for (Int_t ev=0; ev<tree_entries; ev++) {
      tree->GetEntry(ev, 0);

      //General useful variables
      TLorentzVector vPhoton(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);
      if (binVar=="phi_vPh") *boundValue = vPhoton.Phi();
      else if (binVar=="theta_vPh") *boundValue = vPhoton.Theta();

      sort_val.push_back(*boundValue);
    }
    std::sort(sort_val.begin(), sort_val.end() );
    bounds[0] = sort_val.at(0);
    for (Int_t i=0; i<nBounds; i++) {
      bounds[i] = (i) ? sort_val.at( i*sort_val.size()/(1.0*nBounds-1.0) - 1 ) :
	sort_val.at( i*sort_val.size()/(1.0*nBounds-1.0) );
    }
  }
  // }}}
  
  //Kinematic Counting Setup
  ///////////////
  // {{{
  Int_t nBins = nBounds-1;
  
  //Bin filling
  TH1D *hPhi1_LL[nBins], *hPhi1_LO[nBins], *hPhi1_LL_LO[nBins];
  TH1D *hTheta1_LL[nBins], *hTheta1_LO[nBins], *hTheta1_LL_LO[nBins];
  TH1D *hP1_LL[nBins], *hP1_LO[nBins], *hP1_LL_LO[nBins];
  TH1D *hPhi2_LL[nBins], *hPhi2_LO[nBins], *hPhi2_LL_LO[nBins];
  TH1D *hTheta2_LL[nBins], *hTheta2_LO[nBins], *hTheta2_LL_LO[nBins];
  TH1D *hP2_LL[nBins], *hP2_LO[nBins], *hP2_LL_LO[nBins];
  for (Int_t i=0; i<nBins; i++) {

    //Particle 1
    hPhi1_LL[i] = new TH1D(Form("hPhi1_LL_%i",i), Form("hPhi1_LL_%f",bounds[i]),
			   100, -TMath::Pi(), TMath::Pi() );
    hPhi1_LO[i] = new TH1D(Form("hPhi1_LO_%i",i), Form("hPhi1_LO_%f",bounds[i]),
			   100, -TMath::Pi(), TMath::Pi() );
    hPhi1_LL_LO[i] = new TH1D(Form("hPhi1_LL_LO_%i",i),
			      Form("hPhi1_LL_LO_%f",bounds[i]),
			      100, -TMath::Pi(), TMath::Pi() );

    hTheta1_LL[i] = new TH1D(Form("hTheta1_LL_%i",i),
			     Form("hTheta1_LL_%f",bounds[i]),
			     100, 0, 0.21);
    hTheta1_LO[i] = new TH1D(Form("hTheta1_LO_%i",i),
			     Form("hTheta1_LO_%f",bounds[i]),
			     100, 0, 0.21);
    hTheta1_LL_LO[i] = new TH1D(Form("hTheta1_LL_LO_%i",i),
				Form("hTheta1_LL_LO_%f",bounds[i]),
				100, 0, 0.21);

    hP1_LL[i] = new TH1D(Form("hP1_LL_%i",i),
			 Form("hP1_LL_%f",bounds[i]),
			 200, 0, 200);
    hP1_LO[i] = new TH1D(Form("hP1_LO_%i",i),
			 Form("hP1_LO_%f",bounds[i]),
			 200, 0, 200);
    hP1_LL_LO[i] = new TH1D(Form("hP1_LL_LO_%i",i),
			    Form("hP1_LL_LO_%f",bounds[i]),
			    200, 0, 200);

    //Particle 2
    hPhi2_LL[i] = new TH1D(Form("hPhi2_LL_%i",i), Form("hPhi2_LL_%f",bounds[i]),
			   100, -TMath::Pi(), TMath::Pi() );
    hPhi2_LO[i] = new TH1D(Form("hPhi2_LO_%i",i), Form("hPhi2_LO_%f",bounds[i]),
			   100, -TMath::Pi(), TMath::Pi() );
    hPhi2_LL_LO[i] = new TH1D(Form("hPhi2_LL_LO_%i",i),
			      Form("hPhi2_LL_LO_%f",bounds[i]),
			      100, -TMath::Pi(), TMath::Pi() );

    hTheta2_LL[i] = new TH1D(Form("hTheta2_LL_%i",i),
			     Form("hTheta2_LL_%f",bounds[i]),
			     100, 0, 0.21);
    hTheta2_LO[i] = new TH1D(Form("hTheta2_LO_%i",i),
			     Form("hTheta2_LO_%f",bounds[i]),
			     100, 0, 0.21);
    hTheta2_LL_LO[i] = new TH1D(Form("hTheta2_LL_LO_%i",i),
				Form("hTheta2_LL_LO_%f",bounds[i]),
				100, 0, 0.21);

    hP2_LL[i] = new TH1D(Form("hP2_LL_%i",i),
			 Form("hP2_LL_%f",bounds[i]),
			 200, 0, 200);
    hP2_LO[i] = new TH1D(Form("hP2_LO_%i",i),
			 Form("hP2_LO_%f",bounds[i]),
			 200, 0, 200);
    hP2_LL_LO[i] = new TH1D(Form("hP2_LL_LO_%i",i),
			    Form("hP2_LL_LO_%f",bounds[i]),
			    200, 0, 200);

  }
  // }}}

  //Tree loop
  Bool_t first = true;
  cout << "Number of entries in tree: " << tree_entries << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {//Tree Loop 2
    //cout << "Debug mode" << endl; for (Int_t ev=0; ev<1000; ev++) {
    tree->GetEntry(ev, 0);

    if (first || ev==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!!!!!!!" << endl;
      //cout << "Additional Last-Outer cuts" << endl;

      first = false;
    }

    //Basic checks/cuts
    if(binVar=="vxZ_upstream" && targetPosition==1) continue;
    else if(binVar=="vxZ_downstream" && targetPosition==0) continue;
    
    //General useful variables
    TLorentzVector vPhoton(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);
    if (binVar=="phi_vPh") *boundValue = vPhoton.Phi();
    else if (binVar=="theta_vPh") *boundValue = vPhoton.Theta();
        
    if (trigMask == 65792){//Last-Last
      BinDataFill(hPhi1_LL, phi_traj1, nBins, *boundValue, bounds);
      BinDataFill(hPhi2_LL, phi_traj2, nBins, *boundValue, bounds);

      BinDataFill(hTheta1_LL, theta_traj1, nBins, *boundValue, bounds);
      BinDataFill(hTheta2_LL, theta_traj2, nBins, *boundValue, bounds);

      BinDataFill(hP1_LL, qP_traj1, nBins, *boundValue, bounds);
      BinDataFill(hP2_LL, -1.0*qP_traj2, nBins, *boundValue, bounds);
    }
    else if (trigMask == 65540){//Last-Outer
      //additional last-outer cuts
      //if (vx_z < -202.235) continue;
      //if (q_transverse < 0.781648) continue;
      //if (phi_traj1 < 0.7 && phi_traj1 > -0.7) continue;
      //if (phi_traj2 < 0.7 && phi_traj1 > -0.7) continue;
      //if (*boundValue < 0.0356868) continue;
	    
      BinDataFill(hPhi1_LO, phi_traj1, nBins, *boundValue, bounds);
      BinDataFill(hPhi2_LO, phi_traj2, nBins, *boundValue, bounds);

      BinDataFill(hTheta1_LO, theta_traj1, nBins, *boundValue, bounds);
      BinDataFill(hTheta2_LO, theta_traj2, nBins, *boundValue, bounds);

      BinDataFill(hP1_LO, qP_traj1, nBins, *boundValue, bounds);
      BinDataFill(hP2_LO, -1.0*qP_traj2, nBins, *boundValue, bounds);
    }
    else if (trigMask == 65796){//Last-Last and Last-Outer
      BinDataFill(hPhi1_LL_LO, phi_traj1, nBins, *boundValue, bounds);
      BinDataFill(hPhi2_LL_LO, phi_traj2, nBins, *boundValue, bounds);

      BinDataFill(hTheta1_LL_LO, theta_traj1, nBins, *boundValue, bounds);
      BinDataFill(hTheta2_LL_LO, theta_traj2, nBins, *boundValue, bounds);

      BinDataFill(hP1_LL_LO, qP_traj1, nBins, *boundValue, bounds);
      BinDataFill(hP2_LL_LO, -1.0*qP_traj2, nBins, *boundValue, bounds);
    }

    
    //Bin data
    /*if (targetPosition == 0) {//UpStream target
            
      }//End UpStream target
      else if (targetPosition == 1) {//DownStream target
      
      }//End DownStream target*/

  }//End tree loop

  //Draw graphs
  // {{{
  const Int_t nCanvas = 20;
  if (nCanvas < nBins/5.0) {
    cout << "nCanvas is too small..." << endl;
    exit(EXIT_FAILURE); }
  
  TCanvas *cPhi[nCanvas], *cTheta[nCanvas], *cP[nCanvas];
  for (Int_t i=0; i<nBins/5.0; i++) {
    cPhi[i] = new TCanvas(Form("cPhi_%i", i), Form("cPhi_%i", i) );
    (nBins>=5) ? cPhi[i]->Divide(3, 5) : cPhi[i]->Divide(3, nBins); }
  for (Int_t i=0; i<nBins/5.0; i++) {
    cTheta[i] = new TCanvas(Form("cTheta_%i", i), Form("cTheta_%i", i) );
    (nBins>=5) ? cTheta[i]->Divide(3, 5) : cTheta[i]->Divide(3, nBins); }
  for (Int_t i=0; i<nBins/5.0; i++) {
    cP[i] = new TCanvas(Form("cP_%i", i), Form("cP_%i", i) );
    (nBins>=5) ? cP[i]->Divide(3, 5) : cP[i]->Divide(3, nBins); }
  
  for (Int_t ibin=0, ic=0, ipad=1; ibin<nBins; ibin++, ipad++) {
    if (ipad>15) {ic++; ipad=1;}
    cout << "Bin: " << ibin << "    pad: " << ipad << "    canvas: "<<ic<<endl;
    Int_t start_pad = ipad;
    
    cPhi[ic]->cd(ipad);
    hPhi1_LL[ibin]->Draw();
    hPhi2_LL[ibin]->Draw("sames");
    hPhi2_LL[ibin]->SetLineColor(kRed);

    ipad++; cPhi[ic]->cd(ipad);
    hPhi1_LO[ibin]->Draw();
    hPhi2_LO[ibin]->Draw("sames");
    hPhi2_LO[ibin]->SetLineColor(kRed);

    ipad++; cPhi[ic]->cd(ipad);
    hPhi1_LL_LO[ibin]->Draw();
    hPhi2_LL_LO[ibin]->Draw("sames");
    hPhi2_LL_LO[ibin]->SetLineColor(kRed);


    Int_t tmp_pad = start_pad;
    cTheta[ic]->cd(tmp_pad);
    hTheta1_LL[ibin]->Draw();
    hTheta2_LL[ibin]->Draw("sames");
    hTheta2_LL[ibin]->SetLineColor(kRed);

    tmp_pad++; cTheta[ic]->cd(tmp_pad);
    hTheta1_LO[ibin]->Draw();
    hTheta2_LO[ibin]->Draw("sames");
    hTheta2_LO[ibin]->SetLineColor(kRed);

    tmp_pad++; cTheta[ic]->cd(tmp_pad);
    hTheta1_LL_LO[ibin]->Draw();
    hTheta2_LL_LO[ibin]->Draw("sames");
    hTheta2_LL_LO[ibin]->SetLineColor(kRed);

    tmp_pad = start_pad;
    cP[ic]->cd(tmp_pad);
    hP1_LL[ibin]->Draw();
    hP2_LL[ibin]->Draw("sames");
    hP2_LL[ibin]->SetLineColor(kRed);

    tmp_pad++; cP[ic]->cd(tmp_pad);
    hP1_LO[ibin]->Draw();
    hP2_LO[ibin]->Draw("sames");
    hP2_LO[ibin]->SetLineColor(kRed);

    tmp_pad++; cP[ic]->cd(tmp_pad);
    hP1_LL_LO[ibin]->Draw();
    hP2_LL_LO[ibin]->Draw("sames");
    hP2_LL_LO[ibin]->SetLineColor(kRed);
  }//loop of bins
  // }}}

  cout << " " << endl;
  cout << "Bounds for: " << binVar << endl;
  for (Int_t i=0; i<nBounds; i++) {
    cout << bounds[i] << endl;
  }
  cout << " " << endl;

  //Write output
  ///////////////
  // {{{
  if (Qflag || wflag){
    TFile* fout = (Qflag ? new TFile(outFile, "RECREATE")
		   : new TFile("Output.root", "RECREATE") );

    for (Int_t i=0; i<nBins; i++) {
      hPhi1_LL[i]->Write();
      hPhi2_LL[i]->Write();
    }

    fout->Close();

    cout << " " << endl;
    if (Qflag) cout << outFile << " file written" << endl;
    else cout << "Output.root file written" << endl;
  }
  // }}}
  
  theApp.Run();//Needed to make root graphics work on C++
}
