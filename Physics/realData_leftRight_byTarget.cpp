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
    cout << "Option:  -S Left/right asymmetry choice" << endl;
    cout << "	 (True=no spin influence, Spin=spin influence, default=Spin)"
	 << endl;
    cout << "Option:  -P       (Turn off polarization and dilution corrections)"
	 << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t wflag=0, Qflag=0, fflag=0, Sflag=0, Pflag=0;
  Int_t c;
  TString fname = "", outFile = "", leftrightChoice="";
  
  while ((c = getopt (argc, argv, "Pwf:Q:S:")) != -1) {
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
    case 'S':
      Sflag = 1;
      leftrightChoice += optarg;
      break;
    case 'P':
      Pflag = 1;
      break;      
    case '?':
      if (optopt == 'u')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'S')
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
  
  //Opening data files/getting trees
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
  Double_t pT_bounds[nBounds], M_bounds[nBounds];//comment open angle
  //Double_t pT_bounds[nBounds]; //Open angle
  for (Int_t i=0; i<nBounds; i++) {
    xN_bounds[i] = tv_xN_bounds[i];
    xPi_bounds[i] = tv_xPi_bounds[i];
    xF_bounds[i] = tv_xF_bounds[i];
    pT_bounds[i] = tv_pT_bounds[i];
    M_bounds[i] = tv_M_bounds[i];//comment open angle
  }
  //Double_t M_bounds[] = {0.0, 0.06, 0.08, 0.10, 0.11, 0.12, 0.15, 0.35};//Open angle

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

  //Vertex specific
  Double_t vx_z;
  Int_t targetPosition;
  //Drell-Yan Angles
  Double_t PhiS_simple, Theta_CS, vOpenAngle;
  //Event
  Int_t trigMask, MasterTrigMask;
  //Target values
  Double_t Spin_0, Spin_1, Spin_2, Spin_3, Spin_4, Spin_5, Spin_6;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;
  
  //Vertex specific
  tree->SetBranchAddress("vx_z", &vx_z);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("Theta_CS", &Theta_CS);
  tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Event
  tree->SetBranchAddress("trigMask", &trigMask);
  tree->SetBranchAddress("MasterTrigMask", &MasterTrigMask);
  //Target values
  tree->SetBranchAddress("Spin_0", &Spin_0);
  tree->SetBranchAddress("Spin_1", &Spin_1);
  tree->SetBranchAddress("Spin_2", &Spin_2);
  tree->SetBranchAddress("Spin_3", &Spin_3);
  tree->SetBranchAddress("Spin_4", &Spin_4);
  tree->SetBranchAddress("Spin_5", &Spin_5);
  tree->SetBranchAddress("Spin_6", &Spin_6);
  //DY-variables
  tree->SetBranchAddress("x_beam", &x_beam);
  tree->SetBranchAddress("x_target", &x_target);
  tree->SetBranchAddress("x_feynman", &x_feynman);
  tree->SetBranchAddress("q_transverse", &q_transverse);
  tree->SetBranchAddress("Mmumu", &Mmumu);
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

  //By target
  unsigned long long xN_Left_UpStream_Up[nBins], xN_Right_UpStream_Up[nBins];
  unsigned long long xN_Left_UpStream_Down[nBins],xN_Right_UpStream_Down[nBins];
  unsigned long long xN_Left_DownStream_Up[nBins], xN_Right_DownStream_Up[nBins];
  unsigned long long xN_Left_DownStream_Down[nBins],xN_Right_DownStream_Down[nBins];
  unsigned long long xPi_Left_UpStream_Up[nBins], xPi_Right_UpStream_Up[nBins];
  unsigned long long xPi_Left_UpStream_Down[nBins],xPi_Right_UpStream_Down[nBins];
  unsigned long long xPi_Left_DownStream_Up[nBins], xPi_Right_DownStream_Up[nBins];
  unsigned long long xPi_Left_DownStream_Down[nBins],xPi_Right_DownStream_Down[nBins];
  unsigned long long xF_Left_UpStream_Up[nBins], xF_Right_UpStream_Up[nBins];
  unsigned long long xF_Left_UpStream_Down[nBins],xF_Right_UpStream_Down[nBins];
  unsigned long long xF_Left_DownStream_Up[nBins], xF_Right_DownStream_Up[nBins];
  unsigned long long xF_Left_DownStream_Down[nBins],xF_Right_DownStream_Down[nBins];
  unsigned long long pT_Left_UpStream_Up[nBins], pT_Right_UpStream_Up[nBins];
  unsigned long long pT_Left_UpStream_Down[nBins],pT_Right_UpStream_Down[nBins];
  unsigned long long pT_Left_DownStream_Up[nBins], pT_Right_DownStream_Up[nBins];
  unsigned long long pT_Left_DownStream_Down[nBins],pT_Right_DownStream_Down[nBins];
  unsigned long long M_Left_UpStream_Up[nBins], M_Right_UpStream_Up[nBins];
  unsigned long long M_Left_UpStream_Down[nBins],M_Right_UpStream_Down[nBins];
  unsigned long long M_Left_DownStream_Up[nBins], M_Right_DownStream_Up[nBins];
  unsigned long long M_Left_DownStream_Down[nBins],M_Right_DownStream_Down[nBins];
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

    //By target
    xN_Left_UpStream_Up[i]=0; xN_Right_UpStream_Up[i]=0;
    xN_Left_UpStream_Down[i]=0; xN_Right_UpStream_Down[i]=0;
    xN_Left_DownStream_Up[i]=0; xN_Right_DownStream_Up[i]=0;
    xN_Left_DownStream_Down[i]=0; xN_Right_DownStream_Down[i]=0;
    xPi_Left_UpStream_Up[i]=0; xPi_Right_UpStream_Up[i]=0;
    xPi_Left_UpStream_Down[i]=0; xPi_Right_UpStream_Down[i]=0;
    xPi_Left_DownStream_Up[i]=0; xPi_Right_DownStream_Up[i]=0;
    xPi_Left_DownStream_Down[i]=0; xPi_Right_DownStream_Down[i]=0;
    xF_Left_UpStream_Up[i]=0; xF_Right_UpStream_Up[i]=0;
    xF_Left_UpStream_Down[i]=0; xF_Right_UpStream_Down[i]=0;
    xF_Left_DownStream_Up[i]=0; xF_Right_DownStream_Up[i]=0;
    xF_Left_DownStream_Down[i]=0; xF_Right_DownStream_Down[i]=0;
    pT_Left_UpStream_Up[i]=0; pT_Right_UpStream_Up[i]=0;
    pT_Left_UpStream_Down[i]=0; pT_Right_UpStream_Down[i]=0;
    pT_Left_DownStream_Up[i]=0; pT_Right_DownStream_Up[i]=0;
    pT_Left_DownStream_Down[i]=0; pT_Right_DownStream_Down[i]=0;
    M_Left_UpStream_Up[i]=0; M_Right_UpStream_Up[i]=0;
    M_Left_UpStream_Down[i]=0; M_Right_UpStream_Down[i]=0;
    M_Left_DownStream_Up[i]=0; M_Right_DownStream_Up[i]=0;
    M_Left_DownStream_Down[i]=0; M_Right_DownStream_Down[i]=0;
  }
  // }}}
  
  //Phi filling
  ///////////////
  // {{{
  const Int_t phiBins = 16;
  TH1D* h_phi_UpStream = new TH1D("h_phi_UpStream", "h_phi_UpStream", phiBins,
				  -TMath::Pi()/2, 3*TMath::Pi()/2);
  TH1D* h_phi_DownStream = new TH1D("h_phi_DownStream", "h_phi_DownStream",
				    phiBins,
				    -TMath::Pi()/2, 3*TMath::Pi()/2);
  TH1D* h_phi = new TH1D("h_phi", "h_phi", phiBins,
			 -TMath::Pi()/2, 3*TMath::Pi()/2);
  // }}}

  //Tree loop
  Bool_t first = true;
  Int_t tree_entries = tree->GetEntries();//Tree Loop
  cout << "Number of entries in tree: " << tree_entries << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    //cout << "Debug mode" << endl; for (Int_t ev=0; ev<1000; ev++) {
    tree->GetEntry(ev, 0);

    if (first || ev==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!!!!!!!" << endl;
      //cout << "Trig Mask settings:" << endl;
      //cout << "Only Last-Outer" << endl;
      //cout << "Only Last-Last" << endl;
      //cout << " " << endl;
      //cout << "Left = Top, Right = Bottom" << endl;
      //cout << " " << endl;
      //cout << "M bins are binned in vOpenAngle" << endl;
      //cout << " " << endl;

      if(leftrightChoice=="True") {
	cout << "True left/right asymmetry (no spin influence)" << endl;
      }
      else if (leftrightChoice=="Spin" || leftrightChoice==""){
	cout << "Spin influnced left/right asymmetry" << endl;
      }

      first = false;
    }

    //Trig Mask
    //524 == Last-Outer
    //768 == Last-Last
    //780 == Last-Last && Last-Outer
    //if (trigMask == 65792 || trigMask == 65796) continue; //Only Last-Outer
    if (trigMask != 65792) continue; //Only Last-Outer 
    //if (trigMask == 65540| trigMask == 65796) continue; //Only Last-Last

    //Choose Left/Right
    // {{{
    Double_t phi_photon_lab = ShiftPhiSimple(PhiS_simple);
    Bool_t Left=false, Right=false;
    if (leftrightChoice=="True"){
      //True spectrometer left/right (no spin influence)
      if (phi_photon_lab < TMath::Pi()/2 && phi_photon_lab > -TMath::Pi()/2){
	//if (phi_photon_lab < TMath::Pi() && phi_photon_lab > 0.0){//Top 
	Left = true;
      }
      else if (phi_photon_lab < 3*TMath::Pi()/2 && phi_photon_lab>TMath::Pi()/2){
	//else if (phi_photon_lab > TMath::Pi() || phi_photon_lab < 0.0){//Bottom
	Right = true;
      }
      else {
	cout << "No Left or Right choosen" << endl;
	cout << phi_photon_lab << " " << Spin_0 << endl;
	cout << " " << endl;
      }
    }
    else if (leftrightChoice=="Spin" || leftrightChoice==""){
      //Spin influenced left/right
      if (phi_photon_lab < TMath::Pi()/2 && phi_photon_lab > -TMath::Pi()/2
	  && Spin_0 > 0){ // Target spin up
	Left = true;
      }
      else if (phi_photon_lab < 3*TMath::Pi()/2 && phi_photon_lab > TMath::Pi()/2
	       && Spin_0 > 0){// Target spin up
	Right = true;
      }
      else if (phi_photon_lab < 3*TMath::Pi()/2 && phi_photon_lab > TMath::Pi()/2
	       && Spin_0 < 0){ // Target spin down
	Left = true;
      }
      else if (phi_photon_lab < TMath::Pi()/2 && phi_photon_lab > -TMath::Pi()/2
	       && Spin_0 < 0){// Target spin down
	Right = true;
      }
      else {
	cout << "No Left or Right choosen" << endl;
	cout << phi_photon_lab << " " << Spin_0 << endl;
	cout << " " << endl;
      }
    }
    // }}}
    
    //Bin data
    if (targetPosition == 0) {//UpStream target
      
      if (Left){//Left
	// {{{
	BinDataCounts(xN_Left_UpStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Left_UpStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Left_UpStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Left_UpStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Left_UpStream, nBins, Mmumu, M_bounds);//co open angle
	//BinDataCounts(M_Left_UpStream, nBins, vOpenAngle, M_bounds);//openAngle

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Left_UpStream_Up, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_UpStream_Up, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_UpStream_Up, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_UpStream_Up, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_UpStream_Up, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Left_UpStream_Up, nBins, vOpenAngle, M_bounds);//openAngle
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Left_UpStream_Down, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_UpStream_Down, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_UpStream_Down, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_UpStream_Down, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_UpStream_Down, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Left_UpStream_Down, nBins, vOpenAngle, M_bounds);//openAngle
	}
	// }}}
      }//End Left
      else if (Right){//Right
	// {{{
	BinDataCounts(xN_Right_UpStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Right_UpStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Right_UpStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Right_UpStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Right_UpStream, nBins, Mmumu, M_bounds);//openAngle
	//BinDataCounts(M_Right_UpStream, nBins, vOpenAngle, M_bounds);//openAngle

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Right_UpStream_Up, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_UpStream_Up, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_UpStream_Up, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_UpStream_Up, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_UpStream_Up, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Right_UpStream_Up, nBins, vOpenAngle, M_bounds);//openAngle
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Right_UpStream_Down, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_UpStream_Down, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_UpStream_Down, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_UpStream_Down, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_UpStream_Down, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Right_UpStream_Down, nBins, vOpenAngle, M_bounds);//openAngle
	}
	// }}}
      }//End Right

      h_phi_UpStream->Fill(phi_photon_lab);
      h_phi->Fill(phi_photon_lab);
    }//End UpStream target
    else if (targetPosition == 1) {//DownStream target
      
      if (Left){//Left
	// {{{
	BinDataCounts(xN_Left_DownStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Left_DownStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Left_DownStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Left_DownStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Left_DownStream, nBins, Mmumu, M_bounds);//openAngle
	//BinDataCounts(M_Left_DownStream, nBins, vOpenAngle, M_bounds);//openAngle

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Left_DownStream_Up, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_DownStream_Up, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_DownStream_Up, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_DownStream_Up, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_DownStream_Up, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Left_DownStream_Up, nBins, vOpenAngle, M_bounds);//openAngle
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Left_DownStream_Down, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_DownStream_Down, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_DownStream_Down, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_DownStream_Down, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_DownStream_Down, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Left_DownStream_Down, nBins, vOpenAngle, M_bounds);//openAngle
	}
	// }}}
      }//End Left
      else if (Right){//Right
	// {{{
	BinDataCounts(xN_Right_DownStream, nBins, x_target, xN_bounds);
	BinDataCounts(xPi_Right_DownStream, nBins, x_beam, xPi_bounds);
	BinDataCounts(xF_Right_DownStream, nBins, x_feynman, xF_bounds);
	BinDataCounts(pT_Right_DownStream, nBins, q_transverse, pT_bounds);
	BinDataCounts(M_Right_DownStream, nBins, Mmumu, M_bounds);//openAngle
	//BinDataCounts(M_Right_DownStream, nBins, vOpenAngle, M_bounds);//openAngle

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Right_DownStream_Up, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_DownStream_Up, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_DownStream_Up, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_DownStream_Up, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_DownStream_Up, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Right_DownStream_Up, nBins, vOpenAngle, M_bounds);//openAngle
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Right_DownStream_Down, nBins, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_DownStream_Down, nBins, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_DownStream_Down, nBins, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_DownStream_Down, nBins, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_DownStream_Down, nBins, Mmumu, M_bounds);//openAngle
	  //BinDataCounts(M_Right_DownStream_Down, nBins, vOpenAngle, M_bounds);//openAngle
	}
	// }}}
      }//End Right

      h_phi_DownStream->Fill(phi_photon_lab);
      h_phi->Fill(phi_photon_lab);
    }//End DownStream target

  }//End tree loop
  

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
  Double_t xN_Asym_UpStream_Up[nBins], xN_Asym_DownStream_Up[nBins]; //By target
  Double_t xN_Asym_UpStream_Down[nBins], xN_Asym_DownStream_Down[nBins];
  Double_t e_xN_Asym_UpStream_Up[nBins], e_xN_Asym_DownStream_Up[nBins];
  Double_t e_xN_Asym_UpStream_Down[nBins], e_xN_Asym_DownStream_Down[nBins];
  BinnedLeftRight(xN_Left_UpStream_Up, xN_Right_UpStream_Up, 
		  xN_Asym_UpStream_Up, e_xN_Asym_UpStream_Up, nBins);
  BinnedLeftRight(xN_Left_UpStream_Down, xN_Right_UpStream_Down,
		  xN_Asym_UpStream_Down, e_xN_Asym_UpStream_Down, nBins);
  BinnedLeftRight(xN_Left_DownStream_Up, xN_Right_DownStream_Up, 
		  xN_Asym_DownStream_Up, e_xN_Asym_DownStream_Up, nBins);
  BinnedLeftRight(xN_Left_DownStream_Down, xN_Right_DownStream_Down,
		  xN_Asym_DownStream_Down, e_xN_Asym_DownStream_Down, nBins);
  Double_t xPi_Asym_UpStream[nBins], xPi_Asym_DownStream[nBins], xPi_Asym[nBins];
  Double_t e_xPi_Asym_UpStream[nBins], e_xPi_Asym_DownStream[nBins];
  Double_t e_xPi_Asym[nBins];
  BinnedLeftRight(xPi_Left_UpStream, xPi_Right_UpStream, xPi_Left_DownStream,
		  xPi_Right_DownStream, xPi_Asym_UpStream, xPi_Asym_DownStream,
		  xPi_Asym, e_xPi_Asym_UpStream, e_xPi_Asym_DownStream,
		  e_xPi_Asym,
		  nBins);
  Double_t xPi_Asym_UpStream_Up[nBins], xPi_Asym_DownStream_Up[nBins]; //By target
  Double_t xPi_Asym_UpStream_Down[nBins], xPi_Asym_DownStream_Down[nBins];
  Double_t e_xPi_Asym_UpStream_Up[nBins], e_xPi_Asym_DownStream_Up[nBins];
  Double_t e_xPi_Asym_UpStream_Down[nBins], e_xPi_Asym_DownStream_Down[nBins];
  BinnedLeftRight(xPi_Left_UpStream_Up, xPi_Right_UpStream_Up, 
		  xPi_Asym_UpStream_Up, e_xPi_Asym_UpStream_Up, nBins);
  BinnedLeftRight(xPi_Left_UpStream_Down, xPi_Right_UpStream_Down,
		  xPi_Asym_UpStream_Down, e_xPi_Asym_UpStream_Down, nBins);
  BinnedLeftRight(xPi_Left_DownStream_Up, xPi_Right_DownStream_Up, 
		  xPi_Asym_DownStream_Up, e_xPi_Asym_DownStream_Up, nBins);
  BinnedLeftRight(xPi_Left_DownStream_Down, xPi_Right_DownStream_Down,
		  xPi_Asym_DownStream_Down, e_xPi_Asym_DownStream_Down, nBins);
  Double_t xF_Asym_UpStream[nBins], xF_Asym_DownStream[nBins], xF_Asym[nBins];
  Double_t e_xF_Asym_UpStream[nBins], e_xF_Asym_DownStream[nBins];
  Double_t e_xF_Asym[nBins];
  BinnedLeftRight(xF_Left_UpStream, xF_Right_UpStream, xF_Left_DownStream,
		  xF_Right_DownStream, xF_Asym_UpStream, xF_Asym_DownStream,
		  xF_Asym, e_xF_Asym_UpStream, e_xF_Asym_DownStream, e_xF_Asym,
		  nBins);
  Double_t xF_Asym_UpStream_Up[nBins], xF_Asym_DownStream_Up[nBins]; //By target
  Double_t xF_Asym_UpStream_Down[nBins], xF_Asym_DownStream_Down[nBins];
  Double_t e_xF_Asym_UpStream_Up[nBins], e_xF_Asym_DownStream_Up[nBins];
  Double_t e_xF_Asym_UpStream_Down[nBins], e_xF_Asym_DownStream_Down[nBins];
  BinnedLeftRight(xF_Left_UpStream_Up, xF_Right_UpStream_Up, 
		  xF_Asym_UpStream_Up, e_xF_Asym_UpStream_Up, nBins);
  BinnedLeftRight(xF_Left_UpStream_Down, xF_Right_UpStream_Down,
		  xF_Asym_UpStream_Down, e_xF_Asym_UpStream_Down, nBins);
  BinnedLeftRight(xF_Left_DownStream_Up, xF_Right_DownStream_Up, 
		  xF_Asym_DownStream_Up, e_xF_Asym_DownStream_Up, nBins);
  BinnedLeftRight(xF_Left_DownStream_Down, xF_Right_DownStream_Down,
		  xF_Asym_DownStream_Down, e_xF_Asym_DownStream_Down, nBins);
  Double_t pT_Asym_UpStream[nBins], pT_Asym_DownStream[nBins], pT_Asym[nBins];
  Double_t e_pT_Asym_UpStream[nBins], e_pT_Asym_DownStream[nBins];
  Double_t e_pT_Asym[nBins];
  BinnedLeftRight(pT_Left_UpStream, pT_Right_UpStream, pT_Left_DownStream,
  		  pT_Right_DownStream, pT_Asym_UpStream, pT_Asym_DownStream,
  		  pT_Asym, e_pT_Asym_UpStream, e_pT_Asym_DownStream, e_pT_Asym,
  		  nBins);
  Double_t pT_Asym_UpStream_Up[nBins], pT_Asym_DownStream_Up[nBins]; //By target
  Double_t pT_Asym_UpStream_Down[nBins], pT_Asym_DownStream_Down[nBins];
  Double_t e_pT_Asym_UpStream_Up[nBins], e_pT_Asym_DownStream_Up[nBins];
  Double_t e_pT_Asym_UpStream_Down[nBins], e_pT_Asym_DownStream_Down[nBins];
  BinnedLeftRight(pT_Left_UpStream_Up, pT_Right_UpStream_Up, 
		  pT_Asym_UpStream_Up, e_pT_Asym_UpStream_Up, nBins);
  BinnedLeftRight(pT_Left_UpStream_Down, pT_Right_UpStream_Down,
		  pT_Asym_UpStream_Down, e_pT_Asym_UpStream_Down, nBins);
  BinnedLeftRight(pT_Left_DownStream_Up, pT_Right_DownStream_Up, 
		  pT_Asym_DownStream_Up, e_pT_Asym_DownStream_Up, nBins);
  BinnedLeftRight(pT_Left_DownStream_Down, pT_Right_DownStream_Down,
		  pT_Asym_DownStream_Down, e_pT_Asym_DownStream_Down, nBins);
  Double_t M_Asym_UpStream[nBins], M_Asym_DownStream[nBins], M_Asym[nBins];
  Double_t e_M_Asym_UpStream[nBins], e_M_Asym_DownStream[nBins];
  Double_t e_M_Asym[nBins];
  BinnedLeftRight(M_Left_UpStream, M_Right_UpStream, M_Left_DownStream,
		  M_Right_DownStream, M_Asym_UpStream, M_Asym_DownStream,
		  M_Asym, e_M_Asym_UpStream, e_M_Asym_DownStream, e_M_Asym,
		  nBins);
  Double_t M_Asym_UpStream_Up[nBins], M_Asym_DownStream_Up[nBins]; //By target
  Double_t M_Asym_UpStream_Down[nBins], M_Asym_DownStream_Down[nBins];
  Double_t e_M_Asym_UpStream_Up[nBins], e_M_Asym_DownStream_Up[nBins];
  Double_t e_M_Asym_UpStream_Down[nBins], e_M_Asym_DownStream_Down[nBins];
  BinnedLeftRight(M_Left_UpStream_Up, M_Right_UpStream_Up, 
		  M_Asym_UpStream_Up, e_M_Asym_UpStream_Up, nBins);
  BinnedLeftRight(M_Left_UpStream_Down, M_Right_UpStream_Down,
		  M_Asym_UpStream_Down, e_M_Asym_UpStream_Down, nBins);
  BinnedLeftRight(M_Left_DownStream_Up, M_Right_DownStream_Up, 
		  M_Asym_DownStream_Up, e_M_Asym_DownStream_Up, nBins);
  BinnedLeftRight(M_Left_DownStream_Down, M_Right_DownStream_Down,
		  M_Asym_DownStream_Down, e_M_Asym_DownStream_Down, nBins);
  // }}}

  //Correct for dilution factor and polarization
  ///////////////
  // {{{
  if(Pflag) {
    cout << " " << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << "No dilution factor or polarization corrections" << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << " " << endl;
  }
  else {
    cout << " " << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << "Dilution factor and polarization corrections are made" << endl;
    cout << " " << endl;
    
    CorrectDilPol(xN_Asym_UpStream, e_xN_Asym_UpStream,
		  Dil_xN_UpStream, Pol_xN_UpStream, nBins);
    CorrectDilPol(xN_Asym_DownStream, e_xN_Asym_DownStream,
		  Dil_xN_DownStream, Pol_xN_DownStream, nBins);
    CorrectDilPol(xN_Asym, e_xN_Asym,
		  Dil_xN, Pol_xN, nBins);
    CorrectDilPol(xPi_Asym_UpStream, e_xPi_Asym_UpStream,
		  Dil_xPi_UpStream, Pol_xPi_UpStream, nBins);
    CorrectDilPol(xPi_Asym_DownStream, e_xPi_Asym_DownStream,
		  Dil_xPi_DownStream, Pol_xPi_DownStream, nBins);
    CorrectDilPol(xPi_Asym, e_xPi_Asym,
		  Dil_xPi, Pol_xPi, nBins);
    CorrectDilPol(xF_Asym_UpStream, e_xF_Asym_UpStream,
		  Dil_xF_UpStream, Pol_xF_UpStream, nBins);
    CorrectDilPol(xF_Asym_DownStream, e_xF_Asym_DownStream,
		  Dil_xF_DownStream, Pol_xF_DownStream, nBins);
    CorrectDilPol(xF_Asym, e_xF_Asym,
		  Dil_xF, Pol_xF, nBins);
    CorrectDilPol(pT_Asym_UpStream, e_pT_Asym_UpStream,
		  Dil_pT_UpStream, Pol_pT_UpStream, nBins);
    CorrectDilPol(pT_Asym_DownStream, e_pT_Asym_DownStream,
		  Dil_pT_DownStream, Pol_pT_DownStream, nBins);
    CorrectDilPol(pT_Asym, e_pT_Asym,
		  Dil_pT, Pol_pT, nBins);
    CorrectDilPol(M_Asym_UpStream, e_M_Asym_UpStream,
		  Dil_M_UpStream, Pol_M_UpStream, nBins);
    CorrectDilPol(M_Asym_DownStream, e_M_Asym_DownStream,
		  Dil_M_DownStream, Pol_M_DownStream, nBins);
    CorrectDilPol(M_Asym, e_M_Asym,
		  Dil_M, Pol_M, nBins);

    //Correct dilution/polarization by target
    CorrectDilPol(xN_Asym_UpStream_Up, e_xN_Asym_UpStream_Up, Dil_xN_UpStream_Up,
		  Pol_xN_UpStream_Up, nBins);
    CorrectDilPol(xN_Asym_UpStream_Down, e_xN_Asym_UpStream_Down,
		  Dil_xN_UpStream_Down, Pol_xN_UpStream_Down, nBins);
    CorrectDilPol(xN_Asym_DownStream_Up, e_xN_Asym_DownStream_Up,
		  Dil_xN_DownStream_Up, Pol_xN_DownStream_Up, nBins);
    CorrectDilPol(xN_Asym_DownStream_Down, e_xN_Asym_DownStream_Down,
		  Dil_xN_DownStream_Down, Pol_xN_DownStream_Down, nBins);

    CorrectDilPol(xPi_Asym_UpStream_Up, e_xPi_Asym_UpStream_Up, Dil_xPi_UpStream_Up,
		  Pol_xPi_UpStream_Up, nBins);
    CorrectDilPol(xPi_Asym_UpStream_Down, e_xPi_Asym_UpStream_Down,
		  Dil_xPi_UpStream_Down, Pol_xPi_UpStream_Down, nBins);
    CorrectDilPol(xPi_Asym_DownStream_Up, e_xPi_Asym_DownStream_Up,
		  Dil_xPi_DownStream_Up, Pol_xPi_DownStream_Up, nBins);
    CorrectDilPol(xPi_Asym_DownStream_Down, e_xPi_Asym_DownStream_Down,
		  Dil_xPi_DownStream_Down, Pol_xPi_DownStream_Down, nBins);

    CorrectDilPol(xF_Asym_UpStream_Up, e_xF_Asym_UpStream_Up, Dil_xF_UpStream_Up,
		  Pol_xF_UpStream_Up, nBins);
    CorrectDilPol(xF_Asym_UpStream_Down, e_xF_Asym_UpStream_Down,
		  Dil_xF_UpStream_Down, Pol_xF_UpStream_Down, nBins);
    CorrectDilPol(xF_Asym_DownStream_Up, e_xF_Asym_DownStream_Up,
		  Dil_xF_DownStream_Up, Pol_xF_DownStream_Up, nBins);
    CorrectDilPol(xF_Asym_DownStream_Down, e_xF_Asym_DownStream_Down,
		  Dil_xF_DownStream_Down, Pol_xF_DownStream_Down, nBins);

    CorrectDilPol(pT_Asym_UpStream_Up, e_pT_Asym_UpStream_Up, Dil_pT_UpStream_Up,
		  Pol_pT_UpStream_Up, nBins);
    CorrectDilPol(pT_Asym_UpStream_Down, e_pT_Asym_UpStream_Down,
		  Dil_pT_UpStream_Down, Pol_pT_UpStream_Down, nBins);
    CorrectDilPol(pT_Asym_DownStream_Up, e_pT_Asym_DownStream_Up,
		  Dil_pT_DownStream_Up, Pol_pT_DownStream_Up, nBins);
    CorrectDilPol(pT_Asym_DownStream_Down, e_pT_Asym_DownStream_Down,
		  Dil_pT_DownStream_Down, Pol_pT_DownStream_Down, nBins);

    CorrectDilPol(M_Asym_UpStream_Up, e_M_Asym_UpStream_Up, Dil_M_UpStream_Up,
		  Pol_M_UpStream_Up, nBins);
    CorrectDilPol(M_Asym_UpStream_Down, e_M_Asym_UpStream_Down,
		  Dil_M_UpStream_Down, Pol_M_UpStream_Down, nBins);
    CorrectDilPol(M_Asym_DownStream_Up, e_M_Asym_DownStream_Up,
		  Dil_M_DownStream_Up, Pol_M_DownStream_Up, nBins);
    CorrectDilPol(M_Asym_DownStream_Down, e_M_Asym_DownStream_Down,
		  Dil_M_DownStream_Down, Pol_M_DownStream_Down, nBins);
  }
  // }}}
    
  //Phi Left/Right
  ///////////////
  // {{{
  TH1D* h_LR_phi_UpStream = new TH1D("h_LR_phi_UpStream", "h_LR_phi_UpStream",
				     phiBins/2, -TMath::Pi()/2, TMath::Pi()/2 );
  TH1D* h_LR_phiMirror_UpStream = new TH1D("h_LR_phiMirror_UpStream",
					   "h_LR_phiMirror_UpStream",
					   phiBins/2,
					   -TMath::Pi()/2, TMath::Pi()/2 );
  TH1D* h_LR_phi_DownStream = new TH1D("h_LR_phi_DownStream",
				       "h_LR_phi_DownStream",
				       phiBins/2, -TMath::Pi()/2, TMath::Pi()/2 );
  TH1D* h_LR_phiMirror_DownStream = new TH1D("h_LR_phiMirror_DownStream",
					     "h_LR_phiMirror_DownStream",
					     phiBins/2,
					     -TMath::Pi()/2, TMath::Pi()/2 );
  TH1D* h_LR_phi = new TH1D("h_LR_phi", "h_LR_phi",
			    phiBins/2, -TMath::Pi()/2, TMath::Pi()/2 );
  TH1D* h_LR_phiMirror = new TH1D("h_LR_phiMirror", "h_LR_phiMirror", phiBins/2,
				  -TMath::Pi()/2, TMath::Pi()/2 );

  Hist_LRAsym(h_phi_UpStream, h_LR_phi_UpStream, phiBins);
  Hist_LR_MirrorAsym(h_phi_UpStream, h_LR_phiMirror_UpStream, phiBins);
  Hist_LRAsym(h_phi_DownStream, h_LR_phi_DownStream, phiBins);
  Hist_LR_MirrorAsym(h_phi_DownStream, h_LR_phiMirror_DownStream, phiBins);
  Hist_LRAsym(h_phi, h_LR_phi, phiBins);
  Hist_LR_MirrorAsym(h_phi, h_LR_phiMirror, phiBins);
  // }}}
  
  //TGraphs
  // {{{
  TVectorD tv_xN_xval( *( (TVectorD*)fdata->Get("tv_xN_xval") ) );
  TVectorD tv_xPi_xval( *( (TVectorD*)fdata->Get("tv_xPi_xval") ) );
  TVectorD tv_xF_xval( *( (TVectorD*)fdata->Get("tv_xF_xval") ) );
  TVectorD tv_pT_xval( *( (TVectorD*)fdata->Get("tv_pT_xval") ) );
  TVectorD tv_M_xval( *( (TVectorD*)fdata->Get("tv_M_xval") ) );
  Double_t xval_xN[nBins], xval_xPi[nBins], xval_xF[nBins];
  Double_t xval_pT[nBins], xval_M[nBins];//openAngle
  Double_t ex[nBins];
  for (Int_t i=0; i<nBins; i++) {
    xval_xN[i] = tv_xN_xval[i]; xval_xPi[i] = tv_xPi_xval[i];
    xval_xF[i]=tv_xF_xval[i]; xval_pT[i]=tv_pT_xval[i]; xval_M[i]=tv_M_xval[i];//openAngle
    ex[i] = 0.0;
  }

  //Double_t xval_pT[nBins];//openAngle
  //Double_t xval_M[] = {0.03, 0.07, 0.1, 0.115, 0.13, 0.155, 0.2};//openAngle

  //Double_t xval_xN[] = {0.100881, 0.157517, 0.248419}; //High Mass binning
  //Double_t xval_xPi[] = {0.314897, 0.476024, 0.684437};
  //Double_t xval_xF[] = {0.0876318, 0.307122, 0.556132};
  //Double_t xval_pT[] = {0.61, 1.18, 1.98};
  //Double_t xval_M[] = {4.50638, 5.07645, 6.4052};
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
  // }}}

  ////////////////
  //Draw and pretty up graphs
  ////////////////
  //Binned physics values
  ///////////////
  // {{{
  TCanvas* cPhys = new TCanvas();
  cPhys->Divide(3, 5);
  cPhys->cd(1);
  Double_t xmin = gr_xN_LR_UpStream->GetXaxis()->GetXmin();
  Double_t xmax = gr_xN_LR_UpStream->GetXaxis()->GetXmax();
  TLine *l_xN_UpStream=new TLine(xmin,0.0,xmax,0.0);
  gr_xN_LR_UpStream->Draw("AP");
  l_xN_UpStream->Draw("same");
  
  cPhys->cd(2);
  xmin = gr_xN_LR_DownStream->GetXaxis()->GetXmin();
  xmax = gr_xN_LR_DownStream->GetXaxis()->GetXmax();
  TLine *l_xN_DownStream=new TLine(xmin,0.0,xmax,0.0);
  gr_xN_LR_DownStream->Draw("AP");
  l_xN_DownStream->Draw("same");
  cPhys->cd(3);
  
  xmin = gr_xN_LR->GetXaxis()->GetXmin();
  xmax = gr_xN_LR->GetXaxis()->GetXmax();
  TLine *l_xN=new TLine(xmin,0.0,xmax,0.0);
  gr_xN_LR->Draw("AP");
  l_xN->Draw("same");

  SetupTGraph(gr_xN_LR_UpStream, "UpStream", "xN", 1);//xN Setups
  SetupTLine(l_xN_UpStream);
  SetupTGraph(gr_xN_LR_DownStream, "DownStream", "xN", 1);
  SetupTLine(l_xN_DownStream);
  SetupTGraph(gr_xN_LR, "Both Targets", "xN", 1);
  SetupTLine(l_xN);
  
  
  cPhys->cd(4);
  xmin = gr_xPi_LR_UpStream->GetXaxis()->GetXmin();
  xmax = gr_xPi_LR_UpStream->GetXaxis()->GetXmax();
  TLine *l_xPi_UpStream=new TLine(xmin,0.0,xmax,0.0);
  gr_xPi_LR_UpStream->Draw("AP");
  l_xPi_UpStream->Draw("same");
  
  cPhys->cd(5);
  xmin = gr_xPi_LR_DownStream->GetXaxis()->GetXmin();
  xmax = gr_xPi_LR_DownStream->GetXaxis()->GetXmax();
  TLine *l_xPi_DownStream=new TLine(xmin,0.0,xmax,0.0);
  gr_xPi_LR_DownStream->Draw("AP");
  l_xPi_DownStream->Draw("same");
  
  cPhys->cd(6);
  xmin = gr_xPi_LR->GetXaxis()->GetXmin();
  xmax = gr_xPi_LR->GetXaxis()->GetXmax();
  TLine *l_xPi=new TLine(xmin,0.0,xmax,0.0);
  gr_xPi_LR->Draw("AP");
  l_xPi->Draw("same");

  SetupTGraph(gr_xPi_LR_UpStream, "UpStream", "xPi", 1.0);//xPi Setups
  SetupTLine(l_xPi_UpStream);
  SetupTGraph(gr_xPi_LR_DownStream, "DownStream", "xPi", 1.0);
  SetupTLine(l_xPi_DownStream);
  SetupTGraph(gr_xPi_LR, "Both Targets", "xPi", 1.0);
  SetupTLine(l_xPi);
  
  
  cPhys->cd(7);
  xmin = gr_xF_LR_UpStream->GetXaxis()->GetXmin();
  xmax = gr_xF_LR_UpStream->GetXaxis()->GetXmax();
  TLine *l_xF_UpStream=new TLine(xmin,0.0,xmax,0.0);
  gr_xF_LR_UpStream->Draw("AP");
  l_xF_UpStream->Draw("same");
  
  cPhys->cd(8);
  xmin = gr_xF_LR_DownStream->GetXaxis()->GetXmin();
  xmax = gr_xF_LR_DownStream->GetXaxis()->GetXmax();
  TLine *l_xF_DownStream=new TLine(xmin,0.0,xmax,0.0);
  gr_xF_LR_DownStream->Draw("AP");
  l_xF_DownStream->Draw("same");
  
  cPhys->cd(9);
  xmin = gr_xF_LR->GetXaxis()->GetXmin();
  xmax = gr_xF_LR->GetXaxis()->GetXmax();
  TLine *l_xF=new TLine(xmin,0.0,xmax,0.0);
  gr_xF_LR->Draw("AP");
  l_xF->Draw("same");

  SetupTGraph(gr_xF_LR_UpStream, "UpStream", "xF", 1.0);//xF Setups
  SetupTLine(l_xF_UpStream);
  SetupTGraph(gr_xF_LR_DownStream, "DownStream", "xF", 1.0);
  SetupTLine(l_xF_DownStream);
  SetupTGraph(gr_xF_LR, "Both Targets", "xF", 1.0);
  SetupTLine(l_xF);
  

  cPhys->cd(10);
  xmin = gr_pT_LR_UpStream->GetXaxis()->GetXmin();
  xmax = gr_pT_LR_UpStream->GetXaxis()->GetXmax();
  TLine *l_pT_UpStream=new TLine(xmin,0.0,xmax,0.0);
  gr_pT_LR_UpStream->Draw("AP");
  l_pT_UpStream->Draw("same");
  
  cPhys->cd(11);
  xmin = gr_pT_LR_DownStream->GetXaxis()->GetXmin();
  xmax = gr_pT_LR_DownStream->GetXaxis()->GetXmax();
  TLine *l_pT_DownStream=new TLine(xmin,0.0,xmax,0.0);
  gr_pT_LR_DownStream->Draw("AP");
  l_pT_DownStream->Draw("same");
  
  cPhys->cd(12);
  xmin = gr_pT_LR->GetXaxis()->GetXmin();
  xmax = gr_pT_LR->GetXaxis()->GetXmax();
  TLine *l_pT=new TLine(xmin,0.0,xmax,0.0);
  gr_pT_LR->Draw("AP");
  l_pT->Draw("same");

  SetupTGraph(gr_pT_LR_UpStream, "UpStream", "pT", 1.0);//pT Setups
  SetupTLine(l_pT_UpStream);
  SetupTGraph(gr_pT_LR_DownStream, "DownStream", "pT", 1.0);
  SetupTLine(l_pT_DownStream);
  SetupTGraph(gr_pT_LR, "Both Targets", "pT", 1.0);
  SetupTLine(l_pT);
  
  
  cPhys->cd(13);
  xmin = gr_M_LR_UpStream->GetXaxis()->GetXmin();
  xmax = gr_M_LR_UpStream->GetXaxis()->GetXmax();
  TLine *l_M_UpStream=new TLine(xmin,0.0,xmax,0.0);
  gr_M_LR_UpStream->Draw("AP");
  l_M_UpStream->Draw("same");
  
  cPhys->cd(14);
  xmin = gr_M_LR_DownStream->GetXaxis()->GetXmin();
  xmax = gr_M_LR_DownStream->GetXaxis()->GetXmax();
  TLine *l_M_DownStream=new TLine(xmin,0.0,xmax,0.0);
  gr_M_LR_DownStream->Draw("AP");
  l_M_DownStream->Draw("same");
  
  cPhys->cd(15);
  xmin = gr_M_LR->GetXaxis()->GetXmin();
  xmax = gr_M_LR->GetXaxis()->GetXmax();
  TLine *l_M=new TLine(xmin,0.0,xmax,0.0);
  gr_M_LR->Draw("AP");
  l_M->Draw("same");

  SetupTGraph(gr_M_LR_UpStream, "UpStream", "M", 1.0);//M Setups
  SetupTLine(l_M_UpStream);
  SetupTGraph(gr_M_LR_DownStream, "DownStream", "M", 1.0);
  SetupTLine(l_M_DownStream);
  SetupTGraph(gr_M_LR, "Both Targets", "M", 1.0);
  SetupTLine(l_M);
  // }}}
  
  //TGraphs by target
  ///////////////
  // {{{
  TGraphErrors* gr_xN_LR_UpStream_Up = new TGraphErrors(nBins, xval_xN,
						     xN_Asym_UpStream_Up, ex,
						     e_xN_Asym_UpStream_Up);
  TGraphErrors* gr_xN_LR_UpStream_Down = new TGraphErrors(nBins, xval_xN,
						     xN_Asym_UpStream_Down, ex,
						     e_xN_Asym_UpStream_Down);
  TGraphErrors* gr_xN_LR_DownStream_Up = new TGraphErrors(nBins, xval_xN,
						     xN_Asym_DownStream_Up, ex,
						     e_xN_Asym_DownStream_Up);
  TGraphErrors* gr_xN_LR_DownStream_Down = new TGraphErrors(nBins, xval_xN,
						     xN_Asym_DownStream_Down,ex,
						     e_xN_Asym_DownStream_Down);

  Int_t ipad = 1;
  TCanvas *cxN = new TCanvas();
  cxN->Divide(2, 2);
  cxN->cd(ipad); ipad++;
  gr_xN_LR_UpStream_Up->Draw("AP");
  SetupTGraph(gr_xN_LR_UpStream_Up, "UpStream Pol Up", "xN", 1.0);
  xmin = gr_xN_LR_UpStream_Up->GetXaxis()->GetXmin();
  xmax = gr_xN_LR_UpStream_Up->GetXaxis()->GetXmax();
  TLine *l_xN_LR_UpStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_xN_LR_UpStream_Up->Draw("same");
  SetupTLine(l_xN_LR_UpStream_Up);
  cxN->cd(ipad); ipad++;
  gr_xN_LR_UpStream_Down->Draw("AP");
  SetupTGraph(gr_xN_LR_UpStream_Down, "UpStream Pol Down", "xN");
  TLine *l_xN_LR_UpStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_xN_LR_UpStream_Down->Draw("same");
  SetupTLine(l_xN_LR_UpStream_Down);
  cxN->cd(ipad); ipad++;
  gr_xN_LR_DownStream_Up->Draw("AP");
  SetupTGraph(gr_xN_LR_DownStream_Up, "DownStream Pol Up", "xN");
  TLine *l_xN_LR_DownStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_xN_LR_DownStream_Up->Draw("same");
  SetupTLine(l_xN_LR_DownStream_Up);
  cxN->cd(ipad); ipad++;
  gr_xN_LR_DownStream_Down->Draw("AP");
  SetupTGraph(gr_xN_LR_DownStream_Down, "DownStream Pol Down", "xN");
  TLine *l_xN_LR_DownStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_xN_LR_DownStream_Down->Draw("same");
  SetupTLine(l_xN_LR_DownStream_Down);

  
  TGraphErrors* gr_xPi_LR_UpStream_Up = new TGraphErrors(nBins, xval_xPi,
						     xPi_Asym_UpStream_Up, ex,
						     e_xPi_Asym_UpStream_Up);
  TGraphErrors* gr_xPi_LR_UpStream_Down = new TGraphErrors(nBins, xval_xPi,
						     xPi_Asym_UpStream_Down, ex,
						     e_xPi_Asym_UpStream_Down);
  TGraphErrors* gr_xPi_LR_DownStream_Up = new TGraphErrors(nBins, xval_xPi,
						     xPi_Asym_DownStream_Up, ex,
						     e_xPi_Asym_DownStream_Up);
  TGraphErrors* gr_xPi_LR_DownStream_Down = new TGraphErrors(nBins, xval_xPi,
						     xPi_Asym_DownStream_Down,ex,
						     e_xPi_Asym_DownStream_Down);

  ipad = 1;
  TCanvas *cxPi = new TCanvas();
  cxPi->Divide(2, 2);
  cxPi->cd(ipad); ipad++;
  gr_xPi_LR_UpStream_Up->Draw("AP");
  SetupTGraph(gr_xPi_LR_UpStream_Up, "UpStream Pol Up", "xPi", 1.0);
  xmin = gr_xPi_LR_UpStream_Up->GetXaxis()->GetXmin();
  xmax = gr_xPi_LR_UpStream_Up->GetXaxis()->GetXmax();
  TLine *l_xPi_LR_UpStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_xPi_LR_UpStream_Up->Draw("same");
  SetupTLine(l_xPi_LR_UpStream_Up);
  cxPi->cd(ipad); ipad++;
  gr_xPi_LR_UpStream_Down->Draw("AP");
  SetupTGraph(gr_xPi_LR_UpStream_Down, "UpStream Pol Down", "xPi");
  TLine *l_xPi_LR_UpStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_xPi_LR_UpStream_Down->Draw("same");
  SetupTLine(l_xPi_LR_UpStream_Down);
  cxPi->cd(ipad); ipad++;
  gr_xPi_LR_DownStream_Up->Draw("AP");
  SetupTGraph(gr_xPi_LR_DownStream_Up, "DownStream Pol Up", "xPi");
  TLine *l_xPi_LR_DownStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_xPi_LR_DownStream_Up->Draw("same");
  SetupTLine(l_xPi_LR_DownStream_Up);
  cxPi->cd(ipad); ipad++;
  gr_xPi_LR_DownStream_Down->Draw("AP");
  SetupTGraph(gr_xPi_LR_DownStream_Down, "DownStream Pol Down", "xPi");
  TLine *l_xPi_LR_DownStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_xPi_LR_DownStream_Down->Draw("same");
  SetupTLine(l_xPi_LR_DownStream_Down);


  TGraphErrors* gr_xF_LR_UpStream_Up = new TGraphErrors(nBins, xval_xF,
						     xF_Asym_UpStream_Up, ex,
						     e_xF_Asym_UpStream_Up);
  TGraphErrors* gr_xF_LR_UpStream_Down = new TGraphErrors(nBins, xval_xF,
						     xF_Asym_UpStream_Down, ex,
						     e_xF_Asym_UpStream_Down);
  TGraphErrors* gr_xF_LR_DownStream_Up = new TGraphErrors(nBins, xval_xF,
						     xF_Asym_DownStream_Up, ex,
						     e_xF_Asym_DownStream_Up);
  TGraphErrors* gr_xF_LR_DownStream_Down = new TGraphErrors(nBins, xval_xF,
						     xF_Asym_DownStream_Down,ex,
						     e_xF_Asym_DownStream_Down);

  ipad = 1;
  TCanvas *cxF = new TCanvas();
  cxF->Divide(2, 2);
  cxF->cd(ipad); ipad++;
  gr_xF_LR_UpStream_Up->Draw("AP");
  SetupTGraph(gr_xF_LR_UpStream_Up, "UpStream Pol Up", "xF", 1.0);
  xmin = gr_xF_LR_UpStream_Up->GetXaxis()->GetXmin();
  xmax = gr_xF_LR_UpStream_Up->GetXaxis()->GetXmax();
  TLine *l_xF_LR_UpStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_xF_LR_UpStream_Up->Draw("same");
  SetupTLine(l_xF_LR_UpStream_Up);
  cxF->cd(ipad); ipad++;
  gr_xF_LR_UpStream_Down->Draw("AP");
  SetupTGraph(gr_xF_LR_UpStream_Down, "UpStream Pol Down", "xF");
  TLine *l_xF_LR_UpStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_xF_LR_UpStream_Down->Draw("same");
  SetupTLine(l_xF_LR_UpStream_Down);
  cxF->cd(ipad); ipad++;
  gr_xF_LR_DownStream_Up->Draw("AP");
  SetupTGraph(gr_xF_LR_DownStream_Up, "DownStream Pol Up", "xF");
  TLine *l_xF_LR_DownStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_xF_LR_DownStream_Up->Draw("same");
  SetupTLine(l_xF_LR_DownStream_Up);
  cxF->cd(ipad); ipad++;
  gr_xF_LR_DownStream_Down->Draw("AP");
  SetupTGraph(gr_xF_LR_DownStream_Down, "DownStream Pol Down", "xF");
  TLine *l_xF_LR_DownStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_xF_LR_DownStream_Down->Draw("same");
  SetupTLine(l_xF_LR_DownStream_Down);


  TGraphErrors* gr_pT_LR_UpStream_Up = new TGraphErrors(nBins, xval_pT,
						     pT_Asym_UpStream_Up, ex,
						     e_pT_Asym_UpStream_Up);
  TGraphErrors* gr_pT_LR_UpStream_Down = new TGraphErrors(nBins, xval_pT,
						     pT_Asym_UpStream_Down, ex,
						     e_pT_Asym_UpStream_Down);
  TGraphErrors* gr_pT_LR_DownStream_Up = new TGraphErrors(nBins, xval_pT,
						     pT_Asym_DownStream_Up, ex,
						     e_pT_Asym_DownStream_Up);
  TGraphErrors* gr_pT_LR_DownStream_Down = new TGraphErrors(nBins, xval_pT,
						     pT_Asym_DownStream_Down,ex,
						     e_pT_Asym_DownStream_Down);

  ipad = 1;
  TCanvas *cpT = new TCanvas();
  cpT->Divide(2, 2);
  cpT->cd(ipad); ipad++;
  gr_pT_LR_UpStream_Up->Draw("AP");
  SetupTGraph(gr_pT_LR_UpStream_Up, "UpStream Pol Up", "pT", 1.0);
  xmin = gr_pT_LR_UpStream_Up->GetXaxis()->GetXmin();
  xmax = gr_pT_LR_UpStream_Up->GetXaxis()->GetXmax();
  TLine *l_pT_LR_UpStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_pT_LR_UpStream_Up->Draw("same");
  SetupTLine(l_pT_LR_UpStream_Up);
  cpT->cd(ipad); ipad++;
  gr_pT_LR_UpStream_Down->Draw("AP");
  SetupTGraph(gr_pT_LR_UpStream_Down, "UpStream Pol Down", "pT");
  TLine *l_pT_LR_UpStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_pT_LR_UpStream_Down->Draw("same");
  SetupTLine(l_pT_LR_UpStream_Down);
  cpT->cd(ipad); ipad++;
  gr_pT_LR_DownStream_Up->Draw("AP");
  SetupTGraph(gr_pT_LR_DownStream_Up, "DownStream Pol Up", "pT");
  TLine *l_pT_LR_DownStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_pT_LR_DownStream_Up->Draw("same");
  SetupTLine(l_pT_LR_DownStream_Up);
  cpT->cd(ipad); ipad++;
  gr_pT_LR_DownStream_Down->Draw("AP");
  SetupTGraph(gr_pT_LR_DownStream_Down, "DownStream Pol Down", "pT");
  TLine *l_pT_LR_DownStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_pT_LR_DownStream_Down->Draw("same");
  SetupTLine(l_pT_LR_DownStream_Down);


  TGraphErrors* gr_M_LR_UpStream_Up = new TGraphErrors(nBins, xval_M,
						     M_Asym_UpStream_Up, ex,
						     e_M_Asym_UpStream_Up);
  TGraphErrors* gr_M_LR_UpStream_Down = new TGraphErrors(nBins, xval_M,
						     M_Asym_UpStream_Down, ex,
						     e_M_Asym_UpStream_Down);
  TGraphErrors* gr_M_LR_DownStream_Up = new TGraphErrors(nBins, xval_M,
						     M_Asym_DownStream_Up, ex,
						     e_M_Asym_DownStream_Up);
  TGraphErrors* gr_M_LR_DownStream_Down = new TGraphErrors(nBins, xval_M,
						     M_Asym_DownStream_Down,ex,
						     e_M_Asym_DownStream_Down);

  ipad = 1;
  TCanvas *cM = new TCanvas();
  cM->Divide(2, 2);
  cM->cd(ipad); ipad++;
  gr_M_LR_UpStream_Up->Draw("AP");
  SetupTGraph(gr_M_LR_UpStream_Up, "UpStream Pol Up", "M", 1.0);
  xmin = gr_M_LR_UpStream_Up->GetXaxis()->GetXmin();
  xmax = gr_M_LR_UpStream_Up->GetXaxis()->GetXmax();
  TLine *l_M_LR_UpStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_M_LR_UpStream_Up->Draw("same");
  SetupTLine(l_M_LR_UpStream_Up);
  cM->cd(ipad); ipad++;
  gr_M_LR_UpStream_Down->Draw("AP");
  SetupTGraph(gr_M_LR_UpStream_Down, "UpStream Pol Down", "M");
  TLine *l_M_LR_UpStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_M_LR_UpStream_Down->Draw("same");
  SetupTLine(l_M_LR_UpStream_Down);
  cM->cd(ipad); ipad++;
  gr_M_LR_DownStream_Up->Draw("AP");
  SetupTGraph(gr_M_LR_DownStream_Up, "DownStream Pol Up", "M");
  TLine *l_M_LR_DownStream_Up = new TLine(xmin,0.0,xmax,0.0);
  l_M_LR_DownStream_Up->Draw("same");
  SetupTLine(l_M_LR_DownStream_Up);
  cM->cd(ipad); ipad++;
  gr_M_LR_DownStream_Down->Draw("AP");
  SetupTGraph(gr_M_LR_DownStream_Down, "DownStream Pol Down", "M");
  TLine *l_M_LR_DownStream_Down = new TLine(xmin,0.0,xmax,0.0);
  l_M_LR_DownStream_Down->Draw("same");
  SetupTLine(l_M_LR_DownStream_Down);
  // }}}
  
  //Write output
  ///////////////
  // {{{
  if (Qflag || wflag){
    TFile* fout = (Qflag ? new TFile(outFile, "RECREATE")
		   : new TFile("Output.root", "RECREATE") );
    
    gr_xN_LR_UpStream->Write("gr_xN_LR_UpStream");
    gr_xN_LR_DownStream->Write("gr_xN_LR_DownStream");
    gr_xN_LR->Write("gr_xN_LR");

    gr_xN_LR_UpStream_Up->Write("gr_xN_LR_UpStream_Up");//By target
    gr_xN_LR_UpStream_Down->Write("gr_xN_LR_UpStream_Down");
    gr_xN_LR_DownStream_Up->Write("gr_xN_LR_DownStream_Up");
    gr_xN_LR_DownStream_Down->Write("gr_xN_LR_DownStream_Down");

    gr_xPi_LR_UpStream->Write("gr_xPi_LR_UpStream");
    gr_xPi_LR_DownStream->Write("gr_xPi_LR_DownStream");
    gr_xPi_LR->Write("gr_xPi_LR");

    gr_xPi_LR_UpStream_Up->Write("gr_xPi_LR_UpStream_Up");//By target
    gr_xPi_LR_UpStream_Down->Write("gr_xPi_LR_UpStream_Down");
    gr_xPi_LR_DownStream_Up->Write("gr_xPi_LR_DownStream_Up");
    gr_xPi_LR_DownStream_Down->Write("gr_xPi_LR_DownStream_Down");

    gr_xF_LR_UpStream->Write("gr_xF_LR_UpStream");
    gr_xF_LR_DownStream->Write("gr_xF_LR_DownStream");
    gr_xF_LR->Write("gr_xF_LR");

    gr_xF_LR_UpStream_Up->Write("gr_xF_LR_UpStream_Up");//By target
    gr_xF_LR_UpStream_Down->Write("gr_xF_LR_UpStream_Down");
    gr_xF_LR_DownStream_Up->Write("gr_xF_LR_DownStream_Up");
    gr_xF_LR_DownStream_Down->Write("gr_xF_LR_DownStream_Down");

    gr_pT_LR_UpStream->Write("gr_pT_LR_UpStream");
    gr_pT_LR_DownStream->Write("gr_pT_LR_DownStream");
    gr_pT_LR->Write("gr_pT_LR");

    gr_pT_LR_UpStream_Up->Write("gr_pT_LR_UpStream_Up");//By target
    gr_pT_LR_UpStream_Down->Write("gr_pT_LR_UpStream_Down");
    gr_pT_LR_DownStream_Up->Write("gr_pT_LR_DownStream_Up");
    gr_pT_LR_DownStream_Down->Write("gr_pT_LR_DownStream_Down");

    gr_M_LR_UpStream->Write("gr_M_LR_UpStream");
    gr_M_LR_DownStream->Write("gr_M_LR_DownStream");
    gr_M_LR->Write("gr_M_LR");

    gr_M_LR_UpStream_Up->Write("gr_M_LR_UpStream_Up");//By target
    gr_M_LR_UpStream_Down->Write("gr_M_LR_UpStream_Down");
    gr_M_LR_DownStream_Up->Write("gr_M_LR_DownStream_Up");
    gr_M_LR_DownStream_Down->Write("gr_M_LR_DownStream_Down");
  
    h_phi_UpStream->Write();
    h_phi_DownStream->Write();
    h_phi->Write();

    h_LR_phi_UpStream->Write();
    h_LR_phi_DownStream->Write();
    h_LR_phi->Write();

    h_LR_phiMirror_UpStream->Write();
    h_LR_phiMirror_DownStream->Write();
    h_LR_phiMirror->Write();
  
    fout->Close();

    cout << " " << endl;
    if (Qflag) cout << outFile << " file written" << endl;
    else cout << "Output.root file written" << endl;
  }
  // }}}
    
  theApp.Run();//Needed to make root graphics work on C++
}
