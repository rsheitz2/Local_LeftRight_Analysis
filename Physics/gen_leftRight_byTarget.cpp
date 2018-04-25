#include "common.h"
#include "functions.h"
#include "setup.h"

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-bbinfile || -BbinVar] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << " " << endl;
    cout << "----Changing Binning Options----" << endl;
    cout << "Option:  -b binFile with binning information	";
    cout << "(textfile should be made from "
	 << "Macro/Binning/genAvgBinBounds.C)" << endl;
    cout << "Option:  -M (\"HM\", \"JPsi\", \"AMDY\") to specify which mass ";
    cout << "range to use for \"binning information\" " << endl;
    cout << "   (default when only \"-b\" is used=AMDY)" << endl;
    cout << "Option   -B binVar  " 
	 << "(Additional binning variable not in binFile)" << endl;
    cout << "   (Current options: openAngle, phi_pIn, theta_pIn, "
	 << "qP_pIn, phi_photon)" << endl;
    cout << "Option   -N bins" << endl;
    cout << "    (Option can only be used with \"-B\" option for binVar)"<<endl;
    cout << "" << endl;
    cout << "----Writing Options----"<< endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << " " << endl;
    cout << "----Additional Options----" << endl;
    cout << "Option:  -D    (Debug mode, only 1000 events considered)" << endl;
    cout << "Option:  -p                (no cout statements in tree loop";
    cout << "to speed things up)" << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t binFlag=0, Mflag=0, wflag=0, Qflag=0, fflag=0, Dflag=0, pflag=0,Bflag=0;
  Int_t Nflag=0;
  Int_t c;
  TString fname="", outFile="", binFile="", massRange="", binVar="";
  Int_t NVar=0;
  
  while ((c = getopt (argc, argv, "b:M:B:N:wQ:Dpf:")) != -1) {
    switch (c) {
    case 'b':
      binFlag = 1;
      binFile += optarg;
      break;
    case 'M':
      Mflag = 1;
      massRange += optarg;
      break;
    case 'B':
      Bflag = 1;
      binVar += optarg;
      break;
    case 'N':
      Nflag = 1;
      NVar = atoi(optarg);
      break;
    case 'w':
      wflag = 1;
      break;
    case 'Q':
      Qflag = 1;
      outFile += optarg;
      break;
    case 'D':
      Dflag = 1;
      break;
    case 'p':
      pflag = 1;
      break;
    case 'f':
      fflag = 1;
      fname += optarg;
      break;
    case '?':
      if (optopt == 'b')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'M')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'B')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'N')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'Q')
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

  vector<Double_t> xN_bounds, xN_xval; 
  vector<Double_t> xPi_bounds, xPi_xval;
  vector<Double_t> xF_bounds, xF_xval; 
  vector<Double_t> pT_bounds, pT_xval;
  vector<Double_t> M_bounds, M_xval;
  vector<Double_t> rad_bounds, rad_xval;
  vector<Double_t> vxZ_upstream_bounds, vxZ_upstream_xval;
  vector<Double_t> vxZ_downstream_bounds, vxZ_downstream_xval;
  Int_t nBounds;
  if (!binFlag && !Bflag) {
    cout << "No \"-b\" or \"-B\" options used" << endl;
    exit(EXIT_FAILURE);
  }
  else if (binFlag) {
    xN_bounds.push_back(0.0);
    xPi_bounds.push_back(0.0);
    xF_bounds.push_back(-1.0);
    pT_bounds.push_back(0.0);
    rad_bounds.push_back(0.0);
    vxZ_upstream_bounds.push_back(-294.5);
    vxZ_downstream_bounds.push_back(-219.5);

    if (!Mflag || massRange=="AMDY")M_bounds.push_back(0.0);//All Mass DY
    else if (massRange=="HM") M_bounds.push_back(4.3);//High mass
    else if (massRange=="JPsi")M_bounds.push_back(2.5);//JPsi mass
    else {
      cout << "Invalid mass range specified" << endl;
      exit(EXIT_FAILURE);
    }
    
    string line;
    TString dy_type = "";
    Int_t xval = 1;
    ifstream f_bins(binFile);
    if(!f_bins.is_open() ) {
      cout << " " << endl;
      cout << "binFile: " << binFile << " did not open" << endl;
      exit(EXIT_FAILURE); }
    while (!f_bins.eof()) {
      getline(f_bins,line);

      if (line[1] == 'N') {
	if (dy_type == "xN") xval = 1;
	else {
	  dy_type = "xN";
	  xval = 0;
	}			
      }
      else if (line[1] == 'P') {
	if (dy_type == "xPi") xval = 1;
	else {
	  dy_type = "xPi";
	  xval = 0;
	}			
      }
      else if (line[1] == 'F') {
	if (dy_type == "xF") xval = 1;
	else {
	  dy_type = "xF";
	  xval = 0;
	}			
      }
      else if (line[1] == 'T') {
	if (dy_type == "pT") xval = 1;
	else {
	  dy_type = "pT";
	  xval = 0;
	}			
      }
      else if (line[2] == 's') {
	if (dy_type == "M") xval = 1;
	else {
	  dy_type = "M";
	  xval = 0;
	}			
      }
      else if (line[0] == 'r') {
	if (dy_type == "rad") xval = 1;
	else {
	  dy_type = "rad";
	  xval = 0;
	}			
      }
      else if (line[4] == 'u') {
	if (dy_type == "vxZ_upstream") xval = 1;
	else {
	  dy_type = "vxZ_upstream";
	  xval = 0;
	}			
      }
      else if (line[4] == 'd') {
	if (dy_type == "vxZ_downstream") xval = 1;
	else {
	  dy_type = "vxZ_downstream";
	  xval = 0;
	}			
      }

      //Don't read title lines
      if (line[0] == 'x' || line[0] == 'p' || line[1] == 'a' || line[0] == 'v'){
	continue;
      }
      else if (line == "-nan"){
	cout << "\"nan\" in binFile" << endl;
	exit(EXIT_FAILURE);
      }
      
      if (dy_type == "xN"){
	if (xval == 0) xN_bounds.push_back(atof(line.c_str() ) );
	else xN_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "xPi"){
	if (xval == 0) xPi_bounds.push_back(atof(line.c_str() ) );
	else xPi_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "xF"){
	if (xval == 0) xF_bounds.push_back(atof(line.c_str() ) );
	else xF_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "pT"){
	if (xval == 0) pT_bounds.push_back(atof(line.c_str() ) );
	else pT_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "M"){
	if (xval == 0) M_bounds.push_back(atof(line.c_str() ) );
	else M_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "rad"){
	if (xval == 0) rad_bounds.push_back(atof(line.c_str() ) );
	else rad_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "vxZ_upstream"){
	if (xval == 0) vxZ_upstream_bounds.push_back(atof(line.c_str() ) );
	else vxZ_upstream_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "vxZ_downstream"){
	if (xval == 0) vxZ_downstream_bounds.push_back(atof(line.c_str() ) );
	else vxZ_downstream_xval.push_back(atof(line.c_str() ) );
      }
    }//end file loop

    xN_bounds.push_back(1.0);
    xPi_bounds.push_back(1.0);
    xF_bounds.push_back(1.0);
    pT_bounds.push_back(7.0);
    rad_bounds.push_back(1.9);
    vxZ_upstream_bounds.push_back(-239.3);
    vxZ_downstream_bounds.push_back(-164.3);
    if(xN_xval.size()==0 || xPi_xval.size()==0 || xF_xval.size()==0 ||
       pT_xval.size()==0 || M_xval.size()==0 || rad_xval.size()==0 ||
       vxZ_upstream_xval.size()==0 || vxZ_downstream_xval.size()==0){
      cout << "Error:" << endl;
      cout << "Modern xval values not specifed in " << binFile << endl;
      cout << " " << endl;
      exit(EXIT_FAILURE);
    }

    if (!Mflag || massRange=="AMDY")M_bounds.push_back(16.0);//All Mass DY
    else if (massRange=="HM") M_bounds.push_back(8.5);//High mass
    else if (massRange=="JPsi")M_bounds.push_back(4.3);//JPsi mass
    cout << " " << endl;
    cout << "Mass range set to:" << endl;
    (!Mflag) ? cout << "AMDY" << endl : cout << massRange << endl;
    cout << " " << endl;

    if (M_bounds.at(0) > M_bounds.at(1) || M_bounds.back() < M_xval.front() ){
      cout << "Mass Range not setup correct" << endl;
      exit(EXIT_FAILURE);
    }

    nBounds = xN_bounds.size();
  }//end binFlag
  
  Int_t nVar_nBounds = 0;
  if (Bflag && !Nflag) {
    cout << "\"-B\" option used but no \"-N\" option" << endl;
    exit(EXIT_FAILURE);
  }
  else if (Nflag && !NVar){
    cout << " " << endl;
    cout << "Please enter and integer with -N option" << endl;
    exit(EXIT_FAILURE);
  }
  else nVar_nBounds = NVar + 1;
  Double_t upS_bounds[nVar_nBounds], downS_bounds[nVar_nBounds];

  if(Dflag) {
    cout << " " << endl;
    cout << "Debug mode only 1000 events considered" << endl;
    cout << " " << endl;
  }
  // }}}

  //Opening data files/setting up tree
  ///////////////
  // {{{
  TFile *fdata = TFile::Open(fname);
  TChain *tree = new TChain("Events");
  TList *liData = fdata->GetListOfKeys();
  TIter iter( liData->MakeIterator() );
  while(TObject* obj = iter()){
    TKey* theKey = (TKey*)obj;
    if (strncmp (theKey->GetClassName(),"TTree",4) == 0){
      tree->Add( fname+"/"+obj->GetName()+";"+Form("%i",theKey->GetCycle()) );
    }
  }
    
  //Vertex specific
  Double_t vx_z;
  Int_t targetPosition;
  //Drell-Yan Angles
  Double_t PhiS_simple, Theta_CS;
  //Virtual Photon
  Double_t muM_X, muM_Y, muM_Z, muM_E;
  Double_t muP_X, muP_Y, muP_Z, muP_E;
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
  Double_t vOpenAngle;
  //Beam variables
  Double_t phi_pIn, theta_pIn, qP_pIn;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;
  
  //Vertex specific
  tree->SetBranchAddress("vx_z", &vx_z);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("Theta_CS", &Theta_CS);
  //Virtual Photon
  tree->SetBranchAddress("muM_X", &muM_X);
  tree->SetBranchAddress("muM_Y", &muM_Y);
  tree->SetBranchAddress("muM_Z", &muM_Z);
  tree->SetBranchAddress("muM_E", &muM_E);
  tree->SetBranchAddress("muP_X", &muP_X);
  tree->SetBranchAddress("muP_Y", &muP_Y);
  tree->SetBranchAddress("muP_Z", &muP_Z);
  tree->SetBranchAddress("muP_E", &muP_E);
  tree->SetBranchAddress("vPhoton_X", &vPhoton_X);
  tree->SetBranchAddress("vPhoton_Y", &vPhoton_Y);
  tree->SetBranchAddress("vPhoton_Z", &vPhoton_Z);
  tree->SetBranchAddress("vPhoton_E", &vPhoton_E);
  tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Beam variables
  tree->SetBranchAddress("phi_pIn", &phi_pIn);
  tree->SetBranchAddress("theta_pIn", &theta_pIn);
  tree->SetBranchAddress("qP_pIn", &qP_pIn);
  //DY-variables
  tree->SetBranchAddress("x_beam", &x_beam);
  tree->SetBranchAddress("x_target", &x_target);
  tree->SetBranchAddress("x_feynman", &x_feynman);
  tree->SetBranchAddress("q_transverse", &q_transverse);
  tree->SetBranchAddress("Mmumu", &Mmumu);

  //Generic bound value
  Double_t *boundValue;
  if (!Bflag) {}
  else if (binVar=="openAngle") boundValue = &vOpenAngle;
  else if (binVar=="phi_pIn") boundValue = &phi_pIn;
  else if (binVar=="theta_pIn") boundValue = &theta_pIn;
  else if (binVar=="qP_pIn") boundValue = &qP_pIn;
  else if (binVar=="phi_photon") {}//done in loops
  else {
    cout << " " << endl;
    cout << "\"-B\" " << binVar << " option not supported" << endl;
    exit(EXIT_FAILURE);
  }

  // }}}

  //Kinematic Counting Setup
  ///////////////
  // {{{
  Int_t nBins = (binFlag) ? nBounds-1 : 0;
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

  //Generic binning
  unsigned long long Left_UpStream[NVar], Right_UpStream[NVar];
  unsigned long long Left_DownStream[NVar], Right_DownStream[NVar];

  //By target
  unsigned long long Left_UpStream_Up[NVar], Right_UpStream_Up[NVar];
  unsigned long long Left_UpStream_Down[NVar],Right_UpStream_Down[NVar];
  unsigned long long Left_DownStream_Up[NVar], Right_DownStream_Up[NVar];
  unsigned long long Left_DownStream_Down[NVar],Right_DownStream_Down[NVar];
  
  for (Int_t i=0; i<NVar; i++) {
    Left_UpStream[i]=0; Right_UpStream[i]=0;
    Left_DownStream[i]=0; Right_DownStream[i]=0;

    //By target
    Left_UpStream_Up[i]=0; Right_UpStream_Up[i]=0;
    Left_UpStream_Down[i]=0; Right_UpStream_Down[i]=0;
    Left_DownStream_Up[i]=0; Right_DownStream_Up[i]=0;
    Left_DownStream_Down[i]=0; Right_DownStream_Down[i]=0;
  }

  const Int_t num_1Dh = 3;
  TH1D *h1[num_1Dh];
  //h1[0] = new TH1D("h_vPhoton", "h_vPhoton", 100, -TMath::Pi(), 2*TMath::Pi() );
  h1[0] = new TH1D("h_vPhoton", "h_vPhoton", 200, 0, 0.01);
  h1[1] = new TH1D("h_vPhoton_1", "h_vPhoton_1", 200, 0, 0.01);
  h1[2] = new TH1D("h_vPhoton_2", "h_vPhoton_2", 200, 0, 0.01);

  const Int_t num_2Dh = 2;
  TH2D *h2[num_2Dh];
  h2[0] = new TH2D("h_vPmag_vPhoton", "h_vPmag_vPhoton", 100, 0, 5,
		   100, -TMath::Pi(), 2*TMath::Pi() );
  h2[1] = new TH2D("h_openAngle_vPhoton", "h_openAngle_vPhoton", 100, 0, 3,
		   100, -TMath::Pi(), 2*TMath::Pi() );
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

  TH1D* h_phi_UpStream_Up = new TH1D("h_phi_UpStream_Up", "h_phi_UpStream_Up",
				     phiBins, -TMath::Pi()/2, 3*TMath::Pi()/2);
  TH1D* h_phi_UpStream_Down = new TH1D("h_phi_UpStream_Down",
				       "h_phi_UpStream_Down", phiBins,
				       -TMath::Pi()/2, 3*TMath::Pi()/2);
  TH1D* h_phi_DownStream_Up = new TH1D("h_phi_DownStream_Up",
				       "h_phi_DownStream_Up", phiBins,
				       -TMath::Pi()/2, 3*TMath::Pi()/2);
  TH1D* h_phi_DownStream_Down = new TH1D("h_phi_DownStream_Down",
					 "h_phi_DownStream_Down", phiBins,
					 -TMath::Pi()/2, 3*TMath::Pi()/2);
  // }}}
  
  //First tree loop
  //Equal out by target data
  // {{{

  Int_t nUpStream=0, nDownStream=0;
  Int_t tree_entries = (!Dflag) ? tree->GetEntries() : 1000;//Tree Loop

  vector<Double_t> upS_sort_val, downS_sort_val;
  for (Int_t ev=0; ev<tree_entries; ev++) {
      tree->GetEntry(ev, 0);

      //General useful variables
      TLorentzVector vPhoton(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);
      if (binVar=="phi_photon") *boundValue = vPhoton.Phi();
            
      if (vx_z < -229.4) nUpStream++;
      else nDownStream++;

      if (vx_z < -320.0 || vx_z > -135.0) {
	cout << "Target vx_z problems: vx_z= " << vx_z << endl;
	exit(EXIT_FAILURE);
      }

      if (Bflag) {
	if (vx_z < -229.4) upS_sort_val.push_back(*boundValue);
	else downS_sort_val.push_back(*boundValue);
      }
  }
  (nUpStream>nDownStream) ? nUpStream=nDownStream : nDownStream=nUpStream;
  cout << "Number of Entries    UpStream: " << nUpStream << "   DownStream: "
       << nDownStream << endl;

  Double_t upS_xval[NVar], downS_xval[NVar], nVar_ex[NVar];
  if (Bflag){
    std::sort(upS_sort_val.begin(), upS_sort_val.end() );
    std::sort(downS_sort_val.begin(), downS_sort_val.end() );

    Int_t upS_size = upS_sort_val.size();
    Int_t downS_size = downS_sort_val.size();
    for (Int_t i=0; i<nVar_nBounds; i++) {
      upS_bounds[i] = (i) ? upS_sort_val.at( 1.0*i*upS_size/NVar -1) :
	upS_sort_val.at(0);
      downS_bounds[i] = (i) ? downS_sort_val.at( 1.0*i*downS_size/NVar -1) :
	downS_sort_val.at(0);
    }
        
    Int_t upS_counts[NVar]; Int_t downS_counts[NVar]; 
    for (Int_t i=0; i<NVar; i++) {
      upS_counts[i] = 0;
      downS_counts[i] = 0;
    }

    Int_t ibin=0; 
    for (vector<Double_t>::iterator it=upS_sort_val.begin();
	 it!=upS_sort_val.end(); it++, ibin++) {

      if (ibin/(upS_size/NVar) < NVar){
	upS_xval[ibin/(upS_size/NVar)] += *it;
	upS_counts[ibin/(upS_size/NVar)]++;
      }
      else{
	upS_xval[NVar-1] += *it;
	upS_counts[NVar-1]++;
      }
    }

    ibin=0; 
    for (vector<Double_t>::iterator it=downS_sort_val.begin();
	 it!=downS_sort_val.end(); it++, ibin++) {

      if (ibin/(downS_size/NVar) < NVar){
	downS_xval[ibin/(downS_size/NVar)] += *it;
	downS_counts[ibin/(downS_size/NVar)]++;
      }
      else{
	downS_xval[NVar-1] += *it;
	downS_counts[NVar-1]++;
      }
    }

    for (Int_t i=0; i<NVar; i++) {
      upS_xval[i] = upS_xval[i]/(1.0*upS_counts[i]);
      downS_xval[i] = downS_xval[i]/(1.0*downS_counts[i]);
      nVar_ex[i] = 0.0;
    }
  }//Bflag
  // }}}

  //Tree loop
  Bool_t first = true;
  Int_t stopUpStream=0, stopDownStream=0;
  cout << "Number of entries in tree: " << tree->GetEntries() << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    tree->GetEntry(ev, 0);
    
    if (first || ev==tree_entries-1){//first or last
      cout << " " << endl;
      cout << "Setup!!!!!!!!!!!!!!!" << endl;
      cout << "Spin influnced left/right asymmetry" << endl;
      
      first = false;
    }
    else if (stopUpStream>=nUpStream && stopDownStream>=nDownStream ) {
      break;
    }

    ////All data after cuts
    //////////////
    //Choose Left/Right
    // {{{
    Double_t phi_photon_lab = ShiftPhiSimple(PhiS_simple);
    Bool_t Left=false, Right=false;
    //Spin influenced left/right
    if (phi_photon_lab<TMath::Pi()/2 && phi_photon_lab>-TMath::Pi()/2) {
      Left = true;}
    else if (phi_photon_lab<3*TMath::Pi()/2 && phi_photon_lab>TMath::Pi()/2){
      Right = true;}
    else {
      cout << "No Left or Right choosen" << endl;
      cout << phi_photon_lab << "  continuing" << endl;
      cout << " " << endl;
      continue;
    }
    // }}}

    //General useful variables
    TLorentzVector vPhoton(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);
    if (binVar=="phi_photon") *boundValue = vPhoton.Phi();

    //h1[0]->Fill(phi_photon_lab);
    //h1[0]->Fill(*boundValue);
    //h2[0]->Fill( TMath::Sqrt(vPhoton_X*vPhoton_X + vPhoton_Y*vPhoton_Y),
    //phi_photon_lab);
    //h2[1]->Fill(vOpenAngle,
    //phi_photon_lab);
    
    
    //Bin data
    if (vx_z < -229.4) {//UpStream target
      if (stopUpStream>=nUpStream) continue;
      else if (stopUpStream>nUpStream/2.0) {//Switch Left/Right for MC
	Bool_t tmp = Left;
	Left = Right;
	Right = tmp;
      }
      stopUpStream++;

      if (Left){//Left
	// {{{
	if (binFlag){
	  BinDataCounts(xN_Left_UpStream, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_UpStream, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_UpStream, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_UpStream, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_UpStream, Mmumu, M_bounds);
	}
	if (Bflag) BinDataCounts(Left_UpStream, NVar, *boundValue, upS_bounds);

	if (stopUpStream<nUpStream/2.0){//Polarized Up
	  if (binFlag){
	    BinDataCounts(xN_Left_UpStream_Up, x_target, xN_bounds);
	    BinDataCounts(xPi_Left_UpStream_Up, x_beam, xPi_bounds);
	    BinDataCounts(xF_Left_UpStream_Up, x_feynman, xF_bounds);
	    BinDataCounts(pT_Left_UpStream_Up, q_transverse, pT_bounds);
	    BinDataCounts(M_Left_UpStream_Up, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Left_UpStream_Up, NVar, *boundValue,
				   upS_bounds);
	}
	else if (stopUpStream>nUpStream/2.0){//Polarized Down
	  if (binFlag){
	    BinDataCounts(xN_Left_UpStream_Down, x_target, xN_bounds);
	    BinDataCounts(xPi_Left_UpStream_Down, x_beam, xPi_bounds);
	    BinDataCounts(xF_Left_UpStream_Down, x_feynman, xF_bounds);
	    BinDataCounts(pT_Left_UpStream_Down, q_transverse, pT_bounds);
	    BinDataCounts(M_Left_UpStream_Down, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Left_UpStream_Down,NVar,*boundValue,
				   upS_bounds);
	}
	// }}}
      }//End Left
      else if (Right){//Right
	// {{{
	if (binFlag){
	  BinDataCounts(xN_Right_UpStream, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_UpStream, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_UpStream, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_UpStream, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_UpStream, Mmumu, M_bounds);
	}
	if (Bflag) BinDataCounts(Right_UpStream, NVar, *boundValue, upS_bounds);

	if (stopUpStream<nUpStream/2.0){//Polarized Up
	  if (binFlag){
	    BinDataCounts(xN_Right_UpStream_Up, x_target, xN_bounds);
	    BinDataCounts(xPi_Right_UpStream_Up, x_beam, xPi_bounds);
	    BinDataCounts(xF_Right_UpStream_Up, x_feynman, xF_bounds);
	    BinDataCounts(pT_Right_UpStream_Up, q_transverse, pT_bounds);
	    BinDataCounts(M_Right_UpStream_Up, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Right_UpStream_Up, NVar,*boundValue,
				   upS_bounds);
	}
	else if (stopUpStream>nUpStream/2.0){//Polarized Down
	  if (binFlag){
	    BinDataCounts(xN_Right_UpStream_Down, x_target, xN_bounds);
	    BinDataCounts(xPi_Right_UpStream_Down, x_beam, xPi_bounds);
	    BinDataCounts(xF_Right_UpStream_Down, x_feynman, xF_bounds);
	    BinDataCounts(pT_Right_UpStream_Down, q_transverse, pT_bounds);
	    BinDataCounts(M_Right_UpStream_Down, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Right_UpStream_Down,NVar,*boundValue,
				   upS_bounds);
	}
	// }}}
      }//End Right
    }//End UpStream target
    else {//DownStream target
      if (stopDownStream>=nDownStream) continue;
      else if (stopDownStream>nDownStream/2.0) {//Switch Left/Right for MC
	Bool_t tmp = Left;
	Left = Right;
	Right = tmp;
      }
      stopDownStream++;
      
      if (Left){//Left
	// {{{
	if (binFlag){
	  BinDataCounts(xN_Left_DownStream, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_DownStream, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_DownStream, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_DownStream, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_DownStream, Mmumu, M_bounds);
	}
	if (Bflag) BinDataCounts(Left_DownStream, NVar, *boundValue,
				 downS_bounds);

	if (stopDownStream<nDownStream/2.0){//Polarized Up
	  if (binFlag){
	    BinDataCounts(xN_Left_DownStream_Up, x_target, xN_bounds);
	    BinDataCounts(xPi_Left_DownStream_Up, x_beam, xPi_bounds);
	    BinDataCounts(xF_Left_DownStream_Up, x_feynman, xF_bounds);
	    BinDataCounts(pT_Left_DownStream_Up, q_transverse, pT_bounds);
	    BinDataCounts(M_Left_DownStream_Up, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Left_DownStream_Up, NVar,*boundValue,
				   downS_bounds);
	}
	else if (stopDownStream>nDownStream/2.0){//Polarized Down
	  if (binFlag){
	    BinDataCounts(xN_Left_DownStream_Down, x_target, xN_bounds);
	    BinDataCounts(xPi_Left_DownStream_Down, x_beam, xPi_bounds);
	    BinDataCounts(xF_Left_DownStream_Down, x_feynman, xF_bounds);
	    BinDataCounts(pT_Left_DownStream_Down, q_transverse, pT_bounds);
	    BinDataCounts(M_Left_DownStream_Down, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Left_DownStream_Down,NVar,*boundValue,
				   downS_bounds);
	}
	// }}}
      }//End Left
      else if (Right){//Right
	// {{{
	if (binFlag){
	  BinDataCounts(xN_Right_DownStream, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_DownStream, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_DownStream, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_DownStream, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_DownStream, Mmumu, M_bounds);
	}
	if (Bflag) BinDataCounts(Right_DownStream, NVar, *boundValue,
				 downS_bounds);

	if (stopDownStream<nDownStream/2.0){//Polarized Up
	  if (binFlag){
	    BinDataCounts(xN_Right_DownStream_Up, x_target, xN_bounds);
	    BinDataCounts(xPi_Right_DownStream_Up, x_beam, xPi_bounds);
	    BinDataCounts(xF_Right_DownStream_Up, x_feynman, xF_bounds);
	    BinDataCounts(pT_Right_DownStream_Up, q_transverse, pT_bounds);
	    BinDataCounts(M_Right_DownStream_Up, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Right_DownStream_Up, NVar,*boundValue,
				   downS_bounds);
	}
	else if (stopDownStream>nDownStream/2.0){//Polarized Down
	  if (binFlag){
	    BinDataCounts(xN_Right_DownStream_Down, x_target, xN_bounds);
	    BinDataCounts(xPi_Right_DownStream_Down, x_beam, xPi_bounds);
	    BinDataCounts(xF_Right_DownStream_Down, x_feynman, xF_bounds);
	    BinDataCounts(pT_Right_DownStream_Down, q_transverse, pT_bounds);
	    BinDataCounts(M_Right_DownStream_Down, Mmumu, M_bounds);
	  }
	  if (Bflag) BinDataCounts(Right_DownStream_Down, NVar,*boundValue,
				   downS_bounds);
	}
	// }}}
      }//End Right
    }//End DownStream target

  }//End tree loop

  //Draw Basic Dist
  ///////////////
  /*
  // {{{
  TCanvas* cDist = new TCanvas();
  cDist->Divide(3);
  cDist->cd(1);
  h1[0]->Draw();
  cDist->cd(2);
  //h2[0]->Draw("colz");
  h1[1]->Draw("colz");
  cDist->cd(3);
  //h2[1]->Draw("colz");
  h1[2]->Draw("colz");
  // }}}
  */
  
  //Asymmetries
  ///////////////
  // {{{
  Double_t xN_Asym_UpStream[nBins], xN_Asym_DownStream[nBins], xN_Asym[nBins];
  Double_t e_xN_Asym_UpStream[nBins], e_xN_Asym_DownStream[nBins];
  Double_t e_xN_Asym[nBins];
  Double_t xN_Asym_UpStream_Up[nBins], xN_Asym_DownStream_Up[nBins]; //By target
  Double_t xN_Asym_UpStream_Down[nBins], xN_Asym_DownStream_Down[nBins];
  Double_t e_xN_Asym_UpStream_Up[nBins], e_xN_Asym_DownStream_Up[nBins];
  Double_t e_xN_Asym_UpStream_Down[nBins], e_xN_Asym_DownStream_Down[nBins];

  Double_t xPi_Asym_UpStream[nBins], xPi_Asym_DownStream[nBins], xPi_Asym[nBins];
  Double_t e_xPi_Asym_UpStream[nBins], e_xPi_Asym_DownStream[nBins];
  Double_t e_xPi_Asym[nBins];
  Double_t xPi_Asym_UpStream_Up[nBins], xPi_Asym_DownStream_Up[nBins]; //By target
  Double_t xPi_Asym_UpStream_Down[nBins], xPi_Asym_DownStream_Down[nBins];
  Double_t e_xPi_Asym_UpStream_Up[nBins], e_xPi_Asym_DownStream_Up[nBins];
  Double_t e_xPi_Asym_UpStream_Down[nBins], e_xPi_Asym_DownStream_Down[nBins];

  Double_t xF_Asym_UpStream[nBins], xF_Asym_DownStream[nBins], xF_Asym[nBins];
  Double_t e_xF_Asym_UpStream[nBins], e_xF_Asym_DownStream[nBins];
  Double_t e_xF_Asym[nBins];
  Double_t xF_Asym_UpStream_Up[nBins], xF_Asym_DownStream_Up[nBins]; //By target
  Double_t xF_Asym_UpStream_Down[nBins], xF_Asym_DownStream_Down[nBins];
  Double_t e_xF_Asym_UpStream_Up[nBins], e_xF_Asym_DownStream_Up[nBins];
  Double_t e_xF_Asym_UpStream_Down[nBins], e_xF_Asym_DownStream_Down[nBins];

  Double_t pT_Asym_UpStream[nBins], pT_Asym_DownStream[nBins], pT_Asym[nBins];
  Double_t e_pT_Asym_UpStream[nBins], e_pT_Asym_DownStream[nBins];
  Double_t e_pT_Asym[nBins];
  Double_t pT_Asym_UpStream_Up[nBins], pT_Asym_DownStream_Up[nBins]; //By target
  Double_t pT_Asym_UpStream_Down[nBins], pT_Asym_DownStream_Down[nBins];
  Double_t e_pT_Asym_UpStream_Up[nBins], e_pT_Asym_DownStream_Up[nBins];
  Double_t e_pT_Asym_UpStream_Down[nBins], e_pT_Asym_DownStream_Down[nBins];

  Double_t M_Asym_UpStream[nBins], M_Asym_DownStream[nBins], M_Asym[nBins];
  Double_t e_M_Asym_UpStream[nBins], e_M_Asym_DownStream[nBins];
  Double_t e_M_Asym[nBins];
  Double_t M_Asym_UpStream_Up[nBins], M_Asym_DownStream_Up[nBins]; //By target
  Double_t M_Asym_UpStream_Down[nBins], M_Asym_DownStream_Down[nBins];
  Double_t e_M_Asym_UpStream_Up[nBins], e_M_Asym_DownStream_Up[nBins];
  Double_t e_M_Asym_UpStream_Down[nBins], e_M_Asym_DownStream_Down[nBins];

  if (binFlag){
    BinnedLeftRight(xN_Left_UpStream, xN_Right_UpStream, xN_Left_DownStream,
		    xN_Right_DownStream, xN_Asym_UpStream, xN_Asym_DownStream,
		    xN_Asym, e_xN_Asym_UpStream, e_xN_Asym_DownStream, e_xN_Asym,
		    nBins);
    BinnedLeftRight(xN_Left_UpStream_Up, xN_Right_UpStream_Up, 
		    xN_Asym_UpStream_Up, e_xN_Asym_UpStream_Up, nBins);
    BinnedLeftRight(xN_Left_UpStream_Down, xN_Right_UpStream_Down,
		    xN_Asym_UpStream_Down, e_xN_Asym_UpStream_Down, nBins);
    BinnedLeftRight(xN_Left_DownStream_Up, xN_Right_DownStream_Up, 
		    xN_Asym_DownStream_Up, e_xN_Asym_DownStream_Up, nBins);
    BinnedLeftRight(xN_Left_DownStream_Down, xN_Right_DownStream_Down,
		    xN_Asym_DownStream_Down, e_xN_Asym_DownStream_Down, nBins);
  
    BinnedLeftRight(xPi_Left_UpStream, xPi_Right_UpStream, xPi_Left_DownStream,
		    xPi_Right_DownStream, xPi_Asym_UpStream, xPi_Asym_DownStream,
		    xPi_Asym, e_xPi_Asym_UpStream, e_xPi_Asym_DownStream,
		    e_xPi_Asym,
		    nBins);
    BinnedLeftRight(xPi_Left_UpStream_Up, xPi_Right_UpStream_Up, 
		    xPi_Asym_UpStream_Up, e_xPi_Asym_UpStream_Up, nBins);
    BinnedLeftRight(xPi_Left_UpStream_Down, xPi_Right_UpStream_Down,
		    xPi_Asym_UpStream_Down, e_xPi_Asym_UpStream_Down, nBins);
    BinnedLeftRight(xPi_Left_DownStream_Up, xPi_Right_DownStream_Up, 
		    xPi_Asym_DownStream_Up, e_xPi_Asym_DownStream_Up, nBins);
    BinnedLeftRight(xPi_Left_DownStream_Down, xPi_Right_DownStream_Down,
		    xPi_Asym_DownStream_Down, e_xPi_Asym_DownStream_Down, nBins);
  
    BinnedLeftRight(xF_Left_UpStream, xF_Right_UpStream, xF_Left_DownStream,
		    xF_Right_DownStream, xF_Asym_UpStream, xF_Asym_DownStream,
		    xF_Asym, e_xF_Asym_UpStream, e_xF_Asym_DownStream, e_xF_Asym,
		    nBins);
    BinnedLeftRight(xF_Left_UpStream_Up, xF_Right_UpStream_Up, 
		    xF_Asym_UpStream_Up, e_xF_Asym_UpStream_Up, nBins);
    BinnedLeftRight(xF_Left_UpStream_Down, xF_Right_UpStream_Down,
		    xF_Asym_UpStream_Down, e_xF_Asym_UpStream_Down, nBins);
    BinnedLeftRight(xF_Left_DownStream_Up, xF_Right_DownStream_Up, 
		    xF_Asym_DownStream_Up, e_xF_Asym_DownStream_Up, nBins);
    BinnedLeftRight(xF_Left_DownStream_Down, xF_Right_DownStream_Down,
		    xF_Asym_DownStream_Down, e_xF_Asym_DownStream_Down, nBins);

    BinnedLeftRight(pT_Left_UpStream, pT_Right_UpStream, pT_Left_DownStream,
		    pT_Right_DownStream, pT_Asym_UpStream, pT_Asym_DownStream,
		    pT_Asym, e_pT_Asym_UpStream, e_pT_Asym_DownStream, e_pT_Asym,
		    nBins);
    BinnedLeftRight(pT_Left_UpStream_Up, pT_Right_UpStream_Up, 
		    pT_Asym_UpStream_Up, e_pT_Asym_UpStream_Up, nBins);
    BinnedLeftRight(pT_Left_UpStream_Down, pT_Right_UpStream_Down,
		    pT_Asym_UpStream_Down, e_pT_Asym_UpStream_Down, nBins);
    BinnedLeftRight(pT_Left_DownStream_Up, pT_Right_DownStream_Up, 
		    pT_Asym_DownStream_Up, e_pT_Asym_DownStream_Up, nBins);
    BinnedLeftRight(pT_Left_DownStream_Down, pT_Right_DownStream_Down,
		    pT_Asym_DownStream_Down, e_pT_Asym_DownStream_Down, nBins);

    BinnedLeftRight(M_Left_UpStream, M_Right_UpStream, M_Left_DownStream,
		    M_Right_DownStream, M_Asym_UpStream, M_Asym_DownStream,
		    M_Asym, e_M_Asym_UpStream, e_M_Asym_DownStream, e_M_Asym,
		    nBins);
    BinnedLeftRight(M_Left_UpStream_Up, M_Right_UpStream_Up, 
		    M_Asym_UpStream_Up, e_M_Asym_UpStream_Up, nBins);
    BinnedLeftRight(M_Left_UpStream_Down, M_Right_UpStream_Down,
		    M_Asym_UpStream_Down, e_M_Asym_UpStream_Down, nBins);
    BinnedLeftRight(M_Left_DownStream_Up, M_Right_DownStream_Up, 
		    M_Asym_DownStream_Up, e_M_Asym_DownStream_Up, nBins);
    BinnedLeftRight(M_Left_DownStream_Down, M_Right_DownStream_Down,
		    M_Asym_DownStream_Down, e_M_Asym_DownStream_Down, nBins);
  }

  //Generic binning
  Double_t Asym_UpStream[NVar], Asym_DownStream[NVar], Asym[NVar];
  Double_t e_Asym_UpStream[NVar], e_Asym_DownStream[NVar];
  Double_t e_Asym[NVar];
  
  Double_t Asym_UpStream_Up[NVar], Asym_DownStream_Up[NVar]; //By target
  Double_t Asym_UpStream_Down[NVar], Asym_DownStream_Down[NVar];
  Double_t e_Asym_UpStream_Up[NVar], e_Asym_DownStream_Up[NVar];
  Double_t e_Asym_UpStream_Down[NVar], e_Asym_DownStream_Down[NVar];

  if(Bflag){
    BinnedLeftRight(Left_UpStream, Right_UpStream, Left_DownStream,
		    Right_DownStream, Asym_UpStream, Asym_DownStream,
		    Asym, e_Asym_UpStream, e_Asym_DownStream, e_Asym,
		    NVar);
    BinnedLeftRight(Left_UpStream_Up, Right_UpStream_Up, 
		    Asym_UpStream_Up, e_Asym_UpStream_Up, NVar);
    BinnedLeftRight(Left_UpStream_Down, Right_UpStream_Down,
		    Asym_UpStream_Down, e_Asym_UpStream_Down, NVar);
    BinnedLeftRight(Left_DownStream_Up, Right_DownStream_Up, 
		    Asym_DownStream_Up, e_Asym_DownStream_Up, NVar);
    BinnedLeftRight(Left_DownStream_Down, Right_DownStream_Down,
		    Asym_DownStream_Down, e_Asym_DownStream_Down, NVar);
  }
  // }}}

  //TGraphs
  // {{{
  Double_t xval_xN[nBins], xval_xPi[nBins], xval_xF[nBins];
  Double_t xval_pT[nBins], xval_M[nBins];
  Double_t ex[nBins];
  for (Int_t i=0; i<nBins; i++) {
    xval_xN[i] = xN_xval.at(i); xval_xPi[i] = xPi_xval.at(i);
    xval_xF[i]=xF_xval.at(i); xval_pT[i]=pT_xval.at(i);xval_M[i]=M_xval.at(i);
    ex[i] = 0.0;
  }
  
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

  //By Target
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


  //Generic binning
  TGraphErrors* gr_LR_UpStream = new TGraphErrors(NVar, upS_xval,
						     Asym_UpStream, nVar_ex,
						     e_Asym_UpStream);
  TGraphErrors* gr_LR_DownStream = new TGraphErrors(NVar, downS_xval,
						       Asym_DownStream, nVar_ex,
						       e_Asym_DownStream);
  TGraphErrors* gr_LR = new TGraphErrors(NVar, upS_xval,
					    Asym, nVar_ex,
					    e_Asym);

  //By Target
  TGraphErrors* gr_LR_UpStream_Up = new TGraphErrors(NVar, upS_xval,
						     Asym_UpStream_Up, nVar_ex,
						     e_Asym_UpStream_Up);
  TGraphErrors* gr_LR_UpStream_Down = new TGraphErrors(NVar, upS_xval,
						       Asym_UpStream_Down,
						       nVar_ex,
						       e_Asym_UpStream_Down);
  TGraphErrors* gr_LR_DownStream_Up = new TGraphErrors(NVar, downS_xval,
						       Asym_DownStream_Up,
						       nVar_ex,
						       e_Asym_DownStream_Up);
  TGraphErrors* gr_LR_DownStream_Down = new TGraphErrors(NVar, downS_xval,
							 Asym_DownStream_Down,
							 nVar_ex,
							 e_Asym_DownStream_Down);
  // }}}

  ////////////////
  //Draw and pretty up graphs
  ////////////////
  //Binned physics values
  ///////////////
  // {{{
  if (binFlag){
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
  }//has binFlag
  // }}}
  
  //TGraphs by target
  ///////////////
  // {{{

  if (binFlag){
    Int_t ipad = 1;
    TCanvas *cxN = new TCanvas();
    cxN->Divide(2, 2);
    cxN->cd(ipad); ipad++;
    gr_xN_LR_UpStream_Up->Draw("AP");
    SetupTGraph(gr_xN_LR_UpStream_Up, "UpStream Pol Up", "xN", 1.0);
    Double_t xmin = gr_xN_LR_UpStream_Up->GetXaxis()->GetXmin();
    Double_t xmax = gr_xN_LR_UpStream_Up->GetXaxis()->GetXmax();
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
  }//has binFlag

  // }}}

  //Generica binning
  // {{{
  if (Bflag){
    TCanvas *cGen = new TCanvas();
    cGen->Divide(3);
    cGen->cd(1);
    Double_t xmin = gr_LR_UpStream->GetXaxis()->GetXmin();
    Double_t xmax = gr_LR_UpStream->GetXaxis()->GetXmax();
    TLine *l_UpStream=new TLine(xmin,0.0,xmax,0.0);
    gr_LR_UpStream->Draw("AP");
    l_UpStream->Draw("same");
  
    cGen->cd(2);
    xmin = gr_LR_DownStream->GetXaxis()->GetXmin();
    xmax = gr_LR_DownStream->GetXaxis()->GetXmax();
    TLine *l_DownStream=new TLine(xmin,0.0,xmax,0.0);
    gr_LR_DownStream->Draw("AP");
    l_DownStream->Draw("same");
    
    cGen->cd(3);
    xmin = gr_LR->GetXaxis()->GetXmin();
    xmax = gr_LR->GetXaxis()->GetXmax();
    TLine *l=new TLine(xmin,0.0,xmax,0.0);
    gr_LR->Draw("AP");
    l->Draw("same");

    SetupTGraph(gr_LR_UpStream, "UpStream", binVar, 1);//generic setups
    SetupTLine(l_UpStream);
    SetupTGraph(gr_LR_DownStream, "DownStream", binVar, 1);
    SetupTLine(l_DownStream);
    SetupTGraph(gr_LR, "Both Targets", binVar, 1);
    SetupTLine(l);

    //Generic by target
    Int_t ipad = 1;
    TCanvas *cTargGen = new TCanvas();
    cTargGen->Divide(2, 2);
    cTargGen->cd(ipad); ipad++;
    gr_LR_UpStream_Up->Draw("AP");
    SetupTGraph(gr_LR_UpStream_Up, "UpStream Pol Up", "Gen", 1.0);
    xmin = gr_LR_UpStream_Up->GetXaxis()->GetXmin();
    xmax = gr_LR_UpStream_Up->GetXaxis()->GetXmax();
    TLine *l_LR_UpStream_Up = new TLine(xmin,0.0,xmax,0.0);
    l_LR_UpStream_Up->Draw("same");
    SetupTLine(l_LR_UpStream_Up);
    cTargGen->cd(ipad); ipad++;
    gr_LR_UpStream_Down->Draw("AP");
    SetupTGraph(gr_LR_UpStream_Down, "UpStream Pol Down", "Gen");
    TLine *l_LR_UpStream_Down = new TLine(xmin,0.0,xmax,0.0);
    l_LR_UpStream_Down->Draw("same");
    SetupTLine(l_LR_UpStream_Down);
    
    cTargGen->cd(ipad); ipad++;
    xmin = gr_LR_DownStream_Up->GetXaxis()->GetXmin();
    xmax = gr_LR_DownStream_Up->GetXaxis()->GetXmax();
    gr_LR_DownStream_Up->Draw("AP");
    SetupTGraph(gr_LR_DownStream_Up, "DownStream Pol Up", "Gen");
    TLine *l_LR_DownStream_Up = new TLine(xmin,0.0,xmax,0.0);
    l_LR_DownStream_Up->Draw("same");
    SetupTLine(l_LR_DownStream_Up);
    cTargGen->cd(ipad); ipad++;
    gr_LR_DownStream_Down->Draw("AP");
    SetupTGraph(gr_LR_DownStream_Down, "DownStream Pol Down", "Gen");
    TLine *l_LR_DownStream_Down = new TLine(xmin,0.0,xmax,0.0);
    l_LR_DownStream_Down->Draw("same");
    SetupTLine(l_LR_DownStream_Down);
  }
  // }}}
  
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

    if (binFlag){
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
    }

    if (Bflag){
      gr_LR_UpStream->Write("gr_LR_UpStream");
      gr_LR_DownStream->Write("gr_LR_DownStream");
      gr_LR->Write("gr_LR");

      gr_LR_UpStream_Up->Write("gr_LR_UpStream_Up");//By target
      gr_LR_UpStream_Down->Write("gr_LR_UpStream_Down");
      gr_LR_DownStream_Up->Write("gr_LR_DownStream_Up");
      gr_LR_DownStream_Down->Write("gr_LR_DownStream_Down");
    }

    for (Int_t i=0; i<num_1Dh; i++) h1[i]->Write();

    for (Int_t i=0; i<num_2Dh; i++) h2[i]->Write();
  
    fout->Close();
  }
  // }}}
  
  theApp.Run();//Needed to make root graphics work on C++
}
