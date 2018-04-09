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
    cout << "----Writing Options----"<< endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << " " << endl;
    cout << "----Additional Setup Options----"<< endl;
    cout << "Option:  -S Left/right asymmetry choice" << endl;
    cout << "	 (True=no spin influence, Spin=spin influence, default=Spin)"
	 << endl;
    cout << "    (Default=Spin)" << endl;
    cout << "Option:  -P       (Turn off polarization and dilution corrections)"
	 << endl;
    cout << "Option:  -T trig       (Only specific trigger)" << endl;
    cout << "                  (\"LL\"=last-last, \"LO\"=last-outer, " <<
      "\"LL_LO\"=last-last && last-outer)" << endl;
    cout << "Option:  -i min       (Minimum mass to consider)" << endl;
    cout << "    (Default=0.0)" << endl;
    cout << "Option:  -a max       (Maximum mass to consider)" << endl;
    cout << "    (Default=12.0)" << endl;
    cout << " " << endl;
    cout << "----Changing Binning Options----" << endl;
    cout << "Option:  -b textfile with binning information	";
    cout << "(textfile should be made from Macro/Binning/avgBinBounds.C)"<<endl;
    cout << "Option:  -M (\"HM\", \"JPsi\", \"AMDY\") to specify which mass ";
    cout << "range to use for \"binning information\" " << endl;
    cout << "   (Must be used with \"-b\" option)" << endl;
    cout << "   (default when \"-b\" is used=AMDY)" << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t wflag=0, Qflag=0, fflag=0, Sflag=0, Pflag=0, Tflag=0;
  Int_t iflag=0, aflag=0, binFlag=0, Mflag=0;
  Int_t c;
  TString fname = "", outFile = "", leftrightChoice="", trig="";
  TString binFile = "", massRange="";
  Double_t M_min=0.0, M_max=12.0;
  
  while ((c = getopt (argc, argv, "Pwf:Q:S:T:i:a:b:M:")) != -1) {
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
      break;
    case 'S':
      Sflag = 1;
      leftrightChoice += optarg;
      break;
    case 'P':
      Pflag = 1;
      break;
    case 'T':
      Tflag = 1;
      trig += optarg;
      break;
    case 'i':
      iflag = 1;
      M_min = stof(optarg);
      break;
    case 'a':
      aflag = 1;
      M_max = stof(optarg);
      break;
    case 'b':
      binFlag = 1;
      binFile += optarg;
      break;
    case 'M':
      Mflag = 1;
      massRange += optarg;
      break;
    case '?':
      if (optopt == 'u')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'S')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'T')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'i')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'a')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'b')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'M')
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

  Int_t trigChoice=0;
  if (trig=="LL") trigChoice = 65792;
  else if (trig=="LO") trigChoice = 65540;
  else if (trig=="LL_LO") trigChoice = 65796;
  else if (Tflag) {
    cout << " " << endl;
    cout << "Option -T" << trig << " is not a valid choice" << endl;
    exit(EXIT_FAILURE);
  }

  vector<Double_t> xN_bounds, xN_xval; 
  vector<Double_t> xPi_bounds, xPi_xval;
  vector<Double_t> xF_bounds, xF_xval; 
  vector<Double_t> pT_bounds, pT_xval;
  vector<Double_t> M_bounds, M_xval;
  vector<Double_t> rad_bounds, rad_xval;
  vector<Double_t> vxZ_upstream_bounds, vxZ_upstream_xval;
  vector<Double_t> vxZ_downstream_bounds, vxZ_downstream_xval;
  if (binFlag) {
    xN_bounds.push_back(0.0);
    xPi_bounds.push_back(0.0);
    xF_bounds.push_back(-1.0);
    pT_bounds.push_back(0.4);
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
    pT_bounds.push_back(5.0);
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
  }//end binFlag
  else if (Mflag){
    cout << " " << endl;
    cout << "\"-M\" flag specifed but not \"-b\" flag" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
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
  Int_t nBounds;
  if(binFlag) nBounds = xN_bounds.size();
  else {
    nBounds = tv_xN_bounds.GetNoElements();
    for (Int_t i=0; i<nBounds; i++) {
      xN_bounds.push_back(tv_xN_bounds[i]);
      xPi_bounds.push_back(tv_xPi_bounds[i]);
      xF_bounds.push_back(tv_xF_bounds[i]);
      pT_bounds.push_back(tv_pT_bounds[i]);
      M_bounds.push_back(tv_M_bounds[i]);
    }
  }

  //Vertex specific
  Double_t vx_z, vx_x, vx_y;
  Int_t targetPosition;
  //Drell-Yan Angles
  Double_t PhiS_simple, Theta_CS, vOpenAngle;
  //Mu plus
  Double_t theta_traj1;
  //Mu minus
  Double_t theta_traj2;
  //Event
  Int_t trigMask, MasterTrigMask;
  //Target values
  Double_t Spin_0, Spin_1, Spin_2, Spin_3, Spin_4, Spin_5, Spin_6;
  Double_t dilutionFactor, Polarization;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;
  
  //Vertex specific
  tree->SetBranchAddress("vx_z", &vx_z);
  tree->SetBranchAddress("vx_x", &vx_x);
  tree->SetBranchAddress("vx_y", &vx_y);
  tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  tree->SetBranchAddress("Theta_CS", &Theta_CS);
  tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Mu plus
  tree->SetBranchAddress("theta_traj1", &theta_traj1);
  //Mu Minus
  tree->SetBranchAddress("theta_traj2", &theta_traj2);
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
  tree->SetBranchAddress("dilutionFactor", &dilutionFactor);
  tree->SetBranchAddress("Polarization", &Polarization);  
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

  //Averages
  Double_t AvgPolarization=0.0, AvgDilution=0.0, AvgDilution_corrected=0.0;
  Int_t AvgPolarization_count=0, AvgDilution_count=0;

  vector<Double_t> AvgPol_xN(nBins, 0.0), AvgDil_xN(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count(nBins, 0), AvgDil_xN_count(nBins, 0);
  vector<Double_t> AvgPol_xN_UpStream(nBins, 0.0), AvgDil_xN_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_UpStream(nBins, 0), AvgDil_xN_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_xN_DownStream(nBins, 0.0), AvgDil_xN_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_DownStream(nBins, 0), AvgDil_xN_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_xN_UpStream_Up(nBins, 0.0), AvgDil_xN_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_UpStream_Up(nBins, 0), AvgDil_xN_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xN_DownStream_Up(nBins, 0.0), AvgDil_xN_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_DownStream_Up(nBins, 0), AvgDil_xN_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xN_UpStream_Down(nBins, 0.0), AvgDil_xN_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_UpStream_Down(nBins, 0), AvgDil_xN_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_xN_DownStream_Down(nBins, 0.0), AvgDil_xN_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_DownStream_Down(nBins, 0), AvgDil_xN_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_xPi(nBins, 0.0), AvgDil_xPi(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count(nBins, 0), AvgDil_xPi_count(nBins, 0);
  vector<Double_t> AvgPol_xPi_UpStream(nBins, 0.0), AvgDil_xPi_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_UpStream(nBins, 0), AvgDil_xPi_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_xPi_DownStream(nBins, 0.0), AvgDil_xPi_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_DownStream(nBins, 0), AvgDil_xPi_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_xPi_UpStream_Up(nBins, 0.0), AvgDil_xPi_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_UpStream_Up(nBins, 0), AvgDil_xPi_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xPi_DownStream_Up(nBins, 0.0), AvgDil_xPi_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_DownStream_Up(nBins, 0), AvgDil_xPi_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xPi_UpStream_Down(nBins, 0.0), AvgDil_xPi_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_UpStream_Down(nBins, 0), AvgDil_xPi_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_xPi_DownStream_Down(nBins, 0.0), AvgDil_xPi_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_DownStream_Down(nBins, 0), AvgDil_xPi_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_xF(nBins, 0.0), AvgDil_xF(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count(nBins, 0), AvgDil_xF_count(nBins, 0);
  vector<Double_t> AvgPol_xF_UpStream(nBins, 0.0), AvgDil_xF_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_UpStream(nBins, 0), AvgDil_xF_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_xF_DownStream(nBins, 0.0), AvgDil_xF_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_DownStream(nBins, 0), AvgDil_xF_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_xF_UpStream_Up(nBins, 0.0), AvgDil_xF_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_UpStream_Up(nBins, 0), AvgDil_xF_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xF_DownStream_Up(nBins, 0.0), AvgDil_xF_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_DownStream_Up(nBins, 0), AvgDil_xF_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xF_UpStream_Down(nBins, 0.0), AvgDil_xF_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_UpStream_Down(nBins, 0), AvgDil_xF_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_xF_DownStream_Down(nBins, 0.0), AvgDil_xF_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_DownStream_Down(nBins, 0), AvgDil_xF_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_pT(nBins, 0.0), AvgDil_pT(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count(nBins, 0), AvgDil_pT_count(nBins, 0);
  vector<Double_t> AvgPol_pT_UpStream(nBins, 0.0), AvgDil_pT_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_UpStream(nBins, 0), AvgDil_pT_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_pT_DownStream(nBins, 0.0), AvgDil_pT_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_DownStream(nBins, 0), AvgDil_pT_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_pT_UpStream_Up(nBins, 0.0), AvgDil_pT_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_UpStream_Up(nBins, 0), AvgDil_pT_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_pT_DownStream_Up(nBins, 0.0), AvgDil_pT_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_DownStream_Up(nBins, 0), AvgDil_pT_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_pT_UpStream_Down(nBins, 0.0), AvgDil_pT_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_UpStream_Down(nBins, 0), AvgDil_pT_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_pT_DownStream_Down(nBins, 0.0), AvgDil_pT_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_DownStream_Down(nBins, 0), AvgDil_pT_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_M(nBins, 0.0), AvgDil_M(nBins, 0.0);
  vector<Int_t> AvgPol_M_count(nBins, 0), AvgDil_M_count(nBins, 0);
  vector<Double_t> AvgPol_M_UpStream(nBins, 0.0), AvgDil_M_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_UpStream(nBins, 0), AvgDil_M_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_M_DownStream(nBins, 0.0), AvgDil_M_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_DownStream(nBins, 0), AvgDil_M_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_M_UpStream_Up(nBins, 0.0), AvgDil_M_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_UpStream_Up(nBins, 0), AvgDil_M_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_M_DownStream_Up(nBins, 0.0), AvgDil_M_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_DownStream_Up(nBins, 0), AvgDil_M_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_M_UpStream_Down(nBins, 0.0), AvgDil_M_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_UpStream_Down(nBins, 0), AvgDil_M_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_M_DownStream_Down(nBins, 0.0), AvgDil_M_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_DownStream_Down(nBins, 0), AvgDil_M_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_rad(nBins, 0.0), AvgDil_rad(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count(nBins, 0), AvgDil_rad_count(nBins, 0);
  vector<Double_t> AvgPol_rad_UpStream(nBins, 0.0), AvgDil_rad_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_UpStream(nBins, 0), AvgDil_rad_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_rad_DownStream(nBins, 0.0), AvgDil_rad_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_DownStream(nBins, 0), AvgDil_rad_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_rad_UpStream_Up(nBins, 0.0), AvgDil_rad_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_UpStream_Up(nBins, 0), AvgDil_rad_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_rad_DownStream_Up(nBins, 0.0), AvgDil_rad_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_DownStream_Up(nBins, 0), AvgDil_rad_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_rad_UpStream_Down(nBins, 0.0), AvgDil_rad_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_UpStream_Down(nBins, 0), AvgDil_rad_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_rad_DownStream_Down(nBins, 0.0), AvgDil_rad_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_DownStream_Down(nBins, 0), AvgDil_rad_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_vxZ_upstream(nBins, 0.0), AvgDil_vxZ_upstream(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_UpStream(nBins, 0), AvgDil_vxZ_count_UpStream(nBins, 0);

  //Polarization by target
  vector<Double_t> AvgPol_vxZ_upstream_Up(nBins, 0.0), AvgDil_vxZ_upstream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_UpStream_Up(nBins, 0), AvgDil_vxZ_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_vxZ_upstream_Down(nBins, 0.0), AvgDil_vxZ_upstream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_UpStream_Down(nBins, 0), AvgDil_vxZ_count_UpStream_Down(nBins, 0);

  vector<Double_t> AvgPol_vxZ_downstream(nBins, 0.0), AvgDil_vxZ_downstream(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_DownStream(nBins, 0), AvgDil_vxZ_count_DownStream(nBins, 0);

  //Polarization by target
  vector<Double_t> AvgPol_vxZ_downstream_Up(nBins, 0.0), AvgDil_vxZ_downstream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_DownStream_Up(nBins, 0), AvgDil_vxZ_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_vxZ_downstream_Down(nBins, 0.0), AvgDil_vxZ_downstream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_DownStream_Down(nBins, 0), AvgDil_vxZ_count_DownStream_Down(nBins, 0);
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

  const Int_t nBasic=10;
  TH1D* hBasic_low[nBasic], *hBasic_high[nBasic];
  for (Int_t i=0; i<nBasic; i++) {
    hBasic_low[i] = new TH1D(Form("hlow_%i", i), Form("hlow_%i", i),100,2,5);
    hBasic_high[i] = new TH1D(Form("hhigh_%i",i),Form("hhigh_%i", i),100,2,5);
  }
  
  
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
      
      if (iflag || aflag){
	cout << "Additional mass cut" << endl;
	cout << "    Mass range " << M_min << " - " << M_max << endl;
      }
      
      if(leftrightChoice=="True") {
	cout << "True left/right asymmetry (no spin influence)" << endl;
      }
      else if (leftrightChoice=="Spin" || leftrightChoice==""){
	cout << "Spin influnced left/right asymmetry" << endl;
      }
      if (Tflag){
	cout << "Trigger mask set to: " << trig << " only" << endl;
      }

      cout << " " << endl;
      first = false;
    }
    
    //Additional Optional cuts
    if (Tflag && (trigMask != trigChoice)) continue;
    if (Mmumu < M_min || Mmumu > M_max) continue;

    if (theta_traj1 < 0.04 && theta_traj1 > 0.02){
      hBasic_low[0]->Fill(Mmumu);
    }
    else if (theta_traj1 > 0.04){
      hBasic_high[0]->Fill(Mmumu);
    }

    if (theta_traj2 < 0.04 && theta_traj2 > 0.02){
      hBasic_low[1]->Fill(Mmumu);
    }
    else if (theta_traj2 > 0.04){
      hBasic_high[1]->Fill(Mmumu);
    }

    
    ////All data after cuts
    //////////////
    ///////////////General useful quantities
    Double_t radius = TMath::Sqrt(vx_x*vx_x+vx_y*vx_y);

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
	BinDataCounts(xN_Left_UpStream, x_target, xN_bounds);
	BinDataCounts(xPi_Left_UpStream, x_beam, xPi_bounds);
	BinDataCounts(xF_Left_UpStream, x_feynman, xF_bounds);
	BinDataCounts(pT_Left_UpStream, q_transverse, pT_bounds);
	BinDataCounts(M_Left_UpStream, Mmumu, M_bounds);

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Left_UpStream_Up, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_UpStream_Up, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_UpStream_Up, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_UpStream_Up, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_UpStream_Up, Mmumu, M_bounds);
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Left_UpStream_Down, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_UpStream_Down, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_UpStream_Down, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_UpStream_Down, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_UpStream_Down, Mmumu, M_bounds);
	}
	// }}}
      }//End Left
      else if (Right){//Right
	// {{{
	BinDataCounts(xN_Right_UpStream, x_target, xN_bounds);
	BinDataCounts(xPi_Right_UpStream, x_beam, xPi_bounds);
	BinDataCounts(xF_Right_UpStream, x_feynman, xF_bounds);
	BinDataCounts(pT_Right_UpStream, q_transverse, pT_bounds);
	BinDataCounts(M_Right_UpStream, Mmumu, M_bounds);

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Right_UpStream_Up, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_UpStream_Up, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_UpStream_Up, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_UpStream_Up, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_UpStream_Up, Mmumu, M_bounds);
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Right_UpStream_Down, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_UpStream_Down, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_UpStream_Down, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_UpStream_Down, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_UpStream_Down, Mmumu, M_bounds);
	}
	// }}}
      }//End Right


      if (binFlag){
	//Dilution
	// {{{
      Double_t dil = TMath::Abs(dilutionFactor);
      Double_t correct_dil = 0.95*dil;
      AvgDilution += dil;
      AvgDilution_corrected += correct_dil;
      AvgDilution_count++;

      BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
      BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
      BinAvg(AvgDil_M, AvgDil_M_count, Mmumu, M_bounds, correct_dil);
      BinAvg(AvgDil_rad, AvgDil_rad_count, radius, rad_bounds, correct_dil);
      BinAvg(AvgDil_vxZ_upstream, AvgDil_vxZ_count_UpStream, vx_z,
	     vxZ_upstream_bounds, correct_dil);
      
      BinAvg(AvgDil_xN_UpStream, AvgDil_xN_count_UpStream, x_target,
	     xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi_UpStream, AvgDil_xPi_count_UpStream, x_beam,
	     xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF_UpStream, AvgDil_xF_count_UpStream, x_feynman,
	     xF_bounds, correct_dil);
      BinAvg(AvgDil_pT_UpStream, AvgDil_pT_count_UpStream, q_transverse,
	     pT_bounds, correct_dil);
      BinAvg(AvgDil_M_UpStream, AvgDil_M_count_UpStream, Mmumu,
	     M_bounds, correct_dil);
      BinAvg(AvgDil_rad_UpStream, AvgDil_rad_count_UpStream, radius, rad_bounds,
	     correct_dil);
      // }}}

	//Polarization
	// {{{
      Double_t pol = TMath::Abs(Polarization);
      BinAvg(AvgPol_xN_UpStream, AvgPol_xN_count_UpStream, x_target,
	     xN_bounds, pol);
      BinAvg(AvgPol_xPi_UpStream, AvgPol_xPi_count_UpStream, x_beam,
	     xPi_bounds, pol);
      BinAvg(AvgPol_xF_UpStream, AvgPol_xF_count_UpStream, x_feynman,
	     xF_bounds, pol);
      BinAvg(AvgPol_pT_UpStream, AvgPol_pT_count_UpStream, q_transverse,
	     pT_bounds, pol);
      BinAvg(AvgPol_M_UpStream, AvgPol_M_count_UpStream, Mmumu,
	     M_bounds, pol);
      BinAvg(AvgPol_rad_UpStream, AvgPol_rad_count_UpStream, radius,
	     rad_bounds, pol);
      BinAvg(AvgPol_vxZ_upstream,AvgPol_vxZ_count_UpStream, vx_z,
	     vxZ_upstream_bounds, pol);
      // }}}

	if (Spin_0 > 0) {//Polarized Up
	  //Dilution/Polarization
	  // {{{
	BinAvg(AvgDil_xN_UpStream_Up, AvgDil_xN_count_UpStream_Up, x_target,
	       xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_UpStream_Up, AvgDil_xPi_count_UpStream_Up, x_beam,
	       xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_UpStream_Up, AvgDil_xF_count_UpStream_Up, x_feynman,
	       xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_UpStream_Up,AvgDil_pT_count_UpStream_Up,q_transverse,
	       pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_UpStream_Up, AvgDil_M_count_UpStream_Up, Mmumu,
	       M_bounds, correct_dil);
	BinAvg(AvgDil_rad_UpStream_Up, AvgDil_rad_count_UpStream_Up, radius,
	       rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_upstream_Up, AvgDil_vxZ_count_UpStream_Up, vx_z,
	       vxZ_upstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_UpStream_Up, AvgPol_xN_count_UpStream_Up, x_target,
	       xN_bounds, pol);
	BinAvg(AvgPol_xPi_UpStream_Up, AvgPol_xPi_count_UpStream_Up, x_beam,
	       xPi_bounds, pol);
	BinAvg(AvgPol_xF_UpStream_Up, AvgPol_xF_count_UpStream_Up, x_feynman,
	       xF_bounds, pol);
	BinAvg(AvgPol_pT_UpStream_Up,AvgPol_pT_count_UpStream_Up,q_transverse,
	       pT_bounds, pol);  
	BinAvg(AvgPol_M_UpStream_Up, AvgPol_M_count_UpStream_Up, Mmumu,
	       M_bounds, pol);
	BinAvg(AvgPol_rad_UpStream_Up, AvgPol_rad_count_UpStream_Up, radius,
	       rad_bounds, pol);
	BinAvg(AvgPol_vxZ_upstream_Up, AvgPol_vxZ_count_UpStream_Up, vx_z,
	       vxZ_upstream_bounds, pol);
	// }}}
	}
	else if (Spin_0 < 0){//Polarized Down
	  //Dilution/Polarization
	  // {{{
	BinAvg(AvgDil_xN_UpStream_Down, AvgDil_xN_count_UpStream_Down, 
	       x_target, xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_UpStream_Down, AvgDil_xPi_count_UpStream_Down, 
	       x_beam, xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_UpStream_Down, AvgDil_xF_count_UpStream_Down, 
	       x_feynman, xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_UpStream_Down,AvgDil_pT_count_UpStream_Down,
	       q_transverse, pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_UpStream_Down, AvgDil_M_count_UpStream_Down,
	       Mmumu, M_bounds, correct_dil);
	BinAvg(AvgDil_rad_UpStream_Down, AvgDil_rad_count_UpStream_Down, radius,
	       rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_upstream_Down, AvgDil_vxZ_count_UpStream_Down, vx_z,
	       vxZ_upstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_UpStream_Down, AvgPol_xN_count_UpStream_Down, 
	       x_target, xN_bounds, pol);
	BinAvg(AvgPol_xPi_UpStream_Down, AvgPol_xPi_count_UpStream_Down, 
	       x_beam, xPi_bounds, pol);
	BinAvg(AvgPol_xF_UpStream_Down, AvgPol_xF_count_UpStream_Down, 
	       x_feynman, xF_bounds, pol);
	BinAvg(AvgPol_pT_UpStream_Down,AvgPol_pT_count_UpStream_Down,
	       q_transverse, pT_bounds, pol);  
	BinAvg(AvgPol_M_UpStream_Down, AvgPol_M_count_UpStream_Down,
	       Mmumu, M_bounds, pol);
	BinAvg(AvgPol_rad_UpStream_Down, AvgPol_rad_count_UpStream_Down, radius,
	       rad_bounds, pol);
	BinAvg(AvgPol_vxZ_upstream_Down, AvgPol_vxZ_count_UpStream_Down, vx_z,
	       vxZ_upstream_bounds, pol);
	// }}}
	}
      }

      h_phi_UpStream->Fill(phi_photon_lab);
      h_phi->Fill(phi_photon_lab);
    }//End UpStream target
    else if (targetPosition == 1) {//DownStream target
      
      if (Left){//Left
	// {{{
	BinDataCounts(xN_Left_DownStream, x_target, xN_bounds);
	BinDataCounts(xPi_Left_DownStream, x_beam, xPi_bounds);
	BinDataCounts(xF_Left_DownStream, x_feynman, xF_bounds);
	BinDataCounts(pT_Left_DownStream, q_transverse, pT_bounds);
	BinDataCounts(M_Left_DownStream, Mmumu, M_bounds);

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Left_DownStream_Up, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_DownStream_Up, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_DownStream_Up, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_DownStream_Up, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_DownStream_Up, Mmumu, M_bounds);
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Left_DownStream_Down, x_target, xN_bounds);
	  BinDataCounts(xPi_Left_DownStream_Down, x_beam, xPi_bounds);
	  BinDataCounts(xF_Left_DownStream_Down, x_feynman, xF_bounds);
	  BinDataCounts(pT_Left_DownStream_Down, q_transverse, pT_bounds);
	  BinDataCounts(M_Left_DownStream_Down, Mmumu, M_bounds);
	}
	// }}}
      }//End Left
      else if (Right){//Right
	// {{{
	BinDataCounts(xN_Right_DownStream, x_target, xN_bounds);
	BinDataCounts(xPi_Right_DownStream, x_beam, xPi_bounds);
	BinDataCounts(xF_Right_DownStream, x_feynman, xF_bounds);
	BinDataCounts(pT_Right_DownStream, q_transverse, pT_bounds);
	BinDataCounts(M_Right_DownStream, Mmumu, M_bounds);

	if (Spin_0 > 0){//Polarized Up
	  BinDataCounts(xN_Right_DownStream_Up, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_DownStream_Up, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_DownStream_Up, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_DownStream_Up, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_DownStream_Up, Mmumu, M_bounds);
	}
	else if (Spin_0 < 0){//Polarized Down
	  BinDataCounts(xN_Right_DownStream_Down, x_target, xN_bounds);
	  BinDataCounts(xPi_Right_DownStream_Down, x_beam, xPi_bounds);
	  BinDataCounts(xF_Right_DownStream_Down, x_feynman, xF_bounds);
	  BinDataCounts(pT_Right_DownStream_Down, q_transverse, pT_bounds);
	  BinDataCounts(M_Right_DownStream_Down, Mmumu, M_bounds);
	}
	// }}}
      }//End Right


      if (binFlag){
	//Dilution
	// {{{
      Double_t dil = TMath::Abs(dilutionFactor);
      Double_t correct_dil = 0.91*dil;
      AvgDilution += dil;
      AvgDilution_corrected += correct_dil;
      AvgDilution_count++;

      BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
      BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
      BinAvg(AvgDil_M, AvgDil_M_count, Mmumu, M_bounds, correct_dil);
      BinAvg(AvgDil_rad, AvgDil_rad_count, radius, rad_bounds, correct_dil);
      BinAvg(AvgDil_vxZ_downstream, AvgDil_vxZ_count_DownStream, vx_z,
	     vxZ_downstream_bounds, correct_dil);

      BinAvg(AvgDil_xN_DownStream, AvgDil_xN_count_DownStream, x_target,
	     xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi_DownStream, AvgDil_xPi_count_DownStream, x_beam,
	     xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF_DownStream, AvgDil_xF_count_DownStream, x_feynman,
	     xF_bounds, correct_dil);
      BinAvg(AvgDil_pT_DownStream, AvgDil_pT_count_DownStream, q_transverse,
	     pT_bounds, correct_dil);
      BinAvg(AvgDil_M_DownStream, AvgDil_M_count_DownStream, Mmumu,
	     M_bounds, correct_dil);
      BinAvg(AvgDil_rad_DownStream, AvgDil_rad_count_DownStream, radius,
	     rad_bounds, correct_dil);
      BinAvg(AvgDil_vxZ_downstream,AvgDil_vxZ_count_DownStream, vx_z,
	     vxZ_downstream_bounds, correct_dil);
      // }}}

	//Polarization
	// {{{
      Double_t pol = TMath::Abs(Polarization);
      BinAvg(AvgPol_xN_DownStream, AvgPol_xN_count_DownStream, x_target,
	     xN_bounds, pol);
      BinAvg(AvgPol_xPi_DownStream, AvgPol_xPi_count_DownStream, x_beam,
	     xPi_bounds, pol);
      BinAvg(AvgPol_xF_DownStream, AvgPol_xF_count_DownStream, x_feynman,
	     xF_bounds, pol);
      BinAvg(AvgPol_pT_DownStream, AvgPol_pT_count_DownStream, q_transverse,
	     pT_bounds, pol);
      BinAvg(AvgPol_M_DownStream, AvgPol_M_count_DownStream, Mmumu,
	     M_bounds, pol);
      BinAvg(AvgPol_rad_DownStream, AvgPol_rad_count_DownStream, radius,
	     rad_bounds, pol);
      BinAvg(AvgPol_vxZ_downstream,AvgPol_vxZ_count_DownStream, vx_z,
	     vxZ_downstream_bounds, pol);
      // }}}

	if (Spin_0 > 0) {//Polarized Up
	  //Dilution/Polarization
	  // {{{
	BinAvg(AvgDil_xN_DownStream_Up, AvgDil_xN_count_DownStream_Up, 
	       x_target, xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_DownStream_Up, AvgDil_xPi_count_DownStream_Up, 
	       x_beam, xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_DownStream_Up, AvgDil_xF_count_DownStream_Up, 
	       x_feynman, xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_DownStream_Up,AvgDil_pT_count_DownStream_Up,
	       q_transverse, pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_DownStream_Up, AvgDil_M_count_DownStream_Up,
	       Mmumu, M_bounds, correct_dil);
	BinAvg(AvgDil_rad_DownStream_Up, AvgDil_rad_count_DownStream_Up, radius,
	       rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_downstream_Up, AvgDil_vxZ_count_DownStream_Up, vx_z,
	       vxZ_downstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_DownStream_Up, AvgPol_xN_count_DownStream_Up, 
	       x_target, xN_bounds, pol);
	BinAvg(AvgPol_xPi_DownStream_Up, AvgPol_xPi_count_DownStream_Up, 
	       x_beam, xPi_bounds, pol);
	BinAvg(AvgPol_xF_DownStream_Up, AvgPol_xF_count_DownStream_Up, 
	       x_feynman, xF_bounds, pol);
	BinAvg(AvgPol_pT_DownStream_Up,AvgPol_pT_count_DownStream_Up,
	       q_transverse, pT_bounds, pol);  
	BinAvg(AvgPol_M_DownStream_Up, AvgPol_M_count_DownStream_Up,
	       Mmumu, M_bounds, pol);
	BinAvg(AvgPol_rad_DownStream_Up, AvgPol_rad_count_DownStream_Up, radius,
	       rad_bounds, pol);
	BinAvg(AvgPol_vxZ_downstream_Up, AvgPol_vxZ_count_DownStream_Up, vx_z,
	       vxZ_downstream_bounds, pol);
	// }}}
	}
	else if (Spin_0 < 0){//Polarized Down
	  //Dilution/Polarization
	  // {{{
	BinAvg(AvgDil_xN_DownStream_Down, AvgDil_xN_count_DownStream_Down, 
	       x_target, xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_DownStream_Down, AvgDil_xPi_count_DownStream_Down, 
	       x_beam, xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_DownStream_Down, AvgDil_xF_count_DownStream_Down, 
	       x_feynman, xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_DownStream_Down,AvgDil_pT_count_DownStream_Down,
	       q_transverse, pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_DownStream_Down, AvgDil_M_count_DownStream_Down,
	       Mmumu, M_bounds, correct_dil);
	BinAvg(AvgDil_rad_DownStream_Down, AvgDil_rad_count_DownStream_Down,
	       radius, rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_downstream_Down, AvgDil_vxZ_count_DownStream_Down,
	       vx_z, vxZ_downstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_DownStream_Down, AvgPol_xN_count_DownStream_Down, 
	       x_target, xN_bounds, pol);
	BinAvg(AvgPol_xPi_DownStream_Down, AvgPol_xPi_count_DownStream_Down, 
	       x_beam, xPi_bounds, pol);
	BinAvg(AvgPol_xF_DownStream_Down, AvgPol_xF_count_DownStream_Down, 
	       x_feynman, xF_bounds, pol);
	BinAvg(AvgPol_pT_DownStream_Down,AvgPol_pT_count_DownStream_Down,
	       q_transverse, pT_bounds, pol);  
	BinAvg(AvgPol_M_DownStream_Down, AvgPol_M_count_DownStream_Down,
	       Mmumu, M_bounds, pol);
	BinAvg(AvgPol_rad_DownStream_Down, AvgPol_rad_count_DownStream_Down,
	       radius, rad_bounds, pol);
	BinAvg(AvgPol_vxZ_downstream_Down, AvgPol_vxZ_count_DownStream_Down,
	       vx_z, vxZ_downstream_bounds, pol);
	// }}}
	}
      }

      h_phi_DownStream->Fill(phi_photon_lab);
      h_phi->Fill(phi_photon_lab);
    }//End DownStream target

    if (binFlag){
      //Average Polarization
      // {{{
      AvgPolarization += TMath::Abs(Polarization);
      AvgPolarization_count++;
      if(Polarization == 0.0 || dilutionFactor == 0.0){
	cout << "Problems with Polarization or dilution value" << endl;
	cout << "Polarization = " << Polarization << endl;
	cout << "Dilution = " << dilutionFactor << endl;
      }

      Double_t pol = TMath::Abs(Polarization);
      BinAvg(AvgPol_xN, AvgPol_xN_count, x_target, xN_bounds, pol);
      BinAvg(AvgPol_xPi, AvgPol_xPi_count, x_beam, xPi_bounds, pol);
      BinAvg(AvgPol_xF, AvgPol_xF_count, x_feynman, xF_bounds, pol);
      BinAvg(AvgPol_pT, AvgPol_pT_count, q_transverse, pT_bounds, pol);
      BinAvg(AvgPol_M, AvgPol_M_count, Mmumu, M_bounds, pol);
      BinAvg(AvgPol_rad, AvgPol_rad_count, radius, rad_bounds, pol);
      // }}}
    }

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
  //Dilution and polarization corrections
  TVectorD Dil_int(1);
  TVectorD Pol_int(1);

  ////Binned Amplitudes
  //Dilution and Polarization corrections
  TVectorD Dil_xN(nBins);
  TVectorD Dil_xPi(nBins);
  TVectorD Dil_xF(nBins);
  TVectorD Dil_pT(nBins);
  TVectorD Dil_M(nBins);

  TVectorD Dil_xN_UpStream(nBins);
  TVectorD Dil_xPi_UpStream(nBins);
  TVectorD Dil_xF_UpStream(nBins);
  TVectorD Dil_pT_UpStream(nBins);
  TVectorD Dil_M_UpStream(nBins);

  TVectorD Dil_xN_UpStream_Up(nBins);
  TVectorD Dil_xPi_UpStream_Up(nBins);
  TVectorD Dil_xF_UpStream_Up(nBins);
  TVectorD Dil_pT_UpStream_Up(nBins);
  TVectorD Dil_M_UpStream_Up(nBins);

  TVectorD Dil_xN_UpStream_Down(nBins);
  TVectorD Dil_xPi_UpStream_Down(nBins);
  TVectorD Dil_xF_UpStream_Down(nBins);
  TVectorD Dil_pT_UpStream_Down(nBins);
  TVectorD Dil_M_UpStream_Down(nBins);

  TVectorD Dil_xN_DownStream(nBins);
  TVectorD Dil_xPi_DownStream(nBins);
  TVectorD Dil_xF_DownStream(nBins);
  TVectorD Dil_pT_DownStream(nBins);
  TVectorD Dil_M_DownStream(nBins);

  TVectorD Dil_xN_DownStream_Up(nBins);
  TVectorD Dil_xPi_DownStream_Up(nBins);
  TVectorD Dil_xF_DownStream_Up(nBins);
  TVectorD Dil_pT_DownStream_Up(nBins);
  TVectorD Dil_M_DownStream_Up(nBins);

  TVectorD Dil_xN_DownStream_Down(nBins);
  TVectorD Dil_xPi_DownStream_Down(nBins);
  TVectorD Dil_xF_DownStream_Down(nBins);
  TVectorD Dil_pT_DownStream_Down(nBins);
  TVectorD Dil_M_DownStream_Down(nBins);
  
  TVectorD Pol_xN(nBins);
  TVectorD Pol_xPi(nBins);
  TVectorD Pol_xF(nBins);
  TVectorD Pol_pT(nBins);
  TVectorD Pol_M(nBins);

  TVectorD Pol_xN_UpStream(nBins);
  TVectorD Pol_xPi_UpStream(nBins);
  TVectorD Pol_xF_UpStream(nBins);
  TVectorD Pol_pT_UpStream(nBins);
  TVectorD Pol_M_UpStream(nBins);

  TVectorD Pol_xN_UpStream_Up(nBins);
  TVectorD Pol_xPi_UpStream_Up(nBins);
  TVectorD Pol_xF_UpStream_Up(nBins);
  TVectorD Pol_pT_UpStream_Up(nBins);
  TVectorD Pol_M_UpStream_Up(nBins);

  TVectorD Pol_xN_UpStream_Down(nBins);
  TVectorD Pol_xPi_UpStream_Down(nBins);
  TVectorD Pol_xF_UpStream_Down(nBins);
  TVectorD Pol_pT_UpStream_Down(nBins);
  TVectorD Pol_M_UpStream_Down(nBins);

  TVectorD Pol_xN_DownStream(nBins);
  TVectorD Pol_xPi_DownStream(nBins);
  TVectorD Pol_xF_DownStream(nBins);
  TVectorD Pol_pT_DownStream(nBins);
  TVectorD Pol_M_DownStream(nBins);

  TVectorD Pol_xN_DownStream_Up(nBins);
  TVectorD Pol_xPi_DownStream_Up(nBins);
  TVectorD Pol_xF_DownStream_Up(nBins);
  TVectorD Pol_pT_DownStream_Up(nBins);
  TVectorD Pol_M_DownStream_Up(nBins);

  TVectorD Pol_xN_DownStream_Down(nBins);
  TVectorD Pol_xPi_DownStream_Down(nBins);
  TVectorD Pol_xF_DownStream_Down(nBins);
  TVectorD Pol_pT_DownStream_Down(nBins);
  TVectorD Pol_M_DownStream_Down(nBins);

  if (!binFlag && !Pflag){
    Dil_int = *( (TVectorD*)fdata->Get("Dil_int") );
    Pol_int = *( (TVectorD*)fdata->Get("Pol_int") );

    ////Binned Amplitudes
    //Dilution and Polarization corrections
    Dil_xN = *( (TVectorD*)fdata->Get("Dil_xN") );
    Dil_xPi = *( (TVectorD*)fdata->Get("Dil_xPi") );
    Dil_xF = *( (TVectorD*)fdata->Get("Dil_xF") );
    Dil_pT = *( (TVectorD*)fdata->Get("Dil_pT") );
    Dil_M = *( (TVectorD*)fdata->Get("Dil_M") );

    Dil_xN_UpStream = *( (TVectorD*)fdata->Get("Dil_xN_UpStream") );
    Dil_xPi_UpStream = *( (TVectorD*)fdata->Get("Dil_xPi_UpStream") );
    Dil_xF_UpStream = *( (TVectorD*)fdata->Get("Dil_xF_UpStream") );
    Dil_pT_UpStream = *( (TVectorD*)fdata->Get("Dil_pT_UpStream") );
    Dil_M_UpStream = *( (TVectorD*)fdata->Get("Dil_M_UpStream") );

    Dil_xN_UpStream_Up = *( (TVectorD*)fdata->Get("Dil_xN_UpStream_Up"));
    Dil_xPi_UpStream_Up = *((TVectorD*)fdata->Get("Dil_xPi_UpStream_Up"));
    Dil_xF_UpStream_Up = *( (TVectorD*)fdata->Get("Dil_xF_UpStream_Up"));
    Dil_pT_UpStream_Up = *( (TVectorD*)fdata->Get("Dil_pT_UpStream_Up"));
    Dil_M_UpStream_Up = *( (TVectorD*)fdata->Get("Dil_M_UpStream_Up") );

    Dil_xN_UpStream_Down = *( (TVectorD*)fdata->Get("Dil_xN_UpStream_Down"));
    Dil_xPi_UpStream_Down=*((TVectorD*)fdata->Get("Dil_xPi_UpStream_Down"));
    Dil_xF_UpStream_Down = *( (TVectorD*)fdata->Get("Dil_xF_UpStream_Down"));
    Dil_pT_UpStream_Down = *( (TVectorD*)fdata->Get("Dil_pT_UpStream_Down"));
    Dil_M_UpStream_Down = *( (TVectorD*)fdata->Get("Dil_M_UpStream_Down") );

    Dil_xN_DownStream = *( (TVectorD*)fdata->Get("Dil_xN_DownStream") );
    Dil_xPi_DownStream = *( (TVectorD*)fdata->Get("Dil_xPi_DownStream"));
    Dil_xF_DownStream = *( (TVectorD*)fdata->Get("Dil_xF_DownStream") );
    Dil_pT_DownStream = *( (TVectorD*)fdata->Get("Dil_pT_DownStream") );
    Dil_M_DownStream = *( (TVectorD*)fdata->Get("Dil_M_DownStream") );

    Dil_xN_DownStream_Up = *( (TVectorD*)fdata->Get("Dil_xN_DownStream_Up"));
    Dil_xPi_DownStream_Up=*((TVectorD*)fdata->Get("Dil_xPi_DownStream_Up"));
    Dil_xF_DownStream_Up = *( (TVectorD*)fdata->Get("Dil_xF_DownStream_Up"));
    Dil_pT_DownStream_Up = *( (TVectorD*)fdata->Get("Dil_pT_DownStream_Up"));
    Dil_M_DownStream_Up = *( (TVectorD*)fdata->Get("Dil_M_DownStream_Up") );

    Dil_xN_DownStream_Down = *((TVectorD*)fdata->Get("Dil_xN_DownStream_Down"));
    Dil_xPi_DownStream_Down=*((TVectorD*)fdata->Get("Dil_xPi_DownStream_Down"));
    Dil_xF_DownStream_Down = *((TVectorD*)fdata->Get("Dil_xF_DownStream_Down"));
    Dil_pT_DownStream_Down = *((TVectorD*)fdata->Get("Dil_pT_DownStream_Down"));
    Dil_M_DownStream_Down = *( (TVectorD*)fdata->Get("Dil_M_DownStream_Down") );
  
    Pol_xN = *( (TVectorD*)fdata->Get("Pol_xN") );
    Pol_xPi = *( (TVectorD*)fdata->Get("Pol_xPi") );
    Pol_xF = *( (TVectorD*)fdata->Get("Pol_xF") );
    Pol_pT = *( (TVectorD*)fdata->Get("Pol_pT") );
    Pol_M = *( (TVectorD*)fdata->Get("Pol_M") );

    Pol_xN_UpStream = *( (TVectorD*)fdata->Get("Pol_xN_UpStream") );
    Pol_xPi_UpStream = *( (TVectorD*)fdata->Get("Pol_xPi_UpStream") );
    Pol_xF_UpStream = *( (TVectorD*)fdata->Get("Pol_xF_UpStream") );
    Pol_pT_UpStream = *( (TVectorD*)fdata->Get("Pol_pT_UpStream") );
    Pol_M_UpStream = *( (TVectorD*)fdata->Get("Pol_M_UpStream") );

    Pol_xN_UpStream_Up = *( (TVectorD*)fdata->Get("Pol_xN_UpStream_Up") );
    Pol_xPi_UpStream_Up = *( (TVectorD*)fdata->Get("Pol_xPi_UpStream_Up") );
    Pol_xF_UpStream_Up = *( (TVectorD*)fdata->Get("Pol_xF_UpStream_Up") );
    Pol_pT_UpStream_Up = *( (TVectorD*)fdata->Get("Pol_pT_UpStream_Up") );
    Pol_M_UpStream_Up = *( (TVectorD*)fdata->Get("Pol_M_UpStream_Up") );

    Pol_xN_UpStream_Down = *( (TVectorD*)fdata->Get("Pol_xN_UpStream_Down") );
    Pol_xPi_UpStream_Down = *( (TVectorD*)fdata->Get("Pol_xPi_UpStream_Down") );
    Pol_xF_UpStream_Down = *( (TVectorD*)fdata->Get("Pol_xF_UpStream_Down") );
    Pol_pT_UpStream_Down = *( (TVectorD*)fdata->Get("Pol_pT_UpStream_Down") );
    Pol_M_UpStream_Down = *( (TVectorD*)fdata->Get("Pol_M_UpStream_Down") );

    Pol_xN_DownStream = *( (TVectorD*)fdata->Get("Pol_xN_DownStream") );
    Pol_xPi_DownStream = *( (TVectorD*)fdata->Get("Pol_xPi_DownStream"));
    Pol_xF_DownStream = *( (TVectorD*)fdata->Get("Pol_xF_DownStream") );
    Pol_pT_DownStream = *( (TVectorD*)fdata->Get("Pol_pT_DownStream") );
    Pol_M_DownStream = *( (TVectorD*)fdata->Get("Pol_M_DownStream") );

    Pol_xN_DownStream_Up = *( (TVectorD*)fdata->Get("Pol_xN_DownStream_Up") );
    Pol_xPi_DownStream_Up = *( (TVectorD*)fdata->Get("Pol_xPi_DownStream_Up") );
    Pol_xF_DownStream_Up = *( (TVectorD*)fdata->Get("Pol_xF_DownStream_Up") );
    Pol_pT_DownStream_Up = *( (TVectorD*)fdata->Get("Pol_pT_DownStream_Up") );
    Pol_M_DownStream_Up = *( (TVectorD*)fdata->Get("Pol_M_DownStream_Up") );

    Pol_xN_DownStream_Down = *((TVectorD*)fdata->Get("Pol_xN_DownStream_Down"));
    Pol_xPi_DownStream_Down=*((TVectorD*)fdata->Get("Pol_xPi_DownStream_Down"));
    Pol_xF_DownStream_Down= *( (TVectorD*)fdata->Get("Pol_xF_DownStream_Down"));
    Pol_pT_DownStream_Down= *( (TVectorD*)fdata->Get("Pol_pT_DownStream_Down"));
    Pol_M_DownStream_Down = *( (TVectorD*)fdata->Get("Pol_M_DownStream_Down") );
  }//no binFlag
  else if (!Pflag) {
    for (Int_t i=0; i<nBins; i++) {
      //Dilution
      ///////////////
      Dil_xN[i] = AvgDil_xN[i]/AvgDil_xN_count[i];
      Dil_xPi[i] = AvgDil_xPi[i]/AvgDil_xPi_count[i];
      Dil_xF[i] = AvgDil_xF[i]/AvgDil_xF_count[i];
      Dil_pT[i] = AvgDil_pT[i]/AvgDil_pT_count[i];
      Dil_M[i] = AvgDil_M[i]/AvgDil_M_count[i];

      //UpStream
      Dil_xN_UpStream[i] = AvgDil_xN_UpStream[i]/AvgDil_xN_count_UpStream[i];
      Dil_xPi_UpStream[i] = AvgDil_xPi_UpStream[i]/AvgDil_xPi_count_UpStream[i];
      Dil_xF_UpStream[i] = AvgDil_xF_UpStream[i]/AvgDil_xF_count_UpStream[i];
      Dil_pT_UpStream[i] = AvgDil_pT_UpStream[i]/AvgDil_pT_count_UpStream[i];
      Dil_M_UpStream[i] = AvgDil_M_UpStream[i]/AvgDil_M_count_UpStream[i];

      //UpStream Polarized Up
      Dil_xN_UpStream_Up[i] = 
	AvgDil_xN_UpStream_Up[i]/AvgDil_xN_count_UpStream_Up[i];
      Dil_xPi_UpStream_Up[i] = 
	AvgDil_xPi_UpStream_Up[i]/AvgDil_xPi_count_UpStream_Up[i];
      Dil_xF_UpStream_Up[i] = 
	AvgDil_xF_UpStream_Up[i]/AvgDil_xF_count_UpStream_Up[i];
      Dil_pT_UpStream_Up[i] = 
	AvgDil_pT_UpStream_Up[i]/AvgDil_pT_count_UpStream_Up[i];
      Dil_M_UpStream_Up[i] = 
	AvgDil_M_UpStream_Up[i]/AvgDil_M_count_UpStream_Up[i];
      
      //UpStream Polarized Down
      Dil_xN_UpStream_Down[i] = 
	AvgDil_xN_UpStream_Down[i]/AvgDil_xN_count_UpStream_Down[i];
      Dil_xPi_UpStream_Down[i] = 
	AvgDil_xPi_UpStream_Down[i]/AvgDil_xPi_count_UpStream_Down[i];
      Dil_xF_UpStream_Down[i] = 
	AvgDil_xF_UpStream_Down[i]/AvgDil_xF_count_UpStream_Down[i];
      Dil_pT_UpStream_Down[i] =
	AvgDil_pT_UpStream_Down[i]/AvgDil_pT_count_UpStream_Down[i];
      Dil_M_UpStream_Down[i] =
	AvgDil_M_UpStream_Down[i]/AvgDil_M_count_UpStream_Down[i];
      
      //DownStream
      Dil_xN_DownStream[i] =AvgDil_xN_DownStream[i]/AvgDil_xN_count_DownStream[i];
      Dil_xPi_DownStream[i]=
	AvgDil_xPi_DownStream[i]/AvgDil_xPi_count_DownStream[i];
      Dil_xF_DownStream[i] =AvgDil_xF_DownStream[i]/AvgDil_xF_count_DownStream[i];
      Dil_pT_DownStream[i] =AvgDil_pT_DownStream[i]/AvgDil_pT_count_DownStream[i];
      Dil_M_DownStream[i] = AvgDil_M_DownStream[i]/AvgDil_M_count_DownStream[i];
      
      //DownStream Polarized Up
      Dil_xN_DownStream_Up[i] =
	AvgDil_xN_DownStream_Up[i]/AvgDil_xN_count_DownStream_Up[i];
      Dil_xPi_DownStream_Up[i]=
	AvgDil_xPi_DownStream_Up[i]/AvgDil_xPi_count_DownStream_Up[i];
      Dil_xF_DownStream_Up[i] =
	AvgDil_xF_DownStream_Up[i]/AvgDil_xF_count_DownStream_Up[i];
      Dil_pT_DownStream_Up[i] =
	AvgDil_pT_DownStream_Up[i]/AvgDil_pT_count_DownStream_Up[i];
      Dil_M_DownStream_Up[i] =
	AvgDil_M_DownStream_Up[i]/AvgDil_M_count_DownStream_Up[i];
      	
      //DownStream Polarized Down
      Dil_xN_DownStream_Down[i] =
	AvgDil_xN_DownStream_Down[i]/AvgDil_xN_count_DownStream_Down[i];
      Dil_xPi_DownStream_Down[i]=
	AvgDil_xPi_DownStream_Down[i]/AvgDil_xPi_count_DownStream_Down[i];
      Dil_xF_DownStream_Down[i] =
	AvgDil_xF_DownStream_Down[i]/AvgDil_xF_count_DownStream_Down[i];
      Dil_pT_DownStream_Down[i] =
	AvgDil_pT_DownStream_Down[i]/AvgDil_pT_count_DownStream_Down[i];
      Dil_M_DownStream_Down[i] =
	AvgDil_M_DownStream_Down[i]/AvgDil_M_count_DownStream_Down[i];
      
      //Polarization
      ////////////////
      Pol_xN[i] = AvgPol_xN[i]/AvgPol_xN_count[i];
      Pol_xPi[i] = AvgPol_xPi[i]/AvgPol_xPi_count[i];
      Pol_xF[i] = AvgPol_xF[i]/AvgPol_xF_count[i];
      Pol_pT[i] = AvgPol_pT[i]/AvgPol_pT_count[i];
      Pol_M[i] = AvgPol_M[i]/AvgPol_M_count[i];
      
      //UpStream////////////
      Pol_xN_UpStream[i] = AvgPol_xN_UpStream[i]/AvgPol_xN_count_UpStream[i];
      Pol_xPi_UpStream[i] = AvgPol_xPi_UpStream[i]/AvgPol_xPi_count_UpStream[i];
      Pol_xF_UpStream[i] = AvgPol_xF_UpStream[i]/AvgPol_xF_count_UpStream[i];
      Pol_pT_UpStream[i] = AvgPol_pT_UpStream[i]/AvgPol_pT_count_UpStream[i];
      Pol_M_UpStream[i] = AvgPol_M_UpStream[i]/AvgPol_M_count_UpStream[i];
      
      //UpStream Polarized Up
      Pol_xN_UpStream_Up[i] = 
	AvgPol_xN_UpStream_Up[i]/AvgPol_xN_count_UpStream_Up[i];
      Pol_xPi_UpStream_Up[i] = 
	AvgPol_xPi_UpStream_Up[i]/AvgPol_xPi_count_UpStream_Up[i];
      Pol_xF_UpStream_Up[i] = 
	AvgPol_xF_UpStream_Up[i]/AvgPol_xF_count_UpStream_Up[i];
      Pol_pT_UpStream_Up[i] = 
	AvgPol_pT_UpStream_Up[i]/AvgPol_pT_count_UpStream_Up[i];
      Pol_M_UpStream_Up[i] = 
	AvgPol_M_UpStream_Up[i]/AvgPol_M_count_UpStream_Up[i];

      //UpStream Polarized Down
      Pol_xN_UpStream_Down[i] = 
	AvgPol_xN_UpStream_Down[i]/AvgPol_xN_count_UpStream_Down[i];
      Pol_xPi_UpStream_Down[i] = 
	AvgPol_xPi_UpStream_Down[i]/AvgPol_xPi_count_UpStream_Down[i];
      Pol_xF_UpStream_Down[i] = 
	AvgPol_xF_UpStream_Down[i]/AvgPol_xF_count_UpStream_Down[i];
      Pol_pT_UpStream_Down[i] =
	AvgPol_pT_UpStream_Down[i]/AvgPol_pT_count_UpStream_Down[i];
      Pol_M_UpStream_Down[i] =
	AvgPol_M_UpStream_Down[i]/AvgPol_M_count_UpStream_Down[i];
	
      //DownStream////////////
      Pol_xN_DownStream[i] =AvgPol_xN_DownStream[i]/AvgPol_xN_count_DownStream[i];
      Pol_xPi_DownStream[i]=
	AvgPol_xPi_DownStream[i]/AvgPol_xPi_count_DownStream[i];
      Pol_xF_DownStream[i] =AvgPol_xF_DownStream[i]/AvgPol_xF_count_DownStream[i];
      Pol_pT_DownStream[i] =AvgPol_pT_DownStream[i]/AvgPol_pT_count_DownStream[i];
      Pol_M_DownStream[i] = AvgPol_M_DownStream[i]/AvgPol_M_count_DownStream[i];

      //DownStream Polarized Up
      Pol_xN_DownStream_Up[i] =
	AvgPol_xN_DownStream_Up[i]/AvgPol_xN_count_DownStream_Up[i];
      Pol_xPi_DownStream_Up[i]=
	AvgPol_xPi_DownStream_Up[i]/AvgPol_xPi_count_DownStream_Up[i];
      Pol_xF_DownStream_Up[i] =
	AvgPol_xF_DownStream_Up[i]/AvgPol_xF_count_DownStream_Up[i];
      Pol_pT_DownStream_Up[i] =
	AvgPol_pT_DownStream_Up[i]/AvgPol_pT_count_DownStream_Up[i];
      Pol_M_DownStream_Up[i] =
	AvgPol_M_DownStream_Up[i]/AvgPol_M_count_DownStream_Up[i];
	
      //DownStream Polarized Down
      Pol_xN_DownStream_Down[i] =
	AvgPol_xN_DownStream_Down[i]/AvgPol_xN_count_DownStream_Down[i];
      Pol_xPi_DownStream_Down[i]=
	AvgPol_xPi_DownStream_Down[i]/AvgPol_xPi_count_DownStream_Down[i];
      Pol_xF_DownStream_Down[i] =
	AvgPol_xF_DownStream_Down[i]/AvgPol_xF_count_DownStream_Down[i];
      Pol_pT_DownStream_Down[i] =
	AvgPol_pT_DownStream_Down[i]/AvgPol_pT_count_DownStream_Down[i];
      Pol_M_DownStream_Down[i] =
	AvgPol_M_DownStream_Down[i]/AvgPol_M_count_DownStream_Down[i];

    }//Dilution and polariztion setup loop

    Dil_int[0] = AvgDilution_corrected/AvgDilution_count;
    Pol_int[0] = AvgPolarization/AvgPolarization_count;
  }//has binFlag
    
  
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
    if (binFlag){
      xval_xN[i] = xN_xval.at(i); xval_xPi[i] = xPi_xval.at(i);
      xval_xF[i]=xF_xval.at(i); xval_pT[i]=pT_xval.at(i);xval_M[i]=M_xval.at(i);
      ex[i] = 0.0;
    }
    else {xval_xN[i] = tv_xN_xval[i]; xval_xPi[i] = tv_xPi_xval[i];
      xval_xF[i]=tv_xF_xval[i]; xval_pT[i]=tv_pT_xval[i];xval_M[i]=tv_M_xval[i];
      ex[i] = 0.0;
    }
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

  TCanvas* c1 = new TCanvas();
  c1->Divide(2);
  c1->cd(1);
  hBasic_low[1]->Draw();
  hBasic_low[1]->SetLineColor(kRed);
  hBasic_low[0]->Draw("sames");
  c1->cd(2);
  hBasic_high[0]->Draw();
  hBasic_high[1]->Draw("sames");
  hBasic_high[1]->SetLineColor(kRed);
  
  theApp.Run();//Needed to make root graphics work on C++
}
