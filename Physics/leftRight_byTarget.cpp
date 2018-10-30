#include "common.h"
#include "functions.h"
#include "leftright.h"
#include "twotargets.h"
#include "lrSpinCorr.h"
#include "genericBounds.h"
#include "lr_tgraph.h"
#include "binDist.h"

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){//Starting Help information
    cout << "\n\nUsage:" << endl;
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
    cout << "  --Real Data Only Options--"<< endl;
    cout << "Option:  -S Left/right asymmetry choice" << endl;
    cout << "	 (True=no spin influence, Spin=spin influence, default=Spin)"
	 << endl;
    cout << "    (Default=Spin) (don't use this with -C option)" << endl;
    cout << "Option:  -P       (Turn off polarization and dilution corrections)"
	 << endl;
    cout << " " << endl;
    cout << "  --Monte-Carlo Data Only Options--"<< endl;
    cout << "Option:  -C    (To specify reconstructed MC data)" << endl;
    cout << "Option:  -G    (To specify generated MC data)" << endl;
    cout <<"          (Assumed real data if these options are not given)"<<endl;
    cout << " " << endl;
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
    cout << "Option   -V binVar  " 
	 << "(Additional binning variable not in binFile)" << endl;
    cout << "   (Current options: openAngle, phi_pIn, theta_pIn, "
	 << "qP_pIn, phi_photon)" << endl;
    cout << "Option   -N bins" << endl;
    cout << "    (Option can only be used with \"-V\" option for binVar)"<<endl;
    cout << "Option:  -Z change number of bins for histogram distributions/n";
    cout << "    (default nHbins is 150)" << endl;
    cout << " " << endl;
    cout << "----Additional Options----" << endl;
    cout << "Option:  -D    (Debug mode, only 1000 events considered)" << endl;
    cout << "Option:  -H    (to see example inputs for this program)" << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }

  //Read input arguments
  Int_t wflag=0, Qflag=0, fflag=0, Sflag=0, Pflag=0, Tflag=0, Vflag=0, Nflag=0;
  Int_t iflag=0, aflag=0, binFlag=0, Dflag=0, Cflag=0, Hflag=0, Gflag=0,Zflag=0;
  Int_t c;
  TString fname = "", outFile = "", leftrightChoice="", trig="";
  TString binFile = "", binVar="";
  Double_t M_min=0.0, M_max=12.0;
  Int_t NVar, nHbins=150;
  
  while ((c = getopt (argc, argv, "Pwf:Q:S:T:V:N:i:a:b:M:DCGHZ:")) != -1) {
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
    case 'V':
      Vflag = 1;
      binVar += optarg;
      break;
    case 'N':
      Nflag = 1;
      NVar = atoi(optarg);
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
    case 'D':
      Dflag = 1;
      break;
    case 'C':
      Cflag = 1;
      break;
    case 'G':
      Gflag = 1;
      break;
    case 'H':
      Hflag = 1;
      break;
    case 'Z':
      Zflag = 1;
      nHbins = atoi(optarg);
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
      else if (optopt == 'V')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'i')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'a')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'b')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'Z')
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
  
  if (Hflag){
    cout << "MC example:" << endl;
    cout << "./leftRight_byTarget -C -QData/MC_Data/W12Mimic_5bins.root " <<
      "-b../DATA/MC_Data/TriggerMimic/W12Mimic_5bins.txt " <<
      "-f../DATA/MC_Data/TriggerMimic/W12Mimic_5bins.root" << endl;
    cout << " " << endl;
    cout << "Real Data example:" << endl;
    cout << "need to write..." << endl;
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
  
  if(Pflag){
    cout << " " << endl;
    cout << "No polarization or dilution factors performed" << endl;
    cout << " " << endl;
  }

  if(Dflag) {
    cout << " " << endl;
    cout << "Debug mode only 1000 events considered" << endl;
    cout << " " << endl;
  }

  if(Cflag) {
    cout << "" << endl;
    cout << "Using reconstructed Monte-Carlo data!" << endl;
    cout << " " << endl;

    if (Pflag || Sflag){
      cout << "" << endl;
      cout << "-P && -S options are not to be used with -C option" << endl;
      exit(EXIT_FAILURE);
    }

    if(Gflag){
      cout << "" << endl;
      cout << "Cannot be generated and reconstructed MC data" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else if (Gflag){
    cout << "" << endl;
    cout << "Using generated Monte-Carlo data!" << endl;
    cout << " " << endl;

    if (Pflag || Sflag){
      cout << "" << endl;
      cout << "-P && -S options are not to be used with -G option" << endl;
      exit(EXIT_FAILURE);
    }

    if(Cflag){
      cout << "" << endl;
      cout << "Cannot be generated and reconstructed MC data" << endl;
      exit(EXIT_FAILURE);
    }
  }
  else {
    cout << "" << endl;
    cout << "Using REAL data!" << endl;
    cout << " " << endl;
  }

  if (Vflag && !Nflag) {
    cout << "\"-B\" option used but no \"-N\" option" << endl;
    exit(EXIT_FAILURE);
  }
  else if (Nflag && !NVar){
    cout << " " << endl;
    cout << "Please enter and integer with -N option" << endl;
    exit(EXIT_FAILURE);
  }

  if(Zflag){ cout << "Number of histogram bins changed to  " << nHbins << endl;}
  
  //Opening data files/getting trees
  TFile *fdata = TFile::Open(fname);
  TChain *tree = new TChain("Particles");
  Int_t nTree=0;
  TList *li = (TList*)fdata->GetListOfKeys(); TIter iter( li->MakeIterator() );
  while (TObject *obj = iter() ){ TKey *theKey = (TKey*)obj;
    if (strncmp (theKey->GetClassName(),"TTree",4) == 0){
      tree->Add(fname+"/"+obj->GetName()+";"+Form("%i", theKey->GetCycle() ));
      nTree++;
    }
  } tree->LoadTree(0); cout << "Number of trees cycles:    " << nTree << endl; 
  
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
  
  Int_t Errors = 0;
  //Vertex specific
  Errors+=tree->SetBranchAddress("vx_z", &vx_z); 
  Errors+=tree->SetBranchAddress("vx_x", &vx_x);
  Errors+=tree->SetBranchAddress("vx_y", &vx_y);
  Errors+=tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  Errors+=tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  Errors+=tree->SetBranchAddress("Theta_CS", &Theta_CS);
  Errors+=tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Mu plus
  if (!Gflag) Errors+=tree->SetBranchAddress("theta_traj1", &theta_traj1);
  else Errors+=tree->SetBranchAddress("theta_muP", &theta_traj1);
  //Mu Minus
  if (!Gflag) Errors+=tree->SetBranchAddress("theta_traj2", &theta_traj2);
  else Errors+=tree->SetBranchAddress("theta_muM", &theta_traj2);
  //Event
  if (!Gflag) Errors+=tree->SetBranchAddress("trigMask", &trigMask);
  //Target values
  if (!Cflag && !Gflag) {
    Errors+=tree->SetBranchAddress("MasterTrigMask", &MasterTrigMask);
    Errors+=tree->SetBranchAddress("Spin_0", &Spin_0);
    Errors+=tree->SetBranchAddress("Spin_1", &Spin_1);
    Errors+=tree->SetBranchAddress("Spin_2", &Spin_2);
    Errors+=tree->SetBranchAddress("Spin_3", &Spin_3);
    Errors+=tree->SetBranchAddress("Spin_4", &Spin_4);
    Errors+=tree->SetBranchAddress("Spin_5", &Spin_5);
    Errors+=tree->SetBranchAddress("Spin_6", &Spin_6);
    Errors+=tree->SetBranchAddress("dilutionFactor", &dilutionFactor);
    Errors+=tree->SetBranchAddress("Polarization", &Polarization);
  }
  //DY-variables
  Errors+=tree->SetBranchAddress("x_beam", &x_beam);
  Errors+=tree->SetBranchAddress("x_target", &x_target);
  Errors+=tree->SetBranchAddress("x_feynman", &x_feynman);
  Errors+=tree->SetBranchAddress("q_transverse", &q_transverse);
  Errors+=tree->SetBranchAddress("Mmumu", &Mmumu);

  if(Errors) {
    cout << " " << endl;
    cout << "Errors in Setting tree branch addresses" << endl;
    cout << "Probably some variables are not defined in your file" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  Double_t radius;
  Double_t *boundValue=0;
  if (Vflag){
    cout << "Need asymetries written and taken care of..." << endl;
    exit(EXIT_FAILURE);
    if (binVar=="openAngle") boundValue = &vOpenAngle;
    else if (binVar=="radius") boundValue = &radius;
    else if (binVar=="vx_z") boundValue = &vx_z;
    else {
      cout << " " << endl;
      cout << "\"-B\" " << binVar << " option not supported" << endl;
      exit(EXIT_FAILURE);
    }
  }

  //Define variables
  leftright *xN, *xPi, *xF, *pT, *Mass;
  if (binFlag){
    xN = new lr_tgraph(binFile, "xN");
    xPi = new lr_tgraph(binFile, "xPi");
    xF = new lr_tgraph(binFile, "xF");
    pT = new lr_tgraph(binFile, "pT");
    Mass = new lr_tgraph(binFile, "mass", "M");
  }
  else{
    cout << "binFlag must be specifed for now because of new asymmetries"<<endl;
    exit(EXIT_FAILURE);
    
    xN = new lr_tgraph(fdata, "xN");
    xPi = new lr_tgraph(fdata, "xPi");
    xF = new lr_tgraph(fdata, "xF");
    pT = new lr_tgraph(fdata, "pT");
    Mass = new lr_tgraph(fdata, "mass", "M");

    if (!Cflag && !Gflag){
      xN->SetCorr(fdata, "xN");
      xPi->SetCorr(fdata, "xPi");
      xF->SetCorr(fdata, "xF");
      pT->SetCorr(fdata, "pT");
      Mass->SetCorr(fdata, "M");
    }
  }
  const Int_t nBasics = 5;
  leftright *Basics[nBasics] = {xN, xPi, xF, pT, Mass};
  Double_t *basicValues[nBasics] = {&x_target, &x_beam, &x_feynman,
				    &q_transverse, &Mmumu};
  genericBounds *genericPhys = NULL;
  if (Vflag) genericPhys = new lr_tgraph(NVar, binVar);

  binDist *MuMu_b_xN = new binDist("MuMu", "xN", binFile, nHbins, M_min, M_max);
  binDist *MuMu_b_xPi = new binDist("MuMu","xPi",binFile, nHbins, M_min, M_max);
  binDist *MuMu_b_xF = new binDist("MuMu", "xF", binFile, nHbins, M_min, M_max);
  binDist *MuMu_b_pT = new binDist("MuMu", "pT", binFile, nHbins, M_min, M_max);
  binDist *MuMu_b[nBasics-1] = {MuMu_b_xN, MuMu_b_xPi, MuMu_b_xF, MuMu_b_pT};

  //1st tree loop, equal out by target data
  Int_t nUpStream=0, nDownStream=0;
  Int_t tree_entries = (!Dflag) ? tree->GetEntries() : 1000;//Tree Loop
  if (Cflag || Gflag || Vflag) {
    for (Int_t ev=0; ev<tree_entries; ev++) {
      tree->GetEntry(ev, 0);

      //Additional Optional cuts
      if (Tflag && (trigMask != trigChoice)) continue;
      if (Mmumu < M_min || Mmumu > M_max) continue;
      if (Gflag && (targetPosition!=0 && targetPosition!=1) ) continue;
      
      //General useful quantities
      radius = TMath::Sqrt(vx_x*vx_x + vx_y*vx_y);

      if (targetPosition == 0) {
	nUpStream++;
	if (Vflag) genericPhys->SetUpBounds("upstream", *boundValue);
      }
      else if (targetPosition == 1) {
	nDownStream++;
	if (Vflag) genericPhys->SetUpBounds("downstream", *boundValue);
      }
      else {if (!Gflag) cout << "Not in NH3 targets" << endl;}
    }

    if (Cflag || Gflag){
    (nUpStream>nDownStream) ? nUpStream=nDownStream : nDownStream=nUpStream;
    cout << "Number of Entries    UpStream: " << nUpStream << "   DownStream: "
	 << nDownStream << endl;}

    if (Vflag) genericPhys->MakeBounds();
  }

  //Tree loop
  Bool_t first = true;
  Int_t stopUpStream=0, stopDownStream=0;
  cout << "Number of entries in tree: " << tree->GetEntries() << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    tree->GetEntry(ev, 0);
    
    //Tree setup
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
      else if (leftrightChoice=="Spin" || leftrightChoice=="" || Cflag ||Gflag){
	cout << "Spin influnced left/right asymmetry" << endl;
      }
      if (Tflag){
	cout << "Trigger mask set to: " << trig << " only" << endl;
      }
      if (Gflag){ cout << "NH3 target set to true values" << endl;}
      
      cout << " " << endl;
      first = false;
    }
    else if ( (Cflag||Gflag)
	     &&(stopUpStream>=nUpStream && stopDownStream>=nDownStream) ){
      break;}


    //Additional Optional cuts
    if (Tflag && (trigMask != trigChoice)) continue;
    if (Mmumu < M_min || Mmumu > M_max) continue;
    if (Gflag && (targetPosition!=0 && targetPosition!=1) ) continue;

    //General useful quantities
    radius = TMath::Sqrt(vx_x*vx_x + vx_y*vx_y);

    //Choose Left/Right
    Double_t phi_photon_lab = ShiftPhiSimple(PhiS_simple);
    Bool_t Left=false, Right=false;
    if (leftrightChoice=="True"){
      //True spectrometer left/right (no spin influence)
      if (phi_photon_lab < TMath::Pi()/2 && phi_photon_lab > -TMath::Pi()/2){
	Left = true;
      }
      else if (phi_photon_lab <3*TMath::Pi()/2 && phi_photon_lab>TMath::Pi()/2){
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
      if (phi_photon_lab<TMath::Pi()/2 && phi_photon_lab>-TMath::Pi()/2
	  && (Spin_0>0||Cflag||Gflag) ){ // Target spin up
	Left = true;
      }
      else if (phi_photon_lab<3*TMath::Pi()/2 && phi_photon_lab>TMath::Pi()/2
	       && (Spin_0>0||Cflag||Gflag) ){// Target spin up
	Right = true;
      }
      else if (phi_photon_lab< 3*TMath::Pi()/2 && phi_photon_lab > TMath::Pi()/2
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

    if (targetPosition == 0) {//UpStream target
      if ( (Cflag||Gflag) &&(stopUpStream>=nUpStream) ) continue;

      if ( (Cflag||Gflag) &&stopUpStream>nUpStream/2.0) {//Switch Left/Right for MC
	Bool_t tmp = Left;
	Left = Right;
	Right = tmp;
      }
      if (Cflag||Gflag) stopUpStream++;

      if (Left){//Polarized Up
	if ( (!Cflag&&!Gflag&&Spin_0>0)
	     || ( (Cflag||Gflag) &&stopUpStream<nUpStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("left_upstream_up", *basicValues[i]);

	  if (Vflag) genericPhys->BinDataCounts(targetPosition,
						"left_upstream_up",*boundValue);

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("left_upstream_up", Mmumu, *basicValues[i]);
	  
	}//Polarized Down
	else if ( (!Cflag&&!Gflag&&Spin_0<0)
		  || ( (Cflag||Gflag) &&stopUpStream>nUpStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("left_upstream_down", *basicValues[i]);
	  
	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition,
				       "left_upstream_down",*boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("left_upstream_down", Mmumu, *basicValues[i]);
	}
      }//End Left
      else if (Right){//Polarized Up
	if ( (!Cflag&&!Gflag&&Spin_0>0)
	     || ( (Cflag||Gflag) &&stopUpStream<nUpStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("right_upstream_up", *basicValues[i]);
	  
	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition, "right_upstream_up",
				       *boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("right_upstream_up", Mmumu, *basicValues[i]);
	  
	}//Polarized Down
	else if ( (!Cflag&&!Gflag&&Spin_0<0)
		  || ( (Cflag||Gflag) &&stopUpStream>nUpStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("right_upstream_down", *basicValues[i]);
	  
	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition, "right_upstream_down",
				       *boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("right_upstream_down", Mmumu, *basicValues[i]);
	}
      }//End Right

      if (!Cflag &&!Gflag && binFlag){//Dilution/Polarization corrections
	Double_t dil = TMath::Abs(dilutionFactor);
	Double_t correct_dil = 0.95*dil;
	Double_t pol = TMath::Abs(Polarization);

	if (Spin_0 > 0) {//Polarized Up
	  for (Int_t i=0; i<nBasics; i++) {
	    Basics[i]->SetCorr("pol_upstream_up", *basicValues[i], pol, Left);
	    Basics[i]->SetCorr("dil_upstream_up", *basicValues[i], correct_dil,
			       Left);
	  }
	  
	  if (Vflag) {
	    genericPhys->SetCorrBounds("pol_upstream_up", *boundValue, pol);
	    genericPhys->SetCorrBounds("dil_upstream_up", *boundValue,
				       correct_dil);
	  }
	}
	else if (Spin_0 < 0){//Polarized Down
	  for (Int_t i=0; i<nBasics; i++) {
	    Basics[i]->SetCorr("pol_upstream_down", *basicValues[i], pol, Left);
	    Basics[i]->SetCorr("dil_upstream_down", *basicValues[i],
			       correct_dil, Left);
	  }
	  
	  if (Vflag) {
	    genericPhys->SetCorrBounds("pol_upstream_down", *boundValue, pol);
	    genericPhys->SetCorrBounds("dil_upstream_down", *boundValue,
				       correct_dil);
	  }
	}
      }
    }//End UpStream target
    else if (targetPosition == 1) {//DownStream target
      if ( (Cflag||Gflag) &&(stopDownStream>=nDownStream) ) continue;

      if ( (Cflag||Gflag)
	  && stopDownStream>nDownStream/2.0) {//Switch Left/Right for MC
	Bool_t tmp = Left;
	Left = Right;
	Right = tmp;
      }
      if (Cflag||Gflag) stopDownStream++;

      if (Left){//Polarized Up
	if ( (!Cflag&&!Gflag&&Spin_0>0)
	     || ((Cflag||Gflag) &&stopDownStream<nDownStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("left_downstream_up", *basicValues[i]);

	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition, "left_downstream_up",
				       *boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("left_downstream_up", Mmumu, *basicValues[i]);
	  
	}//Polarized Down
	else if ( (!Cflag&&!Gflag&&Spin_0<0)
		  || ((Cflag||Gflag) &&stopDownStream>nDownStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("left_downstream_down", *basicValues[i]);

	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition, "left_downstream_down",
				       *boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("left_downstream_down", Mmumu, *basicValues[i]);
	}
      }//End Left
      else if (Right){//Polarized Up
	if ( (!Cflag&&!Gflag&&Spin_0>0)
	     || ( (Cflag||Gflag) &&stopDownStream<nDownStream/2.0) ){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("right_downstream_up", *basicValues[i]);
	  
	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition, "right_downstream_up",
				       *boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("right_downstream_up", Mmumu, *basicValues[i]);
	  
	}//Polarized Down
	else if ( (!Cflag&&!Gflag&&Spin_0<0)
		  || ( (Cflag||Gflag) &&stopDownStream>nDownStream/2.0)){
	  for (Int_t i=0; i<nBasics; i++) 
	    Basics[i]->BinDataCounts("right_downstream_down", *basicValues[i]);
	  
	  if (Vflag) {
	    genericPhys->BinDataCounts(targetPosition, "right_downstream_down",
				       *boundValue);}

	  for (Int_t i=0; i<nBasics-1; i++) 
	    MuMu_b[i]->BinFill("right_downstream_down", Mmumu, *basicValues[i]);
	}
      }//End Right

      if (!Cflag &&!Gflag && binFlag){//Dilution/Polarization corrections
	Double_t dil = TMath::Abs(dilutionFactor);
	Double_t correct_dil = 0.91*dil;
	Double_t pol = TMath::Abs(Polarization);

	if (Spin_0 > 0) {//Polarized Up
	  for (Int_t i=0; i<nBasics; i++) {
	    Basics[i]->SetCorr("pol_downstream_up", *basicValues[i], pol, Left);
	    Basics[i]->SetCorr("dil_downstream_up", *basicValues[i],
			       correct_dil, Left);
	  }
	  
	  if (Vflag) {
	    genericPhys->SetCorrBounds("pol_downstream_up", *boundValue, pol);
	    genericPhys->SetCorrBounds("dil_downstream_up", *boundValue,
				       correct_dil);
	  }
	}
	else if (Spin_0 < 0){//Polarized Down
	  for (Int_t i=0; i<nBasics; i++) {
	    Basics[i]->SetCorr("pol_downstream_down", *basicValues[i], pol,
			       Left);
	    Basics[i]->SetCorr("dil_downstream_down", *basicValues[i],
			       correct_dil, Left);
	  }
	  
	  if (Vflag) {
	    genericPhys->SetCorrBounds("pol_downstream_down", *boundValue, pol);
	    genericPhys->SetCorrBounds("dil_downstream_down", *boundValue,
				       correct_dil);
	  }
	}
      }//end Dilution/Polarization corrections
    }//End DownStream target

  }//End tree loop


  //Asymmetries
  for (Int_t i=0; i<nBasics; i++) Basics[i]->BinLeftRight();
  if (Vflag) genericPhys->BinLeftRight();

  //Dilution and polarization corrections
  if (!Pflag && !Cflag && !Gflag){
    if(binFlag) {
      for (Int_t i=0; i<nBasics; i++) Basics[i]->AvgCorr();
    }

    cout << "\n!!!!!!!!!!!!!!!" << endl;
    cout << "Dilution factor and polarization corrections are made\n" << endl;

    for (Int_t i=0; i<nBasics; i++) Basics[i]->CorrectDilPol();
    if (Vflag) genericPhys->CorrectDilPol();
  }

  //xN->Print_LR("left_upstream_up");
  //xN->Print_LR("right_upstream_up");
  //xN->Print_Asym("asym_upstream_left");
  //xN->PrintCorr("dil_upstream_left");
  //xN->PrintCorr("pol_upstream_left");

  //Graphs and Drawing
  for (Int_t i=0; i<nBasics; i++) Basics[i]->Fill();
  if (Vflag) genericPhys->Fill("all");
      
  //Write output
  if (Qflag || wflag){
    TFile* fout = (Qflag ? new TFile(outFile, "RECREATE")
		   : new TFile("Output.root", "RECREATE") );

    for (Int_t i=0; i<nBasics; i++) Basics[i]->Write();
    if (Vflag) genericPhys->Write();
    for (Int_t i=0; i<nBasics-1; i++) MuMu_b[i]->Write();
    
    fout->Close();
    
    if (Qflag) cout << outFile << "\n file written" << endl;
    else cout << "\nOutput.root file written" << endl;
  }

  cout << "-------------------------------------------" << endl;
  cout << "---------------code finished---------------" << endl;
  cout << "-------------------------------------------" << endl;

  cout << " " << endl;
  cout << "Notes:" << endl;
  cout << "Do not used    PhiS   (it's not defined well)" << endl;
  cout << " " << endl;
  
}
