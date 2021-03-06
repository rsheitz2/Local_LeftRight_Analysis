#include "common.h"
#include "functions.h"
#include "leftright.h"
#include "twotargets.h"
#include "lrSpinCorr.h"
#include "genericBounds.h"
#include "lr_tgraph.h"
#include "binDist.h"
#include "fourTargets.h"
#include "binDist4Targ.h"

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){//Starting Help information
    cout << "\n\nUsage:" << endl;
    cout << "./falseAsym_leftRight [options] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "\n----Writing Options----"<< endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)";
    cout << "\n\n----Additional Setup Options----"<< endl;
    cout << "  --Real Data Only Options--"<< endl;
    cout << "Option:  -S Left/right asymmetry choice" << endl;
    cout << "	 (True=no spin influence, Spin=spin influence, default=Spin)";
    cout << "\n    (Default=Spin) (don't use this with -C option)" << endl;
    cout << "Option:  -E       (Equal number of runs per sub period)" << endl;
    cout << "Option:  -R       (Even runs sb 1, odd runs sb 2)" << endl;
    cout << "Option:  -T trig       (Only specific trigger)" << endl;
    cout << "                  (\"LL\"=last-last, \"LO\"=last-outer, " <<
      "\"LL_LO\"=last-last && last-outer)" << endl;
    cout << "Option:  -i min       (Minimum mass to consider)" << endl;
    cout << "    (Default=0.0)" << endl;
    cout << "Option:  -a max       (Maximum mass to consider)" << endl;
    cout << "    (Default=12.0)" << endl;
    cout << "\n----Changing Binning Options----" << endl;
    cout << "Option:  -b textfile with binning information	";
    cout << "(textfile should be made from Macro/Binning/avgBinBounds.C)"<<endl;
    cout << "Option:  -Z change number of bins for histogram distributions/n";
    cout << "    (default nHbins is 150)" << endl;
    cout << "\n----Additional Options----" << endl;
    cout << "Option:  -D    (Debug mode, only 1000 events considered)" << endl;
    cout << "Option:  -H    (to see example inputs for this program)\n" << endl;
	
    exit(EXIT_FAILURE);
  }

  //Read input arguments
  Int_t wflag=0, Qflag=0, fflag=0, Sflag=0, Tflag=0, Eflag=0, Rflag=0;
  Int_t iflag=0, aflag=0, binFlag=0, Dflag=0, Hflag=0, Zflag=0;
  Int_t c;
  TString fname = "", outFile = "", leftrightChoice="", trig="";
  TString binFile = "";
  Double_t M_min=0.0, M_max=12.0;
  Int_t nHbins=150;
  
  while ((c = getopt (argc, argv, "wf:Q:S:ERT:i:a:b:M:DHZ:")) != -1) {
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
    case 'E':
      Eflag = 1;
      break;
    case 'R':
      Rflag = 1;
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
    case 'D':
      Dflag = 1;
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
    cout << "\nOption -T" << trig << " is not a valid choice" << endl;
    exit(EXIT_FAILURE);
  }
  if (Eflag) cout << "Runs per sub period evened out" << endl;
  if (Rflag) cout << "Even runs sb 1, Odd runs sb 2" << endl;
  if(Dflag) cout << "\nDebug mode only 1000 events considered\n" << endl;
  if (Sflag) cout <<"Left/right choice set to:   " << leftrightChoice << endl;
  if(Zflag){ cout << "Number of histogram bins changed to  " << nHbins << endl;}
  cout << "\nUsing REAL data!\n" << endl;

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
  Long64_t RunNum, SpillNum, event;
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
  //Drell-Yan Angles
  Errors+=tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  Errors+=tree->SetBranchAddress("Theta_CS", &Theta_CS);
  Errors+=tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //Mu plus
  Errors+=tree->SetBranchAddress("theta_traj1", &theta_traj1);
  //Mu Minus
  Errors+=tree->SetBranchAddress("theta_traj2", &theta_traj2);
  //Event
  Errors+=tree->SetBranchAddress("trigMask", &trigMask);
  Errors+=tree->SetBranchAddress("MasterTrigMask", &MasterTrigMask);
  Errors+=tree->SetBranchAddress("RunNum", &RunNum);
  Errors+=tree->SetBranchAddress("SpillNum", &SpillNum);
  Errors+=tree->SetBranchAddress("event", &event);
  //Target values
  Errors+=tree->SetBranchAddress("Spin_0", &Spin_0);
  Errors+=tree->SetBranchAddress("Spin_1", &Spin_1);
  Errors+=tree->SetBranchAddress("Spin_2", &Spin_2);
  Errors+=tree->SetBranchAddress("Spin_3", &Spin_3);
  Errors+=tree->SetBranchAddress("Spin_4", &Spin_4);
  Errors+=tree->SetBranchAddress("Spin_5", &Spin_5);
  Errors+=tree->SetBranchAddress("Spin_6", &Spin_6);
  Errors+=tree->SetBranchAddress("dilutionFactor", &dilutionFactor);
  Errors+=tree->SetBranchAddress("Polarization", &Polarization);
  //DY-variables
  Errors+=tree->SetBranchAddress("x_beam", &x_beam);
  Errors+=tree->SetBranchAddress("x_target", &x_target);
  Errors+=tree->SetBranchAddress("x_feynman", &x_feynman);
  Errors+=tree->SetBranchAddress("q_transverse", &q_transverse);
  Errors+=tree->SetBranchAddress("Mmumu", &Mmumu);

  if(Errors) {
    cout << "\nErrors in Setting tree branch addresses" << endl;
    cout << "Probably some variables are not defined in your file\n" << endl;
    exit(EXIT_FAILURE);
  }

  //Define variables
  leftright *xN, *xPi, *xF, *pT, *Mass;
  if (binFlag){
    xN = new FourTargets(binFile, "xN");
    xPi = new FourTargets(binFile, "xPi");
    xF = new FourTargets(binFile, "xF");
    pT = new FourTargets(binFile, "pT");
    Mass = new FourTargets(binFile, "mass", "M");
  }
  else{
    cout << "binFlag must be specifed for now because of new asymmetries"<<endl;
    exit(EXIT_FAILURE);
  }
  const Int_t nBasics = 5;
  leftright *Basics[nBasics] = {xN, xPi, xF, pT, Mass};
  Double_t *basicValues[nBasics] = {&x_target, &x_beam, &x_feynman,
				    &q_transverse, &Mmumu};

  binDist *MuMu_b_xN =
    new binDist4Targ("MuMu", "xN", binFile, nHbins, M_min, M_max);
  binDist *MuMu_b_xPi =
    new binDist4Targ("MuMu","xPi",binFile, nHbins, M_min, M_max);
  binDist *MuMu_b_xF =
    new binDist4Targ("MuMu", "xF", binFile, nHbins, M_min, M_max);
  binDist *MuMu_b_pT =
    new binDist4Targ("MuMu", "pT", binFile, nHbins, M_min, M_max);
  binDist *MuMu_b[nBasics-1] = {MuMu_b_xN, MuMu_b_xPi, MuMu_b_xF, MuMu_b_pT};

  Int_t tree_entries = (!Dflag) ? tree->GetEntries() : 1000;//Tree Loop

  //Tree loop
  Bool_t first = true;
  cout << "Number of entries in tree: " << tree->GetEntries() << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    tree->GetEntry(ev, 0);
   
    //Tree setup
    if (first || ev==tree_entries-1){
      cout << "\nSetup!!!!!!!!!!!!!!!" << endl;
      
      if (iflag || aflag){
	cout << "Additional mass cut" << endl;
	cout << "    Mass range " << M_min << " - " << M_max << endl;
      }
      
      if(leftrightChoice=="True") 
	cout << "True left/right asymmetry (no spin influence)" << endl;
      else if (leftrightChoice=="Spin" || leftrightChoice=="")
	cout << "Spin influnced left/right asymmetry" << endl;
	
      if (Tflag)
	cout << "Trigger mask set to: " << trig << " only" << endl;
            
      cout << " " << endl;
      first = false;
    }

    //Additional Optional cuts
    if (Tflag && (trigMask != trigChoice)) continue;
    if (Mmumu < M_min || Mmumu > M_max) continue;

    if (Eflag){//Equal out runs per period (determined for HMDY)
      Int_t per=0;
      if ( (RunNum >= 259363) && (RunNum <= 259677) ){ per = 1; }//W07_sb1
      else if ( (RunNum >= 259744) && (RunNum <= 260016) ){ per = 2; }//W07_sb2
      else if ( (RunNum >= 260074) && (RunNum <= 260264) ){ per = 3; }//W08_sb1
      else if ( (RunNum >= 260317) && (RunNum <= 260565) ){ per = 4; }//W08_sb2
      else if ( (RunNum >= 260627) && (RunNum <= 260852) ){ per = 5; }//W09_sb1
      else if ( (RunNum >= 260895) && (RunNum <= 261496) ){ per = 6; }//W09_sb2
      else if ( (RunNum >= 261515) && (RunNum <= 261761) ){ per = 7; }//W10_sb1
      else if ( (RunNum >= 261970) && (RunNum <= 262221) ){ per = 8; }//W10_sb2
      else if ( (RunNum >= 262370) && (RunNum <= 262772) ){ per = 9; }//W11_sb1
      else if ( (RunNum >= 262831) && (RunNum <= 263090) ){ per = 10; }//W11_sb2
      else if ( (RunNum >= 263143) && (RunNum <= 263347) ){ per = 11; }//W12_sb1
      else if ( (RunNum >= 263386) && (RunNum <= 263603) ){ per = 12; }//W12_sb2
      else if ( (RunNum >= 263655) && (RunNum <= 263853) ){ per = 13; }//W13_sb1
      else if ( (RunNum >= 263926) && (RunNum <= 264134) ){ per = 14; }//W13_sb2
      else if ( (RunNum >= 264170) && (RunNum <= 264330) ){ per = 15; }//W14_sb1
      else if ( (RunNum >= 264429) && (RunNum <= 264562) ){ per = 16; }//W14_sb2
      else if ( (RunNum >= 264619) && (RunNum <= 264672) ){ per = 17; }//W15_sb1
      else if ( (RunNum >= 264736) && (RunNum <= 264857) ){ per = 18; }//W15_sb2
      else { cout << "Run period not specified " << RunNum << endl; }

      if ( (per==1||per==2) && (RunNum>259944) ) continue; //W07
      else if ( (per==3||per==4) && (RunNum>260246) ) continue; //W08
      else if ( (per==5||per==6) && (RunNum>260848) ) continue; //W09
      else if ( (per==7||per==8) && (RunNum>261757) ) continue; //W10
      else if ( (per==9||per==10) && (RunNum>263030) ) continue; //W11
      else if ( (per==11||per==12) && (RunNum>263600) ) continue; //W12
      else if ( (per==13||per==14) && (RunNum>264053) ) continue; //W13
      else if ( (per==15||per==16) && (RunNum>264307) ) continue; //W14
      else if ( (per==17||per==18) && (RunNum>264840) ) continue; //W15
    }
    
    //General useful quantities
    //Double_t radius = TMath::Sqrt(vx_x*vx_x + vx_y*vx_y);
    Double_t dil = TMath::Abs(dilutionFactor), correct_dil;
    Double_t pol = TMath::Abs(Polarization);
    Double_t upS_vxSplit =-266.847;//HMDY
    Double_t downS_vxSplit =-191.829;//HMDY
    if ( (vx_z >= -294.5) && (vx_z <= upS_vxSplit) ) {//upSup
      targetPosition = 0;      
      correct_dil = 0.95*dil;
    }
    else if ( (vx_z > upS_vxSplit) && (vx_z <= -239.3) ) {//upSdown
      targetPosition = 1;
      correct_dil = 0.95*dil;
    }
    else if ( (vx_z >= -219.5) && (vx_z <= downS_vxSplit) ) {//downSup
      targetPosition = 2;
      correct_dil = 0.91*dil;
    }
    else if ( (vx_z > downS_vxSplit) && (vx_z <= -164.3) ) {//downSdown
      targetPosition = 3;
      correct_dil = 0.91*dil;
    }
    else {
      cout << vx_z << "  No target position defined!\n" << endl;
      exit(EXIT_FAILURE);
    }

    //Choose Left/Right
    Double_t phi_photon_lab = ShiftPhiSimple(PhiS_simple);
    if (Rflag){//Even runs upS pol up/downS pol down, opposite for odd runs
      if (RunNum%2){//Odd runs
	if (targetPosition<2) Spin_0 =-1.0;//Upstream
	else Spin_0 = 1.0;//Downstream
      }
      else{//Even runs
	if (targetPosition<2) Spin_0 =1.0;//Upstream
	else Spin_0 = -1.0;//Downstream
      }
    }
    Bool_t Left = ChooseLeftRight(leftrightChoice, phi_photon_lab, Spin_0);
    
    for (Int_t i=0; i<nBasics; i++) {
      Basics[i]->BinDataCounts(targetPosition, Left, Spin_0, *basicValues[i]);
      Basics[i]->SetCorr(*basicValues[i], pol, correct_dil,
			 targetPosition, Left, Spin_0);
    }
	  
    for (Int_t i=0; i<nBasics-1; i++)
      MuMu_b[i]->BinFill(targetPosition, Left, Spin_0, Mmumu, *basicValues[i]);

  }//End tree loop

  //Dilution and polarization averages made
  if(binFlag) {for (Int_t i=0; i<nBasics; i++) {Basics[i]->AvgCorr();}}

  xN->Print_PolDil("upSup");
  //xN->Print_LR("upSdown");
      
  if (Qflag || wflag){//Write output
    TFile* fout = (Qflag ? new TFile(outFile, "RECREATE")
		   : new TFile("Output.root", "RECREATE") );

    for (Int_t i=0; i<nBasics; i++) Basics[i]->WriteAll();
    for (Int_t i=0; i<nBasics-1; i++) MuMu_b[i]->Write();
    
    fout->Close();
    
    if (Qflag) cout << outFile << "\n file written" << endl;
    else cout << "\nOutput.root file written" << endl;
  }

  cout << "-------------------------------------------" << endl;
  cout << "---------------code finished---------------" << endl;
  cout << "-------------------------------------------" << endl;
}
