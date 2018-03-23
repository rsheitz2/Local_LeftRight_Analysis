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
  Int_t nBounds = tv_xN_bounds.GetNoElements();
  Double_t xN_bounds[nBounds], xPi_bounds[nBounds], xF_bounds[nBounds];
  Double_t pT_bounds[nBounds], M_bounds[nBounds];
  Double_t bounds[nBounds];
  for (Int_t i=0; i<nBounds; i++) {
    xN_bounds[i] = tv_xN_bounds[i];
    xPi_bounds[i] = tv_xPi_bounds[i];
    xF_bounds[i] = tv_xF_bounds[i];
    pT_bounds[i] = tv_pT_bounds[i];
    M_bounds[i] = tv_M_bounds[i];

    //Generic bounds
    bounds[i] = M_bounds[i];
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
  Double_t PhiS_simple, Theta_CS, vOpenAngle;
  //Muons
  Double_t phi_traj1, phi_traj2;
  Double_t theta_traj1, theta_traj2;
  Double_t qP_traj1, qP_traj2;
  //Event
  Int_t trigMask, MasterTrigMask;
  //Target values
  Double_t Spin_0, Spin_1, Spin_2, Spin_3, Spin_4, Spin_5, Spin_6;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;
  
  //Vertex specific
  errors += tree->SetBranchAddress("vx_z", &vx_z);
  errors += tree->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  errors += tree->SetBranchAddress("PhiS_simple", &PhiS_simple);
  errors += tree->SetBranchAddress("Theta_CS", &Theta_CS);
  errors += tree->SetBranchAddress("vOpenAngle", &vOpenAngle);
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

  if (errors){
    cout << " " << endl;
    cout << "Errors opening trees variables" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
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
  Int_t tree_entries = tree->GetEntries();//Tree Loop
  cout << "Number of entries in tree: " << tree_entries << endl;
  for (Int_t ev=0; ev<tree_entries; ev++) {
    //cout << "Debug mode" << endl; for (Int_t ev=0; ev<1000; ev++) {
    tree->GetEntry(ev, 0);

    if (first || ev==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!!!!!!!" << endl;

      first = false;
    }

    Double_t boundValue = Mmumu;

    if (trigMask == 65792){
      BinDataFill(hPhi1_LL, phi_traj1, nBins, boundValue, bounds);
      BinDataFill(hPhi2_LL, phi_traj2, nBins, boundValue, bounds);

      BinDataFill(hTheta1_LL, theta_traj1, nBins, boundValue, bounds);
      BinDataFill(hTheta2_LL, theta_traj2, nBins, boundValue, bounds);

      BinDataFill(hP1_LL, qP_traj1, nBins, boundValue, bounds);
      BinDataFill(hP2_LL, -1.0*qP_traj2, nBins, boundValue, bounds);
    }
    else if (trigMask == 65540){
      BinDataFill(hPhi1_LO, phi_traj1, nBins, boundValue, bounds);
      BinDataFill(hPhi2_LO, phi_traj2, nBins, boundValue, bounds);

      BinDataFill(hTheta1_LO, theta_traj1, nBins, boundValue, bounds);
      BinDataFill(hTheta2_LO, theta_traj2, nBins, boundValue, bounds);

      BinDataFill(hP1_LO, qP_traj1, nBins, boundValue, bounds);
      BinDataFill(hP2_LO, -1.0*qP_traj2, nBins, boundValue, bounds);
    }
    else if (trigMask == 65796){
      BinDataFill(hPhi1_LL_LO, phi_traj1, nBins, boundValue, bounds);
      BinDataFill(hPhi2_LL_LO, phi_traj2, nBins, boundValue, bounds);

      BinDataFill(hTheta1_LL_LO, theta_traj1, nBins, boundValue, bounds);
      BinDataFill(hTheta2_LL_LO, theta_traj2, nBins, boundValue, bounds);

      BinDataFill(hP1_LL_LO, qP_traj1, nBins, boundValue, bounds);
      BinDataFill(hP2_LL_LO, -1.0*qP_traj2, nBins, boundValue, bounds);
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
    cPhi[i] = new TCanvas(Form("cPhi_%i", i) );
    cPhi[i]->Divide(3, 5); }
  for (Int_t i=0; i<nBins/5.0; i++) {
    cTheta[i] = new TCanvas(Form("cTheta_%i", i) );
    cTheta[i]->Divide(3, 5); }
  for (Int_t i=0; i<nBins/5.0; i++) {
    cP[i] = new TCanvas(Form("cP_%i", i) );
    cP[i]->Divide(3, 5); }
  
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
