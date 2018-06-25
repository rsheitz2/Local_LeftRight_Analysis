const Int_t nBins=1; Double_t yMax =0.2; TString physType ="xF"; 
//const Int_t nBins=3; Double_t yMax =0.3; 
//const Int_t nBins=5; Double_t yMax =0.2; 
//TString physType ="M"; //xN, xPi, xF, pT, M

//TString massRange ="HM";
//TString massRange ="JPsi3_326";
//TString massRange ="Psi367_386";
TString massRange ="JPsi25_43";

Bool_t toWrite =false; TString period ="WAll";
TString fNameout ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/June26/Data/";

Int_t upS_up=0, upS_down=1, downS_up=2, downS_down=3;

void SetUpTGraph(TGraphErrors* g){
  g->SetMarkerStyle(21);

  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  g->GetYaxis()->SetRangeUser(-yMax, yMax);

  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);
}


Double_t Amp(Double_t NL[][nBins], Double_t NR[][nBins],
	     Double_t P[][nBins], Int_t bi){
  Int_t nTarg=4;

  Double_t Pol=0.0, Nsum=0.0;
  for (Int_t tr=0; tr<nTarg; tr++) {
    Pol += P[tr][bi]*( NL[tr][bi]+NR[tr][bi] );
    Nsum += NL[tr][bi]+NR[tr][bi]; 
  }
  Pol /= Nsum;

  Double_t Lup, Rup;
  Lup = NL[upS_up][bi]*NR[upS_down][bi]; Rup = NR[upS_up][bi]*NL[upS_down][bi];

  Double_t Ldown, Rdown;
  Ldown = NL[downS_down][bi]*NR[downS_up][bi];
  Rdown = NR[downS_down][bi]*NL[downS_up][bi];
  //Ldown = NL[downS_up][bi]*NR[downS_down][bi];
  //Rdown = NR[downS_up][bi]*NL[downS_down][bi];


  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;

  return A;
}


Double_t e_Amp(Double_t NL[][nBins], Double_t NR[][nBins],
	       Double_t P[][nBins], Int_t bi){
  Int_t nTarg=4;

  Double_t Pol=0.0, Nsum=0.0;
  for (Int_t tr=0; tr<nTarg; tr++) {
    Pol += P[tr][bi]*( NL[tr][bi]+NR[tr][bi] );
    Nsum += NL[tr][bi]+NR[tr][bi];
  }
  Pol /= Nsum;

  
  Double_t Lup, Rup;
  Lup = NL[upS_up][bi]*NR[upS_down][bi]; Rup = NR[upS_up][bi]*NL[upS_down][bi];
  Double_t Ldown, Rdown;
  Ldown = NL[downS_down][bi]*NR[downS_up][bi];
  Rdown = NR[downS_down][bi]*NL[downS_up][bi];

  Double_t L=Lup*Ldown, R=Rup*Rdown;
  L = TMath::Power(L, 0.25);
  R = TMath::Power(R, 0.25);

  
  Double_t LinvSum = 1.0/NL[upS_up][bi] + 1.0/NR[upS_down][bi] +
    1.0/NL[downS_down][bi] + 1.0/NR[downS_up][bi];
  Double_t RinvSum = 1.0/NR[upS_up][bi] + 1.0/NL[upS_down][bi] +
    1.0/NR[downS_down][bi] + 1.0/NL[downS_up][bi];
  
  
  Double_t dL = 0.25*L * TMath::Sqrt( LinvSum );
  Double_t dR = 0.25*R * TMath::Sqrt( RinvSum );

  Double_t e = L - R;
  e /= ( L + R );

  Double_t error = TMath::Sqrt( (1-e)*(1-e)*dL*dL + (1+e)*(1+e)*dR*dR ) ;
  error /= ( L + R );
  error /= Pol;

  return error;
}


void allTarg_geoMeanFalseAN(TString fname=""){
  if (fname==""){
    fname += "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/AcceptCorrections/Data/";

    cout << "Using default data from ./AcceptCorrections/Data/" << endl;
    cout << "Data originally made in Presents/May1" << endl;
    cout << " " << endl;
  }
  

  TFile *f_LR
    = TFile::Open(Form("%s/leftRight_byTarget_%s_%s_%i.root", 
		       fname.Data(), period.Data(), massRange.Data(), nBins) );
  TFile *f_LR_noCorr
    = TFile::Open(Form("%s/leftRight_byTarget_%s_%s_%i_noCorr.root",
		       fname.Data(), period.Data(), massRange.Data(), nBins) );

  if ( !(f_LR) || !(f_LR_noCorr) ){//Basic file checks
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nTarg = 4;
  TString targName[nTarg] = {"asym_upstream_up", "asym_upstream_down",
			     "asym_downstream_up", "asym_downstream_down"};
  
  
  TGraphErrors *g_LR[nTarg], *g_LR_noCorr[nTarg];
  for (Int_t i=0; i<nTarg; i++) {
    g_LR[i] = (TGraphErrors*)f_LR->Get(Form("%s_%s", physType.Data(),
					    targName[i].Data() ) );
    g_LR_noCorr[i] = (TGraphErrors*)f_LR_noCorr->Get(Form("%s_%s",
							  physType.Data(),
							  targName[i].Data()) );
  }

  //Get Left/Right counts and Polarization values
  Double_t N_L[nTarg][nBins], N_R[nTarg][nBins], Pol[nTarg][nBins];
  Double_t ex[nBins]= {0.};
  Double_t *xvals = g_LR[0]->GetX();
  for (Int_t tr=0; tr<nTarg; tr++) {
    Double_t *y_lr = g_LR[tr]->GetY();
    Double_t *y_lr_noCorr = g_LR_noCorr[tr]->GetY();
    Double_t *e_y_lr = g_LR[tr]->GetEY();
    Double_t *e_y_lr_noCorr = g_LR_noCorr[tr]->GetEY();

    for (Int_t bi=0; bi<nBins; bi++) {
      Double_t A = y_lr_noCorr[bi];
      Double_t dA = e_y_lr_noCorr[bi];
      
      N_L[tr][bi] = (1-A)*(1+A)*(1+A)/(2*dA*dA);
      N_R[tr][bi] = (1+A)*(1-A)*(1-A)/(2*dA*dA);
      
      Pol[tr][bi] = y_lr_noCorr[bi]/y_lr[bi];

      
      if (std::isnan( Pol[tr][bi] ) ) {//basic check
	cout << " " << endl;
	cout << "Warning Pol is nan for bin:  " << bi <<
	  "  target:  " << targName[tr] << endl;
	cout << "(tr > 0) ? Pol[tr][bi] = Pol[tr-1][bi] : Pol[tr][bi] =" <<
	  " Pol[tr+1][bi];" << endl;
	cout << " " << endl;

	(tr > 0) ? Pol[tr][bi] = Pol[tr-1][bi] : Pol[tr][bi] = Pol[tr+1][bi];
      }
      if (std::isnan( Pol[tr][bi] ) ) {//basic check
	cout << " " << endl;
	cout << "Pol is nan" << endl;
	cout << " " << endl;
	exit(EXIT_FAILURE);
      }
      
    }
  }


  Double_t AN[nBins], e_AN[nBins];
  for (Int_t bi=0; bi<nBins; bi++) {
    AN[bi] = Amp(N_L, N_R, Pol, bi);
    e_AN[bi] = e_Amp(N_L, N_R, Pol, bi);
    
    cout << "AN   " << AN[bi] << "   e_AN  "<< e_AN[bi] << endl;
  }

  
  TCanvas* c1 = new TCanvas(); 
  TGraphErrors *g_AN = new TGraphErrors(nBins, xvals, AN, ex, e_AN);
  SetUpTGraph(g_AN);
  g_AN->Draw("AP");
  
  Double_t xmin =g_AN->GetXaxis()->GetXmin();
  Double_t xmax =g_AN->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0., xmax, 0.);
  li->SetLineColor(kBlue); li->SetLineStyle(8);
  li->Draw("same");

  
  fNameout+="allTarg_geoMeanFalseAN_";
  fNameout+=Form("%i_%s_%s_%s.root", nBins, physType.Data(), period.Data(),
		 massRange.Data() );
  if (toWrite){
    TFile *fOutput = new TFile(fNameout, "RECREATE");
    g_AN->Write("AN");
    fOutput->Close();
  }

  cout << " " << endl;
  cout << "Period is:    "  << period << endl;
  cout << "Physics variable:  " << physType << endl;
  cout << "Number of bins considered: " << nBins << endl;
  cout << "Mass range considered:  " << massRange << endl;
  cout << "File: " << fNameout << "  was writen:  " << toWrite << endl;
  cout << " " << endl;
}
