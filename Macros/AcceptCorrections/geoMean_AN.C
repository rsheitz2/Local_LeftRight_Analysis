//const Int_t nBins=1; Double_t dx =0.02; Double_t yMax =0.4; 
const Int_t nBins=3; Double_t dx =0.005; Double_t yMax =0.5;

TString physType ="xF", period ="WAll";
Bool_t toWrite =true;
TString fNameout ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/June26/Data/";


void SetUpTGraph(TGraphErrors* g, TString name, Int_t ic, Double_t offset){
  g->SetTitle(name);
  
  g->SetMarkerStyle(21);
  g->SetMarkerColor(ic);

  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  g->GetYaxis()->SetRangeUser(-yMax, yMax);

  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  Double_t *xval = g->GetX();
  for (Int_t i=0; i<nBins; i++) xval[i] += offset;
}


Double_t Amp(Double_t L1, Double_t R1, Double_t L2, Double_t R2,
	     Double_t P1, Double_t P2){
  Double_t Pol = ( P1*(L1+R1) + P2*(L2+R2) )/(L1+R1+L2+R2);

  Double_t L = TMath::Sqrt( L1*L2 );
  Double_t R = TMath::Sqrt( R1*R2 );

  Double_t A = L - R;
  A /= ( L + R );
  A /= Pol;

  return A;
}


Double_t e_Amp(Double_t L1, Double_t R1, Double_t L2, Double_t R2,
	     Double_t P1, Double_t P2){
  Double_t Pol = ( P1*(L1+R1) + P2*(L2+R2) )/(L1+R1+L2+R2);
  
  Double_t L = TMath::Sqrt( L1*L2 );
  Double_t R = TMath::Sqrt( R1*R2 );
  Double_t dL = 0.5*L * TMath::Sqrt( 1/L1 + 1/L2 );
  Double_t dR = 0.5*R * TMath::Sqrt( 1/R1 + 1/R2 );

  Double_t e = L - R;
  e /= ( L + R );

  Double_t error = TMath::Sqrt( (1-e)*(1-e)*dL*dL + (1+e)*(1+e)*dR*dR ) ;
  error /= ( L + R );
  error /= Pol;

  return error;
}


void geoMean_AN(TString fname=""){
  if (fname==""){
    fname += "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Presents/May1/Macros/Accept";

    cout << "Using default data from May1/Macros/Accept" << endl;
  }
  

  TFile *f_LR = TFile::Open(Form("%s/%s_%i.root", fname.Data(),
				 period.Data(), nBins) );
  TFile *f_LR_noCorr = TFile::Open(Form("%s/%s_%i_noCorr.root", fname.Data(),
					period.Data(), nBins) );

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
  Double_t N_L[nTarg][nBins], N_R[nTarg][nBins];
  Double_t xvals[nBins], Pol[nTarg][nBins];
  Double_t ex[nBins]= {0.};
  Double_t *x_lr = g_LR[0]->GetX();
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
      
      if (tr==0) xvals[bi] = x_lr[bi];
    }
  }

  const Int_t nAmp = 2;
  Int_t icolor[nAmp] = {3, 4};
  Double_t offsets[nAmp]; for (Int_t i=0; i<nAmp; i++) offsets[i] = i*dx;
  TString ampName[nAmp] = {"AN_upstream", "AN_downstream"};
  
  Double_t AN[nAmp][nBins], e_AN[nAmp][nBins];
  for (Int_t am=0; am<nAmp; am++) {
    for (Int_t bi=0; bi<nBins; bi++) {
      AN[am][bi] = Amp(N_L[am*2][bi], N_R[am*2][bi],
		       N_L[am*2+1][bi], N_R[am*2+1][bi],
		       Pol[am*2][bi], Pol[am*2+1][bi]);
      e_AN[am][bi] = e_Amp(N_L[am*2][bi], N_R[am*2][bi],
			   N_L[am*2+1][bi], N_R[am*2+1][bi],
			   Pol[am*2][bi], Pol[am*2+1][bi]);
    }
  }

  
  TCanvas* c1 = new TCanvas(); 
  TGraphErrors *g_AN[nAmp];
  TLine *li;
  for (Int_t am=0; am<nAmp; am++) {
    g_AN[am] = new TGraphErrors(nBins, xvals, AN[am], ex, e_AN[am]);
    SetUpTGraph(g_AN[am], ampName[am], icolor[am], offsets[am] );

    if (am==0) {
      g_AN[am]->Draw("AP");

      Double_t xmin =g_AN[am]->GetXaxis()->GetXmin();
      Double_t xmax =g_AN[am]->GetXaxis()->GetXmax();
      li = new TLine(xmin, 0., xmax, 0.);
      li->SetLineColor(kBlue);
      li->SetLineStyle(8);
      li->Draw("same");
    }
    else g_AN[am]->Draw("Psame");
  }

  fNameout+="geoMean_";
  fNameout+=Form("%i_%s_%s.root", nBins, physType.Data(), period.Data() );
  if (toWrite){
    TFile *fOutput = new TFile(fNameout, "RECREATE");
    for (Int_t am=0; am<nAmp; am++) {
      g_AN[am]->Write(ampName[am]);
    }
    fOutput->Close();
  }

  cout << " " << endl;
  cout << "Period is:    "  << period << endl;
  cout << "Acceptance from:  " << physType << endl;
  cout << "Number of bins considered: " << nBins << endl;
  cout << "File: " << fNameout << "  was writen:  " << toWrite << endl;
  cout << " " << endl;
}
