const Int_t nBins=1; Double_t dx =0.02; Double_t yMax =0.1; 
//const Int_t nBins=3; Double_t dx =0.005; Double_t yMax =0.5;
//const Int_t nBins=5; Double_t dx =0.005; Double_t yMax =0.5;

Bool_t accCorrected=false;
TString physType ="xF", period ="WAll";
//TString massRange ="HM";
//TString massRange ="JPsi3_326";
TString massRange ="JPsi25_43";
Bool_t toWrite =false;
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


Double_t Amp(Double_t L, Double_t R, Double_t aL, Double_t aR,
	     Double_t P){
  Double_t A = (L*aR - R*aL);
  A /= (L*aR + R*aL );
  A /= P;

  return A;
}


Double_t e_Amp(Double_t L, Double_t R, Double_t aL, Double_t aR,
	     Double_t P){
  Double_t sum = L*aR + R*aL;
  
  Double_t e = 4. * L*aR*aL * R*aL*aR;
  e *= L+R;
  e /= (sum * sum * sum * sum);
  e = TMath::Sqrt( e );
  e /= P;

  return e;
}


void AN_accCorrected(TString fname=""){
  if (fname==""){
    fname += "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/AcceptCorrections/Data";

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
  TFile *f_acc = NULL;
  if (accCorrected){
    f_acc = TFile::Open(Form("%s/Acceptance_%s_%s_%s_%i.root",
			     fname.Data(), physType.Data(), period.Data(),
			     massRange.Data(), nBins) );
  }
    

  //Basic file checks
  if ( !(f_LR) || !(f_LR_noCorr) || (accCorrected && !(f_acc) ) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nTarg = 4;
  TString targName[nTarg] = {"asym_upstream_up", "asym_upstream_down",
			     "asym_downstream_up", "asym_downstream_down"};
  TString accName[nTarg] = {"acc_UpS_Left", "acc_UpS_Right",
			     "acc_DownS_Left", "acc_DownS_Right"};
  Int_t icolor[nTarg] = {3, 4, 8, 9}; 
  Double_t offsets[nTarg];
  for (Int_t i=0; i<nTarg; i++) offsets[i] = i*dx;

  TGraphErrors *g_LR[nTarg], *g_LR_noCorr[nTarg], *g_acc[nTarg];
  for (Int_t i=0; i<nTarg; i++) {
    g_LR[i] = (TGraphErrors*)f_LR->Get(Form("%s_%s", physType.Data(),
					    targName[i].Data() ) );
    g_LR_noCorr[i] = (TGraphErrors*)f_LR_noCorr->Get(Form("%s_%s",
							  physType.Data(),
							  targName[i].Data()) );
    if (accCorrected) g_acc[i] = (TGraphErrors*)f_acc->Get(accName[i]);
  }

  Double_t AN[nTarg][nBins], e_AN[nTarg][nBins];
  Double_t xvals[nBins];
  Double_t ex[nBins]= {0.};
  
  Double_t *y_acc_UpS_L = (accCorrected) ? g_acc[0]->GetY() : NULL;
  Double_t *y_acc_UpS_R = (accCorrected) ? g_acc[1]->GetY() : NULL;
  Double_t *y_acc_DownS_L = (accCorrected) ? g_acc[2]->GetY() : NULL;
  Double_t *y_acc_DownS_R = (accCorrected) ? g_acc[3]->GetY() : NULL;
  Double_t *x_lr = g_LR[0]->GetX();
  for (Int_t tr=0; tr<nTarg; tr++) {
    Double_t *y_lr = g_LR[tr]->GetY();
    Double_t *y_lr_noCorr = g_LR_noCorr[tr]->GetY();
    Double_t *e_y_lr = g_LR[tr]->GetEY();
    Double_t *e_y_lr_noCorr = g_LR_noCorr[tr]->GetEY();

    for (Int_t bi=0; bi<nBins; bi++) {
      Double_t A = y_lr_noCorr[bi];
      Double_t dA = e_y_lr_noCorr[bi];
      
      Double_t N_L = (1-A)*(1+A)*(1+A)/(2*dA*dA);
      Double_t N_R = (1+A)*(1-A)*(1-A)/(2*dA*dA);

      Double_t aL, aR;
      if (accCorrected){
	if (tr==0 || tr==1){//Upstream targets
	  aL = y_acc_UpS_L[bi];
	  aR = y_acc_UpS_R[bi];
	}
	else{
	  aL = y_acc_DownS_L[bi];
	  aR = y_acc_DownS_R[bi];
	}

	if (tr==1 || tr==3){//Switch left/right acceptance for polarization down
	  Double_t tmp;
	  tmp = aL;
	  aL = aR;
	  aR = tmp;
	}
      }
      
      Double_t P = y_lr_noCorr[bi]/y_lr[bi];

      if (accCorrected){
	AN[tr][bi] = Amp(N_L, N_R, aL, aR, P); //acc. correct
	e_AN[tr][bi] = e_Amp(N_L, N_R, aL, aR, P);
      }
      else {
	AN[tr][bi] = y_lr[bi];//no corrections
	e_AN[tr][bi] = e_y_lr[bi];
      }
      
      xvals[bi] = x_lr[bi];
    }
  }


  TCanvas* c4 = new TCanvas(); c4->Divide(2,2);
  TCanvas* c1 = new TCanvas();
  TGraphErrors *g_AN[nTarg];
  TLine *li;
  for (Int_t i=0; i<nTarg; i++) {
    g_AN[i] = new TGraphErrors(nBins, xvals, AN[i], ex, e_AN[i]);
    

    if (i==0) {
      Double_t xmin =g_AN[i]->GetXaxis()->GetXmin();
      Double_t xmax =g_AN[i]->GetXaxis()->GetXmax();
      li = new TLine(xmin, 0., xmax, 0.);
      li->SetLineColor(kBlue);
      li->SetLineStyle(8);
    }

    
    SetUpTGraph(g_AN[i], targName[i], icolor[i], offsets[i] );
    c4->cd(i+1);
    g_AN[i]->Draw("AP");
    li->Draw("same");

    c1->cd();
    (i==0) ? g_AN[i]->Draw("AP") : g_AN[i]->Draw("Psame");
    li->Draw("same");
  }


  if (accCorrected) fNameout+="AN_accCorrected_";
  else fNameout+="AN_accNotCorrected_";
  fNameout+=Form("%i_%s_%s_%s.root", nBins, physType.Data(), period.Data(),
		 massRange.Data() );
  TString AN_name[nTarg] = {"AN_upstream_up", "AN_upstream_down",
			    "AN_downstream_up", "AN_downstream_down"};
  if (toWrite){
    TFile *fOutput = new TFile(fNameout, "RECREATE");
    for (Int_t tr=0; tr<nTarg; tr++) g_AN[tr]->Write( AN_name[tr] );
    fOutput->Close();
  }


  cout << " " << endl;
  cout << "Period is:    "  << period << endl;
  cout << "Acceptance from:  " << physType << endl;
  cout << "Number of bins considered: " << nBins << endl;
  cout << "Acceptance correction?   " << accCorrected << endl;
  cout << "File: " << fNameout << "  was writen:  " << toWrite << endl;
  cout << " " << endl;
}
