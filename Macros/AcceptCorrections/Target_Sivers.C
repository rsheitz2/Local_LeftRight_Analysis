const Int_t nBins=1; Double_t dx =0.02; Double_t yMax =0.4; 
//const Int_t nBins=3; Double_t dx =0.005; Double_t yMax =0.7;

Bool_t aSivers=true;
TString physType ="xF", period ="WAll";


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


Double_t Amp(Double_t a_sum, Double_t a_diff, Double_t R, Double_t P){
  Double_t A = (R*P*a_sum - a_diff)*TMath::Pi();
  A /= (2*P*(a_sum - R*P*a_diff) );

  //Double_t A = R*TMath::Pi()/2.0; //Simple
  return A;
}


Double_t e_Amp(Double_t aL, Double_t aR, Double_t R, Double_t P,
	       Double_t dR){
  Double_t a_sum = aL + aR;
  Double_t a_diff = aL - aR;
  
  Double_t e = 2*TMath::Pi()*aL*aR*dR;
  e /= (a_sum - R*P*a_diff)*(a_sum - R*P*a_diff);

  return e;
}


void Target_Sivers(TString fname=""){
  if (fname==""){
    fname += "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Presents/May1/Macros/Accept";

    cout << "Using default data from May1/Macros/Accept" << endl;
  }

  TFile *f_LR = TFile::Open(Form("%s/%s_%i.root", fname.Data(), period.Data(),
				 nBins) );
  TFile *f_LR_noCorr = TFile::Open(Form("%s/%s_%i_noCorr.root", fname.Data(),
					period.Data(), nBins) );
  TFile *f_acc = TFile::Open(Form("%s/Acceptance_%s_%s_%i.root",
				  fname.Data(), physType.Data(), period.Data(),
				  nBins) );

  if ( !(f_LR) || !(f_LR_noCorr) || !(f_acc) ){//Basic file checks
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
    
  TGraphErrors *g_LR[nTarg], *g_LR_noCorr[nTarg], *g_acc[nTarg];;
  for (Int_t i=0; i<nTarg; i++) {
    g_LR[i] = (TGraphErrors*)f_LR->Get(Form("%s_%s", physType.Data(),
					    targName[i].Data() ) );
    g_LR_noCorr[i] = (TGraphErrors*)f_LR_noCorr->Get(Form("%s_%s",
							  physType.Data(),
							  targName[i].Data()) );
    g_acc[i] = (TGraphErrors*)f_acc->Get(accName[i]);
  }
  
  Double_t Asivers[nTarg][nBins], e_Asivers[nTarg][nBins];
  Double_t xvals[nBins];
  Double_t ex[nBins]= {0.};

  Double_t *y_acc_UpS_L = g_acc[0]->GetY();
  Double_t *y_acc_UpS_R = g_acc[1]->GetY();
  Double_t *y_acc_DownS_L = g_acc[2]->GetY();
  Double_t *y_acc_DownS_R = g_acc[3]->GetY();
  Double_t *x_lr = g_LR[0]->GetX();
  for (Int_t tr=0; tr<nTarg; tr++) {
    Double_t *y_lr = g_LR[tr]->GetY();
    Double_t *y_lr_noCorr = g_LR_noCorr[tr]->GetY();

    Double_t *e_y_lr = g_LR[tr]->GetEY();

    for (Int_t bi=0; bi<nBins; bi++) {
      Double_t P = y_lr_noCorr[bi]/y_lr[bi];
      Double_t R = y_lr[bi];

      Double_t aL, aR;
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
      Double_t a_diff =aL - aR;
      Double_t a_sum =aL + aR;

      //Check output contributions
      /*cout << R*P << " " << R*P*R*P*a_diff/a_sum << " " << a_diff/a_sum <<
      " " << R*P*(a_diff/a_sum)*(a_diff/a_sum) << endl;*/

      Double_t dR = e_y_lr[bi];
      if (aSivers){
	Asivers[tr][bi] = Amp(a_sum, a_diff, R, P);
	e_Asivers[tr][bi] = e_Amp(aL, aR, R, P, dR);
      }
      else{
	Asivers[tr][bi] = R;
	e_Asivers[tr][bi] = dR;
      }
    
      xvals[bi] = x_lr[bi];
    }
  }


  TCanvas* c4 = new TCanvas(); c4->Divide(2,2);
  TCanvas* c1 = new TCanvas();
  TGraphErrors *g_Asivers[nTarg];
  TLine *li;
  for (Int_t i=0; i<nTarg; i++) {
    g_Asivers[i] = new TGraphErrors(nBins, xvals, Asivers[i], ex, e_Asivers[i]);
    

    if (i==0) {
      Double_t xmin =g_Asivers[i]->GetXaxis()->GetXmin();
      Double_t xmax =g_Asivers[i]->GetXaxis()->GetXmax();
      li = new TLine(xmin, 0., xmax, 0.);
      li->SetLineColor(kBlue);
      li->SetLineStyle(8);
    }

    
    SetUpTGraph(g_Asivers[i], targName[i], icolor[i], offsets[i] );
    c4->cd(i+1);
    g_Asivers[i]->Draw("AP");
    li->Draw("same");

    c1->cd();
    (i==0) ? g_Asivers[i]->Draw("AP") : g_Asivers[i]->Draw("Psame");
    li->Draw("same");
  }


  cout << " " << endl;
  cout << "Period is:    "  << period << endl;
  cout << "Acceptance from:  " << physType << endl;
  cout << "Number of bins considered: " << nBins << endl;
  cout << "Output amplitude is Asivers  " << aSivers << endl;
  cout << "Output amplitude is AN  " << !aSivers << endl;
  cout << " " << endl;
}
