void PullDist(){
  const Int_t nPeriod = 8;
    
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/May1/Data/leftRight_byTarget/";
  TFile* f_per[nPeriod];
  TGraphErrors *g_xN[nPeriod], *g_xPi[nPeriod], *g_xF[nPeriod], *g_pT[nPeriod];
  TGraphErrors *g_M[nPeriod];
  
  vector<Double_t> lr, e_lr;
  for (Int_t p=0; p<nPeriod; p++) {
    f_per[p] = (p<3) ? TFile::Open(Form("%sHM_W0%i_3bins.root",path.Data(),p+7))
      : TFile::Open(Form("%sHM_W%i_3bins.root",path.Data(),p+7));
    
    g_xN[p] = (TGraphErrors*)f_per[p]->Get("xN_asym");
    g_xPi[p] = (TGraphErrors*)f_per[p]->Get("xPi_asym");
    g_xF[p] = (TGraphErrors*)f_per[p]->Get("xF_asym");
    g_pT[p] = (TGraphErrors*)f_per[p]->Get("pT_asym");
    g_M[p] = (TGraphErrors*)f_per[p]->Get("M_asym");

    Double_t *yVal_xN=g_xN[p]->GetY(), *e_yVal_xN=g_xN[p]->GetEY();
    Double_t *yVal_xPi=g_xPi[p]->GetY(), *e_yVal_xPi=g_xPi[p]->GetEY();
    Double_t *yVal_xF=g_xF[p]->GetY(), *e_yVal_xF=g_xF[p]->GetEY();
    Double_t *yVal_pT=g_pT[p]->GetY(), *e_yVal_pT=g_pT[p]->GetEY();
    Double_t *yVal_M=g_M[p]->GetY(), *e_yVal_M=g_M[p]->GetEY();
    for (Int_t i=0; i<3; i++) {
      lr.push_back(yVal_xN[i]); e_lr.push_back(e_yVal_xN[i]);
      lr.push_back(yVal_xPi[i]); e_lr.push_back(e_yVal_xPi[i]);
      lr.push_back(yVal_xF[i]); e_lr.push_back(e_yVal_xF[i]);
      lr.push_back(yVal_pT[i]); e_lr.push_back(e_yVal_pT[i]);
      lr.push_back(yVal_M[i]); e_lr.push_back(e_yVal_M[i]);
    }

  }//period loop


  Double_t avg=0., sigma=0.;
  for(vector<Double_t>::iterator it=lr.begin(); it!=lr.end();
      it++){
    avg += *it;
    sigma += (*it)*(*it);
  }
  avg = avg/(1.0*lr.size() ); sigma = sigma/(1.0*lr.size() );
  sigma = sigma - avg*avg;

  //cout << avg << " " << sigma << endl;//cleanup

  TH1D* hPull = new TH1D("hPull", "hPull", 20, -4.0, 4.0);
  for(vector<Double_t>::iterator it=lr.begin(), e_it=e_lr.begin(); it!=lr.end();
      it++, e_it++){
    Double_t pull = *it - avg;
    Double_t denom = ((*e_it)*(*e_it)>sigma) ? TMath::Sqrt( (*e_it)*(*e_it) - sigma ) :
      TMath::Sqrt( sigma - (*e_it)*(*e_it) );

    //cout << *it << " " << pull << " " << denom << endl;//cleanup
    
    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    hPull->Fill(pull);
  }

  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111111);
  TCanvas* c1 = new TCanvas();
  hPull->SetMarkerStyle(20);
  hPull->SetMarkerColor(kBlue);
  hPull->SetLineColor(kBlue);
  hPull->Draw("EP");
}
