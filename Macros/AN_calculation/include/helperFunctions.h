//Conventions
enum TargPol {upS_up, upS_down, downS_up, downS_down};

//Errors Calculations
Double_t FullyCorrelatedRatioError(Double_t A, Double_t B,
			      Double_t eA, Double_t eB){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t r = A/B;
  if (r<0) r *= -1.0;
 
  Double_t e = eA*eA/(A*A) + eB*eB/(B*B) - 2*eA*eB/(A*B);
  e = TMath::Sqrt( e );
  e *= r;

  return e;
}


Double_t CorrelatedRatioError(Double_t A, Double_t B, Double_t eA,
				   Double_t eB, Double_t cov){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t r = A/B;
  if (r<0) r *= -1.0;
 
  Double_t e = eA*eA/(A*A) + eB*eB/(B*B) - 2*cov/(A*B);
  e = TMath::Sqrt( e );
  e *= r;

  return e;
}


Double_t FullyCorrelatedDiffError(Double_t eA, Double_t eB){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }

  Double_t e = eA*eA + eB*eB - 2*eA*eB;
  e = TMath::Sqrt(e);

  return e;
}


Double_t DiffError(Double_t eA, Double_t eB){
  Double_t e2A =eA*eA, e2B =eB*eB;

  return TMath::Sqrt(e2A + e2B);
}


Double_t RatioError(Double_t A, Double_t B,
		    Double_t eA, Double_t eB){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t r = A/B;
  if (r<0) r *= -1.0;
 
  Double_t e = eA*eA/(A*A) + eB*eB/(B*B);
  e = TMath::Sqrt( e );
  e *= r;

  return e;
}


//Additional Calculations
Double_t WeightedAvg(vector<Double_t> &A, vector<Double_t> &eA){

  Double_t avg=0.0, e=0.0;
  for (vector<Double_t>::iterator iA=A.begin(), iE=eA.begin(); iA!=A.end();
       iA++, iE++) {
    avg += *iA/( (*iE) * (*iE) );
    e += 1.0/( (*iE) * (*iE) );
  }
  
  avg /= e;

  return avg;
}


Double_t WeightedAvg(Double_t A, Double_t B,
		     Double_t eA, Double_t eB){
  Double_t e = 1/(eA*eA) + 1/(eB*eB);
  
  Double_t avg = A/(eA*eA) + B/(eB*eB);
  avg /= e;

  return avg;
}


Double_t WeightedErr(vector<Double_t> &eA){

  Double_t e=0.0;
  for (vector<Double_t>::iterator it=eA.begin(); it!=eA.end(); it++) {
    e += 1.0/( (*it) * (*it) );
  }

  e = 1/e;
  e = TMath::Sqrt(e);
  
  return e;
}


Double_t WeightedErr(Double_t eA, Double_t eB){
  Double_t e = 1/(eA*eA) + 1/(eB*eB);
  e = 1/e;
  e = TMath::Sqrt(e);
  
  return e;
}


Double_t GaussInt(Double_t A, Double_t sigma){
  
  return TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
}


Double_t GaussIntError(Double_t A, Double_t sigma,
		       Double_t eA, Double_t eSigma){
  Double_t e = eA*eA/( A*A ) + eSigma*eSigma/( sigma*sigma );
  e = TMath::Sqrt( e );
  e *= TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
  
  return e;
}


Double_t CorrelatedGaussIntError(Double_t A, Double_t sigma,
				 Double_t eA, Double_t eSigma, Double_t cov){
  Double_t e = eA*eA/( A*A ) + eSigma*eSigma/( sigma*sigma ) + 2*cov/(A*sigma);
  e = TMath::Sqrt( e );
  e *= TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
  
  return e;
}


Double_t ExpoInt(Double_t A, Double_t b){
  if (b > 0){
    cout << "Cannot do exponential integral" << endl;
    exit(EXIT_FAILURE);
  }
  
  return -1.0*A/b;
}


//Aesthetics
void SetUpTGraph(TGraphErrors* g, Int_t icolor=1,
		 Double_t offset=0.0, Int_t nBins=0){
  g->SetMarkerStyle(21);
  g->SetMarkerColor(icolor);

  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  if (nBins != 0){
    Double_t *xval = g->GetX();
    for (Int_t i=0; i<nBins; i++) xval[i] += offset;
  }
}


void SetUpTGraph(TGraph* g){
  g->SetMarkerStyle(21);
  
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);
}


void SetUpHist(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);
  
  h->GetXaxis()->SetNdivisions(508);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);
}


void SetUpHist(TH1D* h, Double_t xmin, Double_t xmax){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);
  
  h->GetXaxis()->SetNdivisions(508);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);

  h->GetXaxis()->SetRangeUser(xmin, xmax);
}


void SetUpTF(TF1* f){
  f->GetYaxis()->SetNdivisions(504);
  f->GetYaxis()->SetLabelFont(22);
  f->GetYaxis()->SetLabelSize(0.08);
  
  f->GetXaxis()->SetNdivisions(508);
  f->GetXaxis()->SetLabelFont(22);
  f->GetXaxis()->SetLabelSize(0.08);
}


void DrawLine(TGraphErrors *g, Double_t yval){
  Double_t min_x = g->GetXaxis()->GetXmin();	
  Double_t max_x = g->GetXaxis()->GetXmax();	
  TLine* li = new TLine(min_x, yval, max_x, yval);
  
  li->SetLineColor(1);
  li->SetLineStyle(8);
  li->SetLineWidth(2);
  
  li->Draw("same");
}


void DrawLine(TH1D *h, Double_t yval){
  Double_t min_x = h->GetXaxis()->GetXmin();	
  Double_t max_x = h->GetXaxis()->GetXmax();	
  TLine* li = new TLine(min_x, yval, max_x, yval);
  
  li->SetLineColor(1);
  li->SetLineStyle(8);
  li->SetLineWidth(2);
  
  li->Draw("same");
}

//Usefunctions
void GetPolarization(Double_t *vals_noPcorr, Double_t *vals, Double_t *Pol,
		     Int_t nBins){
  Double_t Pol_mean; Int_t Pol_nan_counts=0;

  for (Int_t bi=0; bi<nBins; bi++) {
    Pol[bi] = vals_noPcorr[bi]/vals[bi];
    if ( isnan(Pol[bi]) ) Pol_nan_counts++;
    else Pol_mean += Pol[bi]; 
  }

  if (Pol_nan_counts == nBins){
    cout << "Error: No polarization values determined" << endl;
    exit(EXIT_FAILURE);
  }
    
  Pol_mean = Pol_mean/(nBins-Pol_nan_counts);
  for (Int_t bi=0; bi<nBins; bi++) {
    if ( isnan(Pol[bi]) ) Pol[bi] = Pol_mean;
  }
}

