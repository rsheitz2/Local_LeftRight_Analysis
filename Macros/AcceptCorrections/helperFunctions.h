//Errors
Double_t RatioError(long long A, long long B){
  Double_t ratio = 1.0*A/(1.0*B);
  Double_t error = TMath::Sqrt( 1.0/(1.0*A) + 1.0/(1.0*B) );
  
  return ratio*error;
}


Double_t RatioError(Double_t A, Double_t B,
		    Double_t dA, Double_t dB){
  Double_t ratio = A/B;
  Double_t error = TMath::Sqrt( dA*dA/(A*A) + dB*dB/(B*B) );

  return ratio*error;
}


//Aesthetics
void SetUpTGraph(TGraphErrors* g, TString name, Int_t ic, Double_t offset,
		 Int_t nBins){
  g->SetTitle(name);
  
  g->SetMarkerStyle(21);
  g->SetMarkerColor(ic);

  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  if (nBins!=0){
    Double_t *xval = g->GetX();
    for (Int_t i=0; i<nBins; i++) xval[i] += offset;
  }
}


void SetUpTGraph(TGraphErrors* g, Double_t offset, Int_t nBins){
  g->SetMarkerStyle(21);
  g->SetMarkerSize(1.4);
  
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);

  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  if (nBins!=0){
    Double_t *xvals = g->GetX();
    for (Int_t i=0; i<nBins; i++) xvals[i] += offset;
  }
}


void SetUpTGraph(TGraphErrors* g){
  g->SetMarkerStyle(21);

  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
 
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);
}


void SetUpTH1(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);

  h->GetXaxis()->SetNdivisions(504);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);

  h->SetLineWidth(3);
}
