#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

//Calculations
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


Double_t WeightedAvg(TGraphErrors *g){
  Double_t *yvals = g->GetY();
  Double_t *e_yvals = g->GetEY();
  
  vector<Double_t> vec_vals;
  vector<Double_t> e_vec_vals;
  
  for (Int_t i=0; i<g->GetN(); i++) {
    vec_vals.push_back(yvals[i]);
    e_vec_vals.push_back(e_yvals[i]);
  }

  return WeightedAvg(vec_vals, e_vec_vals);
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


Double_t WeightedErr(TGraphErrors *g){
  Double_t *e_yvals = g->GetEY();
  
  vector<Double_t> e_vec_vals;
  for (Int_t i=0; i<g->GetN(); i++) e_vec_vals.push_back(e_yvals[i]);

  return WeightedErr(e_vec_vals);
}


Double_t WeightedErr(Double_t eA, Double_t eB){
  Double_t e = 1/(eA*eA) + 1/(eB*eB);
  e = 1/e;
  e = TMath::Sqrt(e);
  
  return e;
}


//Aesthetics
void SetUp(TGraphErrors *g, Int_t icolor, Int_t imarker, Double_t xShift){
  g->SetMarkerColor(icolor);
  g->SetMarkerStyle(imarker);

  Double_t *xvals = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xvals[i] += xShift;
}


void SetUp(TGraphErrors* g, Int_t icolor=1,
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


void SetUp(TGraph* g){
  g->SetMarkerStyle(21);
  
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);
}


void SetUp(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);
  
  h->GetXaxis()->SetNdivisions(504);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);
}


void SetUp(TF1* f){
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


TFile* OpenFile(TString name){
  TFile* f = TFile::Open(name);
  if (!f){
    cout << "File does not exist:  " << name << endl;
    exit(EXIT_FAILURE);
  }

  return f;
}

#endif
