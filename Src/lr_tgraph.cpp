#include "lr_tgraph.h"

lr_tgraph::lr_tgraph() : lrSpinCorr(){
  gr_asym = NULL;
  gr_asym_upstream = NULL;
  gr_asym_downstream = NULL;
  gr_asym_upstream_up = NULL;
  gr_asym_upstream_down = NULL;
  gr_asym_downstream_up = NULL;
  gr_asym_downstream_down = NULL;
  
  li_asym = NULL;
  li_asym_upstream = NULL;
  li_asym_downstream = NULL;
  li_asym_upstream_up = NULL;
  li_asym_upstream_down = NULL;
  li_asym_downstream_up = NULL;
  li_asym_downstream_down = NULL;
}


lr_tgraph::lr_tgraph(int nBins, TString name)
  : lrSpinCorr(nBins, name), genericBounds(nBins, name) {
  gr_asym = NULL;
  gr_asym_upstream = NULL;
  gr_asym_downstream = NULL;
  gr_asym_upstream_up = NULL;
  gr_asym_upstream_down = NULL;
  gr_asym_downstream_up = NULL;
  gr_asym_downstream_down = NULL;
  
  li_asym = NULL;
  li_asym_upstream = NULL;
  li_asym_downstream = NULL;
  li_asym_upstream_up = NULL;
  li_asym_upstream_down = NULL;
  li_asym_downstream_up = NULL;
  li_asym_downstream_down = NULL;

   for (Int_t i=0; i<nBins; i++) {
     ex.push_back(0.0);
   }
}


lr_tgraph::lr_tgraph (int nBins, std::vector<Double_t> &in_bounds,
	      std::vector<Double_t> &in_xval, TString thisName)
    : lrSpinCorr(nBins, in_bounds, in_xval, thisName){
  gr_asym = NULL;
  gr_asym_upstream = NULL;
  gr_asym_downstream = NULL;
  gr_asym_upstream_up = NULL;
  gr_asym_upstream_down = NULL;
  gr_asym_downstream_up = NULL;
  gr_asym_downstream_down = NULL;
  
  li_asym = NULL;
  li_asym_upstream = NULL;
  li_asym_downstream = NULL;
  li_asym_upstream_up = NULL;
  li_asym_upstream_down = NULL;
  li_asym_downstream_up = NULL;
  li_asym_downstream_down = NULL;

  for (Int_t i=0; i<nBins; i++) {
     ex.push_back(0.0);
  }
}


lr_tgraph::lr_tgraph (TString binfile, TString type, TString thisName)
  : lrSpinCorr(binfile, type, thisName){

  gr_asym = NULL;
  gr_asym_upstream = NULL;
  gr_asym_downstream = NULL;
  gr_asym_upstream_up = NULL;
  gr_asym_upstream_down = NULL;
  gr_asym_downstream_up = NULL;
  gr_asym_downstream_down = NULL;
  
  li_asym = NULL;
  li_asym_upstream = NULL;
  li_asym_downstream = NULL;
  li_asym_upstream_up = NULL;
  li_asym_upstream_down = NULL;
  li_asym_downstream_up = NULL;
  li_asym_downstream_down = NULL;

  for (Int_t i=0; i<nBins; i++) {
     ex.push_back(0.0);
  }
}


lr_tgraph::lr_tgraph (TFile *f1, TString type, TString thisName)
  : lrSpinCorr(f1, type, thisName){
  gr_asym = NULL;
  gr_asym_upstream = NULL;
  gr_asym_downstream = NULL;
  gr_asym_upstream_up = NULL;
  gr_asym_upstream_down = NULL;
  gr_asym_downstream_up = NULL;
  gr_asym_downstream_down = NULL;
  
  li_asym = NULL;
  li_asym_upstream = NULL;
  li_asym_downstream = NULL;
  li_asym_upstream_up = NULL;
  li_asym_upstream_down = NULL;
  li_asym_downstream_up = NULL;
  li_asym_downstream_down = NULL;

  for (Int_t i=0; i<nBins; i++) {
     ex.push_back(0.0);
  }
}


void lr_tgraph::Fill(){
  this->gr_asym = new TGraphErrors(this->nBins, &(this->xval[0]),
				   &(this->asym[0]),
				   &(this->ex[0]), &(this->e_asym[0]) );
  this->gr_asym_upstream = new TGraphErrors(this->nBins, &(this->xval[0]),
					    &(this->asym_upstream[0]),
					    &(this->ex[0]),
					    &(this->e_asym_upstream[0]) );
  this->gr_asym_downstream = new TGraphErrors(this->nBins, &(this->xval[0]),
					    &(this->asym_downstream[0]),
					    &(this->ex[0]),
					    &(this->e_asym_downstream[0]) );

  this->gr_asym_upstream_up = new TGraphErrors(this->nBins, &(this->xval[0]),
					       &(this->asym_upstream_up[0]),
					       &(this->ex[0]),
					       &(this->e_asym_upstream_up[0]) );
  this->gr_asym_upstream_down = new TGraphErrors(this->nBins, &(this->xval[0]),
						 &(this->asym_upstream_down[0]),
						 &(this->ex[0]),
						 &(this->e_asym_upstream_down[0]));
  this->gr_asym_downstream_up = new TGraphErrors(this->nBins, &(this->xval[0]),
						 &(this->asym_downstream_up[0]),
						 &(this->ex[0]),
						 &(this->e_asym_downstream_up[0]) );
  this->gr_asym_downstream_down = new TGraphErrors(this->nBins, &(this->xval[0]),
						   &(this->asym_downstream_down[0]),
						   &(this->ex[0]),
						   &(this->e_asym_downstream_down[0]));
  
  double xmin = this->gr_asym->GetXaxis()->GetXmin();
  double xmax = this->gr_asym->GetXaxis()->GetXmax();
  this->li_asym = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_upstream->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream->GetXaxis()->GetXmax();
  this->li_asym_upstream = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream->GetXaxis()->GetXmax();
  this->li_asym_downstream = new TLine(xmin, 0.0, xmax, 0.0);

  xmin = this->gr_asym_upstream_up->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream_up->GetXaxis()->GetXmax();
  this->li_asym_upstream_up = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_upstream_down->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream_down->GetXaxis()->GetXmax();
  this->li_asym_upstream_down = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream_up->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream_up->GetXaxis()->GetXmax();
  this->li_asym_downstream_up = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream_down->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream_down->GetXaxis()->GetXmax();
  this->li_asym_downstream_down = new TLine(xmin, 0.0, xmax, 0.0);
  
  setupTGraph(this->gr_asym); setupTGraph(this->gr_asym_upstream);
  setupTGraph(this->gr_asym_downstream);
  setupTGraph(this->gr_asym_upstream_up);
  setupTGraph(this->gr_asym_upstream_down);
  setupTGraph(this->gr_asym_downstream_up);
  setupTGraph(this->gr_asym_downstream_down);

  setupTLine(this->li_asym); setupTLine(this->li_asym_upstream);
  setupTLine(this->li_asym_downstream);
  setupTLine(this->li_asym_upstream_up);
  setupTLine(this->li_asym_upstream_down);
  setupTLine(this->li_asym_downstream_up);
  setupTLine(this->li_asym_downstream_down);
}


void lr_tgraph::Fill(TString toFill){
  if (toFill!="all") {
    std::cout << "Wrong option: " << toFill
	      << " not supported in lr_tgraph::Fill" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  this->gr_asym = new TGraphErrors(this->nBins, &(this->upstream_xval[0]),
				   &(this->asym[0]),
				   &(this->ex[0]), &(this->e_asym[0]) );
  this->gr_asym_upstream = new TGraphErrors(this->nBins, &(this->upstream_xval[0]),
					    &(this->asym_upstream[0]),
					    &(this->ex[0]),
					    &(this->e_asym_upstream[0]) );
  this->gr_asym_downstream = new TGraphErrors(this->nBins, &(this->downstream_xval[0]),
					    &(this->asym_downstream[0]),
					    &(this->ex[0]),
					    &(this->e_asym_downstream[0]) );

  this->gr_asym_upstream_up = new TGraphErrors(this->nBins, &(this->upstream_xval[0]),
					       &(this->asym_upstream_up[0]),
					       &(this->ex[0]),
					       &(this->e_asym_upstream_up[0]) );
  this->gr_asym_upstream_down = new TGraphErrors(this->nBins, &(this->upstream_xval[0]),
						 &(this->asym_upstream_down[0]),
						 &(this->ex[0]),
						 &(this->e_asym_upstream_down[0]));
  this->gr_asym_downstream_up = new TGraphErrors(this->nBins, &(this->downstream_xval[0]),
						 &(this->asym_downstream_up[0]),
						 &(this->ex[0]),
						 &(this->e_asym_downstream_up[0]) );
  this->gr_asym_downstream_down = new TGraphErrors(this->nBins, &(this->downstream_xval[0]),
						   &(this->asym_downstream_down[0]),
						   &(this->ex[0]),
						   &(this->e_asym_downstream_down[0]));
  
  double xmin = this->gr_asym->GetXaxis()->GetXmin();
  double xmax = this->gr_asym->GetXaxis()->GetXmax();
  this->li_asym = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_upstream->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream->GetXaxis()->GetXmax();
  this->li_asym_upstream = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream->GetXaxis()->GetXmax();
  this->li_asym_downstream = new TLine(xmin, 0.0, xmax, 0.0);

  xmin = this->gr_asym_upstream_up->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream_up->GetXaxis()->GetXmax();
  this->li_asym_upstream_up = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_upstream_down->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream_down->GetXaxis()->GetXmax();
  this->li_asym_upstream_down = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream_up->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream_up->GetXaxis()->GetXmax();
  this->li_asym_downstream_up = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream_down->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream_down->GetXaxis()->GetXmax();
  this->li_asym_downstream_down = new TLine(xmin, 0.0, xmax, 0.0);
  
  setupTGraph(this->gr_asym); setupTGraph(this->gr_asym_upstream);
  setupTGraph(this->gr_asym_downstream);
  setupTGraph(this->gr_asym_upstream_up);
  setupTGraph(this->gr_asym_upstream_down);
  setupTGraph(this->gr_asym_downstream_up);
  setupTGraph(this->gr_asym_downstream_down);

  setupTLine(this->li_asym); setupTLine(this->li_asym_upstream);
  setupTLine(this->li_asym_downstream);
  setupTLine(this->li_asym_upstream_up);
  setupTLine(this->li_asym_upstream_down);
  setupTLine(this->li_asym_downstream_up);
  setupTLine(this->li_asym_downstream_down);
}


void lr_tgraph::Draw(TString name, TString opt){

  if (name == "asym"){
    this->gr_asym->Draw(opt);
    this->li_asym->Draw("same");
  }
  else if (name == "asym_upstream")
    this->gr_asym_upstream->Draw(opt);
    this->li_asym_upstream->Draw("same");
}


void lr_tgraph::Write(){
  this->gr_asym->Write(Form("%s_asym", this->thisName.Data() ) );
  this->gr_asym_upstream->Write(Form("%s_asym_upstream",
				     this->thisName.Data() ) );
  this->gr_asym_downstream->Write(Form("%s_asym_downstream",
				       this->thisName.Data() ) );

  this->gr_asym_upstream_up->Write(Form("%s_asym_upstream_up",
					this->thisName.Data() ) );
  this->gr_asym_upstream_down->Write(Form("%s_asym_upstream_down",
					  this->thisName.Data() ) );
  this->gr_asym_downstream_up->Write(Form("%s_asym_downstream_up",
					  this->thisName.Data() ) );
  this->gr_asym_downstream_down->Write(Form("%s_asym_downstream_down",
					    this->thisName.Data() ) );
}


//Helper functions
void lr_tgraph::setupTLine(TLine *l){
  l->SetLineColor(kBlue);
  l->SetLineStyle(8);
}


void lr_tgraph::setupTGraph(TGraphErrors* gr){
  gr->SetMarkerStyle(21);
}
