#include "lr_tgraph.h"

lr_tgraph::lr_tgraph(int nBins, TString name)
  : lrSpinCorr(nBins, name), genericBounds(nBins, name) {

  lr_tgraph::setZero(nBins);
}


lr_tgraph::lr_tgraph (int nBins, std::vector<Double_t> &in_bounds,
		      std::vector<Double_t> &in_xval, TString thisName)
  : lrSpinCorr(nBins, in_bounds, in_xval, thisName){

  lr_tgraph::setZero(nBins);
}


lr_tgraph::lr_tgraph (TString binfile, TString type, TString thisName)
  : lrSpinCorr(binfile, type, thisName){

  lr_tgraph::setZero(nBins);
}


lr_tgraph::lr_tgraph (TFile *f1, TString type, TString thisName)
  : lrSpinCorr(f1, type, thisName){

  lr_tgraph::setZero(nBins);
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

  this->gr_asym_upstream_left
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_upstream_left[0]), &(this->ex[0]),
		       &(this->e_asym_upstream_left[0]));
  this->gr_asym_upstream_right
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_upstream_right[0]), &(this->ex[0]),
		       &(this->e_asym_upstream_right[0]));
  this->gr_asym_downstream_left
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_downstream_left[0]), &(this->ex[0]),
		       &(this->e_asym_downstream_left[0]));
  this->gr_asym_downstream_right
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_downstream_right[0]), &(this->ex[0]),
		       &(this->e_asym_downstream_right[0]));
  this->gr_asym_updown_left
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_updown_left[0]), &(this->ex[0]),
		       &(this->e_asym_updown_left[0]));
  this->gr_asym_updown_right
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_updown_right[0]), &(this->ex[0]),
		       &(this->e_asym_updown_right[0]));
  this->gr_asym_downup_left
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_downup_left[0]), &(this->ex[0]),
		       &(this->e_asym_downup_left[0]));
  this->gr_asym_downup_right
    = new TGraphErrors(this->nBins, &(this->xval[0]),
		       &(this->asym_downup_right[0]), &(this->ex[0]),
		       &(this->e_asym_downup_right[0]));
  
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

  xmin = this->gr_asym_upstream_left->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream_left->GetXaxis()->GetXmax();
  this->li_asym_upstream_left = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_upstream_right->GetXaxis()->GetXmin();
  xmax = this->gr_asym_upstream_right->GetXaxis()->GetXmax();
  this->li_asym_upstream_right = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream_left->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream_left->GetXaxis()->GetXmax();
  this->li_asym_downstream_left = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downstream_right->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downstream_right->GetXaxis()->GetXmax();
  this->li_asym_downstream_right = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_updown_left->GetXaxis()->GetXmin();
  xmax = this->gr_asym_updown_left->GetXaxis()->GetXmax();
  this->li_asym_updown_left = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_updown_right->GetXaxis()->GetXmin();
  xmax = this->gr_asym_updown_right->GetXaxis()->GetXmax();
  this->li_asym_updown_right = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downup_left->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downup_left->GetXaxis()->GetXmax();
  this->li_asym_downup_left = new TLine(xmin, 0.0, xmax, 0.0);
  xmin = this->gr_asym_downup_right->GetXaxis()->GetXmin();
  xmax = this->gr_asym_downup_right->GetXaxis()->GetXmax();
  this->li_asym_downup_right = new TLine(xmin, 0.0, xmax, 0.0);

  this->gr_Pol =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol[0]) );
  this->gr_Dil =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil[0]) );

  this->gr_Pol_upstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_upstream[0]) );
  this->gr_Dil_upstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_upstream[0]) );

  this->gr_Pol_upstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_upstream_up[0]) );
  this->gr_Dil_upstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_upstream_up[0]) );

  this->gr_Pol_upstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_upstream_down[0]));
  this->gr_Dil_upstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_upstream_down[0]));

    this->gr_Pol_downstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_downstream[0]) );
  this->gr_Dil_downstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_downstream[0]) );

  this->gr_Pol_downstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_downstream_up[0]));
  this->gr_Dil_downstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_downstream_up[0]));

  this->gr_Pol_downstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_downstream_down[0]));
  this->gr_Dil_downstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_downstream_down[0]));
  
  setupTGraph(this->gr_asym); setupTGraph(this->gr_asym_upstream);
  setupTGraph(this->gr_asym_downstream);
  setupTGraph(this->gr_asym_upstream_up);
  setupTGraph(this->gr_asym_upstream_down);
  setupTGraph(this->gr_asym_downstream_up);
  setupTGraph(this->gr_asym_downstream_down);
  setupTGraph(this->gr_asym_upstream_left);
  setupTGraph(this->gr_asym_upstream_right);
  setupTGraph(this->gr_asym_downstream_left);
  setupTGraph(this->gr_asym_downstream_right);
  setupTGraph(this->gr_asym_updown_left);
  setupTGraph(this->gr_asym_updown_right);
  setupTGraph(this->gr_asym_downup_left);
  setupTGraph(this->gr_asym_downup_right);

  setupTLine(this->li_asym); setupTLine(this->li_asym_upstream);
  setupTLine(this->li_asym_downstream);
  setupTLine(this->li_asym_upstream_up);
  setupTLine(this->li_asym_upstream_down);
  setupTLine(this->li_asym_downstream_up);
  setupTLine(this->li_asym_downstream_down);
  setupTLine(this->li_asym_upstream_left);
  setupTLine(this->li_asym_upstream_right);
  setupTLine(this->li_asym_downstream_left);
  setupTLine(this->li_asym_downstream_right);
  setupTLine(this->li_asym_updown_left);
  setupTLine(this->li_asym_updown_right);
  setupTLine(this->li_asym_downup_left);
  setupTLine(this->li_asym_downup_right);

  setupTGraph(this->gr_Pol);
  setupTGraph(this->gr_Dil);

  setupTGraph(this->gr_Pol_upstream);
  setupTGraph(this->gr_Dil_upstream);

  setupTGraph(this->gr_Pol_upstream_up);
  setupTGraph(this->gr_Dil_upstream_up);

  setupTGraph(this->gr_Pol_upstream_down);
  setupTGraph(this->gr_Dil_upstream_down);

  setupTGraph(this->gr_Pol_downstream);
  setupTGraph(this->gr_Dil_downstream);

  setupTGraph(this->gr_Pol_downstream_up);
  setupTGraph(this->gr_Dil_downstream_up);

  setupTGraph(this->gr_Pol_downstream_down);
  setupTGraph(this->gr_Dil_downstream_down);
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

  this->gr_Pol =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol[0]) );
  this->gr_Dil =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil[0]) );

  this->gr_Pol_upstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_upstream[0]) );
  this->gr_Dil_upstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_upstream[0]) );

  this->gr_Pol_upstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_upstream_up[0]) );
  this->gr_Dil_upstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_upstream_up[0]) );

  this->gr_Pol_upstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_upstream_down[0]));
  this->gr_Dil_upstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_upstream_down[0]));

    this->gr_Pol_downstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_downstream[0]) );
  this->gr_Dil_downstream =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_downstream[0]) );

  this->gr_Pol_downstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_downstream_up[0]));
  this->gr_Dil_downstream_up =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_downstream_up[0]));

  this->gr_Pol_downstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgPol_downstream_down[0]));
  this->gr_Dil_downstream_down =
    new TGraph(this->nBins, &(this->xval[0]), &(this->avgDil_downstream_down[0]));
  
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

  setupTGraph(this->gr_Pol);
  setupTGraph(this->gr_Dil);

  setupTGraph(this->gr_Pol_upstream);
  setupTGraph(this->gr_Dil_upstream);

  setupTGraph(this->gr_Pol_upstream_up);
  setupTGraph(this->gr_Dil_upstream_up);

  setupTGraph(this->gr_Pol_upstream_down);
  setupTGraph(this->gr_Dil_upstream_down);

  setupTGraph(this->gr_Pol_downstream);
  setupTGraph(this->gr_Dil_downstream);

  setupTGraph(this->gr_Pol_downstream_up);
  setupTGraph(this->gr_Dil_downstream_up);

  setupTGraph(this->gr_Pol_downstream_down);
  setupTGraph(this->gr_Dil_downstream_down);
}


void lr_tgraph::Draw(TString name, TString opt){

  if (name == "asym"){
    this->gr_asym->Draw(opt);
    this->li_asym->Draw("same");
  }
  else if (name == "asym_upstream") {
    this->gr_asym_upstream->Draw(opt);
    this->li_asym_upstream->Draw("same");
  }
  //not completed (probably doesn't need to be)

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

  this->gr_asym_upstream_left->Write(Form("%s_asym_upstream_left",
					  this->thisName.Data() ) );
  this->gr_asym_upstream_right->Write(Form("%s_asym_upstream_right",
					   this->thisName.Data() ) );
  this->gr_asym_downstream_left->Write(Form("%s_asym_downstream_left",
					    this->thisName.Data() ) );
  this->gr_asym_downstream_right->Write(Form("%s_asym_downstream_right",
					     this->thisName.Data() ) );

  this->gr_asym_updown_left->Write(Form("%s_asym_updown_left",
					this->thisName.Data() ) );
  this->gr_asym_updown_right->Write(Form("%s_asym_updown_right",
					 this->thisName.Data() ) );
  this->gr_asym_downup_left->Write(Form("%s_asym_downup_left",
					this->thisName.Data() ) );
  this->gr_asym_downup_right->Write(Form("%s_asym_downup_right",
					 this->thisName.Data() ) );

  this->gr_Pol->Write(Form("%s_Pol", this->thisName.Data()));
  this->gr_Dil->Write(Form("%s_Dil", this->thisName.Data()));

  this->gr_Pol_upstream->Write(Form("%s_Pol_upstream", this->thisName.Data()));
  this->gr_Dil_upstream->Write(Form("%s_Dil_upstream", this->thisName.Data()));

  this->gr_Pol_upstream_up->Write(Form("%s_Pol_upstream_up",
				       this->thisName.Data()));
  this->gr_Dil_upstream_up->Write(Form("%s_Dil_upstream_up",
				       this->thisName.Data()));

  this->gr_Pol_upstream_down->Write(Form("%s_Pol_upstream_down",
				       this->thisName.Data()));
  this->gr_Dil_upstream_down->Write(Form("%s_Dil_upstream_down",
				       this->thisName.Data()));

  this->gr_Pol_downstream->Write(Form("%s_Pol_downstream",
				      this->thisName.Data()));
  this->gr_Dil_downstream->Write(Form("%s_Dil_downstream",
				      this->thisName.Data()));

  this->gr_Pol_downstream_up->Write(Form("%s_Pol_downstream_up",
				       this->thisName.Data()));
  this->gr_Dil_downstream_up->Write(Form("%s_Dil_downstream_up",
				       this->thisName.Data()));

  this->gr_Pol_downstream_down->Write(Form("%s_Pol_downstream_down",
				       this->thisName.Data()));
  this->gr_Dil_downstream_down->Write(Form("%s_Dil_downstream_down",
				       this->thisName.Data()));

  //Counts
  std::vector<Double_t> left_upS
    (this->left_upstream.begin(), this->left_upstream.end());
  TGraph *g_left_upstream
    = new TGraph(this->nBins, &(this->xval[0]), &(left_upS[0]));

  std::vector<Double_t> left_upS_up
    (this->left_upstream_up.begin(), this->left_upstream_up.end());
  TGraph *g_left_upstream_up
    = new TGraph(this->nBins, &(this->xval[0]), &(left_upS_up[0]));
  std::vector<Double_t> left_upS_down
    (this->left_upstream_down.begin(), this->left_upstream_down.end());
  TGraph *g_left_upstream_down
    = new TGraph(this->nBins, &(this->xval[0]), &(left_upS_down[0]));

  std::vector<Double_t> right_upS
    (this->right_upstream.begin(), this->right_upstream.end());
  TGraph *g_right_upstream
    = new TGraph(this->nBins, &(this->xval[0]), &(right_upS[0]));

  std::vector<Double_t> right_upS_up
    (this->right_upstream_up.begin(), this->right_upstream_up.end());
  TGraph *g_right_upstream_up
    = new TGraph(this->nBins, &(this->xval[0]), &(right_upS_up[0]));
  std::vector<Double_t> right_upS_down
    (this->right_upstream_down.begin(), this->right_upstream_down.end());
  TGraph *g_right_upstream_down
    = new TGraph(this->nBins, &(this->xval[0]), &(right_upS_down[0]));

  std::vector<Double_t> left_downS
    (this->left_downstream.begin(), this->left_downstream.end());
  TGraph *g_left_downstream
    = new TGraph(this->nBins, &(this->xval[0]), &(left_downS[0]));

  std::vector<Double_t> left_downS_up
    (this->left_downstream_up.begin(), this->left_downstream_up.end());
  TGraph *g_left_downstream_up
    = new TGraph(this->nBins, &(this->xval[0]), &(left_downS_up[0]));
  std::vector<Double_t> left_downS_down
    (this->left_downstream_down.begin(), this->left_downstream_down.end());
  TGraph *g_left_downstream_down
    = new TGraph(this->nBins, &(this->xval[0]), &(left_downS_down[0]));

  std::vector<Double_t> right_downS
    (this->right_downstream.begin(), this->right_downstream.end());
  TGraph *g_right_downstream
    = new TGraph(this->nBins, &(this->xval[0]), &(right_downS[0]));

  std::vector<Double_t> right_downS_up
    (this->right_downstream_up.begin(), this->right_downstream_up.end());
  TGraph *g_right_downstream_up
    = new TGraph(this->nBins, &(this->xval[0]), &(right_downS_up[0]));
  std::vector<Double_t> right_downS_down
    (this->right_downstream_down.begin(), this->right_downstream_down.end());
  TGraph *g_right_downstream_down
    = new TGraph(this->nBins, &(this->xval[0]), &(right_downS_down[0]));

  setupTGraph(g_left_upstream);
  setupTGraph(g_left_upstream_up); setupTGraph(g_left_upstream_down);
  
  setupTGraph(g_right_upstream);
  setupTGraph(g_right_upstream_up); setupTGraph(g_right_upstream_down);
  
  setupTGraph(g_left_downstream);
  setupTGraph(g_left_downstream_up); setupTGraph(g_left_downstream_down);

  setupTGraph(g_right_downstream);
  setupTGraph(g_right_downstream_up); setupTGraph(g_right_downstream_down);

  g_left_upstream->Write(Form("%s_left_upstream", this->thisName.Data()) );
  g_left_upstream_up->Write(Form("%s_left_upstream_up", this->thisName.Data()));
  g_left_upstream_down->Write(Form("%s_left_upstream_down",
				   this->thisName.Data()));

  g_right_upstream->Write(Form("%s_right_upstream", this->thisName.Data()) );
  g_right_upstream_up->Write(Form("%s_right_upstream_up",
				  this->thisName.Data()));
  g_right_upstream_down->Write(Form("%s_right_upstream_down",
				    this->thisName.Data()));

  g_left_downstream->Write(Form("%s_left_downstream", this->thisName.Data()) );
  g_left_downstream_up->Write(Form("%s_left_downstream_up",
				   this->thisName.Data()));
  g_left_downstream_down->Write(Form("%s_left_downstream_down",
				     this->thisName.Data()));

  g_right_downstream->Write(Form("%s_right_downstream", this->thisName.Data()));
  g_right_downstream_up->Write(Form("%s_right_downstream_up",
				    this->thisName.Data()));
  g_right_downstream_down->Write(Form("%s_right_downstream_down",
				      this->thisName.Data()));
}


//Helper functions
void lr_tgraph::setupTLine(TLine *l){
  l->SetLineColor(kBlue);
  l->SetLineStyle(8);
}


void lr_tgraph::setupTGraph(TGraphErrors* gr){
  gr->SetMarkerStyle(21);
}

void lr_tgraph::setupTGraph(TGraph* gr){
  gr->SetMarkerStyle(21);
}


void lr_tgraph::setZero(){
  
  this->gr_asym = NULL;
  this->gr_asym_upstream = NULL;
  this->gr_asym_downstream = NULL;
  this->gr_asym_upstream_up = NULL;
  this->gr_asym_upstream_down = NULL;
  this->gr_asym_downstream_up = NULL;
  this->gr_asym_downstream_down = NULL;
  this->gr_asym_upstream_left = NULL;
  this->gr_asym_upstream_right = NULL;
  this->gr_asym_downstream_left = NULL;
  this->gr_asym_downstream_right = NULL;
  this->gr_asym_updown_left = NULL;
  this->gr_asym_updown_right = NULL;
  this->gr_asym_downup_left = NULL;
  this->gr_asym_downup_right = NULL;
  
  this->li_asym = NULL;
  this->li_asym_upstream = NULL;
  this->li_asym_downstream = NULL;
  this->li_asym_upstream_up = NULL;
  this->li_asym_upstream_down = NULL;
  this->li_asym_downstream_up = NULL;
  this->li_asym_downstream_down = NULL;
  this->li_asym_upstream_left = NULL;
  this->li_asym_upstream_right = NULL;
  this->li_asym_downstream_left = NULL;
  this->li_asym_downstream_right = NULL;
  this->li_asym_updown_left = NULL;
  this->li_asym_updown_right = NULL;
  this->li_asym_downup_left = NULL;
  this->li_asym_downup_right = NULL;

  this->gr_Pol = NULL;
  this->gr_Dil = NULL;

  this->gr_Pol_upstream = NULL;
  this->gr_Dil_upstream = NULL;

  this->gr_Pol_upstream_up = NULL;
  this->gr_Dil_upstream_up = NULL;

  this->gr_Pol_upstream_down = NULL;
  this->gr_Dil_upstream_down = NULL;

  this->gr_Pol_downstream = NULL;
  this->gr_Dil_downstream = NULL;

  this->gr_Pol_downstream_up = NULL;
  this->gr_Dil_downstream_up = NULL;

  this->gr_Pol_downstream_down = NULL;
  this->gr_Dil_downstream_down = NULL;
}


void lr_tgraph::setZero(Int_t nBins){
  lr_tgraph::setZero();

  for (Int_t i=0; i<nBins; i++) this->ex.push_back(0.0);
}
