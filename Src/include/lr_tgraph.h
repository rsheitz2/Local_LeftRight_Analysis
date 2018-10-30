#ifndef LR_TGRAPH_H
#define LR_TGRAPH_H

#include "lrSpinCorr.h"
#include "genericBounds.h"
#include "common.h"

class lr_tgraph: virtual public lrSpinCorr, public genericBounds {
private:
  TGraphErrors *gr_asym, *gr_asym_upstream, *gr_asym_downstream;
  TGraphErrors *gr_asym_upstream_up, *gr_asym_upstream_down;
  TGraphErrors *gr_asym_downstream_up, *gr_asym_downstream_down;
  TGraphErrors *gr_asym_upstream_left, *gr_asym_upstream_right;
  TGraphErrors *gr_asym_downstream_left, *gr_asym_downstream_right;
  TGraphErrors *gr_asym_updown_left, *gr_asym_updown_right;
  TGraphErrors *gr_asym_downup_left, *gr_asym_downup_right;
  
  TLine *li_asym, *li_asym_upstream, *li_asym_downstream;
  TLine *li_asym_upstream_up, *li_asym_upstream_down;
  TLine *li_asym_downstream_up, *li_asym_downstream_down;
  TLine *li_asym_upstream_left, *li_asym_upstream_right;
  TLine *li_asym_downstream_left, *li_asym_downstream_right;
  TLine *li_asym_updown_left, *li_asym_updown_right;
  TLine *li_asym_downup_left, *li_asym_downup_right;

  TGraph *gr_Pol, *gr_Dil;

  std::vector<double> ex;

  //Helper functions
  void setupTLine(TLine *l);
  void setupTGraph(TGraphErrors* gr);
  void setupTGraph(TGraph* gr);

private:
  void setZero();
  void setZero(Int_t nBins);
    
public:
  lr_tgraph() : lrSpinCorr(){lr_tgraph::setZero();}
  lr_tgraph(int nBins, TString name="noName");
  lr_tgraph (int nBins, std::vector<Double_t> &in_bounds,
	     std::vector<Double_t> &in_xval, TString thisName="noName");
  lr_tgraph(TString binfile, TString type, TString thisName="noName");
  lr_tgraph(TFile *f1, TString type, TString thisName="noName");
  
  void Fill();
  void Fill(TString toFill);
  void Draw(TString name, TString opt="AP");
  void Write();
  
};

#endif
