#ifndef BINDIST_H
#define BINDIST_H

#include "common.h"

class binDist{
private:
  std::vector<double> bounds;
  
  TString thisName;
  Int_t hbins;
  Double_t xMin, xMax;

  Int_t nBins;
  std::vector<TH1D*> left_upstream_up, right_upstream_up;
  std::vector<TH1D*> left_upstream_down, right_upstream_down;

  std::vector<TH1D*> left_downstream_up, right_downstream_up;
  std::vector<TH1D*> left_downstream_down, right_downstream_down;

  //Helper functions
  void DefineOneHist(TString hName, std::vector<TH1D*> &vect);
  
public:
  binDist(){std::cout<< "No distribution made" << std::endl;}
  binDist(TString name, TString binVar, TString binfile,
	       Int_t h_bins, Double_t xmn, Double_t xmx);

  Bool_t BinFill(TString target, Double_t distVar, Double_t binVar);
  
  void Write();
};

#endif
