#ifndef BINDIST_H
#define BINDIST_H

#include "common.h"

class binDist{
protected:
  std::vector<double> bounds;
  
  TString thisName;
  Int_t hbins;
  Double_t xMin, xMax;
  Bool_t varBinning;

  Int_t nBins;
  std::vector<TH1D*> left_upstream_up, right_upstream_up;
  std::vector<TH1D*> left_upstream_down, right_upstream_down;

  std::vector<TH1D*> left_downstream_up, right_downstream_up;
  std::vector<TH1D*> left_downstream_down, right_downstream_down;

private:
  //Helper functions
  void GetBounds(TString binVar, TString binfile);
  void DefineOneHist(TString hName, std::vector<TH1D*> &vect);
  void DefineOneHist(TString hName, std::vector<TH1D*> &vect, Double_t *bins);
  
public:
  binDist(){std::cout<< "No distribution made" << std::endl;}
  binDist(TString name, TString binVar, TString binfile,
	       Int_t h_bins, Double_t xmn, Double_t xmx);
  binDist(TString name, TString binVar, TString binfile,
	  Int_t h_bins, Double_t *bins);
    
  virtual Bool_t BinFill(Int_t target, Bool_t Left, Double_t spin,
			 Double_t distVar, Double_t binVar){return true;}
  virtual Bool_t BinFill(TString target, Double_t distVar, Double_t binVar);
  Bool_t BinFill(TString target, Double_t distVar, Int_t bin);

  virtual void Print();
  virtual void Write();
};

#endif
