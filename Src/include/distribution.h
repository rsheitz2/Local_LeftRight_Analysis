#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "common.h"

class distribution{
private:
  TString thisName;
  Int_t nBins;
  Double_t xMin, xMax;
  
public:
  TH1D* dist = NULL;
  
  distribution(){std::cout<< "No distribution made" << std::endl;}
  distribution(TString name, Int_t bins, Double_t xmn, Double_t xmx);
  
  virtual void Write(){dist->Write(this->thisName);}
};

#endif
