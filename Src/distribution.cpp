#include "distribution.h"

distribution::distribution(TString name, Int_t bins, Double_t xmn,
			   Double_t xmx){
  this->thisName = name;
  this->nBins = bins;
  this->xMin = xmn;
  this->xMax = xmx;
  
  this->dist = new TH1D(this->thisName, this->thisName, this->nBins,
			this->xMin, this->xMax);
}
