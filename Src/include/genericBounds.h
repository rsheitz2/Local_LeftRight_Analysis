#ifndef GENERICBOUNDS_H
#define GENERICBOUNDS_H

#include "lrSpinCorr.h"
#include "common.h"

class genericBounds: virtual public lrSpinCorr {
protected:
  std::vector<double> upstream_bounds, upstream_xval;
  std::vector<double> downstream_bounds, downstream_xval;

  //Helper functions
  void avgOneBound(std::vector<double> &bounds, std::vector<double> &xval);
  bool binOneAvgBounds(std::vector<double> &bounds, std::vector<double> &Avg,
		       std::vector<int> &count, double binVal, double avgVal);
    
public:
  genericBounds() : lrSpinCorr() {}
  genericBounds(int nBins, TString name="noName");
  genericBounds (int nBins, std::vector<Double_t> &in_bounds,
		 std::vector<Double_t> &in_xval, TString thisName="noName");
  genericBounds (TString binfile, TString type, TString thisName="noName");
  genericBounds (TFile *f1, TString type, TString thisName="noName");
  
  bool SetCorrBounds(TString target, double binVal,
		     double avgVal);
    
  void Print_Bounds(TString target);
  void Print_xVal(TString target);

  virtual void Fill(){}
  virtual void Fill(TString toFill){}
  virtual void Draw(TString name, TString opt="AP"){}
  virtual void Write(){}
    
  void SetUpBounds(TString target, double boundVal);
  void MakeBounds();

  bool BinDataCounts(int whichTarget, TString target, double binVal,
		     int noprint=0);

};

#endif
