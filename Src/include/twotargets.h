#ifndef TWOTARGETS_H
#define TWOTARGETS_H

#include "leftright.h"
#include "common.h"

class twotargets: public leftright {
protected:
  std::vector<unsigned long long> left_upstream, right_upstream;
  std::vector<unsigned long long> left_downstream, right_downstream;

  std::vector<unsigned long long> left_upstream_up, right_upstream_up;
  std::vector<unsigned long long> left_upstream_down, right_upstream_down;

  std::vector<unsigned long long> left_downstream_up, right_downstream_up;
  std::vector<unsigned long long> left_downstream_down, right_downstream_down;

  std::vector<double> asym, asym_upstream, asym_downstream;
  std::vector<double> e_asym, e_asym_upstream, e_asym_downstream;

  std::vector<double> asym_upstream_up, asym_upstream_down;
  std::vector<double> e_asym_upstream_up, e_asym_upstream_down;

  std::vector<double> asym_downstream_up, asym_downstream_down;
  std::vector<double> e_asym_downstream_up, e_asym_downstream_down;

  //1 Period asymmetries
  std::vector<double> asym_upstream_left, asym_upstream_right;
  std::vector<double> e_asym_upstream_left, e_asym_upstream_right;
  std::vector<double> asym_downstream_left, asym_downstream_right;
  std::vector<double> e_asym_downstream_left, e_asym_downstream_right;

  //sub Period asymmetries
  std::vector<double> asym_updown_left, asym_downup_right;
  std::vector<double> e_asym_updown_left, e_asym_downup_right;
  std::vector<double> asym_downup_left, asym_updown_right;
  std::vector<double> e_asym_downup_left, e_asym_updown_right;

private:  
  //Helper functions
  void setZero(Int_t nBins);

public:
  twotargets() : leftright(){}
  twotargets(int nBins, TString name="noName");
  twotargets (int nBins, std::vector<Double_t> &in_bounds,
	      std::vector<Double_t> &in_xval, TString thisName="noName");
  twotargets (TString binfile, TString type, TString thisName="noName");
  twotargets (TFile *f1, TString type, TString thisName="noName");
    
  bool BinDataCounts(TString target, double binVal,
		     std::vector<double> &boundVals, int noprint=0);
  virtual bool BinDataCounts(TString target, double binVal, int noprint=0);
  virtual bool BinLeftRight();//make all left/right asymmetries

  void Print_LR(TString target);
  void Print_Asym(TString asymName);
  
};

#endif
