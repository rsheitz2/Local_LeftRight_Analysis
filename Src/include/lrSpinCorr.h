#ifndef LRSPINCORR_H
#define LRSPINCORR_H

#include "twotargets.h"
#include "common.h"

class lrSpinCorr: public twotargets {
protected:
  double pol_int, dil_int;
  int pol_int_count, dil_int_count;
  
  std::vector<double> avgPol, avgDil;
  std::vector<int> avgPol_count, avgDil_count;
  
  std::vector<double> avgPol_upstream, avgDil_upstream;
  std::vector<int> avgPol_count_upstream, avgDil_count_upstream;
  std::vector<double> avgPol_downstream, avgDil_downstream;
  std::vector<int> avgPol_count_downstream, avgDil_count_downstream;

  std::vector<double> avgPol_upstream_up, avgDil_upstream_up;
  std::vector<int> avgPol_count_upstream_up, avgDil_count_upstream_up;
  std::vector<double> avgPol_upstream_down, avgDil_upstream_down;
  std::vector<int> avgPol_count_upstream_down, avgDil_count_upstream_down;
  std::vector<double> avgPol_downstream_up, avgDil_downstream_up;
  std::vector<int> avgPol_count_downstream_up, avgDil_count_downstream_up;
  std::vector<double> avgPol_downstream_down, avgDil_downstream_down;
  std::vector<int> avgPol_count_downstream_down, avgDil_count_downstream_down;

  //Helper functions
  void setOneCorr(TFile *f1, std::vector<double> &corr, TString corrType);
  bool binOneAvg(std::vector<double> &Avg, std::vector<int> &count, 
		 double binVal, double avgVal);
  bool makeOneAvg(std::vector<double> &Avg, std::vector<int> &count,
		  TString target);
  bool correctOneDilPol(std::vector<double> &asym, std::vector<double> &e_asym,
			std::vector<double> &pol, std::vector<double> &dil);

public:
  lrSpinCorr() : twotargets(){}
  lrSpinCorr(int nBins, TString name="noName");
  lrSpinCorr (int nBins, std::vector<Double_t> &in_bounds,
	      std::vector<Double_t> &in_xval, TString thisName="noName");
  lrSpinCorr (TString binfile, TString type, TString thisName="noName");
  lrSpinCorr (TFile *f1, TString type, TString thisName="noName");

  bool SetCorr(TFile *f1, TString type);
  bool SetCorr(TString target, double binVal, double avgVal);
  void AvgCorr();
  void CorrectDilPol();

  void PrintCorr(TString name);
};

#endif
