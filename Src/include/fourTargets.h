#ifndef FOURTARGETS_H
#define FOURTARGETS_H

#include "twotargets.h"
#include "common.h"

class FourTargets: public lrSpinCorr {
protected:
  std::vector<Double_t> left_upSup_upP, right_upSup_upP;
  std::vector<Double_t> left_upSdown_upP, right_upSdown_upP;
  std::vector<Double_t> left_upSup_downP, right_upSup_downP;
  std::vector<Double_t> left_upSdown_downP, right_upSdown_downP;

  std::vector<Double_t> left_downSup_upP, right_downSup_upP;
  std::vector<Double_t> left_downSdown_upP, right_downSdown_upP;
  std::vector<Double_t> left_downSup_downP, right_downSup_downP;
  std::vector<Double_t> left_downSdown_downP, right_downSdown_downP;

  std::vector<double> avgPolDil_upSup, avgPolDil_upSdown;
  std::vector<int> avgPolDil_count_upSup, avgPolDil_count_upSdown;
  
  std::vector<double> avgPolDil_downSup, avgPolDil_downSdown;
  std::vector<int> avgPolDil_count_downSup, avgPolDil_count_downSdown;

private:  
  //Helper functions
  void setZero(Int_t nBins);
  bool binOneAvg(std::vector<double> &Avg, std::vector<int> &count, 
		 double binVal, double avgVal);
  bool makeOneAvg(std::vector<double> &Avg, std::vector<int> &count,
		  TString target);
  void setup(TGraph *gr) {
    gr->SetMarkerStyle(21);
  }

public:
  FourTargets() : lrSpinCorr(){}
  FourTargets (TString binfile, TString type, TString thisName="noName");
  
  bool BinDataCounts(Int_t target, Bool_t Left, Double_t Spin, double binVal,
		     int noprint=0) override;
  bool SetCorr(double binVal, double avgPol, double avgDil,
	       Int_t target, bool left, Double_t spin) override;
  void AvgCorr() override;

  void Print_LR(TString target) override;
  void Print_PolDil(TString target) override;

  void WriteAll() override;
};

#endif
