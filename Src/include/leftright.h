#ifndef LEFTRIGHT_H
#define LEFTRIGHT_H

#include "common.h"

class leftright {
protected:
  std::vector<unsigned long long> left, right;
  std::vector<double> bounds, xval;
    
  int nBins;
  TString thisName;

  //Helper functions
  void CoutLoop(std::vector<unsigned long long> &counts, TString name);
  void CoutLoop(std::vector<double> &counts, TString name);

  double AsymmetryError(long long A, long long B);

  bool BinOne_LeftRight(std::vector<unsigned long long> &left,
			std::vector<unsigned long long> &right,
			std::vector<double> &asym, std::vector<double> &e_asym,
			TString name);

  bool MakeAvg(std::vector<double_t> &avg, std::vector<int> &count);

  void setOneBounds(TFile *f1, TString bound, TString xval);
    
public:
  leftright();
  leftright(int nBins, TString name="noName");
  leftright (int nBins, std::vector<Double_t> &in_bounds,
	     std::vector<Double_t> &in_xval, TString thisName="noName");
  leftright (TString binfile, TString type, TString thisName="noName");
  leftright (TFile *f1, TString type, TString thisName="noName");

  virtual bool BinDataCounts(TString target, double binVal,
			     std::vector<double> &boundVals, int noprint=0) = 0;
  virtual bool BinDataCounts(TString target, double binVal, int noprint=0) = 0;
  virtual bool BinDataCounts(int whichTarget, TString target, double binVal,
			     int noprint=0){}
  virtual bool BinLeftRight() = 0;//make all left/right asymmetries

  virtual void Print_LR(TString target) = 0;
  virtual void Print_Asym(TString asymName) = 0;
  virtual void Print_Bounds() { CoutLoop(this->bounds, "bounds"); }
  virtual void Print_Bounds(TString target) { CoutLoop(this->bounds, "bounds");}
  virtual void Print_xVal() { CoutLoop(this->xval, "xval"); }
  virtual void Print_xVal(TString target) { CoutLoop(this->xval, "xval"); }

  virtual bool SetCorr(TFile *f1, TString type){return true;}
  virtual bool SetCorr(TString target,double binVal,double avgVal){return true;}
  virtual void AvgCorr(){}
  virtual void CorrectDilPol(){}
  virtual void PrintCorr(TString name){}

  virtual void Fill(){}
  virtual void Fill(TString toFill){}
  virtual void Draw(TString name, TString opt="AP"){}
  virtual void Write(){}
};

#endif
