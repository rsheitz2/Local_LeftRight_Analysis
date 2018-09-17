#ifndef BINDIST4TARG_H
#define BINDIST4TARG_H

#include "binDist.h"
#include "common.h"

class binDist4Targ: public binDist{
protected:
  std::vector<TH1D*> left_upSup_upP, right_upSup_upP;
  std::vector<TH1D*> left_upSup_downP, right_upSup_downP;
  std::vector<TH1D*> left_upSdown_upP, right_upSdown_upP;
  std::vector<TH1D*> left_upSdown_downP, right_upSdown_downP;

  std::vector<TH1D*> left_downSup_upP, right_downSup_upP;
  std::vector<TH1D*> left_downSup_downP, right_downSup_downP;
  std::vector<TH1D*> left_downSdown_upP, right_downSdown_upP;
  std::vector<TH1D*> left_downSdown_downP, right_downSdown_downP;


private:
  //Helper functions
  void DefineOneHist(TString hName, std::vector<TH1D*> &vect);
  
public:
  binDist4Targ() : binDist(){};
  binDist4Targ(TString name, TString binVar, TString binfile,
	       Int_t h_bins, Double_t xmn, Double_t xmx);
  
  Bool_t BinFill(Int_t target, Bool_t Left, Double_t spin,
		 Double_t distVar, Double_t binVar) override;

  void Write() override;
};

#endif
