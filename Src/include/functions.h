Double_t AsymmetryError(long long A, long long B);

Bool_t BinDataCounts(unsigned long long *counts, Int_t nBins, Double_t binVal,
		     Double_t *binValBounds);

Bool_t BinDataCounts(unsigned long long *counts, Int_t nBins, Double_t binVal,
		     Double_t *binValBounds, Int_t noprint);

Bool_t BinDataFill(TH1D** h1, Double_t fillVal, Int_t nBins, Double_t binVal,
		   Double_t *binValBounds);

Double_t ShiftPhiSimple(Double_t);

Bool_t BinnedLeftRight(unsigned long long*, unsigned long long*,
		       unsigned long long*, unsigned long long*,
		       Double_t*, Double_t*, Double_t*,
		       Double_t*, Double_t*, Double_t*, Int_t);

Bool_t BinnedLeftRight(unsigned long long *Left, unsigned long long *Right,
		       Double_t* Asym, Double_t *e_Asym, Int_t nBins);

void Hist_LRAsym(TH1D* h_input, TH1D* h_output, Int_t nBins);

void Hist_LR_MirrorAsym(TH1D* h_input, TH1D* h_output, Int_t nBins);

void CorrectDilPol(Double_t* Asym_UpStream, Double_t* Asym_DownStream,
		   Double_t* Asym, Double_t* e_Asym_UpStream,
		   Double_t* e_Asym_DownStream, Double_t* e_Asym,
		   TVectorD Dil, TVectorD Pol, Int_t nBins);

void CorrectDilPol(Double_t* Asym, Double_t* e_Asym,
		   TVectorD Dil, TVectorD Pol, Int_t nBins);
