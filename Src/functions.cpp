#include "common.h"

Double_t AsymmetryError(long long A, long long B){
  return TMath::Sqrt( 4.0*A*B/( (A+B)*(A+B)*(A+B) ) );
}


Bool_t BinDataCounts(unsigned long long *counts, Int_t nBins, Double_t binVal,
		     Double_t *binValBounds){

  if(binVal < binValBounds[0] ) {
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "bin value too low!!!!" << std::endl;
    std::cout << binValBounds[0] << " " << binVal << std::endl;
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " " << std::endl;
    return false;
  }
  for (Int_t i=0; i<nBins; i++) {
    if(binVal <= binValBounds[i+1] ) {
      counts[i]++;
      return true;
    }
  }
  
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << binValBounds[nBins] << " " << binVal << std::endl;
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << " " << std::endl;
  return false;
}


Bool_t BinDataCounts(unsigned long long *counts, Int_t nBins, Double_t binVal,
		     Double_t *binValBounds, Int_t noprint){

  if(binVal < binValBounds[0] ) {
    if (noprint) return false;//Don't cout anything
    
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "bin value too low!!!!" << std::endl;
    std::cout << binValBounds[0] << " " << binVal << std::endl;
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " " << std::endl;
    return false;
  }
  for (Int_t i=0; i<nBins; i++) {
    if(binVal <= binValBounds[i+1] ) {
      counts[i]++;
      return true;
    }
  }

  if (noprint) return false;//Don't cout anything
  
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << binValBounds[nBins] << " " << binVal << std::endl;
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << " " << std::endl;
  return false;
}


Bool_t BinDataFill(TH1D** h1, Double_t fillVal, Int_t nBins, Double_t binVal,
		     Double_t *binValBounds){

  if(binVal < binValBounds[0] ) {
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << "bin value too low!!!!" << std::endl;
    std::cout << binValBounds[0] << " " << binVal << std::endl;
    std::cout << "!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " " << std::endl;
    return false;
  }
  for (Int_t i=0; i<nBins; i++) {
    if(binVal <= binValBounds[i+1] ) {
      h1[i]->Fill(fillVal);
      return true;
    }
  }
  
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << binValBounds[nBins] << " " << binVal << std::endl;
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << " " << std::endl;
  return false;
}//BinDataFill


Double_t  ShiftPhiSimple (Double_t PhiS_simple) {
  //Lab reference frame:
  //    -Pi/2 ->   P/i/2 = left
  //    Pi/2  -> 3*P/i/2 = right
  Double_t phi = TMath::Pi()/2 - PhiS_simple;
  
  return phi;
}//ShiftPhiSimple


Bool_t BinnedLeftRight(unsigned long long *Left_UpStream,
		       unsigned long long *Right_UpStream,
		       unsigned long long *Left_DownStream,
		       unsigned long long *Right_DownStream,
		       Double_t* Asym_UpStream, Double_t *Asym_DownStream,
		       Double_t *Asym, Double_t *e_Asym_UpStream,
		       Double_t *e_Asym_DownStream, Double_t *e_Asym,
		       Int_t nBins){
  
  for (Int_t i=0; i<nBins; i++) {
    //Asymmetries
    long long UpStream_diff = Left_UpStream[i] - Right_UpStream[i];
    long long UpStream_sum = Left_UpStream[i] + Right_UpStream[i];
    Asym_UpStream[i] = 1.0*UpStream_diff/(1.0*UpStream_sum);

    long long DownStream_diff = Left_DownStream[i] - Right_DownStream[i];
    long long DownStream_sum = Left_DownStream[i] + Right_DownStream[i];
    Asym_DownStream[i] = 1.0*DownStream_diff/(1.0*DownStream_sum);

    Asym[i] = (UpStream_diff +
	       DownStream_diff)/(1.0*(UpStream_sum + DownStream_sum));

    //Error cals
    e_Asym_UpStream[i] = AsymmetryError(Left_UpStream[i], Right_UpStream[i]);
    e_Asym_DownStream[i] = AsymmetryError(Left_DownStream[i],
					  Right_DownStream[i]);
    
    e_Asym[i] = AsymmetryError(UpStream_sum, DownStream_sum);
			       
  }

  return true;
}//BinnedLeftRight


Bool_t BinnedLeftRight(unsigned long long *Left, unsigned long long *Right,
		       Double_t* Asym, Double_t *e_Asym, Int_t nBins){
  
  for (Int_t i=0; i<nBins; i++) {
    //Asymmetries
    long long diff = Left[i] - Right[i];
    long long sum = Left[i] + Right[i];
    Asym[i] = 1.0*diff/(1.0*sum);

    //Error cals
    e_Asym[i] = AsymmetryError(Left[i], Right[i]);
  }

  return true;
}//BinnedLeftRight


void Hist_LRAsym(TH1D* h_input, TH1D* h_output, Int_t nBins){
  if (nBins % 2 != 0){
    std::cout << "Histogram bin number should be even" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  for (Int_t i=0; i<nBins/2; i++) {
    long long Left = h_input->GetBinContent(i+1);
    long long Right = h_input->GetBinContent(nBins/2+i+1);
    Double_t LR = 1.0*(Left - Right)/( 1.0*(Left + Right) );
    h_output->SetBinContent(i+1, LR);
    Double_t error = AsymmetryError(Left, Right);
    h_output->SetBinError(i+1, error);
  }
}//Hist_LRAsym


void Hist_LR_MirrorAsym(TH1D* h_input, TH1D* h_output, Int_t nBins){
  if (nBins % 2 != 0){
    std::cout << "Histogram bin number should be even" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  for (Int_t i=0; i<nBins/2; i++) {
    long long Left = h_input->GetBinContent(i+1);
    long long RightMirror = h_input->GetBinContent(nBins-i);
    Double_t LR_Mirror = 1.0*(Left - RightMirror)/( 1.0*(Left + RightMirror) );
    h_output->SetBinContent(i+1, LR_Mirror);
    Double_t error = AsymmetryError(Left, RightMirror);
    h_output->SetBinError(i+1, error);
  }
}//Hist_LR_MirrorAsym


void CorrectDilPol(Double_t* Asym_UpStream, Double_t* Asym_DownStream,
		   Double_t* Asym, Double_t* e_Asym_UpStream,
		   Double_t* e_Asym_DownStream, Double_t* e_Asym,
		   TVectorD Dil, TVectorD Pol, Int_t nBins){

  for (Int_t i=0; i<nBins; i++) {
    Asym_UpStream[i] = Asym_UpStream[i]/(Dil[i]*Pol[i]);
    Asym_DownStream[i] = Asym_DownStream[i]/(Dil[i]*Pol[i]);
    Asym[i] = Asym[i]/(Dil[i]*Pol[i]);

    e_Asym_UpStream[i] = e_Asym_UpStream[i]/(Dil[i]*Pol[i]);
    e_Asym_DownStream[i] = e_Asym_DownStream[i]/(Dil[i]*Pol[i]);
    e_Asym[i] = e_Asym[i]/(Dil[i]*Pol[i]);

  }
}//End CorrectDilPol


void CorrectDilPol(Double_t* Asym, Double_t* e_Asym,
		   TVectorD Dil, TVectorD Pol, Int_t nBins){

  for (Int_t i=0; i<nBins; i++) {
    Asym[i] = Asym[i]/(Dil[i]*Pol[i]);

    e_Asym[i] = e_Asym[i]/(Dil[i]*Pol[i]);
  }
}//End CorrectDilPol
