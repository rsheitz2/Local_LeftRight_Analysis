#include "common.h"

Double_t AsymmetryError(long long A, long long B){
  //if (TMath::Abs(4.0*A*B/( (A+B)*(A+B)*(A+B) ) ) < 10e-6 ){
  //  std::cout << "Errors approximated as Sqrt(1/N)" << std::endl;
  //  return TMath::Sqrt( 1.0/(A+B) );
  //}//might be needed for generated errors

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


Bool_t BinDataCounts(unsigned long long *counts, Double_t binVal,
		     std::vector<Double_t> &binValBounds, Int_t noprint){

  Int_t iter=-1;
  for (std::vector<Double_t>::iterator it=binValBounds.begin();
       it!=binValBounds.end(); it++, iter++){
    if(binVal <= *it ) {
      if(iter==-1){
	if (noprint) return false;//Don't cout anything
	
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << "bin value too low!!!!" << std::endl;
	std::cout << *it << " " << binVal << std::endl;
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << " " << std::endl;
	return false;
      }
      
      counts[iter]++;
      return true;
    }
  }

  if (noprint) return false;//Don't cout anything
  
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << binValBounds.back() << " " << binVal << std::endl;
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

    if (!UpStream_sum) {
      Asym_UpStream[i] = -10.0;//No data
      e_Asym_UpStream[i] = 0.0;

      std::cout<<"!!!!!! No Data Upstream for L/R Asymmetry !!!!!"<<std::endl;
    }
    else {
      Asym_UpStream[i] = 1.0*UpStream_diff/(1.0*UpStream_sum);
      e_Asym_UpStream[i] = AsymmetryError(Left_UpStream[i], Right_UpStream[i]);
    }

    
    long long DownStream_diff = Left_DownStream[i] - Right_DownStream[i];
    long long DownStream_sum = Left_DownStream[i] + Right_DownStream[i];

    if (!DownStream_sum) {
      Asym_DownStream[i] = -10.0;//No data
      e_Asym_DownStream[i] = 0.0;

      std::cout<<"!!!!!! No Data Downstream for L/R Asymmetry !!!!!"<<std::endl;
    }
    else {
      Asym_DownStream[i] = 1.0*DownStream_diff/(1.0*DownStream_sum);
      e_Asym_DownStream[i] = AsymmetryError(Left_DownStream[i],
					    Right_DownStream[i]);
    }

    
    if (!UpStream_sum && !DownStream_sum) {
      Asym[i] = -10.0;//No data
      e_Asym[i] = 0.0;

      std::cout<<"!!!!!! No Data   UP/DOWN   stream for L/R Asymmetry !!!!!"
	       <<std::endl;
    }
    else {
      Asym[i] = (UpStream_diff +
		 DownStream_diff)/(1.0*(UpStream_sum + DownStream_sum));
      e_Asym[i] = AsymmetryError(UpStream_sum, DownStream_sum);
    }
			       
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


Bool_t BinAvg(std::vector<Double_t> &Avg, std::vector<Int_t> &count, 
	      Double_t binVal, std::vector<Double_t> &binValBounds, 
	      Double_t avgVal){
  Int_t iter = 0;
  for (std::vector<Double_t>::iterator it = binValBounds.begin();
       it != binValBounds.end(); it++){
    if(iter == 0 && binVal < *it ) {
      std::cout << "!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "bin value too low!!!!" << std::endl;
      std::cout << "Lower bound: " << *it << "    val: "
		<< binVal << std::endl;
      std::cout << " " << std::endl;

      return false;
    }
    else if (binVal < *(it+1)){
      Avg.at(iter) += avgVal;
      count.at(iter)++;

      return true;
    }

    iter++;
  }
  
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << "Upper bound: " << binValBounds.back() << "    val:"
	    << binVal << std::endl;
  std::cout << " " << std::endl;
  return false;
}


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


Bool_t ChooseLeftRight(TString leftrightChoice, Double_t phi_photon,
		       Double_t spin){
  
  if (leftrightChoice=="Spin" || leftrightChoice==""){//Spin influenced left/right
    if (phi_photon<TMath::Pi()/2 && phi_photon>-TMath::Pi()/2 && spin>0 )// Target spin up
      return true;
    else if (phi_photon<3*TMath::Pi()/2 && phi_photon>TMath::Pi()/2 && spin>0 )// Target spin up
      return false;
    else if (phi_photon< 3*TMath::Pi()/2 && phi_photon > TMath::Pi()/2 && spin < 0) // Target spin down
      return true;
    else if (phi_photon < TMath::Pi()/2 && phi_photon > -TMath::Pi()/2 && spin < 0)// Target spin down
      return false;
    else {
      std::cout << "No Left or Right choosen" << std::endl;
      std::cout << phi_photon << " " << spin << "\n" << std::endl;
      exit(EXIT_FAILURE);
    }
  }
  else if (leftrightChoice=="True"){//True spectrometer left/right (no spin influence)
    if (phi_photon < TMath::Pi()/2 && phi_photon > -TMath::Pi()/2)
      return true;
    else if (phi_photon <3*TMath::Pi()/2 && phi_photon>TMath::Pi()/2)
      return false;
    else {
      std::cout << "No Left or Right choosen" << std::endl;
      std::cout << phi_photon << " " << spin << "\n" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  return false;
}//End ChooseLeftRight
