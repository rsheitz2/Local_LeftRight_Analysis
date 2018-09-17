#include "twotargets.h"
#include "lrSpinCorr.h"

lrSpinCorr::lrSpinCorr(int nBins, TString name)
  : twotargets(nBins, name){
  this->pol_int = 0.0, this->dil_int = 0.0;
  this->pol_int_count = 0, this->dil_int_count = 0;

  lrSpinCorr::setZero(nBins);
}


lrSpinCorr::lrSpinCorr (int nBins, std::vector<Double_t> &in_bounds,
			std::vector<Double_t> &in_xval, TString thisName)
  : twotargets(nBins, in_bounds, in_xval, thisName){
  this->pol_int = 0.0, this->dil_int = 0.0;
  this->pol_int_count = 0, this->dil_int_count = 0;
  
  lrSpinCorr::setZero(nBins);
}


lrSpinCorr::lrSpinCorr(TString binfile, TString type, TString thisName)
  : twotargets(binfile, type, thisName){
  this->pol_int = 0.0, this->dil_int = 0.0;
  this->pol_int_count = 0, this->dil_int_count = 0;
  
  lrSpinCorr::setZero(nBins);
}


lrSpinCorr::lrSpinCorr(TFile *f1, TString type, TString thisName)
  : twotargets(f1, type, thisName){
  this->pol_int = 0.0, this->dil_int = 0.0;
  this->pol_int_count = 0, this->dil_int_count = 0;
  
  lrSpinCorr::setZero(nBins);
}


bool lrSpinCorr::SetCorr(TFile *f1, TString type){
  setOneCorr(f1, this->avgPol, Form("Pol_%s",type.Data()));
  
  setOneCorr(f1, this->avgPol_upstream, Form("Pol_%s_UpStream",type.Data()));
  setOneCorr(f1, this->avgPol_downstream,Form("Pol_%s_DownStream",type.Data()));
  setOneCorr(f1, this->avgPol_upstream_up,
	     Form("Pol_%s_UpStream_Up",type.Data()));
  setOneCorr(f1, this->avgPol_upstream_down,
	     Form("Pol_%s_UpStream_Down",type.Data()));
  setOneCorr(f1, this->avgPol_downstream_up,
	     Form("Pol_%s_DownStream_Up",type.Data()));
  setOneCorr(f1, this->avgPol_downstream_down,
	     Form("Pol_%s_DownStream_Down",type.Data()));
  
  setOneCorr(f1, this->avgDil, Form("Dil_%s",type.Data()));
  setOneCorr(f1, this->avgDil_upstream, Form("Dil_%s_UpStream",type.Data()));
  setOneCorr(f1, this->avgDil_downstream,Form("Dil_%s_DownStream",type.Data()));
  setOneCorr(f1, this->avgDil_upstream_up,
	     Form("Dil_%s_UpStream_Up",type.Data()));
  setOneCorr(f1, this->avgDil_upstream_down,
	     Form("Dil_%s_UpStream_Down",type.Data()));
  setOneCorr(f1, this->avgDil_downstream_up,
	     Form("Dil_%s_DownStream_Up",type.Data()));
  setOneCorr(f1, this->avgDil_downstream_down,
	     Form("Dil_%s_DownStream_Down",type.Data()));
  
  TVectorD pol_tmp = *( (TVectorD*)f1->Get("Pol_int") );
  this->pol_int = pol_tmp[0];
  TVectorD dil_tmp = *( (TVectorD*)f1->Get("Dil_int") );
  this->dil_int = dil_tmp[0];

  return true;
}


bool lrSpinCorr::SetCorr(TString target, double binVal, double avgVal,
			 bool left){
  
  if (target=="pol_upstream_up") {
    binOneAvg(this->avgPol_upstream_up, this->avgPol_count_upstream_up,
	      binVal, avgVal);
    binOneAvg(this->avgPol_upstream, this->avgPol_count_upstream,
	      binVal, avgVal);
    binOneAvg(this->avgPol, this->avgPol_count, binVal, avgVal);
    
    if (left){
      binOneAvg(this->avgPol_upstream_left, this->avgPol_count_upstream_left,
		binVal, avgVal);
      binOneAvg(this->avgPol_updown_left, this->avgPol_count_updown_left,
		binVal, avgVal);
    }
    else{
      binOneAvg(this->avgPol_upstream_right, this->avgPol_count_upstream_right,
		binVal, avgVal);
      binOneAvg(this->avgPol_updown_right, this->avgPol_count_updown_right,
		binVal, avgVal);
    }
    
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="pol_upstream_down") {
    binOneAvg(this->avgPol_upstream_down, this->avgPol_count_upstream_down,
	      binVal, avgVal);
    binOneAvg(this->avgPol_upstream, this->avgPol_count_upstream,
	      binVal, avgVal);
    binOneAvg(this->avgPol, this->avgPol_count, binVal, avgVal);
    
    if (left){
      binOneAvg(this->avgPol_upstream_right, this->avgPol_count_upstream_right,
		binVal, avgVal);
      binOneAvg(this->avgPol_downup_right, this->avgPol_count_downup_right,
		binVal, avgVal);
      
    }
    else {
      binOneAvg(this->avgPol_upstream_left, this->avgPol_count_upstream_left,
		binVal, avgVal);
      binOneAvg(this->avgPol_downup_left, this->avgPol_count_downup_left,
		binVal, avgVal);
    }

    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="pol_downstream_up") {
    binOneAvg(this->avgPol_downstream_up, this->avgPol_count_downstream_up,
	      binVal, avgVal);
    binOneAvg(this->avgPol_downstream, this->avgPol_count_downstream,
	      binVal, avgVal);
    binOneAvg(this->avgPol, this->avgPol_count, binVal, avgVal);

    if (left){
      binOneAvg(this->avgPol_downstream_left,this->avgPol_count_downstream_left,
		binVal, avgVal);
      binOneAvg(this->avgPol_downup_left, this->avgPol_count_downup_left,
		binVal, avgVal);
    }
    else {
      binOneAvg(this->avgPol_downstream_right,
		this->avgPol_count_downstream_right, binVal, avgVal);
      binOneAvg(this->avgPol_downup_right, this->avgPol_count_downup_right,
		binVal, avgVal);
    }
    
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="pol_downstream_down") {
    binOneAvg(this->avgPol_downstream_down, this->avgPol_count_downstream_down,
	      binVal, avgVal);
    binOneAvg(this->avgPol_downstream, this->avgPol_count_downstream,
	      binVal, avgVal);
    binOneAvg(this->avgPol, this->avgPol_count, binVal, avgVal);

    if (left){
      binOneAvg(this->avgPol_downstream_right,
		this->avgPol_count_downstream_right, binVal, avgVal);
      binOneAvg(this->avgPol_updown_right, this->avgPol_count_updown_right,
		binVal, avgVal);
    }
    else{
      binOneAvg(this->avgPol_downstream_left,this->avgPol_count_downstream_left,
		binVal, avgVal);
      binOneAvg(this->avgPol_updown_left, this->avgPol_count_updown_left,
		binVal, avgVal);
    }
    
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="dil_upstream_up") {
    binOneAvg(this->avgDil_upstream_up, this->avgDil_count_upstream_up,
	      binVal, avgVal);
    binOneAvg(this->avgDil_upstream, this->avgDil_count_upstream,
	      binVal, avgVal);
    binOneAvg(this->avgDil, this->avgDil_count, binVal, avgVal);

    if (left){
      binOneAvg(this->avgDil_upstream_left, this->avgDil_count_upstream_left,
		binVal, avgVal);
      binOneAvg(this->avgDil_updown_left, this->avgDil_count_updown_left,
		binVal, avgVal);
    }
    else{
      binOneAvg(this->avgDil_upstream_right, this->avgDil_count_upstream_right,
		binVal, avgVal);
      binOneAvg(this->avgDil_updown_right, this->avgDil_count_updown_right,
		binVal, avgVal);
    }
    
    dil_int += avgVal; dil_int_count++;
  }
  else if (target=="dil_upstream_down") {
    binOneAvg(this->avgDil_upstream_down, this->avgDil_count_upstream_down,
	      binVal, avgVal);
    binOneAvg(this->avgDil_upstream, this->avgDil_count_upstream,
	      binVal, avgVal);
    binOneAvg(this->avgDil, this->avgDil_count, binVal, avgVal);

    if (left){
      binOneAvg(this->avgDil_upstream_right, this->avgDil_count_upstream_right,
		binVal, avgVal);
      binOneAvg(this->avgDil_downup_right, this->avgDil_count_downup_right,
		binVal, avgVal);
      
    }
    else {
      binOneAvg(this->avgDil_upstream_left, this->avgDil_count_upstream_left,
		binVal, avgVal);
      binOneAvg(this->avgDil_downup_left, this->avgDil_count_downup_left,
		binVal, avgVal);
    }
    
    dil_int += avgVal; dil_int_count++;
  }
  else if (target=="dil_downstream_up") {
    binOneAvg(this->avgDil_downstream_up, this->avgDil_count_downstream_up,
	      binVal, avgVal);
    binOneAvg(this->avgDil_downstream, this->avgDil_count_downstream,
	      binVal, avgVal);
    binOneAvg(this->avgDil, this->avgDil_count, binVal, avgVal);

    if (left){
      binOneAvg(this->avgDil_downstream_left,this->avgDil_count_downstream_left,
		binVal, avgVal);
      binOneAvg(this->avgDil_downup_left, this->avgDil_count_downup_left,
		binVal, avgVal);
    }
    else {
      binOneAvg(this->avgDil_downstream_right,
		this->avgDil_count_downstream_right, binVal, avgVal);
      binOneAvg(this->avgDil_downup_right, this->avgDil_count_downup_right,
		binVal, avgVal);
    }
    
    dil_int += avgVal; dil_int_count++;
  }
  else if (target=="dil_downstream_down") {
    binOneAvg(this->avgDil_downstream_down, this->avgDil_count_downstream_down,
	      binVal, avgVal);
    binOneAvg(this->avgDil_downstream, this->avgDil_count_downstream,
	      binVal, avgVal);
    binOneAvg(this->avgDil, this->avgDil_count, binVal, avgVal);

    if (left){
      binOneAvg(this->avgDil_downstream_right,
		this->avgDil_count_downstream_right, binVal, avgVal);
      binOneAvg(this->avgDil_updown_right, this->avgDil_count_updown_right,
		binVal, avgVal);
    }
    else{
      binOneAvg(this->avgDil_downstream_left,this->avgDil_count_downstream_left,
		binVal, avgVal);
      binOneAvg(this->avgDil_updown_left, this->avgDil_count_updown_left,
		binVal, avgVal);
    }
    
    dil_int += avgVal; dil_int_count++;
  }
  else {
    std::cout << this->thisName << " Wrong target name to " <<
      "lrSpinCorr::SetCorr   " << target << std::endl;
    return false;
  }

  return true;
}


void lrSpinCorr::AvgCorr(){
  this->pol_int = this->pol_int/(1.0*this->pol_int_count);
  makeOneAvg(this->avgPol, this->avgPol_count, "avgPol");
  makeOneAvg(this->avgPol_upstream, this->avgPol_count_upstream,
	     "avgPol_upstream");
  makeOneAvg(this->avgPol_downstream, this->avgPol_count_downstream,
	     "avgPol_downstream");
  makeOneAvg(this->avgPol_upstream_up, this->avgPol_count_upstream_up,
	     "avgPol_upstream_up");
  makeOneAvg(this->avgPol_upstream_down, this->avgPol_count_upstream_down,
	     "avgPol_upstream_down");
  makeOneAvg(this->avgPol_downstream_up, this->avgPol_count_downstream_up,
	     "avgPol_downstream_up");
  makeOneAvg(this->avgPol_downstream_down, this->avgPol_count_downstream_down,
	     "avgPol_downstream_down");
  makeOneAvg(this->avgPol_upstream_left, this->avgPol_count_upstream_left,
	     "avgPol_upstream_left");
  makeOneAvg(this->avgPol_upstream_right, this->avgPol_count_upstream_right,
	     "avgPol_upstream_right");
  makeOneAvg(this->avgPol_downstream_left, this->avgPol_count_downstream_left,
	     "avgPol_downstream_left");
  makeOneAvg(this->avgPol_downstream_right, this->avgPol_count_downstream_right,
	     "avgPol_downstream_right");
  makeOneAvg(this->avgPol_updown_left, this->avgPol_count_updown_left,
	     "avgPol_updown_left");
  makeOneAvg(this->avgPol_updown_right, this->avgPol_count_updown_right,
	     "avgPol_updown_right");
  makeOneAvg(this->avgPol_downup_left, this->avgPol_count_downup_left,
	     "avgPol_downup_left");
  makeOneAvg(this->avgPol_downup_right, this->avgPol_count_downup_right,
	     "avgPol_downup_right");
  
  this->dil_int = this->dil_int/(1.0*this->dil_int_count);
  makeOneAvg(this->avgDil, this->avgDil_count, "avgDil");
  makeOneAvg(this->avgDil_upstream, this->avgDil_count_upstream,
	     "avgDil_upstream");
  makeOneAvg(this->avgDil_downstream, this->avgDil_count_downstream,
	     "avgDil_downstream");
  makeOneAvg(this->avgDil_upstream_up, this->avgDil_count_upstream_up,
	     "avgDil_upstream_up");
  makeOneAvg(this->avgDil_upstream_down, this->avgDil_count_upstream_down,
	     "avgDil_upstream_down");
  makeOneAvg(this->avgDil_downstream_up, this->avgDil_count_downstream_up,
	     "avgDil_downstream_up");
  makeOneAvg(this->avgDil_downstream_down, this->avgDil_count_downstream_down,
	     "avgDil_downstream_down");
  makeOneAvg(this->avgDil_upstream_left, this->avgDil_count_upstream_left,
	     "avgDil_upstream_left");
  makeOneAvg(this->avgDil_upstream_right, this->avgDil_count_upstream_right,
	     "avgDil_upstream_right");
  makeOneAvg(this->avgDil_downstream_left, this->avgDil_count_downstream_left,
	     "avgDil_downstream_left");
  makeOneAvg(this->avgDil_downstream_right, this->avgDil_count_downstream_right,
	     "avgDil_downstream_right");
  makeOneAvg(this->avgDil_updown_left, this->avgDil_count_updown_left,
	     "avgDil_updown_left");
  makeOneAvg(this->avgDil_updown_right, this->avgDil_count_updown_right,
	     "avgDil_updown_right");
  makeOneAvg(this->avgDil_downup_left, this->avgDil_count_downup_left,
	     "avgDil_downup_left");
  makeOneAvg(this->avgDil_downup_right, this->avgDil_count_downup_right,
	     "avgDil_downup_right");
}


void lrSpinCorr::CorrectDilPol(){
  correctOneDilPol(this->asym, this->e_asym, this->avgPol, this->avgDil);

  correctOneDilPol(this->asym_upstream, this->e_asym_upstream,
		   this->avgPol_upstream, this->avgDil_upstream);
  correctOneDilPol(this->asym_downstream, this->e_asym_downstream,
		   this->avgPol_downstream, this->avgDil_downstream);

  correctOneDilPol(this->asym_upstream_up, this->e_asym_upstream_up,
		   this->avgPol_upstream_up, this->avgDil_upstream_up);
  correctOneDilPol(this->asym_upstream_down, this->e_asym_upstream_down,
		   this->avgPol_upstream_down, this->avgDil_upstream_down);

  correctOneDilPol(this->asym_downstream_up, this->e_asym_downstream_up,
		   this->avgPol_downstream_up, this->avgDil_downstream_up);
  correctOneDilPol(this->asym_downstream_down, this->e_asym_downstream_down,
		   this->avgPol_downstream_down, this->avgDil_downstream_down);

  correctOneDilPol(this->asym_upstream_left, this->e_asym_upstream_left,
		   this->avgPol_upstream_left, this->avgDil_upstream_left);
  correctOneDilPol(this->asym_upstream_right, this->e_asym_upstream_right,
		   this->avgPol_upstream_right, this->avgDil_upstream_right);
  correctOneDilPol(this->asym_downstream_left, this->e_asym_downstream_left,
		   this->avgPol_downstream_left, this->avgDil_downstream_left);
  correctOneDilPol(this->asym_downstream_right, this->e_asym_downstream_right,
		   this->avgPol_downstream_right, this->avgDil_downstream_right);

  correctOneDilPol(this->asym_updown_left, this->e_asym_updown_left,
		   this->avgPol_updown_left, this->avgDil_updown_left);
  correctOneDilPol(this->asym_updown_right, this->e_asym_updown_right,
		   this->avgPol_updown_right, this->avgDil_updown_right);
  correctOneDilPol(this->asym_downup_left, this->e_asym_downup_left,
		   this->avgPol_downup_left, this->avgDil_downup_left);
  correctOneDilPol(this->asym_downup_right, this->e_asym_downup_right,
		   this->avgPol_downup_right, this->avgDil_downup_right);
}


void lrSpinCorr::PrintCorr(TString name){

  if (name=="pol_int") {
    std::cout << "Printing content for: " << name 
	      << "  from: " << this->thisName << std::endl;
    std::cout << this->pol_int << std::endl;
    std::cout << " " << std::endl;}
  else if (name=="pol") CoutLoop(this->avgPol, name);
  else if (name=="pol_upstream") CoutLoop(this->avgPol_upstream, name);
  else if (name=="pol_downstream") CoutLoop(this->avgPol_downstream, name);
  else if (name=="pol_upstream_up") CoutLoop(this->avgPol_upstream_up, name);
  else if (name=="pol_upstream_down") CoutLoop(this->avgPol_upstream_down,name);
  else if (name=="pol_downstream_up") CoutLoop(this->avgPol_downstream_up,name);
  else if (name=="pol_downstream_down") {
    CoutLoop(this->avgPol_downstream_down, name);}
  else if (name=="pol_upstream_left") {
    CoutLoop(this->avgPol_upstream_left, name);}
  else if (name=="dil_int") {
    std::cout << "Printing content for: " << name 
	      << "  from: " << this->thisName << std::endl;
    std::cout << this->dil_int << std::endl;
    std::cout << " " << std::endl;}
  else if (name=="dil") CoutLoop(this->avgDil, name);
  else if (name=="dil_upstream") CoutLoop(this->avgDil_upstream, name);
  else if (name=="dil_downstream") CoutLoop(this->avgDil_downstream, name);
  else if (name=="dil_upstream_up") CoutLoop(this->avgDil_upstream_up, name);
  else if (name=="dil_upstream_down") CoutLoop(this->avgDil_upstream_down,name);
  else if (name=="dil_downstream_up") CoutLoop(this->avgDil_downstream_up,name);
  else if (name=="dil_downstream_down") {
    CoutLoop(this->avgDil_downstream_down, name);}
  else if (name=="dil_upstream_left") {
    CoutLoop(this->avgDil_upstream_left, name);}
  else std::cout << "Wrong lrSpinCorr::PrintCorr input name" << std::endl;

}


//Helper functions
void lrSpinCorr::setZero(Int_t nBins){
  for (Int_t i=0; i<nBins; i++) {
    this->avgPol.push_back(0.0); this->avgDil.push_back(0.0);
    this->avgPol_count.push_back(0); this->avgDil_count.push_back(0);

    this->avgPol_upstream.push_back(0.0);
    this->avgDil_upstream.push_back(0.0);
    this->avgPol_count_upstream.push_back(0);
    this->avgDil_count_upstream.push_back(0);
    this->avgPol_downstream.push_back(0.0);
    this->avgDil_downstream.push_back(0.0);
    this->avgPol_count_downstream.push_back(0);
    this->avgDil_count_downstream.push_back(0);

    this->avgPol_upstream_up.push_back(0.0);
    this->avgDil_upstream_up.push_back(0.0);
    this->avgPol_count_upstream_up.push_back(0);
    this->avgDil_count_upstream_up.push_back(0);
    this->avgPol_upstream_down.push_back(0.0);
    this->avgDil_upstream_down.push_back(0.0);
    this->avgPol_count_upstream_down.push_back(0);
    this->avgDil_count_upstream_down.push_back(0);
    this->avgPol_downstream_up.push_back(0.0);
    this->avgDil_downstream_up.push_back(0.0);
    this->avgPol_count_downstream_up.push_back(0);
    this->avgDil_count_downstream_up.push_back(0);
    this->avgPol_downstream_down.push_back(0.0);
    this->avgDil_downstream_down.push_back(0.0);
    this->avgPol_count_downstream_down.push_back(0);
    this->avgDil_count_downstream_down.push_back(0);

    this->avgPol_upstream_left.push_back(0.0);
    this->avgDil_upstream_left.push_back(0.0);
    this->avgPol_count_upstream_left.push_back(0);
    this->avgDil_count_upstream_left.push_back(0);
    this->avgPol_upstream_right.push_back(0.0);
    this->avgDil_upstream_right.push_back(0.0);
    this->avgPol_count_upstream_right.push_back(0);
    this->avgDil_count_upstream_right.push_back(0);
    this->avgPol_downstream_left.push_back(0.0);
    this->avgDil_downstream_left.push_back(0.0);
    this->avgPol_count_downstream_left.push_back(0);
    this->avgDil_count_downstream_left.push_back(0);
    this->avgPol_downstream_right.push_back(0.0);
    this->avgDil_downstream_right.push_back(0.0);
    this->avgPol_count_downstream_right.push_back(0);
    this->avgDil_count_downstream_right.push_back(0);

    this->avgPol_updown_left.push_back(0.0);
    this->avgDil_updown_left.push_back(0.0);
    this->avgPol_count_updown_left.push_back(0);
    this->avgDil_count_updown_left.push_back(0);
    this->avgPol_updown_right.push_back(0.0);
    this->avgDil_updown_right.push_back(0.0);
    this->avgPol_count_updown_right.push_back(0);
    this->avgDil_count_updown_right.push_back(0);
    this->avgPol_downup_left.push_back(0.0);
    this->avgDil_downup_left.push_back(0.0);
    this->avgPol_count_downup_left.push_back(0);
    this->avgDil_count_downup_left.push_back(0);
    this->avgPol_downup_right.push_back(0.0);
    this->avgDil_downup_right.push_back(0.0);
    this->avgPol_count_downup_right.push_back(0);
    this->avgDil_count_downup_right.push_back(0);
  }
}//setZero


void lrSpinCorr::setOneCorr(TFile *f1, std::vector<double> &corr,
			    TString corrType){
  corr.clear();
  TVectorD tmp = *( (TVectorD*)f1->Get(corrType) );
  if (tmp.GetNrows() == 0) {
    std::cout << corrType << " did not open from lrSpinCorr::setOneCorr"
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  for (Int_t i=0; i<this->nBins; i++) corr.push_back(tmp[i]);
};


bool lrSpinCorr::binOneAvg(std::vector<double> &Avg, std::vector<int> &count, 
			   double binVal, double avgVal){
  
  int iter = -1;
  for (std::vector<double>::iterator it=this->bounds.begin();
       it!=this->bounds.end(); it++, iter++){
   
    if(iter == -1 && binVal < *it ) {
      std::cout << "bin value too low!!!!" << std::endl;
      std::cout << "Lower bound: " << *it << "    val: "
		<< binVal << "   from: " << this->thisName << std::endl;
      std::cout << " " << std::endl;

      return false;
    }
    else if (binVal < *(it)){
      Avg.at(iter) += avgVal;
      count.at(iter)++;

      return true;
    }
  }
  
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << "Upper bound: " << bounds.back() << "    val:"
	    << binVal << "   from: " << this->thisName << std::endl;
  std::cout << " " << std::endl;
  return false;
}


bool lrSpinCorr::makeOneAvg(std::vector<double> &Avg, std::vector<int> &count,
			    TString target){
  if (Avg.size() != count.size() ){
    std::cout << "lrSpinCorr::makeOneAvg size mismatch in: " << target <<
      "     from: " << this->thisName << std::endl;
    return false;
  }
  
  int iter = 0;
  for (std::vector<double>::iterator it=Avg.begin();it!=Avg.end(); it++, iter++)
    {
      *it = *it/( 1.0*count.at(iter) );
    }

  return true;
}


bool lrSpinCorr::correctOneDilPol(std::vector<double> &asym,
				  std::vector<double> &e_asym,
				  std::vector<double> &pol,
				  std::vector<double> &dil){

  if (asym.size() != pol.size() || asym.size() != dil.size() ){
    std::cout << "lrSpinCorr::correctOneDilPol size mismatch from: " <<
      this->thisName << std::endl;

    return false;
  }

  int iter = 0;
  for (std::vector<double>::iterator it=asym.begin(); it!=asym.end();
       it++, iter++){
    *it = *it/( pol.at(iter)*dil.at(iter) );
  }

  iter = 0;
  for (std::vector<double>::iterator it=e_asym.begin(); it!=e_asym.end();
       it++, iter++){
    *it = *it/( pol.at(iter)*dil.at(iter) );
  }

  return true;
}
