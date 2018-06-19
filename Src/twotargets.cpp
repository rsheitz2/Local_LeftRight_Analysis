#include "leftright.h"
#include "twotargets.h"

twotargets::twotargets (int nBins, TString thisName)
  : leftright (nBins, thisName){

  twotargets::setZero(nBins);
}


twotargets::twotargets (int nBins, std::vector<Double_t> &in_bounds,
			std::vector<Double_t> &in_xval, TString thisName)
  : leftright(nBins, in_bounds, in_xval, thisName){

  twotargets::setZero(nBins);
}


twotargets::twotargets (TString binfile, TString type, TString thisName)
  : leftright(binfile, type, thisName){

  twotargets::setZero(nBins);
}


twotargets::twotargets (TFile *f1, TString type, TString thisName)
  : leftright(f1, type, thisName){

  twotargets::setZero(nBins);
}


bool twotargets::BinDataCounts(TString target, Double_t binVal, 
			       std::vector<Double_t> &binValBounds,int noprint){
  
  Int_t iter=-1;
  for (std::vector<Double_t>::iterator it=binValBounds.begin();
       it!=binValBounds.end(); it++, iter++){
    if(binVal <= *it ) {
      if(iter==-1){
	if (noprint) return false;//Don't cout anything
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << "bin value too low!!!!" << std::endl;
	std::cout << this->thisName << " has " <<
	  *it << " < " << binVal <<std::endl;
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << " " << std::endl;
	return false;
      }
      
      if (target=="left_upstream_up") {
	this->left_upstream_up.at(iter)++;
	this->left_upstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_upstream_up") {
	this->right_upstream_up.at(iter)++;
	this->right_upstream.at(iter)++;
	this->right.at(iter)++;}
      else if (target=="left_upstream_down") {
	this->left_upstream_down.at(iter)++;
	this->left_upstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_upstream_down") {
	this->right_upstream_down.at(iter)++;
	this->right_upstream.at(iter)++;
	this->right.at(iter)++;}
      else if (target=="left_downstream_up") {
	this->left_downstream_up.at(iter)++;
	this->left_downstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_downstream_up") {
	this->right_downstream_up.at(iter)++;
	this->right_downstream.at(iter)++;
	this->right.at(iter)++;}
      else if (target=="left_downstream_down") {
	this->left_downstream_down.at(iter)++;
	this->left_downstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_downstream_down") {
	this->right_downstream_down.at(iter)++;
	this->right_downstream.at(iter)++;
	this->right.at(iter)++;}
      else {
	std::cout << this->thisName << " Wrong target name to " <<
	  "twotargets::BinDataCounts" << std::endl;
	return false;
      }
      
      return true;
    }
  }

  if (noprint) return false;//Don't cout anything
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << this->thisName << " has " <<
    binValBounds.back() << " > " << binVal << std::endl;
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << " " << std::endl;
  return false;
}//BinDataCounts


bool twotargets::BinDataCounts(TString target, Double_t binVal, int noprint){
  
  Int_t iter=-1;
  for (std::vector<Double_t>::iterator it=this->bounds.begin();
       it!=this->bounds.end(); it++, iter++){

    if(binVal <= *it ) {
      if(iter==-1){
	if (noprint) return false;//Don't cout anything
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << "bin value too low!!!!" << std::endl;
	std::cout << this->thisName << " has " <<
	  *it << " < " << binVal <<std::endl;
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << " " << std::endl;
	return false;
      }
      
      if (target=="left_upstream_up") {
	this->left_upstream_up.at(iter)++;
	this->left_upstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_upstream_up") {
	this->right_upstream_up.at(iter)++;
	this->right_upstream.at(iter)++;
	this->right.at(iter)++;}
      else if (target=="left_upstream_down") {
	this->left_upstream_down.at(iter)++;
	this->left_upstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_upstream_down") {
	this->right_upstream_down.at(iter)++;
	this->right_upstream.at(iter)++;
	this->right.at(iter)++;}
      else if (target=="left_downstream_up") {
	this->left_downstream_up.at(iter)++;
	this->left_downstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_downstream_up") {
	this->right_downstream_up.at(iter)++;
	this->right_downstream.at(iter)++;
	this->right.at(iter)++;}
      else if (target=="left_downstream_down") {
	this->left_downstream_down.at(iter)++;
	this->left_downstream.at(iter)++;
	this->left.at(iter)++;}
      else if (target=="right_downstream_down") {
	this->right_downstream_down.at(iter)++;
	this->right_downstream.at(iter)++;
	this->right.at(iter)++;}
      else {
	std::cout << this->thisName << " Wrong target name to " <<
	  "twotargets::BinDataCounts" << std::endl;
	return false;
      }
      
      return true;
    }
  }

  if (noprint) return false;//Don't cout anything
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << this->thisName << " has " <<
    this->bounds.back() << " < " << binVal << std::endl;
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << " " << std::endl;
  return false;
}//BinDataCounts


Bool_t twotargets::BinLeftRight(){
  
  BinOne_LeftRight(this->left, this->right, this->asym, this->e_asym, "asym");
  
  BinOne_LeftRight(this->left_upstream, this->right_upstream,
		   this->asym_upstream, this->e_asym_upstream, "upstream");
  BinOne_LeftRight(this->left_downstream, this->right_downstream,
		   this->asym_downstream, this->e_asym_downstream,"downstream");
  
  BinOne_LeftRight(this->left_upstream_up, this->right_upstream_up,
		   this->asym_upstream_up, this->e_asym_upstream_up,
		   "upstream_up");
  BinOne_LeftRight(this->left_upstream_down, this->right_upstream_down,
		   this->asym_upstream_down, this->e_asym_upstream_down,
		   "upstream_down");
  BinOne_LeftRight(this->left_downstream_up, this->right_downstream_up,
		   this->asym_downstream_up, this->e_asym_downstream_up,
		   "downstream_up");
  BinOne_LeftRight(this->left_downstream_down, this->right_downstream_down,
		   this->asym_downstream_down, this->e_asym_downstream_down,
		   "downstream_down");

  //1 Period asymmetries
  BinOne_LeftRight(this->left_upstream_up, this->right_upstream_down,
		   this->asym_upstream_left, this->e_asym_upstream_left,
		   "upstream_left");
  BinOne_LeftRight(this->left_upstream_down, this->right_upstream_up,
		   this->asym_upstream_right, this->e_asym_upstream_right,
		   "upstream_right");
  BinOne_LeftRight(this->left_downstream_up, this->right_downstream_down,
		   this->asym_downstream_left, this->e_asym_downstream_left,
		   "downstream_left");
  BinOne_LeftRight(this->left_downstream_down, this->right_downstream_up,
		   this->asym_downstream_right, this->e_asym_downstream_right,
		   "downstream_right");

  //sub Period asymmetries
  BinOne_LeftRight(this->left_upstream_up, this->right_downstream_down,
		   this->asym_updown_left, this->e_asym_updown_left,
		   "updown_left");
  BinOne_LeftRight(this->left_upstream_down, this->right_downstream_up,
		   this->asym_downup_right, this->e_asym_downup_right,
		   "downup_right");
  BinOne_LeftRight(this->left_downstream_up, this->right_upstream_down,
		   this->asym_downup_left, this->e_asym_downup_left,
		   "downup_left");
  BinOne_LeftRight(this->left_downstream_down, this->right_upstream_up,
		   this->asym_updown_right, this->e_asym_updown_right,
		   "updown_right");

  return true;
}//BinLeftRight


void twotargets::Print_LR(TString target){

  if (target=="left") CoutLoop(this->left, target);
  else if (target=="right") CoutLoop(this->right, target);
  else if (target=="left_upstream") CoutLoop(this->left_upstream, target);
  else if (target=="right_upstream") CoutLoop(this->right_upstream, target);
  else if (target=="left_downstream") CoutLoop(this->left_downstream, target);
  else if (target=="right_downstream") CoutLoop(this->right_downstream, target);
  else if (target=="left_upstream_up") CoutLoop(this->left_upstream_up, target);
  else if (target=="left_upstream_down") CoutLoop(this->left_upstream_down,
						  target);
  else if (target=="right_upstream_up") CoutLoop(this->right_upstream_up,
						 target);
  else if (target=="right_upstream_down") CoutLoop(this->right_upstream_down,
						   target);
  else if (target=="left_downstream_up") CoutLoop(this->left_downstream_up,
						  target);
  else if (target=="left_downstream_down") CoutLoop(this->left_downstream_down,
						    target);
  else if (target=="right_downstream_up") CoutLoop(this->right_downstream_up,
						   target);
  else if (target=="right_downstream_down")CoutLoop(this->right_downstream_down,
						    target);
  else {
    std::cout << "Wrong option " << target << " to twotargets::Print_LR"
	      << std::endl;
  }
}


void twotargets::Print_Asym(TString asymName){

  if (asymName=="asym") CoutLoop(this->asym, asymName);
  else if (asymName=="asym_upstream") CoutLoop(this->asym_upstream, asymName);
  else if (asymName=="asym_downstream") CoutLoop(this->asym_downstream,
						 asymName);
  else if (asymName=="asym_upstream_up") CoutLoop(this->asym_upstream_up,
						  asymName);
  else if (asymName=="asym_upstream_down") CoutLoop(this->asym_upstream_down,
						    asymName);
  else if (asymName=="asym_downstream_up") CoutLoop(this->asym_downstream_up,
						    asymName);
  else if (asymName=="asym_downstream_down")CoutLoop(this->asym_downstream_down,
						     asymName);
  else if (asymName=="asym_upstream_left")CoutLoop(this->asym_upstream_left,
						     asymName);
  else if (asymName=="e_asym") CoutLoop(this->e_asym, asymName);
  else if (asymName=="e_asym_upstream") CoutLoop(this->e_asym_upstream, asymName);
  else if (asymName=="e_asym_downstream") CoutLoop(this->e_asym_downstream,
						 asymName);
  else if (asymName=="e_asym_upstream_up") CoutLoop(this->e_asym_upstream_up,
						  asymName);
  else if (asymName=="e_asym_upstream_down") CoutLoop(this->e_asym_upstream_down,
						    asymName);
  else if (asymName=="e_asym_downstream_up") CoutLoop(this->e_asym_downstream_up,
						    asymName);
  else if (asymName=="e_asym_downstream_down")CoutLoop(this->e_asym_downstream_down,
						     asymName);
  else if (asymName=="e_asym_upstream_left")CoutLoop(this->e_asym_upstream_left,
						     asymName);
  else {
    std::cout << "Wrong option " << asymName << " to twotargets::Print_Asym"
	      << std::endl;
  }
}


//Helper functions
void twotargets::setZero(Int_t nBins){
  for (Int_t i=0; i<nBins; i++) {
    this->left_upstream.push_back(0);
    this->right_upstream.push_back(0);
    this->left_downstream.push_back(0);
    this->right_downstream.push_back(0);

    this->left_upstream_up.push_back(0);
    this->right_upstream_up.push_back(0);
    this->left_upstream_down.push_back(0);
    this->right_upstream_down.push_back(0);

    this->left_downstream_up.push_back(0);
    this->right_downstream_up.push_back(0);
    this->left_downstream_down.push_back(0);
    this->right_downstream_down.push_back(0);
  }
}
