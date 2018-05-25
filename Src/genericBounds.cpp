#include "genericBounds.h"

genericBounds::genericBounds(int nBins, TString name) :
  lrSpinCorr(nBins, name){

  for (int i=0; i<nBins; i++) {
    this->upstream_xval.push_back(0.0);
    this->downstream_xval.push_back(0.0);
  }
}


genericBounds::genericBounds (int nBins, std::vector<Double_t> &in_bounds,
			      std::vector<Double_t> &in_xval,
			      TString thisName) :
  lrSpinCorr(nBins, in_bounds, in_xval, thisName){

  for (std::vector<double>::iterator it=in_xval.begin(); it!=in_xval.end();
       it++){
    this->upstream_xval.push_back(*it);
    this->downstream_xval.push_back(*it);
  }

  for (std::vector<double>::iterator it=in_bounds.begin(); it!=in_bounds.end();
       it++){
    this->upstream_bounds.push_back(*it);
    this->downstream_bounds.push_back(*it);
  }
}


genericBounds::genericBounds (TString binfile, TString type,
			      TString thisName) :
  lrSpinCorr(binfile, type, thisName){

  for (std::vector<double>::iterator it=this->xval.begin();
       it!=this->xval.end(); it++){
    this->upstream_xval.push_back(*it);
    this->downstream_xval.push_back(*it);
  }

  for (std::vector<double>::iterator it=this->bounds.begin();
       it!=this->bounds.end(); it++){
    this->upstream_bounds.push_back(*it);
    this->downstream_bounds.push_back(*it);
  }
}


genericBounds::genericBounds (TFile *f1, TString type,TString thisName)
  : lrSpinCorr(f1, type, thisName){

    for (std::vector<double>::iterator it=this->xval.begin();
       it!=this->xval.end(); it++){
    this->upstream_xval.push_back(*it);
    this->downstream_xval.push_back(*it);
  }

  for (std::vector<double>::iterator it=this->bounds.begin();
       it!=this->bounds.end(); it++){
    this->upstream_bounds.push_back(*it);
    this->downstream_bounds.push_back(*it);
  }
}


bool genericBounds::SetCorrBounds(TString target,
				  double binVal, double avgVal){

  if (target=="pol_upstream_up") {
    binOneAvgBounds(this->upstream_bounds, this->avgPol_upstream_up,
		    this->avgPol_count_upstream_up, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgPol_upstream,
		    this->avgPol_count_upstream, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgPol, this->avgPol_count,
		    binVal, avgVal);
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="pol_upstream_down") {
    binOneAvgBounds(this->upstream_bounds, this->avgPol_upstream_down,
		    this->avgPol_count_upstream_down, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgPol_upstream,
		    this->avgPol_count_upstream, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgPol, this->avgPol_count,
		    binVal, avgVal);
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="pol_downstream_up") {
    binOneAvgBounds(this->downstream_bounds, this->avgPol_downstream_up,
		    this->avgPol_count_downstream_up, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgPol_downstream,
		    this->avgPol_count_downstream, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgPol, this->avgPol_count,
		    binVal, avgVal);
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="pol_downstream_down") {
    binOneAvgBounds(this->downstream_bounds, this->avgPol_downstream_down,
		    this->avgPol_count_downstream_down, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgPol_downstream,
		    this->avgPol_count_downstream, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgPol, this->avgPol_count,
		    binVal, avgVal);
    pol_int += avgVal; pol_int_count++;
  }
  else if (target=="dil_upstream_up") {
    binOneAvgBounds(this->upstream_bounds, this->avgDil_upstream_up,
		    this->avgDil_count_upstream_up, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgDil_upstream,
		    this->avgDil_count_upstream, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgDil, this->avgDil_count,
		    binVal, avgVal);
    dil_int += avgVal; dil_int_count++;
  }
  else if (target=="dil_upstream_down") {
    binOneAvgBounds(this->upstream_bounds, this->avgDil_upstream_down,
		    this->avgDil_count_upstream_down, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgDil_upstream,
		    this->avgDil_count_upstream, binVal, avgVal);
    binOneAvgBounds(this->upstream_bounds, this->avgDil, this->avgDil_count,
		    binVal, avgVal);
    dil_int += avgVal; dil_int_count++;
  }
  else if (target=="dil_downstream_up") {
    binOneAvgBounds(this->downstream_bounds, this->avgDil_downstream_up,
		    this->avgDil_count_downstream_up, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgDil_downstream,
		    this->avgDil_count_downstream, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgDil, this->avgDil_count,
		    binVal, avgVal);
    dil_int += avgVal; dil_int_count++;
  }
  else if (target=="dil_downstream_down") {
    binOneAvgBounds(this->downstream_bounds, this->avgDil_downstream_down,
		    this->avgDil_count_downstream_down, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgDil_downstream,
		    this->avgDil_count_downstream, binVal, avgVal);
    binOneAvgBounds(this->downstream_bounds, this->avgDil, this->avgDil_count,
		    binVal, avgVal);
    dil_int += avgVal; dil_int_count++;
  }
  else {
	std::cout << this->thisName << " Wrong target name to " <<
	  "lrSpinCorr::SetCorr" << std::endl;
	return false;
  }

  return true;
}


void genericBounds::Print_Bounds(TString target){
  if (target=="upstream") CoutLoop(this->upstream_bounds, target);
  else if (target=="downstream") CoutLoop(this->downstream_bounds, target);
  else{
    std::cout << "Invald target: " << target <<
      "  to genericBounds::Print_Bounds" << std::endl;
    exit(EXIT_FAILURE);
  }
}


void genericBounds::Print_xVal(TString target){
  if (target=="upstream") CoutLoop(this->upstream_xval, target);
  else if (target=="downstream") CoutLoop(this->downstream_xval, target);
  else{
    std::cout << "Invald target: " << target <<
      "  to genericBounds::Print_xVal" << std::endl;
    exit(EXIT_FAILURE);
  }
}


void genericBounds::SetUpBounds(TString target, double boundVal){
  if (target=="upstream") this->upstream_bounds.push_back(boundVal);
  else if (target=="downstream") this->downstream_bounds.push_back(boundVal);
  else {
    std::cout << "Wrong target: " << target
	      << "   in genericBounds::SetUpBounds from: " << this->thisName
	      << std::endl;
    exit(EXIT_FAILURE);
  }
}


void genericBounds::MakeBounds(){
  if (this->upstream_bounds.size()==0 || this->downstream_bounds.size()==0 ){
    std::cout <<"Nothing filled in bound vectors in genericBounds"<< std::endl;
    exit(EXIT_FAILURE);
  }

  std::sort(this->upstream_bounds.begin(), upstream_bounds.end() );
  std::sort(this->downstream_bounds.begin(), downstream_bounds.end() );

  std::vector<double> tmp_upstream, tmp_downstream;
  Int_t upstream_size = upstream_bounds.size();
  Int_t downstream_size = downstream_bounds.size();
  for (Int_t i=0; i<this->nBins+1; i++) {
    if (i) {
      tmp_upstream.push_back(upstream_bounds.at( 1.0*i*upstream_size/this->nBins-1) );
      tmp_downstream.push_back(downstream_bounds.at( 1.0*i*downstream_size/this->nBins-1) );
    }
    else{
      tmp_upstream.push_back(upstream_bounds.at(0) );
      tmp_downstream.push_back(downstream_bounds.at(0) );
    }
  }

  avgOneBound(this->upstream_bounds, this->upstream_xval);
  avgOneBound(this->downstream_bounds, this->downstream_xval);

  this->upstream_bounds.clear(); this->bounds.clear();
  this->downstream_bounds.clear();

  for (std::vector<double>::iterator it=tmp_upstream.begin();
       it!=tmp_upstream.end(); it++) {
    this->upstream_bounds.push_back(*it);
    this->bounds.push_back(*it);
  }
  for (std::vector<double>::iterator it=tmp_downstream.begin();
       it!=tmp_downstream.end(); it++) this->downstream_bounds.push_back(*it);
}


bool genericBounds::BinDataCounts(int whichTarget, TString target,
				  Double_t binVal, int noprint){

  if ((this->upstream_bounds.size() != this->nBins+1)
      || (this->upstream_bounds.size() != this->downstream_bounds.size()) ){
    std::cout << "Size missmatch in genericBounds::BinDataCounts" << std::endl;
    exit(EXIT_FAILURE);
  }
  
  std::vector<Double_t>::iterator it;
  if (whichTarget) it = this->downstream_bounds.begin();
  else it = this->upstream_bounds.begin();
  for (int iter=0; iter<this->nBins+1; iter++, it++) {

    if(binVal <= *it ) {
      if(iter==0 && binVal<*it){
	if (noprint) return false;//Don't cout anything
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << "bin value too low!!!!" << std::endl;
	std::cout << this->thisName << " has " <<
	  *it << " < " << binVal <<std::endl;
	std::cout << "!!!!!!!!!!!!!!!" << std::endl;
	std::cout << " " << std::endl;
	return false;
      }
      else if (iter==0) iter = 1;
      
      if (target=="left_upstream_up") {
	this->left_upstream_up.at(iter-1)++;
	this->left_upstream.at(iter-1)++;
	this->left.at(iter-1)++;
      }
      else if (target=="right_upstream_up") {
	this->right_upstream_up.at(iter-1)++;
	this->right_upstream.at(iter-1)++;
	this->right.at(iter-1)++;}
      else if (target=="left_upstream_down") {
	this->left_upstream_down.at(iter-1)++;
	this->left_upstream.at(iter-1)++;
	this->left.at(iter-1)++;}
      else if (target=="right_upstream_down") {
	this->right_upstream_down.at(iter-1)++;
	this->right_upstream.at(iter-1)++;
	this->right.at(iter-1)++;}
      else if (target=="left_downstream_up") {
	this->left_downstream_up.at(iter-1)++;
	this->left_downstream.at(iter-1)++;
	this->left.at(iter-1)++;}
      else if (target=="right_downstream_up") {
	this->right_downstream_up.at(iter-1)++;
	this->right_downstream.at(iter-1)++;
	this->right.at(iter-1)++;}
      else if (target=="left_downstream_down") {
	this->left_downstream_down.at(iter-1)++;
	this->left_downstream.at(iter-1)++;
	this->left.at(iter-1)++;}
      else if (target=="right_downstream_down") {
	this->right_downstream_down.at(iter-1)++;
	this->right_downstream.at(iter-1)++;
	this->right.at(iter-1)++;}
      else {
	std::cout << this->thisName << " Wrong target name to " <<
	  "genericBounds::BinDataCounts" << std::endl;
	return false;
      }
      
      return true;
    }
  }

  if (noprint) return false;//Don't cout anything
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << this->thisName << " has binVal: " << binVal << std::endl;
  std::cout << " " << std::endl;
  return false;
}//BinDataCounts


//Helper functions
void genericBounds::avgOneBound(std::vector<double> &bounds,
				std::vector<double> &xval){

  std::vector<int> tmp_counts;
  for (Int_t i=0; i<this->nBins; i++) tmp_counts.push_back(0);

  int size = bounds.size();
  int ibin=0;
  for (std::vector<Double_t>::iterator it=bounds.begin(); it!=bounds.end();
       it++, ibin++) {

    if (ibin/(size/this->nBins) < this->nBins){
      xval[ibin/(size/this->nBins)] += *it;
      tmp_counts[ibin/(size/this->nBins)]++;
    }
    else{
      xval[this->nBins-1] += *it;
      tmp_counts[this->nBins-1]++;
    }
  }

  for (int i=0; i<this->nBins; i++) xval[i] = xval[i]/(1.0*tmp_counts[i]);
}


bool genericBounds::binOneAvgBounds(std::vector<double> &bounds,
				    std::vector<double> &Avg,
				    std::vector<int> &count, double binVal,
				    double avgVal){

  if ( (Avg.size()!=bounds.size()-1) || (Avg.size()!=count.size()) ) {
    std::cout << "Size missmatch in genericBounds::binOneAvgBounds from: " <<
      this->thisName << std::endl;
    exit(EXIT_FAILURE);
  }


  int iter = -1;
  for (std::vector<double>::iterator it=bounds.begin(); it!=bounds.end();
       it++, iter++){

    if(iter == -1 && binVal < *it ) {
      std::cout << "bin value too low!!!!" << std::endl;
      std::cout << "Lower bound: " << *it << "    val: "
		<< binVal << "   from: " << this->thisName << std::endl;
      std::cout << " " << std::endl;

      return false;
    }
    else if (binVal <= *(it)){
      if (iter==-1) iter=0;
      
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
