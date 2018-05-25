#include "leftright.h"

leftright::leftright (){
  this->nBins = 0;
  this->thisName = "noName";
}


leftright::leftright (int nBins, TString thisName){
  for (Int_t i=0; i<nBins; i++) {
    left.push_back(0); right.push_back(0);
    
    bounds.push_back(0.0);
    xval.push_back(0.0);
  }
  bounds.push_back(0.0);

  this->nBins = nBins;
  this->thisName = thisName;
}


leftright::leftright (int nBins, std::vector<Double_t> &in_bounds,
		      std::vector<Double_t> &in_xval, TString thisName){
  for (Int_t i=0; i<nBins; i++) {
    left.push_back(0); right.push_back(0);
  }

  for (std::vector<Double_t>::iterator it=in_bounds.begin();it!=in_bounds.end();
       it++) this->bounds.push_back(*it);
  for (std::vector<Double_t>::iterator it=in_xval.begin();it!=in_xval.end();
       it++) this->xval.push_back(*it);

  this->nBins = nBins;
  this->thisName = thisName;
}


leftright::leftright (TString binfile, TString type, TString thisName){
  if (type=="xN"||type=="xPi"||type=="pT"||type=="mass"||type=="rad") {
    bounds.push_back(0.0);}
  else if (type=="xF") bounds.push_back(-1.0);
  else if (type=="vxZ_upstream") bounds.push_back(-294.5);
  else if (type=="vxZ_downstream") bounds.push_back(-219.5);
  else {
    std::cout << "Invalid type: " << type << " in leftright::leftright"
	      << std::endl;
    exit(EXIT_FAILURE);
  }

  TString boundaries = type; boundaries += " bin boundaries";
  TString averages = type; averages += " bin averages";
  
  std::string line;
  bool first = false, found = false;
  std::ifstream f_bins(binfile);
  if(!f_bins.is_open() ) {
    std::cout << " " << std::endl;
    std::cout << "binFile did not open from: " << thisName << std::endl;
    exit(EXIT_FAILURE); }
  while (!f_bins.eof()) {
    getline(f_bins,line);
    TString tline (line);

    if (tline == boundaries) {
      found = true;
      first = true;
      continue;
    }
    else if (tline == averages) {
      first = false;
      continue;
    }
    else if (!found) continue;

    if (first) this->bounds.push_back(atof(line.c_str() ) );
    else {
      if (this->xval.size() >= this->bounds.size()) break;
      this->xval.push_back(atof(line.c_str() ) );
    }
  }

  if (type=="xN"||type=="xPi"||type=="xF") bounds.push_back(1.0);
  else if (type=="pT") bounds.push_back(5.0);
  else if (type=="mass") bounds.push_back(12.0);
  else if (type=="rad") bounds.push_back(2.0);
  else if (type=="vxZ_upstream") bounds.push_back(-239.3);
  else if (type=="vxZ_downstream") bounds.push_back(-164.3);

  if ( (this->xval.size()==0)|| (this->xval.size()!=this->bounds.size()-1) ){
    std::cout << "xval did not fill well in leftright::leftright   from: " <<
	      type << std::endl;
    exit(EXIT_FAILURE);
  }

  this->nBins = this->xval.size();
  if (thisName == "noName") this->thisName = type;
  else this->thisName = thisName;

  for (Int_t i=0; i<this->nBins; i++) {
    left.push_back(0); right.push_back(0);
  }
}


leftright::leftright (TFile *f1, TString type, TString thisName){
  
  if (type=="xN") setOneBounds(f1, "tv_xN_bounds", "tv_xN_xval");
  else if (type=="xPi") setOneBounds(f1, "tv_xPi_bounds", "tv_xPi_xval");
  else if (type=="xF") setOneBounds(f1, "tv_xF_bounds", "tv_xF_xval");
  else if (type=="pT") setOneBounds(f1, "tv_pT_bounds", "tv_pT_xval");
  else if (type=="mass") setOneBounds(f1, "tv_M_bounds", "tv_M_xval");
  else{
    std::cout << "Option: " << type << " is not valid for leftright::leftright"
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  this->nBins = this->xval.size();
  if (thisName == "noName") this->thisName = type;
  else this->thisName = thisName;

  for (Int_t i=0; i<nBins; i++) {
    left.push_back(0); right.push_back(0);
  }
}


//Helper functions
void leftright::CoutLoop(std::vector<unsigned long long> &counts, TString name){
  std::cout << "Printing content for: " << name
	    << "  from: " << this->thisName << std::endl;
  
  for (std::vector<unsigned long long>::iterator it=counts.begin();
       it!=counts.end(); it++){
    std::cout << *it << std::endl;
  }

  std::cout << " " << std::endl;
}


void leftright::CoutLoop(std::vector<double> &counts, TString name){
  std::cout << "Printing content for: " << name 
	    << "  from: " << this->thisName << std::endl;
  
  for (std::vector<double>::iterator it=counts.begin();
       it!=counts.end(); it++){
    std::cout << *it << std::endl;
  }

  std::cout << " " << std::endl;
}


double leftright::AsymmetryError(long long A, long long B){
  //if (TMath::Abs(4.0*A*B/( (A+B)*(A+B)*(A+B) ) ) < 10e-6 ){
  //  std::cout << "Errors approximated as Sqrt(1/N)" << std::endl;
  //  return TMath::Sqrt( 1.0/(A+B) );
  //}//might be needed for generated errors

  return TMath::Sqrt( 4.0*A*B/( (A+B)*(A+B)*(A+B) ) );
}


bool leftright::BinOne_LeftRight(std::vector<unsigned long long> &left,
				 std::vector<unsigned long long> &right,
				 std::vector<double> &asym,
				 std::vector<double> &e_asym, TString name){

  for (Int_t i=0; i<this->nBins; i++) {
    unsigned long long l = left.at(i);
    unsigned long long r = right.at(i);

    long long diff = l - r;
    long long sum = l + r;
    if (!sum) {
      asym.push_back(-10.0);//no data
      e_asym.push_back(0.0);

      std::cout<<"leftright::BinOne_LeftRight: no data for l/r asymmetry in: "
	       << name << "    from: " << this->thisName
	       << "     bin: " << i << std::endl;
    }
    else{
      asym.push_back( 1.0*diff/(1.0*sum) );
      e_asym.push_back( AsymmetryError(l, r) );
    }
  }
  
  return true;
}//BinnedLeftRight


bool leftright::MakeAvg(std::vector<double_t> &avg, std::vector<int> &count){
  if (avg.size() != count.size()  || avg.size()==0){
    std::cout << this->thisName << ": leftright::MakeAvg size mismatch"
	      << " or zero vector size" << std::endl;
    return false;
  }
  for (Int_t i=0; i<this->nBins; i++) avg.at(i) = avg.at(i)/count.at(i);
  
  return true;
}


void leftright::setOneBounds(TFile *f1, TString bound, TString xval){

  TVectorD tmp_bound = *( (TVectorD*)f1->Get(bound) );
  for (Int_t i=0; i<tmp_bound.GetNrows(); i++) {
    this->bounds.push_back(tmp_bound[i]);}

  TVectorD tmp_xval = *( (TVectorD*)f1->Get(xval) );
  for (Int_t i=0; i<tmp_xval.GetNrows(); i++) {
    this->xval.push_back(tmp_xval[i]);}
};
