#include "binDist.h"

binDist::binDist(TString name, TString binVar, TString binfile,
		 Int_t hbins, Double_t xmn, Double_t xmx){
  this->thisName = name + "_" + binVar;
  this->hbins = hbins;
  this->xMin = xmn;
  this->xMax = xmx;
  this->varBinning =false;

  GetBounds(binVar, binfile);//Get bounds from file
  
  this->nBins = this->bounds.size()-1;
  for (Int_t i=0; i<this->nBins; i++) {
    DefineOneHist(Form("%s_left_upstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upstream_up);
    DefineOneHist(Form("%s_right_upstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upstream_up);
    DefineOneHist(Form("%s_left_upstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upstream_down);
    DefineOneHist(Form("%s_right_upstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upstream_down);

    DefineOneHist(Form("%s_left_downstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downstream_up);
    DefineOneHist(Form("%s_right_downstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downstream_up);
    DefineOneHist(Form("%s_left_downstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downstream_down);
    DefineOneHist(Form("%s_right_downstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downstream_down);
  }
}//binDist


binDist::binDist(TString name, TString binVar, TString binfile,
		 Int_t hbins, Double_t *bins){
  this->thisName = name + "_" + binVar;
  this->hbins = hbins;
  this->xMin = bins[0];
  this->xMax = bins[hbins];
  this->varBinning =true;

  GetBounds(binVar, binfile);//Get bounds from file
  
  this->nBins = this->bounds.size()-1;
  for (Int_t i=0; i<this->nBins; i++) {
    DefineOneHist(Form("%s_left_upstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upstream_up, bins);
    DefineOneHist(Form("%s_right_upstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upstream_up, bins);
    DefineOneHist(Form("%s_left_upstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upstream_down, bins);
    DefineOneHist(Form("%s_right_upstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upstream_down, bins);

    DefineOneHist(Form("%s_left_downstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downstream_up, bins);
    DefineOneHist(Form("%s_right_downstream_up_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downstream_up, bins);
    DefineOneHist(Form("%s_left_downstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downstream_down, bins);
    DefineOneHist(Form("%s_right_downstream_down_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downstream_down, bins);
  }
}//binDist


Bool_t binDist::BinFill(TString target, Double_t distVar, Double_t binVar){
  Int_t iter=-1;
  for (std::vector<Double_t>::iterator it=this->bounds.begin();
       it!=this->bounds.end(); it++, iter++){

    if(binVar <= *it ) {
      if(iter==-1){
	std::cout << "!!!!!!!!!!!!!!!\n" << "bin value too low!!!!\n";
	std::cout << this->thisName << " has " <<
	  *it << " < " << binVar <<std::endl;
	std::cout << "!!!!!!!!!!!!!!!\n\n";
	return false;
      }
      
      if (target=="left_upstream_up") 
	this->left_upstream_up.at(iter)->Fill(distVar);
      else if (target=="right_upstream_up") 
	this->right_upstream_up.at(iter)->Fill(distVar);
      else if (target=="left_upstream_down")
	this->left_upstream_down.at(iter)->Fill(distVar);
      else if (target=="right_upstream_down") 
	this->right_upstream_down.at(iter)->Fill(distVar);
      else if (target=="left_downstream_up")
	this->left_downstream_up.at(iter)->Fill(distVar);
      else if (target=="right_downstream_up")
	this->right_downstream_up.at(iter)->Fill(distVar);
      else if (target=="left_downstream_down")
	this->left_downstream_down.at(iter)->Fill(distVar);
      else if (target=="right_downstream_down")
	this->right_downstream_down.at(iter)->Fill(distVar);
      else {
	std::cout << this->thisName << " Wrong target name to " <<
	  "twotargets::BinDataCounts" << std::endl;
	return false;
      }
      
      return true;
    }
  }

  std::cout << "!!!!!!!!!!!!!!!\n" <<  "bin value too high!!!!\n" ;
  std::cout << this->thisName << " has " <<
    this->bounds.back() << " < " << binVar << std::endl;
  std::cout << "!!!!!!!!!!!!!!!\n\n";
  return false;
  
}//BinFill



Bool_t binDist::BinFill(TString target, Double_t distVar, Int_t bin){

  if (bin < 0 || bin > this->nBins){
    std::cout << "binDist::BinFill bin too large or too small" << std::endl;
    exit(EXIT_FAILURE);
    return false;
  }
  
  if (target=="left_upstream_up") 
    this->left_upstream_up.at(bin)->Fill(distVar);
  else if (target=="right_upstream_up") 
    this->right_upstream_up.at(bin)->Fill(distVar);
  else if (target=="left_upstream_down")
    this->left_upstream_down.at(bin)->Fill(distVar);
  else if (target=="right_upstream_down") 
    this->right_upstream_down.at(bin)->Fill(distVar);
  else if (target=="left_downstream_up")
    this->left_downstream_up.at(bin)->Fill(distVar);
  else if (target=="right_downstream_up")
    this->right_downstream_up.at(bin)->Fill(distVar);
  else if (target=="left_downstream_down")
    this->left_downstream_down.at(bin)->Fill(distVar);
  else if (target=="right_downstream_down")
    this->right_downstream_down.at(bin)->Fill(distVar);
  else {
    std::cout << this->thisName << " Wrong target name to " <<
      "twotargets::BinDataCounts" << std::endl;
    return false;
  }
      
  return true;
    
}//BinFill


void binDist::Print(){
  std::cout << "\nPrinting information for: " << this->thisName;
  std::cout << "\nhbins = " << this->hbins << "  xMin = " << this->xMin
	    << "   xMax = " << this->xMax << std::endl;
  std::cout << "Binned in physics nBins = " << nBins << std::endl;
  std::cout << "Has variable histogram binning:  " << this->varBinning <<"\n\n";
}


void binDist::Write(){
  for (Int_t i=0; i<this->nBins; i++) {
    this->left_upstream_up.at(i)->Write();
    this->right_upstream_up.at(i)->Write();
    this->left_upstream_down.at(i)->Write();
    this->right_upstream_down.at(i)->Write();

    this->left_downstream_up.at(i)->Write();
    this->right_downstream_up.at(i)->Write();
    this->left_downstream_down.at(i)->Write();
    this->right_downstream_down.at(i)->Write();
  }
}//Write


//Helper functions
void binDist::GetBounds(TString binVar, TString binfile){
  if (binVar=="xN"||binVar=="xPi"||binVar=="pT"||binVar=="mass"||binVar=="rad"){
    this->bounds.push_back(0.0);}
  else if (binVar=="xF") this->bounds.push_back(-1.0);
  else if (binVar=="vxZ_upstream") this->bounds.push_back(-294.5);
  else if (binVar=="vxZ_downstream") this->bounds.push_back(-219.5);
  else {
    std::cout << "Invalid binVar: " << binVar << " in leftright::leftright"
	      << std::endl;
    exit(EXIT_FAILURE);
  }

  TString boundaries = binVar; boundaries += " bin boundaries";
  TString averages = binVar; averages += " bin averages";
  
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
      break;
    }
    else if (!found) continue;

    this->bounds.push_back(atof(line.c_str() ) );
  }

  if (binVar=="xN"||binVar=="xPi"||binVar=="xF") this->bounds.push_back(1.0);
  else if (binVar=="pT") this->bounds.push_back(5.0);
  else if (binVar=="mass") this->bounds.push_back(12.0);
  else if (binVar=="rad") this->bounds.push_back(2.0);
  else if (binVar=="vxZ_upstream") this->bounds.push_back(-239.3);
  else if (binVar=="vxZ_downstream") this->bounds.push_back(-164.3);

  if ( this->bounds.size()==0 ){
    std::cout << "bounds did not fill well in binDist::binDist   from: " <<
      binVar << std::endl;
    exit(EXIT_FAILURE);
  }  
}


void binDist::DefineOneHist(TString hName, std::vector<TH1D*> &vect){
  TH1D *h = new TH1D(hName, hName, this->hbins, this->xMin, this->xMax);
  
  vect.push_back(h);
}


void binDist::DefineOneHist(TString hName, std::vector<TH1D*> &vect,
			    Double_t *bins){
  TH1D *h = new TH1D(hName, hName, this->hbins, bins);
  
  vect.push_back(h);
}
