#include "binDist4Targ.h"

binDist4Targ::binDist4Targ(TString name, TString binVar, TString binfile,
		 Int_t hbins, Double_t xmn, Double_t xmx) :
  binDist(name, binVar, binfile, hbins, xmn, xmx){
  
  this->nBins = this->bounds.size()-1;
  for (Int_t i=0; i<this->nBins; i++) {
    DefineOneHist(Form("%s_left_upSup_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upSup_upP);
    DefineOneHist(Form("%s_left_upSdown_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upSdown_upP);

    DefineOneHist(Form("%s_left_upSup_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upSup_downP);
    DefineOneHist(Form("%s_left_upSdown_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_upSdown_downP);

    DefineOneHist(Form("%s_right_upSup_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upSup_upP);
    DefineOneHist(Form("%s_right_upSdown_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upSdown_upP);

    DefineOneHist(Form("%s_right_upSup_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upSup_downP);
    DefineOneHist(Form("%s_right_upSdown_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_upSdown_downP);

    DefineOneHist(Form("%s_left_downSup_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downSup_upP);
    DefineOneHist(Form("%s_left_downSdown_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downSdown_upP);

    DefineOneHist(Form("%s_left_downSup_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downSup_downP);
    DefineOneHist(Form("%s_left_downSdown_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->left_downSdown_downP);

    DefineOneHist(Form("%s_right_downSup_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downSup_upP);
    DefineOneHist(Form("%s_right_downSdown_upP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downSdown_upP);

    DefineOneHist(Form("%s_right_downSup_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downSup_downP);
    DefineOneHist(Form("%s_right_downSdown_downP_%s%i", name.Data(),
		       binVar.Data(), i), this->right_downSdown_downP);
  }
}//binDist4Targ


Bool_t binDist4Targ::BinFill(Int_t target, Bool_t Left, Double_t Spin,
			     Double_t distVar, Double_t binVar){
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

      TString targetName;
      if (target==0){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_upstream_up";
	    left_upSup_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "left_upstream_down";
	    left_upSup_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_upstream_up";
	    right_upSup_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "right_upstream_down";
	    right_upSup_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Right
      }//Target 0
      else if (target==1){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_upstream_up";
	    left_upSdown_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "left_upstream_down";
	    left_upSdown_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_upstream_up";
	    right_upSdown_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "right_upstream_down";
	    right_upSdown_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Right
      }//Target 1
      else if (target==2){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_downstream_up";
	    left_downSup_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "left_downstream_down";
	    left_downSup_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_downstream_up";
	    right_downSup_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "right_downstream_down";
	    right_downSup_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Right
      }//Target 2
      else if (target==3){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_downstream_up";
	    left_downSdown_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "left_downstream_down";
	    left_downSdown_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_downstream_up";
	    right_downSdown_upP.at(iter)->Fill(distVar);
	  }//Spin up
	  else{//Spin down
	    targetName = "right_downstream_down";
	    right_downSdown_downP.at(iter)->Fill(distVar);
	  }//Spin down
	}//Right
      }//Target 3

      binDist::BinFill(targetName, distVar, iter);
      
      return true;
    }
  }

  std::cout << "!!!!!!!!!!!!!!!\n" <<  "bin value too high!!!!\n" ;
  std::cout << this->thisName << " has " <<
    this->bounds.back() << " < " << binVar << std::endl;
  std::cout << "!!!!!!!!!!!!!!!\n\n";
  return false;
}//BinFill


void binDist4Targ::Write(){
  
  binDist::Write();
  
  for (Int_t i=0; i<this->nBins; i++) {
    this->left_upSup_upP.at(i)->Write();
    this->right_upSup_upP.at(i)->Write();
    this->left_upSup_downP.at(i)->Write();
    this->right_upSup_downP.at(i)->Write();

    this->left_upSdown_upP.at(i)->Write();
    this->right_upSdown_upP.at(i)->Write();
    this->left_upSdown_downP.at(i)->Write();
    this->right_upSdown_downP.at(i)->Write();

    this->left_downSup_upP.at(i)->Write();
    this->right_downSup_upP.at(i)->Write();
    this->left_downSup_downP.at(i)->Write();
    this->right_downSup_downP.at(i)->Write();

    this->left_downSdown_upP.at(i)->Write();
    this->right_downSdown_upP.at(i)->Write();
    this->left_downSdown_downP.at(i)->Write();
    this->right_downSdown_downP.at(i)->Write();
  }
}//Write


//Helper functions
void binDist4Targ::DefineOneHist(TString hName, std::vector<TH1D*> &vect){
  TH1D *h = new TH1D(hName, hName, this->hbins, this->xMin, this->xMax);
  
  vect.push_back(h);
}
