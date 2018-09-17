#include "lrSpinCorr.h"
#include "fourTargets.h"


FourTargets::FourTargets (TString binfile, TString type, TString thisName)
  : lrSpinCorr(binfile, type, thisName){

  FourTargets::setZero(nBins);
}


bool FourTargets::BinDataCounts(Int_t target, Bool_t Left, Double_t Spin,
				Double_t binVal, int noprint){
  
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
      
      TString targetName;
      if (target==0){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_upstream_up";
	    left_upSup_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "left_upstream_down";
	    left_upSup_downP.at(iter)++;
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_upstream_up";
	    right_upSup_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "right_upstream_down";
	    right_upSup_downP.at(iter)++;
	  }//Spin down
	}//Right
      }//Target 0
      else if (target==1){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_upstream_up";
	    left_upSdown_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "left_upstream_down";
	    left_upSdown_downP.at(iter)++;
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_upstream_up";
	    right_upSdown_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "right_upstream_down";
	    right_upSdown_downP.at(iter)++;
	  }//Spin down
	}//Right
      }//Target 1
      else if (target==2){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_downstream_up";
	    left_downSup_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "left_downstream_down";
	    left_downSup_downP.at(iter)++;
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_downstream_up";
	    right_downSup_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "right_downstream_down";
	    right_downSup_downP.at(iter)++;
	  }//Spin down
	}//Right
      }//Target 2
      else if (target==3){
	if (Left){
	  if (Spin > 0){
	    targetName = "left_downstream_up";
	    left_downSdown_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "left_downstream_down";
	    left_downSdown_downP.at(iter)++;
	  }//Spin down
	}//Left
	else{//Right
	  if (Spin > 0){
	    targetName = "right_downstream_up";
	    right_downSdown_upP.at(iter)++;
	  }//Spin up
	  else{//Spin down
	    targetName = "right_downstream_down";
	    right_downSdown_downP.at(iter)++;
	  }//Spin down
	}//Right
      }//Target 3

      lrSpinCorr::BinDataCounts(targetName, iter);

      return true;
    }
  }

  if (noprint) return false;//Don't cout anything
  std::cout << "!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << this->thisName << " has " <<
    this->bounds.back() << " < " << binVal << std::endl;
  std::cout << "!!!!!!!!!!!!!!!\n" << std::endl;
  return false;
}//BinDataCounts


bool FourTargets::SetCorr(double binVal, double avgPol, double avgDil,
			  Int_t target, Bool_t Left, Double_t spin){

  TString pol_targetName, dil_targetName;
  if (target==0){
    binOneAvg(this->avgPolDil_upSup, this->avgPolDil_count_upSup, binVal,
	      avgPol*avgDil);

    if (spin > 0){
      pol_targetName = "pol_upstream_up";
      dil_targetName = "dil_upstream_up";
    }
    else{//Spin down
      pol_targetName = "pol_upstream_down";
      dil_targetName = "dil_upstream_down";
    }
  }//Target 0
  else if (target==1){
    binOneAvg(this->avgPolDil_upSdown, this->avgPolDil_count_upSdown, binVal,
	      avgPol*avgDil);
    
    if (spin > 0){
      pol_targetName = "pol_upstream_up";
      dil_targetName = "dil_upstream_up";
    }
    else{//Spin down
      pol_targetName = "pol_upstream_down";
      dil_targetName = "dil_upstream_down";
    }
  }//Target 1
  else if (target==2){
    binOneAvg(this->avgPolDil_downSup, this->avgPolDil_count_downSup, binVal,
	      avgPol*avgDil);
    
    if (spin > 0){
      pol_targetName = "pol_downstream_up";
      dil_targetName = "dil_downstream_up";
    }
    else{//Spin down
      pol_targetName = "pol_downstream_down";
      dil_targetName = "dil_downstream_down";
    }
  }//Target 2
  else if (target==3){
    binOneAvg(this->avgPolDil_downSdown, this->avgPolDil_count_downSdown,
	      binVal, avgPol*avgDil);
    
    if (spin > 0){
      pol_targetName = "pol_downstream_up";
      dil_targetName = "dil_downstream_up";
    }
    else{//Spin down
      pol_targetName = "pol_downstream_down";
      dil_targetName = "dil_downstream_down";
    }
  }//Target 3

  lrSpinCorr::SetCorr(pol_targetName, binVal, avgPol, Left);
  lrSpinCorr::SetCorr(dil_targetName, binVal, avgDil, Left);
      
  return true;
}//SetCorr


void FourTargets::AvgCorr(){
  makeOneAvg(this->avgPolDil_upSup, this->avgPolDil_count_upSup, "avgPolDil_upSup");
  makeOneAvg(this->avgPolDil_upSdown, this->avgPolDil_count_upSdown,
	     "avgPolDil_upSdown");
  makeOneAvg(this->avgPolDil_downSup, this->avgPolDil_count_downSup,"avgPolDil_downSup");
  makeOneAvg(this->avgPolDil_downSdown, this->avgPolDil_count_downSdown,
	     "avgPolDil_downSdown");
  
  lrSpinCorr::AvgCorr();
}


void FourTargets::Print_LR(TString target){
  if (target=="upSup") {
    CoutLoop(this->left_upSup_upP, "left_"+target+"_upP");
    std::cout << "\n";
    CoutLoop(this->left_upSup_downP, "left_"+target+"_downP");
    std::cout << "\n";
    CoutLoop(this->right_upSup_upP, "right_"+target+"_upP");
    std::cout << "\n";
    CoutLoop(this->right_upSup_downP, "right_"+target+"_downP");
  }
  if (target=="upSdown") {
    CoutLoop(this->left_upSdown_upP, "left_"+target+"_upP");
    std::cout << "\n";
    CoutLoop(this->left_upSdown_downP, "left_"+target+"_downP");
    std::cout << "\n";
    CoutLoop(this->right_upSdown_upP, "right_"+target+"_upP");
    std::cout << "\n";
    CoutLoop(this->right_upSdown_downP, "right_"+target+"_downP");
  }
}


void FourTargets::Print_PolDil(TString target){
  if (target=="upSup") CoutLoop(this->avgPolDil_upSup, "avgPolDil_"+target);
  if (target=="upSdown") CoutLoop(this->avgPolDil_upSdown, "avgPolDil_"+target);
}


void FourTargets::WriteAll(){
  TGraph *g_left_upSup_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_upSup_upP[0]));
  TGraph *g_left_upSup_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_upSup_downP[0]));
  TGraph *g_left_upSdown_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_upSdown_upP[0]));
  TGraph *g_left_upSdown_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_upSdown_downP[0]));

  TGraph *g_right_upSup_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_upSup_upP[0]));
  TGraph *g_right_upSup_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_upSup_downP[0]));
  TGraph *g_right_upSdown_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_upSdown_upP[0]));
  TGraph *g_right_upSdown_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_upSdown_downP[0]));

  TGraph *g_left_downSup_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_downSup_upP[0]));
  TGraph *g_left_downSup_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_downSup_downP[0]));
  TGraph *g_left_downSdown_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_downSdown_upP[0]));
  TGraph *g_left_downSdown_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->left_downSdown_downP[0]));

  TGraph *g_right_downSup_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_downSup_upP[0]));
  TGraph *g_right_downSup_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_downSup_downP[0]));
  TGraph *g_right_downSdown_upP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_downSdown_upP[0]));
  TGraph *g_right_downSdown_downP
    = new TGraph(this->nBins, &(this->xval[0]), &(this->right_downSdown_downP[0]));

  setup(g_left_upSup_upP);
  setup(g_left_upSup_downP);
  setup(g_left_upSdown_upP);
  setup(g_left_upSdown_downP);

  setup(g_right_upSup_upP);
  setup(g_right_upSup_downP);
  setup(g_right_upSdown_upP);
  setup(g_right_upSdown_downP);

  setup(g_left_downSup_upP);
  setup(g_left_downSup_downP);
  setup(g_left_downSdown_upP);
  setup(g_left_downSdown_downP);

  setup(g_right_downSup_upP);
  setup(g_right_downSup_downP);
  setup(g_right_downSdown_upP);
  setup(g_right_downSdown_downP);
  
  g_left_upSup_upP->Write(Form("%s_left_upSup_upP", this->thisName.Data()) );
  g_left_upSup_downP->Write(Form("%s_left_upSup_downP", this->thisName.Data()));
  g_left_upSdown_upP->Write(Form("%s_left_upSdown_upP", this->thisName.Data()) );
  g_left_upSdown_downP->Write(Form("%s_left_upSdown_downP", this->thisName.Data()));
  g_right_upSup_upP->Write(Form("%s_right_upSup_upP", this->thisName.Data()) );
  g_right_upSup_downP->Write(Form("%s_right_upSup_downP", this->thisName.Data()));
  g_right_upSdown_upP->Write(Form("%s_right_upSdown_upP", this->thisName.Data()) );
  g_right_upSdown_downP->Write(Form("%s_right_upSdown_downP", this->thisName.Data()));

  g_left_downSup_upP->Write(Form("%s_left_downSup_upP", this->thisName.Data()) );
  g_left_downSup_downP->Write(Form("%s_left_downSup_downP", this->thisName.Data()));
  g_left_downSdown_upP->Write(Form("%s_left_downSdown_upP", this->thisName.Data()) );
  g_left_downSdown_downP->Write(Form("%s_left_downSdown_downP", this->thisName.Data()));
  g_right_downSup_upP->Write(Form("%s_right_downSup_upP", this->thisName.Data()) );
  g_right_downSup_downP->Write(Form("%s_right_downSup_downP", this->thisName.Data()));
  g_right_downSdown_upP->Write(Form("%s_right_downSdown_upP", this->thisName.Data()) );
  g_right_downSdown_downP->Write(Form("%s_right_downSdown_downP", this->thisName.Data()));

  
  TGraph *g_avgPolDil_upSup
    = new TGraph(this->nBins, &(this->xval[0]), &(this->avgPolDil_upSup[0]));
  TGraph *g_avgPolDil_upSdown
    = new TGraph(this->nBins, &(this->xval[0]), &(this->avgPolDil_upSdown[0]));
  TGraph *g_avgPolDil_downSup
    = new TGraph(this->nBins, &(this->xval[0]), &(this->avgPolDil_downSup[0]));
  TGraph *g_avgPolDil_downSdown
    = new TGraph(this->nBins,&(this->xval[0]), &(this->avgPolDil_downSdown[0]));
  setup(g_avgPolDil_upSup); setup(g_avgPolDil_upSdown);
  setup(g_avgPolDil_downSup); setup(g_avgPolDil_downSdown);

  g_avgPolDil_upSup->Write(Form("%s_avgPolDil_upSup", this->thisName.Data()) );
  g_avgPolDil_upSdown->Write(Form("%s_avgPolDil_upSdown",
				  this->thisName.Data()) );
  g_avgPolDil_downSup->Write(Form("%s_avgPolDil_downSup",
				  this->thisName.Data()) );
  g_avgPolDil_downSdown->Write(Form("%s_avgPolDil_downSdown",
				    this->thisName.Data()));
}

//Helper functions
void FourTargets::setZero(Int_t nBins){
  for (Int_t i=0; i<nBins; i++) {
    this->left_upSup_upP.push_back(0.0);
    this->left_upSdown_upP.push_back(0.0);
    this->left_upSup_downP.push_back(0.0);
    this->left_upSdown_downP.push_back(0.0);

    this->left_downSup_upP.push_back(0.0);
    this->left_downSdown_upP.push_back(0.0);
    this->left_downSup_downP.push_back(0.0);
    this->left_downSdown_downP.push_back(0.0);

    this->right_upSup_upP.push_back(0.0);
    this->right_upSdown_upP.push_back(0.0);
    this->right_upSup_downP.push_back(0.0);
    this->right_upSdown_downP.push_back(0.0);

    this->right_downSup_upP.push_back(0.0);
    this->right_downSdown_upP.push_back(0.0);
    this->right_downSup_downP.push_back(0.0);
    this->right_downSdown_downP.push_back(0.0);

    this->avgPolDil_upSup.push_back(0.0);
    this->avgPolDil_upSdown.push_back(0.0);
    this->avgPolDil_downSup.push_back(0.0);
    this->avgPolDil_downSdown.push_back(0.0);

    this->avgPolDil_count_upSup.push_back(0);
    this->avgPolDil_count_upSdown.push_back(0);
    this->avgPolDil_count_downSup.push_back(0);
    this->avgPolDil_count_downSdown.push_back(0);
  }
}


bool FourTargets::binOneAvg(std::vector<double> &Avg, std::vector<int> &count, 
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


bool FourTargets::makeOneAvg(std::vector<double> &Avg, std::vector<int> &count,
			    TString target){
  if (Avg.size() != count.size() ){
    std::cout << "FourTargets::makeOneAvg size mismatch in: " << target <<
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
