
#include "/afs/cern.ch/work/r/rheitz/Analysis/TGeant/Main/src/common.hxx"

void OneSampleMakeTableOfCuts(TString phastFile="noFile",
			      TString mainFile="noFile",
			      TString outFile="OneSampleCuts",
			      TString userEvent="UserEvent411"){
  //Used to determine cut impacts per period and put outputs in
  //OneSampleCuts.txt
  //Only need to change data locations to use for another data set

  if (phastFile=="noFile" && mainFile=="noFile"){
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << "Macro used to create an impact of cuts table" << endl;
    cout << "Can be used with impact after phast or after main or any";
    cout << " combination" << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << " " << endl;
    cout << "Usage root:" << endl;
    cout<< "\'OneSampleMakeTableOfCuts.C(\"phastFile.root\",\"mainFile.root\",";
    cout << "\"outputFile.txt\")\'" << endl;
    cout << " " << endl;
    cout << "Example usage without phastFile:";
    cout << "\'OneSampleMakeTableOfCuts.C(\"noFile\", \"mainFile.root\", ";
    cout << "\"outputFile.txt\",\"UserEvent411\")\'" << endl;
    cout << " " << endl;
    cout << "Typical UserEvents:" << endl;
    cout << "UserEvent411 for Monte Carlo data" << endl;
    cout << "UserEvent404 for real data" << endl;
    cout << " " << endl;
    exit();
  }

  TFile *fPhast, *fMain;
  TH1D *hPhastCuts, *hMainCuts;
  ofstream cuts("OneSampleCuts.txt");
  
  Bool_t hasPhast=false;
  if (phastFile!="noFile") {
    hasPhast = true;
    fPhast = TFile::Open(phastFile);
    userEvent+="/DiMuonCuts";
    hPhastCuts = (TH1D*)fPhast->Get(userEvent);

    for (Int_t bi=1; bi<hPhastCuts->GetNbinsX()+1; bi++) {
      if (hPhastCuts->GetBinContent(bi) > 0){
	cuts << hPhastCuts->GetBinContent(bi) << " ";
	
	TString stringCut(hPhastCuts->GetXaxis()->GetBinLabel(bi));
	stringCut.ReplaceAll(" ", "");
	cuts << stringCut.Data() << " ";
      }
    }
  }
  if (mainFile!="noFile") {
    fMain = TFile::Open(mainFile);
    hMainCuts = (TH1D*)fMain->Get("hCuts");

    Bool_t first=true;
    for (Int_t bi=1; bi<hMainCuts->GetNbinsX()+1; bi++) {
      if (hMainCuts->GetBinContent(bi) > 0){
	
	if (first && hasPhast){//Skip dublicate cuts
	  first = false;
	  continue;
	}
	
	cuts << hMainCuts->GetBinContent(bi) << " ";
	
	TString stringCut(hMainCuts->GetXaxis()->GetBinLabel(bi+1));
	stringCut.ReplaceAll(" ", "");
	cuts << stringCut.Data() << " ";
      }
    }
  }
  
  cuts.close();
}
