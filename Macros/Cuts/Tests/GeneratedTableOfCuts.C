
void GeneratedTableOfCuts(TString fname="noFile",
			      TString outFile="GeneratedCuts.txt"){
  //Used to determine cut impacts per period and put outputs in
  //OneSampleCuts.txt
  //Only need to change data locations to use for another data set

  if (fname=="noFile"){
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << "Macro used to create an impact of cuts table from ";
    cout << "generated events" << endl;
    cout << "!!!!!!!!!!!!!!!" << endl;
    cout << " " << endl;
    cout << "Usage root:" << endl;
    cout<< "\'OneSampleMakeTableOfCuts.C(\"fname.root\",";
    cout << "\"outputFile.txt\")\'" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  ofstream cuts(outFile);
  
  if (fname!="noFile") {
    TFile *f = TFile::Open(fname);
    TH1D *hCuts = (TH1D*)f->Get("h_Cuts");

    for (Int_t bi=1; bi<hCuts->GetNbinsX()+1; bi++) {
      if (hCuts->GetBinContent(bi) > 0){
	cuts << hCuts->GetBinContent(bi) << " ";
	
	TString stringCut(hCuts->GetXaxis()->GetBinLabel(bi-1));
	stringCut.ReplaceAll(" ", "");
	cuts << stringCut.Data() << " ";
      }
    }
  }
  
  cuts.close();
}
