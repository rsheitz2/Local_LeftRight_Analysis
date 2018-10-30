void tableOfCuts(TString start=""){
  //Setup_______________
  TString period ="W07";
  TString Mtype ="LowM_AMDY";
  TString production ="t3";//"t3"=t3, "slot1"=t5
  TString whichCuts ="FinalCuts";//"PhastCuts", "FinalCuts"
    
  Bool_t toWrite =true;
  //Setup_______________

  TString dataPath = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Presents/DATA/RealData/";
  TString inputFile = Form("%s/%s_%s%s.root", Mtype.Data(), period.Data(),
			   Mtype.Data(), production.Data());  ;
  if (production=="t3"){//No production name for t3 input
    inputFile=Form("%s/%s_%s.root", Mtype.Data(), period.Data(),
		   Mtype.Data());  
  }
  
  
  if (start==""){//Basic info
    cout << "\nScript goes through cuts histogram in input file and outputs";
    cout << "cut and numbers in a file" << endl;
    cout << "Only works for real data, FinalCuts at the moment....." << endl;
    cout << "\n\nUtilization:" << endl;
    cout << "root \'tableOfCuts.C(1)\'" <<endl;
    cout << "\n\nCurrent Settings________" <<endl;
    cout << "Period considered:        " << period << endl;
    cout << "Mass type considered:     " << Mtype << endl;
    cout << "Production considered:    " << production << endl;
    cout << "Which cuts considered:    " << whichCuts << endl;
    cout << "\nTo write output file:     " << toWrite << endl;
    cout << "\nData coming from    :     "  << dataPath + inputFile<<"\n"<<endl;
    exit(EXIT_FAILURE);
  }

  //Open data file
  TFile *fData = TFile::Open(dataPath + inputFile);
  if (!fData){
    cout << "fData file does not exist" << endl;
    exit(EXIT_FAILURE);
  }

  //Basic Checks
  if ( (production != "t3") && (production != "slot1")) {
    cout << "Invalid production setting   " << production << endl;
    exit(EXIT_FAILURE);
  }
  if ( (whichCuts!="FinalCuts") && (whichCuts!="PhastCuts") ){
    cout << "Invalid cuts setting   " << whichCuts << endl;
    exit(EXIT_FAILURE);
  }

  //Output file setup
  TString thisDirPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Cuts";
  TString outputName =Form("%s/Data/TableOfCuts/tableOfCuts_%s_%s_%s_%s.csv",
			   thisDirPath.Data(), period.Data(), Mtype.Data(),
			   production.Data(), whichCuts.Data());
  ofstream fCuts;
  if (toWrite) fCuts.open(outputName);
  
  //Get data from file and write to file
  //Integer for which cut/lable and has this integer been written yet
  Int_t iCut =-1, iLabel =-1, iWrite=false;
  Int_t numCuts;
  TString cutLabel;
  TH1D *hCuts = (TH1D*)fData->Get("hCuts");
  if (!hCuts){
    cout << "Histogram of cuts:  hCuts   does not exist in this file" << endl;
    exit(EXIT_FAILURE);
  }
  for (Int_t bi=1; bi<hCuts->GetNbinsX()+1; bi++) {
    if (hCuts->GetBinContent(bi) > 0){
      numCuts =hCuts->GetBinContent(bi);
      iCut++;
      iWrite =false;
    }
    if (strncmp(hCuts->GetXaxis()->GetBinLabel(bi), "", 1) != 0){
      cutLabel = hCuts->GetXaxis()->GetBinLabel(bi);
      iLabel++;
      iWrite =false;
    }

    if ((iCut == iLabel) && (iCut >= 0) && (!iWrite)){
      if (toWrite) {fCuts << cutLabel.Data() << ",   " << numCuts << "\n"; }
      iWrite=true;
    }
  }

  if (toWrite) fCuts.close();

  //final settings
  cout << "\n\nCurrent Settings________" <<endl;
  cout << "Period considered:        " << period << endl;
  cout << "Mass type considered:     " << Mtype << endl;
  cout << "Production considered:    " << production << endl;
  cout << "Which cuts considered:    " << whichCuts << endl;
  cout << "\nData         from   :    "  << dataPath + inputFile << endl;
  if (toWrite){
    cout << "\nFile:  " << outputName << "   was written" << endl;
  }
  else cout << "\nFile: " << outputName << " was NOT written" << endl;
}
