#include "include/helperFunctions.h"


void pullDist(TString start=""){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/\
Data/physBinned";
  const Int_t nPeriod = 9;
  const Int_t nPhysBinned =4;
  
  //Setup_______________
  const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50"; //Must == fitMrange for true fits
  TString fitMrange ="4.30_8.50"; 

  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  //TString whichFit[nPhysBinned] ={"seven", "seven", "six", "six"};
  TString whichFit[nPhysBinned] ={"true", "true", "true", "true"};

  Bool_t toWrite =false;
  //Setup_______________  
    
  if (start==""){
    cout <<"\nThis Script computes the pull distributions for an L/R asymmetry";
    cout <<"\nThe asymmetry is determined ";
    cout << "(# periods)*(# physics kinematics)*(# physics kinematic bins)\n";
    cout << "Input data should be a TGraphErrors for each period and";
    cout << "physical kinematic";
    cout << "\n\nCurrent Setup:" << endl;
    cout << "Data comes from:                 " << path << endl;
    cout << "Number of periods considered:    " << nPeriod << endl;
    cout << "Number of physics kinematic bins " << nBins << endl;
    cout << "Mass type considered:            " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << "Fit physbinning:    " << physBinned[i] << "   fit considered:  ";
      cout << whichFit[i] << endl;
    }
    cout << "\nOutput is to be written:     " << toWrite << endl;
    cout << "\nUsage:" << endl;
    cout <<"root \'PullDist(1)\'\n" << endl;
    exit(EXIT_FAILURE);
  }

  //Get Data by physBinned
  TFile* f_physBinned[nPhysBinned];
  TString docName=""; TString physBinnedNames =""; TString fitNames="";
  vector<Double_t> lr, e_lr;
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString fInput;
    if (whichFit[phys] == "true" ){
      if (lrMrange != fitMrange){
	cout << "Error: lrMrange != fitMrange with whichFit==true" << endl;
	exit(EXIT_FAILURE);
      }
      
      fInput =Form("%s/physBinnedPeriod_%s_%s_%s%s_%s%i.root", path.Data(),
		   whichFit[phys].Data(), fitMrangeType.Data(), process.Data(),
		   fitMrange.Data(), physBinned[phys].Data(), nBins);
    }
    else{
      fInput =Form("%s/physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin.root",
		   path.Data(), whichFit[phys].Data(), fitMrange.Data(),
		   fitMrangeType.Data(), process.Data(), lrMrange.Data(),
		   physBinned[phys].Data(), nBins, hbins);
    }
    f_physBinned[phys] = TFile::Open(fInput);

    if ( !f_physBinned[phys]){
      cout << "File does not exist:  # " << phys << endl;
      exit(EXIT_FAILURE);
    }
    docName += fInput+"\n";
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit[phys]+" ";

    TGraphErrors *g_AN[nPeriod];
    for (Int_t p=0; p<nPeriod; p++) {
      g_AN[p] = (p+7<10)
	? (TGraphErrors*)f_physBinned[phys]->Get(Form("AN_0%i", p+7))
	: (TGraphErrors*)f_physBinned[phys]->Get(Form("AN_%i", p+7));

      Double_t *yVal=g_AN[p]->GetY();
      Double_t *e_yVal=g_AN[p]->GetEY();

      for (Int_t bi=0; bi<nBins; bi++) {
	lr.push_back(yVal[bi]);
	e_lr.push_back(e_yVal[bi]);
      }

    }//period loop
  }//nPhysBinned loop

  //Determine Avg/Std
  Double_t avg = WeightedAvg(lr, e_lr);
  Double_t sig = WeightedErr(e_lr);
  Double_t sigma2 = sig*sig;
  
  //Make and fill pull dist
  gStyle->SetOptFit(1111);
  TH1D* hPull = new TH1D("hPull", "hPull", 20, -4, 4);
  SetUp(hPull);
  for(vector<Double_t>::iterator it=lr.begin(), e_it=e_lr.begin(); it!=lr.end();
      it++, e_it++){
    Double_t pull = *it - avg;
    Double_t denom = (*e_it)*(*e_it) - sigma2;
    
    if (denom < 0.0) denom *= -1.0;
    denom = TMath::Sqrt(denom);

    pull = pull/denom;

    if ( isnan(pull) ) {
      cout << "Pull is nan" << endl;
      exit(EXIT_FAILURE);
    }

    hPull->Fill(pull);
  }

  //Draw Pull
  TCanvas* c1 = new TCanvas();
  hPull->SetMarkerStyle(20);
  hPull->SetMarkerColor(kBlue);
  hPull->SetLineColor(kBlue);
  hPull->Draw("EP");
  hPull->Fit("gaus");

  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility";
  TString fOutput;
  if  (lrMrange == fitMrange){
    fOutput =Form("%s/Data/pullDist/pullDist_%s_%s%s_%ibins.root",
		  thisDirPath.Data(), fitMrangeType.Data(), process.Data(),
		  lrMrange.Data(), nBins);
  }
  else {
    fOutput =Form("%s/Data/pullDist/pullDist_%s_%s_%s%s_%ibins_%ihbin.root",
		  thisDirPath.Data(), fitMrange.Data(), fitMrangeType.Data(),
		  process.Data(), lrMrange.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TNamed docs ("InputData", docName.Data());
    TNamed pBinNam ("physBinned", physBinnedNames.Data());
    TNamed fitNam ("fitNames", fitNames.Data());
    docs.Write();
    pBinNam.Write();
    fitNam.Write();
    
    hPull->Write();
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "pullDist nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    for (Int_t i=0; i<nPhysBinned; i++) {
      cout << "Fit physbinning:    " << physBinned[i] << "   fit considered:  ";
      cout << whichFit[i] << endl;
    }
    cout << "\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
  }
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
