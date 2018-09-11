#include "include/helperFunctions.h"


void OffSet(TGraphErrors *g, Double_t offset){
  Double_t *xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xval[i] += offset;

}


void allPhysBinned(TString start=""){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/\
Data/physBinned";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  
  //Setup_______________
  const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit[nPhysBinned] = {"true", "true", "true", "true"};
  //TString whichFit[nPhysBinned] = {"six", "six", "seven", "seven"};

  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout << "Script draws AN per period and physics binning on a nice plot";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C  ->  GeoMean4Targ.C  ->";
    cout << "  physBinnedPeriod.C  ->  allPhysBinned.C" << endl;
    cout << "or Scripts/physBinnedPipeline.sh  -> allPhysBinned.C" << endl;
    cout << "\nUsage:" << endl;
    cout << "root \'allPhysBinned.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) 
      cout << whichFit[i] << " ";
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  //Period name setups
  const Int_t nPeriods =9;
  TString periods[nPeriods]
    = {"07", "08", "09", "10", "11", "12", "13", "14", "15"};
  
  //Aesthetics setup
  TCanvas* c1 = new TCanvas(); c1->Divide(4, 1, 0, 0.01);
  Double_t offsets[nPhysBinned] = {0.002, 0.004, 0.004, 0.02};
  Double_t yMax =0.9;
  
  //Get Data file/Get graphs and plot
  TString docName=""; TString physBinnedNames =""; TString fitNames="";
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString fname;
    if (whichFit[phys] == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      fname =Form("%s/physBinnedPeriod_%s_%s_%s%s_%s%i.root", path.Data(),
		  whichFit[phys].Data(), fitMrangeType.Data(), process.Data(),
		  lrMrange.Data(), physBinned[phys].Data(), nBins);
    }
    else {
      fname =Form("%s/physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin.root",
		  path.Data(), whichFit[phys].Data(), fitMrange.Data(),
		  fitMrangeType.Data(), process.Data(),lrMrange.Data(),
		  physBinned[phys].Data(), nBins, hbins);
    }
    
    TFile *f_in = TFile::Open(fname);
    if (!f_in){
      cout << "RD or RD_noCorr file does not exist " << endl;
      exit(EXIT_FAILURE);
    }
    docName += fname+"\n";
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit[phys]+" ";

    c1->cd(phys+1);
    for (Int_t p=0; p<nPeriods; p++) {
      TGraphErrors *g_AN
	= (TGraphErrors*) f_in->Get(Form("AN_%s", periods[p].Data() ));
      OffSet(g_AN, offsets[phys]*p);
      
      if (p==0) {
	g_AN->Draw("AP");
	g_AN->GetYaxis()->SetRangeUser(-yMax, yMax);
	
	if (physBinned[phys] == "xN")
	  g_AN->GetXaxis()->SetLimits(0.05, 0.3);
	else if (physBinned[phys] == "xPi")
	  g_AN->GetXaxis()->SetLimits(0.25, 0.8);
	else if (physBinned[phys] == "xF")
	  g_AN->GetXaxis()->SetLimits(0.05, 0.7);
	else if (physBinned[phys] == "pT")
	  g_AN->GetXaxis()->SetLimits(0.5, 2.5);

	DrawLine(g_AN, 0.0);
      }
      else g_AN->Draw("Psame");
    }//end period loop
    
  }//phys binned loop
  
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility";
  TString fOutput;
  if (whichFit[0] == "true"){
    fOutput =
      Form("%s/Data/allPhysBinned/allPhysBinned_true_%s_%s%s_%ibins.root",
	   thisDirPath.Data(), fitMrangeType.Data(), process.Data(),
	   lrMrange.Data(), nBins);
  }
  else {
    fOutput =
      Form("%s/Data/allPhysBinned/allPhysBinned_%s_%s_%s%s_%ibins_%ihbin.root",
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
    
    c1->Write();
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << fitMrangeType << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fits considered:       " << endl;
    for (Int_t i=0; i<nPhysBinned; i++) 
      cout << whichFit[i] << " ";
  }
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;
}
