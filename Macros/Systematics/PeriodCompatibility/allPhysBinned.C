#include "include/helperFunctions.h"


void OffSet(TGraphErrors *g, Double_t offset){
  Double_t *xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xval[i] += offset;

}


void allPhysBinned(TString start=""){
  //Setup_______________
  const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.53";//*/

  /*const Int_t nBins =5;
  TString fitMrangeType ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.00_5.00";
  TString fitMrange ="2.00_8.50";
  TString binRange ="25_43";
  TString whichFit ="thirteen";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.53";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/\
Data/physBinned";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  
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
    cout << "Which fits considered:       " <<  whichFit << endl;
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
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      fname =Form("%s/physBinnedPeriod_%s_%s_%s%s_%s%i_%s_%s.root", path.Data(),
		  whichFit.Data(), fitMrangeType.Data(), process.Data(),
		  lrMrange.Data(), physBinned[phys].Data(), nBins,
		  production.Data(), additionalCuts.Data());
    }
    else {
      fname =Form("%s/physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
		  path.Data(), whichFit.Data(), fitMrange.Data(),
		  fitMrangeType.Data(), process.Data(),lrMrange.Data(),
		  physBinned[phys].Data(), nBins, hbins, production.Data(),
		  additionalCuts.Data() );
    }
    
    TFile *f_in = TFile::Open(fname);
    if (!f_in){
      cout << "RD or RD_noCorr file does not exist " << endl;
      exit(EXIT_FAILURE);
    }
    docName += fname+"\n";
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit+" ";

    c1->cd(phys+1);
    for (Int_t p=0; p<nPeriods; p++) {
      TGraphErrors *g_AN
	= (TGraphErrors*) f_in->Get(Form("AN_W%s", periods[p].Data() ));
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
  if (whichFit == "true"){
    fOutput =
      Form("%s/Data/allPhysBinned/allPhysBinned_true_%s_%s%s_%ibins_%s_%s.root",
	   thisDirPath.Data(), fitMrangeType.Data(), process.Data(),
	   lrMrange.Data(), nBins, production.Data(), additionalCuts.Data());
  }
  else {
    fOutput =
      Form("%s/Data/allPhysBinned/\
allPhysBinned_%s_%s_%s%s_%ibins_%ihbin_%s_%s.root",
	   thisDirPath.Data(), fitMrange.Data(), fitMrangeType.Data(),
	   process.Data(), lrMrange.Data(), nBins, hbins, production.Data(),
	   additionalCuts.Data());
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
