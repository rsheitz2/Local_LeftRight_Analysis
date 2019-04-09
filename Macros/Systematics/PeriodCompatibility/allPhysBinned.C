#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void allPhysBinned(TString start=""){
  //Setup_______________
  /*const Int_t nBins =3;
  TString fitMrangeType ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  const Int_t nBins =4;
  TString fitMrangeType ="LowM_AMDY";
  Int_t hbins =150;
  TString process ="JPsi";//JPsi, psi, DY
  TString lrMrange ="2.87_3.38";
  TString fitMrange ="2.87_3.38";
  TString binRange ="29_34";
  TString whichFit ="true";
  TString production ="slot1";//"t3", "slot1"
  TString additionalCuts ="phiS0.0";//*/

  Bool_t toWrite =false;
  //Setup_______________

  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/PeriodCompatibility/\
Data/physBinned";
  const Int_t nPhysBinned =6;
  TString physBinned[nPhysBinned] ={"xN", "xN", "xPi", "xF", "pT", "M"};
  //1st xN is for integrated
  
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
  = {"07", "08", "09", "10", "11", "12", "13", "14", "15"};//*/
  /*const Int_t nPeriods =3;
  TString periods[nPeriods] = {"07", "11", "15"};//*/
  
  //Aesthetics setup
  gStyle->SetLineWidth(2);
  TCanvas* c1 = new TCanvas();
  c1->Divide(nPhysBinned, 1, 0, 0); c1->SetLeftMargin(0.2);
  TLegend *legend = new TLegend(0.12, 0.8, 0.8, 0.99);
  legend->SetNColumns(3);
  Double_t offsets[nPhysBinned] = {0.0, 0.0023, 0.006, 0.006, 0.02, 0.02};
  Double_t yMax = (fitMrangeType=="HMDY") ? 0.9 : 0.12;
  TString xNames[nPhysBinned] =
    {"", "x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};
  
  //Get Data file/Get graphs and plot
  for (Int_t phys=0; phys<nPhysBinned; phys++) {

    Int_t nBinsName =nBins, integrated =0;
    if (phys == 0) {//used for integrated
      nBinsName =1;
      integrated =1;
    }
    TString fname;
    if (whichFit == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      fname=Form("%s/physBinnedPeriod_%s_%s_%s%s_%s%s%i_%s_%s.root",path.Data(),
		 whichFit.Data(), fitMrangeType.Data(), process.Data(),
		 lrMrange.Data(), binRange.Data(), physBinned[phys].Data(),
		 nBinsName, production.Data(), additionalCuts.Data());
    }
    else {
      fname =Form("%s/physBinnedPeriod_%s%s_%s_%s%s_%s%i_%ihbin_%s_%s.root",
		  path.Data(), whichFit.Data(), fitMrange.Data(),
		  fitMrangeType.Data(), process.Data(),lrMrange.Data(),
		  physBinned[phys].Data(), nBinsName, hbins, production.Data(),
		  additionalCuts.Data() );
    }
    TFile *f_in = OpenFile(fname);

    c1->cd(phys+1);
    for (Int_t p=0; p<nPeriods; p++) {
      TGraphErrors *g_AN
	= (TGraphErrors*) f_in->Get(Form("AN_W%s", periods[p].Data() ));
      OffSet(g_AN, offsets[phys]*p);
      
      if (p==0) {
	gPad->SetFrameLineWidth(2);
	FinalSetup(g_AN); FinalClearTitles(g_AN);
	g_AN->Draw("AP");
	g_AN->GetYaxis()->SetRangeUser(-yMax, yMax);
	SetTitleName(g_AN, xNames[phys], "x");

	if (fitMrangeType=="HMDY"){//Set xaxis range
	  if (integrated)
	    g_AN->GetXaxis()->SetLimits(0.165, 0.179);
	  else if (physBinned[phys] == "xN")
	    g_AN->GetXaxis()->SetLimits(0.05, 0.3);
	  else if (physBinned[phys] == "xPi")
	    g_AN->GetXaxis()->SetLimits(0.22, 0.82);
	  else if (physBinned[phys] == "xF")
	    g_AN->GetXaxis()->SetLimits(0.05, 0.7);
	  else if (physBinned[phys] == "pT")
	    g_AN->GetXaxis()->SetLimits(0.5, 2.4);
	  else if (physBinned[phys] == "M")
	    g_AN->GetXaxis()->SetLimits(4.3, 7.2);
	}
	else if (fitMrangeType=="LowM_AMDY"){//Set xaxis range
	  if (nBins==3){
	    if (integrated)
	      g_AN->GetXaxis()->SetLimits(0.162, 0.18);
	    else if (physBinned[phys] == "xN")
	      g_AN->GetXaxis()->SetLimits(0.04, 0.18);
	    else if (physBinned[phys] == "xPi")
	      g_AN->GetXaxis()->SetLimits(0.16, 0.55);
	    else if (physBinned[phys] == "xF")
	      g_AN->GetXaxis()->SetLimits(0.0, 0.5);
	    else if (physBinned[phys] == "pT")
	      g_AN->GetXaxis()->SetLimits(0.5, 2.2);
	    else if (physBinned[phys] == "M")
	      g_AN->GetXaxis()->SetLimits(2.5, 4.3);
	  }
	  else if (nBins==4){
	    if (integrated)
	      g_AN->GetXaxis()->SetLimits(0.162, 0.18);
	    else if (physBinned[phys] == "xN")
	      g_AN->GetXaxis()->SetLimits(0.04, 0.18);
	    else if (physBinned[phys] == "xPi")
	      g_AN->GetXaxis()->SetLimits(0.16, 0.55);
	    else if (physBinned[phys] == "xF")
	      g_AN->GetXaxis()->SetLimits(0.0, 0.5);
	    else if (physBinned[phys] == "pT")
	      g_AN->GetXaxis()->SetLimits(0.5, 2.2);
	    else if (physBinned[phys] == "M")
	      g_AN->GetXaxis()->SetLimits(2.5, 4.3);
	  }
	}

	DrawLine(g_AN, 0.0);
      }
      else g_AN->Draw("Psame");

      if (integrated){
	g_AN->GetXaxis()->SetLabelSize(0.0);
	g_AN->GetXaxis()->SetTickSize(0.0);
	legend->AddEntry(g_AN, Form("W%s", periods[p].Data()), "p");
      }

      g_AN->SetMarkerSize(1.8);
      
    }//end period loop
  }//phys binned loop

  //Legend
  c1->cd(1);
  legend->SetBorderSize(0); legend->SetTextFont(132); legend->SetTextSize(0.08);
  legend->Draw("same");
  
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
