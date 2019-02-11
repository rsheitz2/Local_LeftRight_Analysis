#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"

void comparison(TString start=""){
  //Setup_______________
  const Int_t nBins =3;//HMDY
  TString Mtype ="HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString binRange ="43_85";
  TString physBinned ="xN";
  TString whichFit ="true";
  TString production ="slot1";
  TString additionalCuts ="phiS0.0";//*/

  TString xChecker ="michael";
  TString corrected ="uncorr"; //"uncorr"=NOT corrected or ""=corrected

  Bool_t toWrite =false;
  //Setup_______________

  //Basic Setup
  TString localPath ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Comparisons/XCheck";
  TString pathAN="/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/AN_calculation/Data/GeoMean4Targ";

  /*if (start==""){
    cout <<"Script draws Final AN asymmetry with systematics binned in physics";
    cout << "\nSystematics are input as a TGraph for each physics bin";
    cout << "\n\n\nUsage:" << endl;
    cout << "root \'binned4TargGeoMean.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Asym data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " <<  whichFit << endl;
    cout << "\n\nOutput is to be written:     " << toWrite << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
    }//*/

  //Basic setup
  const Int_t nPer =3;
  TString period[nPer] =
    {"W07", "W08", "W09"};
  
  Double_t *xvals, ex[nPer] ={0.0};
  TCanvas* cComp = new TCanvas(); cComp->Divide(nPer);

  //Period loop
  for (Int_t p=0; p<nPer; p++) {
    cComp->cd(p+1);
    
    //My calculation
    TString n_Wper;
    if (whichFit=="true"){
      n_Wper =
	Form("%s/GeoMean4Targ_%s_%s_%s_%s%s_%s%s%i_%s_%s.root",pathAN.Data(),
	     whichFit.Data(), period[p].Data(), Mtype.Data(),
	     process.Data(), lrMrange.Data(), binRange.Data(),physBinned.Data(),
	     nBins, production.Data(), additionalCuts.Data());
    }
    else{
      n_Wper =
	Form("%s/GeoMean4Targ_%s%s_%s_%s_%s%s_%s%s%i_%ihbin_%s_%s.root",
	     pathAN.Data(), whichFit.Data(), fitMrange.Data(), period[p].Data(),
	     Mtype.Data(), process.Data(), lrMrange.Data(), binRange.Data(),
	     physBinned.Data(), nBins, hbins, production.Data(),
	     additionalCuts.Data());
    }
    TFile *f_Wper = OpenFile(n_Wper);
    TGraphErrors *g_Wper =
      (TGraphErrors*) f_Wper->Get(Form("AN%s", corrected.Data()));

    if (p==0) xvals = g_Wper->GetX();
    
    //XCheck calculation
    TString xCheckFname =
      Form("%s/Data/CompValues/%s%s%s.txt", localPath.Data(), xChecker.Data(),
	   period[p].Data(), corrected.Data());

    ifstream xCheckFile(xCheckFname);
    if (!xCheckFile){
      cout << "File did not open:  " << xCheckFname << endl;
      exit(EXIT_FAILURE);
    }

    Double_t An[nBins], e_An[nBins];
    Int_t ibin =0;
    while( xCheckFile >> An[ibin] >> e_An[ibin]){
      ibin++;
    }

    TGraphErrors* gXcheck =
      new TGraphErrors(nPer, xvals, An, ex, e_An);
    SetUp(gXcheck);

    //Draw by period
    g_Wper->Draw("AP"); DrawLine(g_Wper, 0.0);
    g_Wper->GetYaxis()->SetRangeUser(-0.08, 0.08);
    g_Wper->SetTitle(Form("%s", period[p].Data()));
    g_Wper->SetMarkerStyle(20);
    
    gXcheck->Draw("Psame"); gXcheck->SetMarkerColor(kBlue);
    gXcheck->SetMarkerStyle(20);
    OffSet(gXcheck, 0.005);
  }
  
  /*//Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/FinalAsym/Data/PhysBinnedData";
  TString fOutput;
  if (whichFit=="true"){
    fOutput = Form("%s/physBinnedData_%s_%s_%s%s_%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), Mtype.Data(),
		   process.Data(), lrMrange.Data(), nBins,
		   hbins, production.Data(), additionalCuts.Data());
  }
  else{
    fOutput = Form("%s/physBinnedData_%s%s_%s_%s%s_%i_%ihbin_%s_%s.root",
		   thisDirPath.Data(), whichFit.Data(), fitMrange.Data(),
		   Mtype.Data(), process.Data(), lrMrange.Data(),
		   nBins, hbins, production.Data(), 
		   additionalCuts.Data());
  }

  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    for (Int_t i=0; i<nPhysBinned; i++) {
      g_AN[i]->Write(Form("AN_%s", physBinned[i].Data() ));
    }

    cAsym->Write();
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Asymmetry data coming from:            " << pathAN << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << Mtype << endl;
    cout << "AN physical process:        " << process << endl;
    cout << "LR integral mass range:     " << lrMrange << endl;
    cout << "Fit mass range:     " << fitMrange << endl;
    cout << "Which fit considered:       " <<  whichFit << endl;
  }
  cout << " " << endl;
  if (toWrite) cout << "File:  " << fOutput << "   was written" << endl;
  else cout << "File: " << fOutput << " was NOT written" << endl;//*/
}

