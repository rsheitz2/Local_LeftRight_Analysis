#include "include/helperFunctions.h"


void OffSet(TGraphErrors *g, Double_t offset){
  Double_t *xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) xval[i] += offset;

}


void FillValues(TH1D *h, TGraphErrors *g){

  Double_t *yVals = g->GetY();
  Double_t *e_yVals = g->GetEY();
  
  for (Int_t i=0; i<g->GetN(); i++){
    Double_t deviation = yVals[i]/e_yVals[i];
    if (deviation < 0) deviation *= -1.0;
    
    h->Fill(deviation);
  }
}


void allSysError(TString start=""){
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/\
TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/\
Data";
  const Int_t nPhysBinned =4;
  TString physBinned[nPhysBinned] ={"xN", "xPi", "xF", "pT"};
  
  //Setup_______________
  const Int_t nBins =3;
  TString period_Mtype ="WAll_HMDY";
  Int_t hbins =150;
  TString process ="DY";//JPsi, psi, DY
  TString lrMrange ="4.30_8.50";
  TString fitMrange ="4.30_8.50";
  TString whichFit[nPhysBinned] = {"true", "true", "true", "true"};
  //TString whichFit[nPhysBinned] = {"six", "six", "seven", "seven"};

  Bool_t toWrite =false;
  //Setup_______________  
  
  if (start==""){
    cout << "Script draws false asymmetries and systematic error from false asymmetries";
    cout << " per period and physics binning on a nice plot";
    cout << "\n\nTotal local pipeline needed for this script" << endl;
    cout << "leftRight_byTarget  ->  functMFit.C  ->  GeoMean4Targ.C  ->";
    cout << "falseGeoMean4Targ_targFlip.C  ->  sysErrorFA.C  ->  allSysError.C";
    cout << "\n\nUsage:" << endl;
    cout << "root \'allSysError.C(1)\'" << endl;
    cout << "\nCurrent settings:" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << period_Mtype << endl;
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
  
  //Aesthetics setup
  TCanvas* cFA = new TCanvas(); cFA->Divide(4, 2, 0, 0.01);
  Double_t offsets[nPhysBinned] = {0.006, 0.01, 0.01, 0.05};
  Double_t yMax =0.9;
  
  //Get Data file/Get graphs and plot
  TH1D* hSys = new TH1D("hSys", "hSys", 10, 0, 1.5);
  TString physBinnedNames ="", fitNames="";
  for (Int_t phys=0; phys<nPhysBinned; phys++) {
    TString FAname, sysName;
    if (whichFit[phys] == "true"){
      if (fitMrange != lrMrange){
	cout << "fit Mass range != left/right mass range with true fit" << endl;
	exit(EXIT_FAILURE);
      }
    
      FAname =
	Form("%s/TargFlip/falseGeoMean4Targ_%s_%s_%s%s_%s%i.root", path.Data(),
	     whichFit[phys].Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);

      sysName =
	Form("%s/sysError/sysErrorFA_%s_%s_%s%s_%s%i.root", path.Data(),
	     whichFit[phys].Data(), period_Mtype.Data(), process.Data(),
	     lrMrange.Data(), physBinned[phys].Data(), nBins);
    }
    else {
      FAname =
	Form("%s/TargFlip/falseGeoMean4Targ_%s%s_%s_%s%s_%s%i_%ihbin.root",
	     path.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);

      sysName =
	Form("%s/sysError/sysErrorFA_%s%s_%s_%s%s_%s%i_%ihbin.root",
	     path.Data(), whichFit[phys].Data(), fitMrange.Data(),
	     period_Mtype.Data(), process.Data(),lrMrange.Data(),
	     physBinned[phys].Data(), nBins, hbins);
    }
    
    TFile *f_FA = TFile::Open(FAname); TFile *f_sys = TFile::Open(sysName);
    if (!f_FA || !f_sys){
      cout << "False asymmetries or systematic error file does not exist"<<endl;
      exit(EXIT_FAILURE);
    }
    physBinnedNames += physBinned[phys]+" ";
    fitNames += whichFit[phys]+" ";

    cFA->cd(phys+1);
    TGraphErrors *g_FA_pol =(TGraphErrors*)f_FA->Get("falseAN_pol");
    TGraphErrors *g_FA_subper =(TGraphErrors*)f_FA->Get("falseAN_subper");
    OffSet(g_FA_subper, offsets[phys]);
    
    g_FA_pol->Draw("AP");
    g_FA_subper->Draw("Psame");
    g_FA_pol->SetTitle("");
    DrawLine(g_FA_pol, 0.0);

    cFA->cd(nPhysBinned+phys+1);
    TGraphErrors *g_sys =(TGraphErrors*)f_sys->Get("gSys");
    g_sys->Draw("AP");
    g_sys->SetTitle("");
    g_sys->GetYaxis()->SetRangeUser(0, 1.5);

    FillValues(hSys, g_FA_pol);
    FillValues(hSys, g_FA_subper);
  }//phys binned loop

  TCanvas* cDist = new TCanvas();
  hSys->Draw("E");
  SetUp(hSys);
  
  //Write Output/Final Settings
  TString thisDirPath="/Users/robertheitz/Documents/Research/DrellYan/Analysis\
/TGeant/Local_LeftRight_Analysis/Macros/Systematics/FalseAsym/Data/allSysError";
  TString fOutput;
  if (whichFit[0] == "true"){
    fOutput =
      Form("%s/allSysError_true_%s_%s%s_%ibins.root", thisDirPath.Data(),
	   period_Mtype.Data(), process.Data(), lrMrange.Data(), nBins);
  }
  else {
    fOutput =
      Form("%s/allSysError_%s_%s_%s%s_%ibins_%ihbin.root", thisDirPath.Data(),
	   fitMrange.Data(), period_Mtype.Data(), process.Data(),
	   lrMrange.Data(), nBins, hbins);
  }
  if(toWrite){
    TFile *fResults = new TFile(fOutput, "RECREATE");
    TNamed pBinNam ("physBinned", physBinnedNames.Data());
    TNamed fitNam ("fitNames", fitNames.Data());
    pBinNam.Write();
    fitNam.Write();
    
    cFA->Write();
    hSys->Write();
  }

  if (start!=1){
    cout << " " << endl;
    cout << "Settings______" << endl;
    cout << "Data coming from:            " << path << endl;
    cout << "physBinned nBins times:     " << nBins << endl;
    cout << "Mass type considered:   " << period_Mtype << endl;
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

