#include "Include/helperFunctions.h"

void writeFile_Mmumu(){
	//Setup_____
	TString whichProcess ="OC"; //"Jpsi", "Psi", "OC", "t3", "slot1"

	Bool_t toWrite =false;
	//Setup_____

	cout << "\nMacro combines writes normalized invariant mass distribution to file ";
	cout << "\nWhich MC component considered: " << whichProcess << endl;
	cout << "\nOutput is written to a file:   " << toWrite << endl;
	cout << "" << endl;

	TString path =
		"/u/sciteam/heitz/Analysis/TGeant/Main/Data/";
	TString nameFile;
	if (whichProcess =="t3" || whichProcess =="slot1"){
		if (whichProcess =="slot1") nameFile ="slot1";
		nameFile +="WAll_LowM_AMDY.root";
		path +="RealDataOutput/LowM_AMDY/";
	}
	else {
		nameFile =Form("Charles_W12_%s.root", whichProcess.Data());
		path +="MainOutput/Charles_Official/";
	}
	TFile *f = TFile::Open(path+nameFile);
	TTree *t = (TTree*) f->Get("pT_Weighted");

	TH1D *h = new TH1D(Form("h_%s", whichProcess.Data()), Form("h_%s", whichProcess.Data()),
										 200, 2.0, 8.5);
	t->Draw(Form("Mmumu>>h_%s", whichProcess.Data()), "Mmumu>2.0&&Mmumu<8.5", "0");

	h->Sumw2();
	if (whichProcess =="t3" || whichProcess =="slot1"){}
	else{	h->Scale(1.0/(h->Integral() ) ); }

	TCanvas *c1 = new TCanvas(); 
	gPad->SetLogy();
	h->Draw();

	TString outName =Form("M_%s.root", whichProcess.Data());
	if (toWrite){
		TString thisPath =
			"/u/sciteam/heitz/Analysis/TGeant/Main/Data/MainOutput/Charles_Official/";
		thisPath += "Macros/Data/MassDist/";
		cout << "     Full file path and name:  ";
		cout << thisPath + outName << endl;
		TFile *fOutput = 
			new TFile(thisPath+outName, "RECREATE");
		h->Write();
		fOutput->Close();
	}

	//Final output
	cout << "\nWhich MC component considered: " << whichProcess << endl;
	cout << "\nFile  " << outName << "   was written:   " << toWrite << endl;
	cout << "" << endl;
}
