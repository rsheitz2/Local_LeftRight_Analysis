void ByTarg_ByProcess_Integrated(){

  const Int_t nGr = 7;
  TString gr_Titles[] = {"xN_LR", "xN_LR_UpStream", "xN_LR_DownStream",
			 "xN_LR_UpStream_Up", "xN_LR_UpStream_Down",
			 "xN_LR_DownStream_Up", "xN_LR_DownStream_Down"};

  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/April17/Data/MC_Data/";
  TFile* f_W12Mimic_int =
    TFile::Open(path+"W12Mimic_1bin.root");
  TGraphErrors* g_W12Mimic_int[nGr];

  TFile* f_StandardMimic_int =
    //TFile::Open(path+"StandardMimic_1bin.root");
    TFile::Open(path+"StandardMimic_1bin_switch.root");
  TGraphErrors* g_StandardMimic_int[nGr];

  //Basic file checks
  if ( !(f_W12Mimic_int) || !(f_StandardMimic_int) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  const Int_t nPro = 2;
  TGraphErrors **g_Process[nPro] = {g_W12Mimic_int, g_StandardMimic_int};
  
  Double_t offset = 0.04;
  gStyle->SetOptFit(1111);
  TCanvas* c1 = new TCanvas();
  for (Int_t i=0, ipad=1; i<nGr; i++) {

    //Get Graphs
    g_W12Mimic_int[i] = (TGraphErrors*)f_W12Mimic_int->Get("gr_"+gr_Titles[i]);
    g_W12Mimic_int[i]->SetMarkerSize(1.5);
    
    g_StandardMimic_int[i] = (TGraphErrors*)f_StandardMimic_int->Get("gr_"+gr_Titles[i]);
    g_StandardMimic_int[i]->SetMarkerColor(kBlue);
    g_StandardMimic_int[i]->SetMarkerSize(1.5);
    
    //Add x-val offset
    if (i<3){
      Double_t *xvals = g_W12Mimic_int[i]->GetX();
      xvals[0] += i*offset;
      
      xvals = g_StandardMimic_int[i]->GetX();
      xvals[0] += i*offset+offset/5.0;
    }
    else{
      Double_t *xvals = g_W12Mimic_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;
      
      xvals = g_StandardMimic_int[i]->GetX();
      xvals[0] += i*offset + 2*offset+offset/5.0;
    }
  }//Loop over nGr

  for (Int_t pr=0; pr<nPro; pr++) {
    for (Int_t gr=0; gr<nGr; gr++) {
      (gr==0 && pr==0) ? g_Process[pr][gr]->Draw("AP") :
	g_Process[pr][gr]->Draw("Psame");
    }
  }

  Double_t yMax = 0.01;
  g_Process[0][0]->GetXaxis()->SetLimits(0.1, 0.5);
  g_Process[0][0]->GetYaxis()->SetRangeUser(-yMax, yMax);
  Double_t xmin = g_Process[0][0]->GetXaxis()->GetXmin();
  Double_t xmax = g_Process[0][0]->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
  li->SetLineColor(kBlack);
  li->SetLineStyle(8);
  li->Draw("same");

}
