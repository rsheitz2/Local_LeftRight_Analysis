void ByTarg_Mimic(){

  const Int_t nGr = 15;
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/April17/Data/MC_Data/";
  TFile* f_W12Mimic =
    TFile::Open(path+"W12Mimic_5bins.root");
  TGraphErrors* g_W12Mimic[nGr];

  TFile* f_W12Mimic_int =
    TFile::Open(path+"W12Mimic_1bin.root");
  TGraphErrors* g_W12Mimic_int[nGr];

  
  TFile* f_StandardMimic =
    TFile::Open(path+"StandardMimic_5bins.root");
  TGraphErrors* g_StandardMimic[nGr];

  TFile* f_StandardMimic_int =
    TFile::Open(path+"StandardMimic_1bin.root");
  TGraphErrors* g_StandardMimic_int[nGr];


  //Basic file checks
  if (!(f_W12Mimic) || !(f_W12Mimic_int)
      || !(f_StandardMimic) || !(f_StandardMimic_int) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  TString gr_Titles[] = {"xN_LR_UpStream", "xN_LR_DownStream", "xN_LR",
			 "xPi_LR_UpStream", "xPi_LR_DownStream", "xPi_LR",
			 "xF_LR_UpStream", "xF_LR_DownStream", "xF_LR",
			 "pT_LR_UpStream", "pT_LR_DownStream", "pT_LR",
			 "M_LR_UpStream", "M_LR_DownStream", "M_LR"};
  
  Double_t offset = 0.004;
  Double_t offset2 = 0.03;

  gStyle->SetOptFit(1111);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3, 2);
  for (Int_t i=0, ipad=1; i<nGr; i++) {

    //Get Graphs
    g_W12Mimic[i] = (TGraphErrors*)f_W12Mimic->Get("gr_"+gr_Titles[i]);
    g_W12Mimic_int[i] = (TGraphErrors*)f_W12Mimic_int->Get("gr_"+gr_Titles[i]);
    

    g_StandardMimic[i] = (TGraphErrors*)f_StandardMimic->Get("gr_"+gr_Titles[i]);
    g_StandardMimic[i]->SetMarkerColor(kBlue);
    g_StandardMimic_int[i] = (TGraphErrors*)f_StandardMimic_int->Get("gr_"+gr_Titles[i]);
    g_StandardMimic_int[i]->SetMarkerColor(kBlue);

    //Add x-val offset
    if (ipad<4){
      Double_t *xvals = g_StandardMimic[i]->GetX();
      for (Int_t nb=0; nb<g_StandardMimic[i]->GetN(); nb++) xvals[nb] += offset;
      xvals = g_StandardMimic_int[i]->GetX();
      xvals[0] += offset;
    }
    else{
      Double_t *xvals = g_StandardMimic[i]->GetX();
      for (Int_t nb=0; nb<g_StandardMimic[i]->GetN(); nb++) xvals[nb] += offset2;
      xvals = g_StandardMimic_int[i]->GetX();
      xvals[0] += offset2;
    }
    

    //Draw Graphs
    //if ( !( (i)%3 ) ) {cout << "UpStream" << endl;//Upstream
    if ( !( (i-1)%3 ) ) {cout << "DownStream" <<endl;//Downstream
    //if ( !( (i+1)%3 ) ) {cout << "Combined targets" << endl;//Combined targets
      cout << "Drawing graph #: " << i << "    on pad #: " << ipad << endl;
      c1->cd(ipad); ipad++;
      
      g_W12Mimic[i]->Draw("AP");
      Double_t yMax = 0.01;
      g_W12Mimic[i]->GetYaxis()->SetRangeUser(-yMax, yMax);
      
      g_StandardMimic[i]->Draw("Psame");
 
      Double_t xmin = g_W12Mimic[i]->GetXaxis()->GetXmin();
      Double_t xmax = g_W12Mimic[i]->GetXaxis()->GetXmax();
      TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
      li->SetLineColor(kBlack);
      li->SetLineStyle(8);
      li->Draw("same");
    }
    
  }//Loop over nGr

  c1->cd(6);
  Int_t iDraw = 2;//Combined targets
  g_W12Mimic_int[iDraw]->Draw("AP");

  g_StandardMimic_int[iDraw]->Draw("Psame");

  Double_t yMax = 0.01;
  g_W12Mimic_int[iDraw]->GetYaxis()->SetRangeUser(-yMax, yMax);
  Double_t xmin = g_W12Mimic_int[iDraw]->GetXaxis()->GetXmin();
  Double_t xmax = g_W12Mimic_int[iDraw]->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
  li->SetLineColor(kBlack);
  li->SetLineStyle(8);
  li->Draw("same");

}
