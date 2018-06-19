void ByTarg_realIntegrated(){

  const Int_t nGr = 7;
  TString gr_Titles[] = {"xN_LR", "xN_LR_UpStream", "xN_LR_DownStream",
			 "xN_LR_UpStream_Up", "xN_LR_UpStream_Down",
			 "xN_LR_DownStream_Up", "xN_LR_DownStream_Down"};
			 
			 
		  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/";
  TFile* f_HM_int =
    TFile::Open(path+"April3/Data/RealData/ByTarget/HM/Wall_1bin.root");
  TGraphErrors* g_HM_int[nGr];
  
  TFile* f_JPsi_int =
    TFile::Open(path+"March27/Data/RealData/ByTarget/Wall_JPsi_1bin_SpinLR_withCorr.root");
    //TFile::Open(path+"April17/Data/RealData/Wall_JPsi_ByTarget_1bin.root");
  TGraphErrors* g_JPsi_int[nGr];

  //Basic file checks
  if ( !(f_HM_int) || !(f_JPsi_int) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  TGraphErrors **g_Process[] = {g_HM_int, g_JPsi_int};
  
  Double_t offset = 0.04;
  gStyle->SetOptFit(1111);
  TCanvas* c1 = new TCanvas();
  for (Int_t i=0, ipad=1; i<nGr; i++) {

    //Get Graphs
    g_HM_int[i] = (TGraphErrors*)f_HM_int->Get("gr_"+gr_Titles[i]);
    g_HM_int[i]->SetMarkerSize(1.5);
    
    g_JPsi_int[i] = (TGraphErrors*)f_JPsi_int->Get("gr_"+gr_Titles[i]);
    g_JPsi_int[i]->SetMarkerColor(kBlue);
    g_JPsi_int[i]->SetMarkerSize(1.5);

    //Add x-val offset
    if (i<3){
      Double_t *xvals = g_HM_int[i]->GetX();
      xvals[0] += i*offset;
      
      xvals = g_JPsi_int[i]->GetX();
      xvals[0] += i*offset;
    }
    else{
      Double_t *xvals = g_HM_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;
      
      xvals = g_JPsi_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;
    }
  }//Loop over nGr

  //Int_t pDraw = 0;//HM
  Int_t pDraw = 1;//JPsi
  for (Int_t i=0; i<nGr; i++) {
    (i==0) ? g_Process[pDraw][i]->Draw("AP") : g_Process[pDraw][i]->Draw("Psame");
  }

  Double_t yMax = 0.1;
  g_Process[pDraw][0]->GetXaxis()->SetLimits(0, 0.45);
  g_Process[pDraw][0]->GetYaxis()->SetRangeUser(-yMax, yMax);
  Double_t xmin = g_Process[pDraw][0]->GetXaxis()->GetXmin();
  Double_t xmax = g_Process[pDraw][0]->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
  li->SetLineColor(kBlack);
  li->SetLineStyle(8);
  li->Draw("same");
}
  
