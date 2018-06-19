void ByTarg_ByProcess_Integrated(){

  const Int_t nGr = 7;
  TString gr_Titles[] = {"xN_LR", "xN_LR_UpStream", "xN_LR_DownStream",
			 "xN_LR_UpStream_Up", "xN_LR_UpStream_Down",
			 "xN_LR_DownStream_Up", "xN_LR_DownStream_Down"};

  
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/April17/Data/MC_Data/";
  TFile* f_AMDY_int =
    TFile::Open(path+"AMDY_byTarget_1bin.root");
  TGraphErrors* g_AMDY_int[nGr];
  
  TFile* f_JPsi_int =
    TFile::Open(path+"JPsi_byTarget_1bin.root");
  TGraphErrors* g_JPsi_int[nGr];

  TFile* f_psi_int =
    TFile::Open(path+"psi_byTarget_1bin.root");
  TGraphErrors* g_psi_int[nGr];
  
  TFile* f_OC_int =
    TFile::Open(path+"OC_byTarget_1bin.root");
  TGraphErrors* g_OC_int[nGr];

  TFile* f_AMDY_OC_JPsi_psi_int =
    TFile::Open(path+"AMDY_OC_JPsi_psi_1bin.root");
  TGraphErrors* g_AMDY_OC_JPsi_psi_int[nGr];

  //Basic file checks
  if ( !(f_AMDY_int) || !(f_JPsi_int)
       || !(f_psi_int) || !(f_OC_int) || !(f_AMDY_OC_JPsi_psi_int) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  const Int_t nPro = 5;
  TGraphErrors **g_Process[nPro] = {g_AMDY_int, g_JPsi_int, g_psi_int, g_OC_int,
				g_AMDY_OC_JPsi_psi_int};
  
  Double_t offset = 0.04;
  gStyle->SetOptFit(1111);
  TCanvas* c1 = new TCanvas();
  for (Int_t i=0, ipad=1; i<nGr; i++) {

    //Get Graphs
    g_AMDY_int[i] = (TGraphErrors*)f_AMDY_int->Get("gr_"+gr_Titles[i]);
    g_AMDY_int[i]->SetMarkerSize(1.5);
    
    g_JPsi_int[i] = (TGraphErrors*)f_JPsi_int->Get("gr_"+gr_Titles[i]);
    g_JPsi_int[i]->SetMarkerColor(kBlue);
    g_JPsi_int[i]->SetMarkerSize(1.5);

    g_psi_int[i] = (TGraphErrors*)f_psi_int->Get("gr_"+gr_Titles[i]);
    g_psi_int[i]->SetMarkerColor(kPink);
    g_psi_int[i]->SetMarkerSize(1.5);

    g_OC_int[i] = (TGraphErrors*)f_OC_int->Get("gr_"+gr_Titles[i]);
    g_OC_int[i]->SetMarkerColor(kGreen);
    g_OC_int[i]->SetMarkerSize(1.5);

    g_AMDY_OC_JPsi_psi_int[i] = (TGraphErrors*)f_AMDY_OC_JPsi_psi_int->Get("gr_"+gr_Titles[i]);
    g_AMDY_OC_JPsi_psi_int[i]->SetMarkerColor(9);//purple
    g_AMDY_OC_JPsi_psi_int[i]->SetMarkerSize(1.5);
    
    //Add x-val offset
    if (i<3){
      Double_t *xvals = g_AMDY_int[i]->GetX();
      xvals[0] += i*offset;
      
      xvals = g_JPsi_int[i]->GetX();
      xvals[0] += i*offset;

      xvals = g_psi_int[i]->GetX();
      xvals[0] += i*offset;

      xvals = g_OC_int[i]->GetX();
      xvals[0] += i*offset;

      xvals = g_AMDY_OC_JPsi_psi_int[i]->GetX();
      xvals[0] += i*offset;
    }
    else{
      Double_t *xvals = g_AMDY_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;
      
      xvals = g_JPsi_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;

      xvals = g_psi_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;

      xvals = g_OC_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;

      xvals = g_AMDY_OC_JPsi_psi_int[i]->GetX();
      xvals[0] += i*offset + 2*offset;
    }
  }//Loop over nGr

  //Int_t pDraw = 0;//AMDY
  //Int_t pDraw = 1;//JPsi
  //Int_t pDraw = 2;//psi
  //Int_t pDraw = 3;//OC
  Int_t pDraw = 4;//AMDY_OC_JPsi_psi
  
  for (Int_t i=0; i<nGr; i++) {
    (i==0) ? g_Process[pDraw][i]->Draw("AP") : g_Process[pDraw][i]->Draw("Psame");
  }

  Double_t yMax = 0.04;
  g_Process[pDraw][0]->GetXaxis()->SetLimits(0, 0.45);
  g_Process[pDraw][0]->GetYaxis()->SetRangeUser(-yMax, yMax);
  Double_t xmin = g_Process[pDraw][0]->GetXaxis()->GetXmin();
  Double_t xmax = g_Process[pDraw][0]->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
  li->SetLineColor(kBlack);
  li->SetLineStyle(8);
  li->Draw("same");

}
