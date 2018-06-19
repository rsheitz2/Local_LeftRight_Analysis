void ByTargPol_PhysBinned(){

  const Int_t nGr = 20;
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/April17/Data/MC_Data/";
  
  TFile* f_AMDY =
    TFile::Open(path+"AMDY_5bins.root");
  TGraphErrors* g_AMDY[nGr];
  
  TFile* f_JPsi =
    TFile::Open(path+"JPsi_5bins.root");
  TGraphErrors* g_JPsi[nGr];

  TFile* f_psi =
    TFile::Open(path+"psi_5bins.root");
  TGraphErrors* g_psi[nGr];

  TFile* f_OC =
    TFile::Open(path+"OC_5bins.root");
  TGraphErrors* g_OC[nGr];

  //Basic file checks
  if (!(f_AMDY) || !(f_JPsi) || !(f_psi) || !(f_OC) ){
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  const Int_t nPro = 4;
  TGraphErrors **g_Process[nPro] = {g_AMDY, g_JPsi, g_psi, g_OC};
  TString gr_Titles[nGr] = {"xN_LR_UpStream_Up", "xN_LR_UpStream_Down",
			    "xN_LR_DownStream_Up", "xN_LR_DownStream_Down",
			    "xPi_LR_UpStream_Up", "xPi_LR_UpStream_Down",
			    "xPi_LR_DownStream_Up", "xPi_LR_DownStream_Down",
			    "xF_LR_UpStream_Up", "xF_LR_UpStream_Down",
			    "xF_LR_DownStream_Up", "xF_LR_DownStream_Down",
			    "pT_LR_UpStream_Up", "pT_LR_UpStream_Down",
			    "pT_LR_DownStream_Up", "pT_LR_DownStream_Down",
			    "M_LR_UpStream_Up", "M_LR_UpStream_Down",
			    "M_LR_DownStream_Up", "M_LR_DownStream_Down"};
  
  //Int_t pDraw = 0;//AMDY
  //Int_t pDraw = 1;//JPsi
  //Int_t pDraw = 2;//psi
  Int_t pDraw = 3;//OC

  Double_t offset = 0.0015;
  gStyle->SetOptFit(1111);
  TCanvas* c1 = new TCanvas();
  c1->Divide(3, 2);
  for (Int_t i=0, ipad=1, icolor=1; i<nGr; i++, icolor++) {
    if (icolor==5) {
      icolor=1;
      ipad++;
    }
    
    //Get Graphs
    g_AMDY[i] = (TGraphErrors*)f_AMDY->Get("gr_"+gr_Titles[i]);
    g_AMDY[i]->SetMarkerColor(icolor);
    
    g_JPsi[i] = (TGraphErrors*)f_JPsi->Get("gr_"+gr_Titles[i]);
    g_JPsi[i]->SetMarkerColor(icolor);
    
    g_psi[i] = (TGraphErrors*)f_psi->Get("gr_"+gr_Titles[i]);
    g_psi[i]->SetMarkerColor(icolor);
    
    g_OC[i] = (TGraphErrors*)f_OC->Get("gr_"+gr_Titles[i]);
    g_OC[i]->SetMarkerColor(icolor);
        
    //Add x-val offset
    for (Int_t pr=0; pr<nPro; pr++) {
      Double_t *xvals = g_Process[pr][i]->GetX();
      for (Int_t nb=0; nb<g_Process[pr][i]->GetN(); nb++) {
	  xvals[nb] = (ipad<4) ?
	    xvals[nb] + icolor*offset : xvals[nb]+icolor*5*offset;
      }
    }
    
    c1->cd(ipad);
    if ( !(i%4) ) {
      g_Process[pDraw][i]->Draw("AP");
      Double_t yMax = 0.08;
      g_Process[pDraw][i]->GetYaxis()->SetRangeUser(-yMax, yMax);
      g_Process[pDraw][i]->GetXaxis()->SetNdivisions(7);
      
      Double_t xmin = g_Process[pDraw][i]->GetXaxis()->GetXmin();
      Double_t xmax = g_Process[pDraw][i]->GetXaxis()->GetXmax();
      TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
      li->SetLineColor(kBlack);
      li->SetLineStyle(8);
      li->Draw("same");
    }
    else {
      g_Process[pDraw][i]->Draw("Psame");
    }
    
  }//Loop over nGr

}
