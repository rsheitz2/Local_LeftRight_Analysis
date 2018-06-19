void ByTarg_ByProcess(){

  const Int_t nGr = 15;
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/April17/Data/MC_Data/";
  TFile* f_AMDY =
    TFile::Open(path+"AMDY_byTarget_10bins.root");
  TGraphErrors* g_AMDY[nGr];

  TFile* f_AMDY_int =
    TFile::Open(path+"AMDY_byTarget_1bin.root");
  TGraphErrors* g_AMDY_int[nGr];

  
  TFile* f_JPsi =
    TFile::Open(path+"JPsi_byTarget_10bins.root");
  TGraphErrors* g_JPsi[nGr];

  TFile* f_JPsi_int =
    TFile::Open(path+"JPsi_byTarget_1bin.root");
  TGraphErrors* g_JPsi_int[nGr];

  
  TFile* f_psi =
    TFile::Open(path+"psi_byTarget_10bins.root");
  TGraphErrors* g_psi[nGr];

  TFile* f_psi_int =
    TFile::Open(path+"psi_byTarget_1bin.root");
  TGraphErrors* g_psi_int[nGr];
  

  TFile* f_OC =
    TFile::Open(path+"OC_byTarget_10bins.root");
  TGraphErrors* g_OC[nGr];

  TFile* f_OC_int =
    TFile::Open(path+"OC_byTarget_1bin.root");
  TGraphErrors* g_OC_int[nGr];

  //Basic file checks
  if (!(f_AMDY) || !(f_AMDY_int) || !(f_JPsi) || !(f_JPsi_int)
      || !(f_psi) || !(f_psi_int) || !(f_OC) || !(f_OC_int) ){
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
    g_AMDY[i] = (TGraphErrors*)f_AMDY->Get("gr_"+gr_Titles[i]);
    g_AMDY_int[i] = (TGraphErrors*)f_AMDY_int->Get("gr_"+gr_Titles[i]);
    

    g_JPsi[i] = (TGraphErrors*)f_JPsi->Get("gr_"+gr_Titles[i]);
    g_JPsi[i]->SetMarkerColor(kBlue);
    g_JPsi_int[i] = (TGraphErrors*)f_JPsi_int->Get("gr_"+gr_Titles[i]);
    g_JPsi_int[i]->SetMarkerColor(kBlue);

    g_psi[i] = (TGraphErrors*)f_psi->Get("gr_"+gr_Titles[i]);
    g_psi[i]->SetMarkerColor(kPink);
    g_psi_int[i] = (TGraphErrors*)f_psi_int->Get("gr_"+gr_Titles[i]);
    g_psi_int[i]->SetMarkerColor(kPink);

    g_OC[i] = (TGraphErrors*)f_OC->Get("gr_"+gr_Titles[i]);
    g_OC[i]->SetMarkerColor(kGreen);
    g_OC_int[i] = (TGraphErrors*)f_OC_int->Get("gr_"+gr_Titles[i]);
    g_OC_int[i]->SetMarkerColor(kGreen);
    
    //Add x-val offset
    if (ipad<4){
      Double_t *xvals = g_JPsi[i]->GetX();
      for (Int_t nb=0; nb<g_JPsi[i]->GetN(); nb++) xvals[nb] += offset;
      xvals = g_JPsi_int[i]->GetX();
      xvals[0] += offset;

      xvals = g_psi[i]->GetX();
      for (Int_t nb=0; nb<g_psi[i]->GetN(); nb++) xvals[nb] += 2*offset;
      xvals = g_psi_int[i]->GetX();
      xvals[0] += offset;

      xvals = g_OC[i]->GetX();
      for (Int_t nb=0; nb<g_OC[i]->GetN(); nb++) xvals[nb] += 2*offset;
      xvals = g_OC_int[i]->GetX();
      xvals[0] += offset;
    }
    else{
      Double_t *xvals = g_JPsi[i]->GetX();
      for (Int_t nb=0; nb<g_JPsi[i]->GetN(); nb++) xvals[nb] += offset2;
      xvals = g_JPsi_int[i]->GetX();
      xvals[0] += offset2;

      xvals = g_psi[i]->GetX();
      for (Int_t nb=0; nb<g_psi[i]->GetN(); nb++) xvals[nb] += 2*offset2;
      xvals = g_psi_int[i]->GetX();
      xvals[0] += offset2;

      xvals = g_OC[i]->GetX();
      for (Int_t nb=0; nb<g_OC[i]->GetN(); nb++) xvals[nb] += 2*offset2;
      xvals = g_OC_int[i]->GetX();
      xvals[0] += offset2;
    }
    

    //Draw Graphs
    //if ( !( (i)%3 ) ) {cout << "UpStream" << endl;//Upstream
    //if ( !( (i-1)%3 ) ) {cout << "DownStream" <<endl;//Downstream
    if ( !( (i+1)%3 ) ) {cout << "Combined targets" << endl;//Combined targets
      cout << "Drawing graph #: " << i << "    on pad #: " << ipad << endl;
      c1->cd(ipad); ipad++;
      
      g_AMDY[i]->Draw("AP");
      Double_t yMax = 0.04;
      g_AMDY[i]->GetYaxis()->SetRangeUser(-yMax, yMax);
      
      g_JPsi[i]->Draw("Psame");
      g_psi[i]->Draw("Psame");
      g_OC[i]->Draw("Psame");
 
      Double_t xmin = g_AMDY[i]->GetXaxis()->GetXmin();
      Double_t xmax = g_AMDY[i]->GetXaxis()->GetXmax();
      TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
      li->SetLineColor(kBlack);
      li->SetLineStyle(8);
      li->Draw("same");
    }
    
  }//Loop over nGr

  c1->cd(6);
  Int_t iDraw = 2;//Combined targets
  g_AMDY_int[iDraw]->Draw("AP");

  g_JPsi_int[iDraw]->Draw("Psame");
  g_psi_int[iDraw]->Draw("Psame");
  g_OC_int[iDraw]->Draw("Psame");

  Double_t yMax = 0.015;
  g_AMDY_int[iDraw]->GetYaxis()->SetRangeUser(-yMax, yMax);
  Double_t xmin = g_AMDY_int[iDraw]->GetXaxis()->GetXmin();
  Double_t xmax = g_AMDY_int[iDraw]->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
  li->SetLineColor(kBlack);
  li->SetLineStyle(8);
  li->Draw("same");

}
