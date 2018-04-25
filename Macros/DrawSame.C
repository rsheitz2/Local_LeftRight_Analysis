void DrawSame(){

  const Int_t nGr = 15;
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/";
  TFile* f_new =
    TFile::Open(path+"April17/Yu_Wall_AMDY_gt4_3bins.root");
  TGraphErrors* g_new[nGr];

  TFile* f_downFlip =
    TFile::Open(path+"April17/Yu_Wall_AMDY_gt3_3bins_downSflip.root");
  TGraphErrors* g_downFlip[nGr];

  TFile* f_upFlip =
    TFile::Open(path+"April17/Yu_Wall_AMDY_gt3_3bins_upSflip.root");
  TGraphErrors* g_upFlip[nGr];

  TFile* f_old =
    TFile::Open(path+"April3/Data/MC_Data/Yu_BW/Wall_AMDY_gt4GeV_3bins_2.root");
  TGraphErrors* g_old[nGr];

  //Basic file checks
  if (!(f_new) || !(f_downFlip) || !(f_upFlip) || !(f_old) ){
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
    g_new[i] = (TGraphErrors*)f_new->Get("gr_"+gr_Titles[i]);
    //g_new[i]->SetMarkerColor(kGreen);

    g_old[i] = (TGraphErrors*)f_old->Get("gr_"+gr_Titles[i]);
    g_old[i]->SetMarkerColor(kGreen);

    g_downFlip[i] = (TGraphErrors*)f_downFlip->Get("gr_"+gr_Titles[i]);
    g_downFlip[i]->SetMarkerColor(kRed);

    g_upFlip[i] = (TGraphErrors*)f_upFlip->Get("gr_"+gr_Titles[i]);
    g_upFlip[i]->SetMarkerColor(kBlue);
    
    //Add x-val offset
    if (ipad<4){
      Double_t *xvals = g_downFlip[i]->GetX();
      for (Int_t nb=0; nb<g_downFlip[i]->GetN(); nb++) xvals[nb] += offset;

      xvals = g_upFlip[i]->GetX();
      for (Int_t nb=0; nb<g_upFlip[i]->GetN(); nb++) xvals[nb] += 2*offset;

      xvals = g_old[i]->GetX();
      for (Int_t nb=0; nb<g_old[i]->GetN(); nb++) xvals[nb] += 2*offset;
    
      //xvals = g_new[i]->GetX();
      //for (Int_t nb=0; nb<g_new[i]->GetN(); nb++) xvals[nb] += 3*offset;
    }
    else{
      Double_t *xvals = g_downFlip[i]->GetX();
      for (Int_t nb=0; nb<g_downFlip[i]->GetN(); nb++) xvals[nb] += offset2;

      xvals = g_upFlip[i]->GetX();
      for (Int_t nb=0; nb<g_upFlip[i]->GetN(); nb++) xvals[nb] += 2*offset2;

      xvals = g_old[i]->GetX();
      for (Int_t nb=0; nb<g_old[i]->GetN(); nb++) xvals[nb] += 2*offset2;
    
      //xvals = g_new[i]->GetX();
      //for (Int_t nb=0; nb<g_new[i]->GetN(); nb++) xvals[nb] += 3*offset2;
    }
    

    //Draw Graphs
    //if ( !( (i)%3 ) ) {cout << "UpStream" << endl;//Upstream
    //if ( !( (i-1)%3 ) ) {cout << "DownStream" <<endl;//Downstream
    if ( !( (i+1)%3 ) ) {cout << "Combined targets" << endl;//Combined targets
      cout << "Drawing graph #: " << i << "    on pad #: " << ipad << endl;
      
      c1->cd(ipad); ipad++;
      
      g_new[i]->Draw("AP");
      Double_t yMax = 0.05;
      g_new[i]->GetYaxis()->SetRangeUser(-yMax, yMax);
      
      //g_downFlip[i]->Draw("Psame");
      //g_upFlip[i]->Draw("Psame");
      g_old[i]->Draw("Psame");
 
      Double_t xmin = g_new[i]->GetXaxis()->GetXmin();
      Double_t xmax = g_new[i]->GetXaxis()->GetXmax();
      TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
      li->SetLineColor(kBlack);
      li->SetLineStyle(8);
      li->Draw("same");
    }
    
  }//Loop over nGr

}
