void DrawSame(){

  const Int_t nGr = 15;
  TString path = "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/Presents/";
  TFile* f_rec =
    TFile::Open(path+"March27/Data/ByTarget/YuShiangMC_HM_run1.1_7bins.root");
  TGraphErrors* g_rec[nGr];

  TFile* f_GenWrec =
    TFile::Open(path+"March27/Data/ByTarget/YuShiangMC_HM_run1.1_GenPhiS_7bins.root");
  TGraphErrors* g_GenWrec[nGr];

  TFile* f_Gen = //Not drawn
    TFile::Open(path+"March13/Data/Generated/ByTarget/gen_Charles_HMDY_byTarget_PhiS_3bins.root");
  TGraphErrors* g_Gen[nGr];

  const Int_t nGr_real = 35;
  TFile* f_Real = //Not drawn
    //TFile::Open(path+"March13/Data/RealData/ByTarget/Wall_HM_SpinInfluenced_noCorrections_LR_3bins.root");
    TFile::Open(path+"March13/Data/RealData/ByTarget/Wall_HM_SpinInfluenced_LR_3bins.root");
  TGraphErrors* g_Real[nGr_real];

  TFile *f_Real_int =
    TFile::Open(path+"March27/Data/RealData/ByTarget/Wall_HM_1bin_SpinLR_withCorr.root");
  TGraphErrors* g_Real_int[nGr_real];

  //Basic file checks
  if (!(f_rec) || !(f_GenWrec) || !(f_Gen) || !(f_Real) ){
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
    g_rec[i] = (TGraphErrors*)f_rec->Get("gr_"+gr_Titles[i]);
    g_rec[i]->SetMarkerColor(kGreen);

    g_GenWrec[i] = (TGraphErrors*)f_GenWrec->Get("gr_"+gr_Titles[i]);
    g_GenWrec[i]->SetMarkerColor(kRed);
    
    g_Gen[i] = (TGraphErrors*)f_Gen->Get("gr_"+gr_Titles[i]);
    g_Gen[i]->SetMarkerColor(kBlue);

    g_Real[i] = (TGraphErrors*)f_Real->Get("gr_"+gr_Titles[i]);
    //g_Real[i]->SetMarkerColor(kBlue);

    g_Real_int[i] = (TGraphErrors*)f_Real_int->Get("gr_"+gr_Titles[i]);
    
    //Add x-val offset
    if (ipad<4){
      Double_t *xvals = g_GenWrec[i]->GetX();
      for (Int_t nb=0; nb<g_GenWrec[i]->GetN(); nb++) xvals[nb] += offset;
    
      xvals = g_Gen[i]->GetX();
      for (Int_t nb=0; nb<g_Gen[i]->GetN(); nb++) xvals[nb] += 2*offset;

      //xvals = g_rec[i]->GetX();
      //for (Int_t nb=0; nb<g_rec[i]->GetN(); nb++) xvals[nb] += 3*offset;
    }
    else{
      Double_t *xvals = g_GenWrec[i]->GetX();
      for (Int_t nb=0; nb<g_GenWrec[i]->GetN(); nb++) xvals[nb] += offset2;
    
      xvals = g_Gen[i]->GetX();
      for (Int_t nb=0; nb<g_Gen[i]->GetN(); nb++) xvals[nb] += 2*offset2;

      //xvals = g_rec[i]->GetX();
      //for (Int_t nb=0; nb<g_rec[i]->GetN(); nb++) xvals[nb] += 3*offset2;
    }
    

    //Draw Graphs
    //if ( !( (i)%3 ) ) {cout << "UpStream" << endl;//Upstream
    //if ( !( (i-1)%3 ) ) {cout << "DownStream" <<endl;//Downstream
    if ( !( (i+1)%3 ) ) {cout << "Combined targets" << endl;//Combined targets
      cout << "Drawing graph #: " << i << "    on pad #: " << ipad << endl;
      
      c1->cd(ipad); ipad++;
      
      g_rec[i]->Draw("AP");
      g_rec[i]->GetYaxis()->SetRangeUser(-0.2, 0.2);
      //g_rec[i]->GetYaxis()->SetRangeUser(-0.035, 0.035);
      //g_rec[i]->GetYaxis()->SetRangeUser(-0.045, 0.045);
      
      
      g_GenWrec[i]->Draw("Psame");
      g_Real[i]->Draw("Psame");
      //g_Gen[i]->Draw("Psame");
      //g_Real[i]->Draw("Psame");
      
      //Draw Gen only and fit 
      //g_Gen[i]->Draw("AP");
      //g_Gen[i]->GetYaxis()->SetRangeUser(-0.03, 0.03);
      //g_Gen[i]->Fit("pol0");
 
      Double_t xmin = g_rec[i]->GetXaxis()->GetXmin();
      Double_t xmax = g_rec[i]->GetXaxis()->GetXmax();
      TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
      li->SetLineColor(kBlack);
      li->SetLineStyle(8);
      li->Draw("same");
    }
    
  }//Loop over nGr

  c1->cd(6);
  g_Real_int[3]->Draw("AP");
  g_Real_int[3]->GetYaxis()->SetRangeUser(-0.2, 0.2);
  Double_t xmin = g_Real_int[3]->GetXaxis()->GetXmin();
  Double_t xmax = g_Real_int[3]->GetXaxis()->GetXmax();
  TLine *li = new TLine(xmin, 0.0, xmax, 0.0);
  li->SetLineColor(kBlack);
  li->SetLineStyle(8);
  li->Draw("same");

}
