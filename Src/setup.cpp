#include "common.h"
#include "setup.h"
#include "functions.h"

void SetupTGraph(TGraphErrors* gr, TString title, TString xTitle,
			  Double_t ymax){
  gr->SetMarkerStyle(21);
  gr->SetTitle(title);
  gr->GetXaxis()->SetTitle(xTitle);
    
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetLabelSize(0.06);
  gr->GetYaxis()->SetLabelSize(0.06);
  gr->GetXaxis()->SetTitleOffset(-0.4);
  //gr->GetXaxis()->SetTitleOffset(0.3);
  
  //gr->GetYaxis()->SetRangeUser(-1.0*ymax, ymax);
}


void SetupTLine(TLine* l){
  l->SetLineStyle(8);
  l->SetLineColor(kBlue);
}


void SetupHist(TH1D* h){
  h->GetXaxis()->SetTitleSize(0.08);
  h->GetXaxis()->SetLabelSize(0.08);
  h->GetYaxis()->SetLabelSize(0.08);
  h->SetLineColor(kBlack);
  h->GetYaxis()->SetTitleSize(0.08);
  h->GetXaxis()->SetTitleOffset(-0.4);
}
