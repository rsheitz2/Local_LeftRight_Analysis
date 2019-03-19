#ifndef FINALSETUP_H
#define FINALSETUP_H

void FinalSetup(TGraphErrors *g){
  g->GetXaxis()->SetLabelFont(133);
  g->GetXaxis()->SetLabelSize(25);
  g->GetXaxis()->SetLabelOffset(0.0);
    
  g->GetYaxis()->SetLabelFont(133);
  g->GetYaxis()->SetLabelSize(25);
}


void FinalSetup(TGraph *g){
  g->GetXaxis()->SetLabelFont(133);
  g->GetXaxis()->SetLabelSize(25);
  g->GetXaxis()->SetLabelOffset(0.0);
    
  g->GetYaxis()->SetLabelFont(133);
  g->GetYaxis()->SetLabelSize(25);
}


void FinalSetup(TH1D *h){
  gPad->Update();
  TPaveStats *ps = (TPaveStats*)h->FindObject("stats");
  ps->SetTextFont(132);

  h->GetXaxis()->SetLabelFont(133);
  h->GetXaxis()->SetLabelSize(25);
  h->GetXaxis()->SetLabelOffset(0.0);
    
  h->GetYaxis()->SetLabelFont(133);
  h->GetYaxis()->SetLabelSize(25);
}


void FinalClearTitles(TGraphErrors *g){
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");  
}


void FinalClearTitles(TGraph *g){
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");  
}


void FinalClearTitles(TH1D *h){
  h->SetTitle("");
  h->GetXaxis()->SetTitle("");
  h->GetYaxis()->SetTitle("");  
}


void SetTitleName(TGraphErrors *g, TString name, TString whichAxis="x"){
  if (whichAxis=="x"){
    g->GetXaxis()->SetTitle(name);
    g->GetXaxis()->SetTitleFont(133);
    g->GetXaxis()->SetTitleSize(25);
    g->GetXaxis()->CenterTitle(true);
  }
  else if (whichAxis=="y"){
    TLatex *latex = new TLatex();
    latex->SetTextSize(25);
    latex->SetTextFont(133);
    //latex->DrawLatexNDC(0.1, 0.6, name);
    latex->DrawLatexNDC(0.01, 0.6, name);
  }
}


void SetTitleName(TGraphAsymmErrors *g, TString name, TString whichAxis="x"){
  if (whichAxis=="x"){
    g->GetXaxis()->SetTitle(name);
    g->GetXaxis()->SetTitleFont(133);
    g->GetXaxis()->SetTitleSize(25);
    g->GetXaxis()->CenterTitle(true);
  }
  else if (whichAxis=="y"){
    TLatex *latex = new TLatex();
    latex->SetTextSize(25);
    latex->SetTextFont(133);
    //latex->DrawLatexNDC(0.1, 0.6, name);
    latex->DrawLatexNDC(0.01, 0.6, name);
  }
}


void SetTitleName(TGraph *g, TString name, TString whichAxis="x"){
  if (whichAxis=="x"){
    g->GetXaxis()->SetTitle(name);
    g->GetXaxis()->SetTitleFont(133);
    g->GetXaxis()->SetTitleSize(25);
    g->GetXaxis()->CenterTitle(true);
  }
  else if (whichAxis=="y"){
    TLatex *latex = new TLatex();
    latex->SetTextSize(25);
    latex->SetTextFont(133);
    //latex->DrawLatexNDC(0.1, 0.6, name);
    latex->DrawLatexNDC(0.01, 0.6, name);
  }
}


void SetTitleName(TH1D *h, TString name, TString whichAxis="x"){
  if (whichAxis=="x"){
    h->GetXaxis()->SetTitle(name);
    h->GetXaxis()->SetTitleFont(133);
    h->GetXaxis()->SetTitleSize(25);
    h->GetXaxis()->CenterTitle(true);
  }
  else if (whichAxis=="y"){
    TLatex *latex = new TLatex();
    latex->SetTextSize(25);
    latex->SetTextFont(133);
    latex->DrawLatexNDC(0.1, 0.6, name);
  }
}

#endif
