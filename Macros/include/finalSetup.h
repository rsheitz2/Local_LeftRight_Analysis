void FinalSetup(TGraphErrors *g){
  g->GetXaxis()->SetLabelFont(132);
  g->GetYaxis()->SetLabelFont(132);
}


void FinalClearTitles(TGraphErrors *g){
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");  
}


void FinalSetup(TGraph *g){
  g->GetXaxis()->SetLabelFont(132);
  g->GetYaxis()->SetLabelFont(132);
}


void FinalClearTitles(TGraph *g){
  g->SetTitle("");
  g->GetXaxis()->SetTitle("");
  g->GetYaxis()->SetTitle("");  
}
