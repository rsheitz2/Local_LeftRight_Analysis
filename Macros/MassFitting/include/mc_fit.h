#ifndef MC_FIT_H
#define MC_FIT_H
//Reconstructed Monte Carlo with:  JPsi, Psi', OC, Drell-Yan
TH1D *hist[4] = {NULL, NULL, NULL, NULL}; //{JPsi, psi, OC, AMDY}

Double_t Fit_MC(Double_t *x, Double_t *par){
  Double_t nJPsi = hist[0]->GetBinContent(hist[0]->FindBin(x[0] ) );
  Double_t npsi = hist[1]->GetBinContent(hist[1]->FindBin(x[0] ) );
  Double_t nOC = hist[2]->GetBinContent(hist[2]->FindBin(x[0] ) );
  Double_t nAMDY = hist[3]->GetBinContent(hist[3]->FindBin(x[0] ) );
    
  return par[0]*nJPsi +par[1]*npsi +par[2]*nOC +par[3]*nAMDY;
}

#endif
