//Aestheitics
void SetUp(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);

  h->GetXaxis()->SetNdivisions(504);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);
}


void SetUp(TGraphErrors* g, Int_t istyle, Int_t icolor){
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);

  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  g->SetMarkerSize(1.3);
  g->SetMarkerStyle(istyle);
  g->SetMarkerColor(icolor);
}


void SetUp(TGraph* g){
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);

  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  g->SetMarkerSize(1.3);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(4);
}


void SetUp(TGraph2D* g){
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);

  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  g->SetMarkerSize(1.3);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(4);
}


void SetUpTGraph(TGraphErrors* g, Int_t icolor=1,
		 Double_t offset=0.0, Int_t nBins=0){
  g->SetMarkerStyle(21);
  g->SetMarkerColor(icolor);


  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);

  if (nBins != 0){
    Double_t *xval = g->GetX();
    for (Int_t i=0; i<nBins; i++) xval[i] += offset;
  }
}


void SetUpTGraph(TGraph* g){
  g->SetMarkerStyle(21);
  
  g->GetYaxis()->SetNdivisions(504);
  g->GetYaxis()->SetLabelFont(22);
  g->GetYaxis()->SetLabelSize(0.08);
  
  g->GetXaxis()->SetNdivisions(504);
  g->GetXaxis()->SetLabelFont(22);
  g->GetXaxis()->SetLabelSize(0.08);
}


void SetUpHist(TH1D* h){
  h->GetYaxis()->SetNdivisions(504);
  h->GetYaxis()->SetLabelFont(22);
  h->GetYaxis()->SetLabelSize(0.08);
  
  h->GetXaxis()->SetNdivisions(508);
  h->GetXaxis()->SetLabelFont(22);
  h->GetXaxis()->SetLabelSize(0.08);
}


void SetUpTF(TF1* f){
  f->GetYaxis()->SetNdivisions(504);
  f->GetYaxis()->SetLabelFont(22);
  f->GetYaxis()->SetLabelSize(0.08);
  
  f->GetXaxis()->SetNdivisions(508);
  f->GetXaxis()->SetLabelFont(22);
  f->GetXaxis()->SetLabelSize(0.08);
}

//Functions
Double_t f_CrystalBall(Double_t *x, Double_t *par){
  //par[3] =alpha, par[4] =n
  Double_t A = TMath::Power(par[4]/par[3], par[4])*TMath::Exp(-par[3]*par[3]/2.0);
  Double_t B = par[4]/par[3] - par[3];
  Double_t C = (par[4]/par[3])*(1.0/(par[4]-1.0))*TMath::Exp(-par[3]*par[3]/2.0);
  Double_t D = TMath::Pi()*(1+TMath::Erf( par[3]/TMath::Sqrt(2)) );
  
  Double_t Norm = 1.0/(par[2]*(C+D) );

  Double_t arg = (x[0] - par[1])/par[2];
  if (arg > -par[3] ) return par[0]*Norm*TMath::Exp(-0.5*arg*arg);
  else return par[0]*Norm*A*TMath::Power((B - arg), -par[1]);
}


Double_t f_Gauss(Double_t *x, Double_t *par){
  
  Double_t arg = (x[0] - par[1])/par[2];
  return par[0]*TMath::Exp(-0.5*arg*arg);
}


//Drawing
void DrawLine(TH1D *h, Double_t yval){
  Double_t min_x = h->GetXaxis()->GetXmin();	
  Double_t max_x = h->GetXaxis()->GetXmax();	
  TLine* li = new TLine(min_x, yval, max_x, yval);
  
  li->SetLineColor(1);
  li->SetLineStyle(8);
  li->SetLineWidth(2);
  
  li->Draw("same");
}


void DrawLineY(TH1D *h, Double_t xval, Double_t min_y, Double_t max_y,
	       Int_t icolor){
  TLine* li = new TLine(xval, min_y, xval, max_y);
  
  li->SetLineColor(icolor);
  li->SetLineStyle(8);
  li->SetLineWidth(2);
  
  li->Draw("same");
}


void Draw_2Gauss_3Expo(Double_t A_JPsi, Double_t M_JPsi, Double_t W_JPsi,
		       Double_t A_psi, Double_t M_psi, Double_t W_psi,
		       Double_t A_e1, Double_t c_e1,
		       Double_t A_e2, Double_t c_e2,
		       Double_t A_e3, Double_t c_e3,
		       Double_t Mmin, Double_t Mmax){

  TF1 *f_JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_JPsi->SetParameters(A_JPsi, M_JPsi, W_JPsi);
  f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_psi->SetParameters(A_psi, M_psi, W_psi);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");
  
  TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-[2]))", Mmin,Mmax);
  f_CombBg->SetParameters(A_e1, c_e1, Mmin); 
  f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

  TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-[2]))", Mmin, Mmax);
  f_OC->SetParameters(A_e2, c_e2, Mmin);
  f_OC->SetLineColor(6); f_OC->Draw("same");

  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  f_DY->SetParameters(A_e3, c_e3, Mmin);
  f_DY->SetLineColor(kBlue); f_DY->Draw("same");
}


void Draw_2Gauss_3Expo(Double_t *pars, Double_t Mmin, Double_t Mmax){

  TF1 *f_JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_JPsi->SetParameters(pars[0], pars[1], pars[2]);
  f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_psi->SetParameters(pars[3], pars[4], pars[5]);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");
  
  TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-[2]))", Mmin,Mmax);
  f_CombBg->SetParameters(pars[6], pars[7], Mmin); 
  f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

  TF1 *f_OC = new TF1("f_OC", "[0]*TMath::Exp([1]*(x-[2]))", Mmin, Mmax);
  f_OC->SetParameters(pars[8], pars[9], Mmin);
  f_OC->SetLineColor(6); f_OC->Draw("same");

  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  f_DY->SetParameters(pars[10], pars[11], Mmin);
  f_DY->SetLineColor(kBlue); f_DY->Draw("same");
}


void Draw_2Gauss_2Expo(Double_t A_JPsi, Double_t M_JPsi, Double_t W_JPsi,
		       Double_t A_psi, Double_t M_psi, Double_t W_psi,
		       Double_t A_e1, Double_t c_e1,
		       Double_t A_e2, Double_t c_e2,
		       Double_t Mmin, Double_t Mmax){

  TF1 *f_JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_JPsi->SetParameters(A_JPsi, M_JPsi, W_JPsi);
  f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_psi->SetParameters(A_psi, M_psi, W_psi);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");
  
  TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-[2]))", Mmin,Mmax);
  f_CombBg->SetParameters(A_e1, c_e1, Mmin); 
  f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  f_DY->SetParameters(A_e2, c_e2, Mmin);
  f_DY->SetLineColor(kBlue); f_DY->Draw("same");
}


void Draw_2Gauss_2Expo(Double_t *pars, Double_t Mmin, Double_t Mmax){

  TF1 *f_JPsi =
    new TF1("f_JPsi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_JPsi->SetParameters(pars[0], pars[1], pars[2]);
  f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

  TF1 *f_psi =
    new TF1("f_psi", "[0]*TMath::Exp(  -0.5*( ([1]-x)*([1]-x) )/( [2]*[2]))",
	    0, Mmax);
  f_psi->SetParameters(pars[3], pars[4], pars[5]);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");
  
  TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-[2]))", Mmin,Mmax);
  f_CombBg->SetParameters(pars[6], pars[7], Mmin); 
  f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  f_DY->SetParameters(pars[8], pars[9], Mmin);
  f_DY->SetLineColor(kBlue); f_DY->Draw("same");
}


void Draw_2Crystal_2Expo(Double_t Norm_JPsi, Double_t M_JPsi, Double_t W_JPsi,
			 Double_t alpha_JPsi, Double_t n_JPsi,
			 Double_t Norm_psi, Double_t M_psi, Double_t W_psi,
			 Double_t alpha_psi, Double_t n_psi,
			 Double_t A_e1, Double_t c_e1,
			 Double_t A_e2, Double_t c_e2,
			 Double_t Mmin, Double_t Mmax){

  TF1 *f_JPsi =
    new TF1("f_JPsi", f_CrystalBall, 0, Mmax, 5);
  f_JPsi->SetParameters(Norm_JPsi, M_JPsi, W_JPsi, alpha_JPsi, n_JPsi);
  f_JPsi->SetLineColor(7); f_JPsi->Draw("same");

  TF1 *f_psi = new TF1("f_psi", f_CrystalBall, 0, Mmax, 5);
  f_psi->SetParameters(Norm_psi, M_psi, W_psi, alpha_psi, n_psi);
  f_psi->SetLineColor(kGreen); f_psi->Draw("same");
  
  TF1 *f_CombBg = new TF1("f_CombBg", "[0]*TMath::Exp([1]*(x-[2]))",0, Mmax);
  f_CombBg->SetParameters(A_e1, c_e1, Mmin); 
  f_CombBg->SetLineColor(28); f_CombBg->Draw("same");

  TF1 *f_DY = new TF1("f_DY", "[0]*TMath::Exp([1]*(x-[2]) )", 0, Mmax);
  f_DY->SetParameters(A_e2, c_e2, Mmin);
  f_DY->SetLineColor(kBlue); f_DY->Draw("same");
}


//Error Calculations
Double_t CorrelatedRatioError(Double_t A, Double_t B,
			      Double_t eA, Double_t eB){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t r = A/B;
  if (r<0) r *= -1.0;
 
  Double_t e = eA*eA/(A*A) + eB*eB/(B*B) - 2*eA*eB/(A*B);
  e = TMath::Sqrt( e );
  e *= r;

  return e;
}


Double_t CorrelatedDiffError(Double_t eA, Double_t eB){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }

  Double_t e = eA*eA + eB*eB - 2*eA*eB;
  e = TMath::Sqrt(e);

  return e;
}


Double_t DiffError(Double_t eA, Double_t eB){
  Double_t e2A =eA*eA, e2B =eB*eB;

  return TMath::Sqrt(e2A + e2B);
}


Double_t RatioError(Double_t A, Double_t B,
		    Double_t eA, Double_t eB){
  if (eA<0 || eB<0){
    cout << "negative errors" << endl;
    exit(EXIT_FAILURE);
  }
  
  Double_t r = A/B;
  if (r<0) r *= -1.0;
  if (A == 0.0) return eB;

  Double_t e = eA*eA/(A*A) + eB*eB/(B*B);
  e = TMath::Sqrt( e );
  
  e *= r;
  
  return e;
}


//Additional Calculations
Double_t WeightedAvg(vector<Double_t> &A, vector<Double_t> &eA){

  Double_t avg=0.0, e=0.0;
  for (vector<Double_t>::iterator iA=A.begin(), iE=eA.begin(); iA!=A.end();
       iA++, iE++) {
    avg += *iA/( (*iE) * (*iE) );
    e += 1.0/( (*iE) * (*iE) );
  }
  
  avg /= e;

  return avg;
}


Double_t WeightedAvg(Double_t A, Double_t B,
		     Double_t eA, Double_t eB){
  Double_t e = 1/(eA*eA) + 1/(eB*eB);
  
  Double_t avg = A/(eA*eA) + B/(eB*eB);
  avg /= e;

  return avg;
}


Double_t WeightedErr(vector<Double_t> &eA){

  Double_t e=0.0;
  for (vector<Double_t>::iterator it=eA.begin(); it!=eA.end(); it++) {
    e += 1.0/( (*it) * (*it) );
  }

  e = 1/e;
  e = TMath::Sqrt(e);
  
  return e;
}


Double_t WeightedErr(Double_t eA, Double_t eB){
  Double_t e = 1/(eA*eA) + 1/(eB*eB);
  e = 1/e;
  e = TMath::Sqrt(e);
  
  return e;
}


Double_t GaussInt(Double_t A, Double_t sigma){
  
  return TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
}


Double_t GaussIntError(Double_t A, Double_t sigma,
		       Double_t eA, Double_t eSigma){
  Double_t e = eA*eA/( A*A ) + eSigma*eSigma/( sigma*sigma );
  e = TMath::Sqrt( e );
  e *= TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
  
  return e;
}


Double_t CorrelatedGaussIntError(Double_t A, Double_t sigma,
				 Double_t eA, Double_t eSigma, Double_t cov){
  Double_t e = eA*eA/( A*A ) + eSigma*eSigma/( sigma*sigma ) + 2*cov/(A*sigma);
  e = TMath::Sqrt( e );
  e *= TMath::Sqrt( 2*TMath::Pi() )*A*sigma;
  
  return e;
}


Double_t ExpoInt(Double_t A, Double_t b){
  if (b > 0){
    cout << "Cannot do exponential integral" << endl;
    exit(EXIT_FAILURE);
  }
  
  return -1.0*A/b;
}


//Functions
Double_t ChiSquareDistr(Double_t *x,Double_t *par)
{
  // Chisquare density distribution for nrFree degrees of freedom
  Double_t nrFree = par[0];
  Double_t A = par[1];
  Double_t chi2 = x[0];

  if (chi2 > 0) {
    Double_t lambda = nrFree/2.;
    Double_t norm = TMath::Gamma(lambda)*TMath::Power(2.,lambda)/A;
    //Double_t norm = TMath::Gamma(lambda)*TMath::Power(2.,lambda);
    return TMath::Power(chi2,lambda-1)*TMath::Exp(-0.5*chi2)/norm;
  }
  else return 0.0;
}


Double_t RedChiSquareDistr(Double_t *x,Double_t *par)
{
  // Reduced Chisquare density distribution for nrFree degrees of freedom
  Double_t nrFree = par[0];
  x[0] = x[0]*nrFree;
  
  return ChiSquareDistr(x, par);
}
