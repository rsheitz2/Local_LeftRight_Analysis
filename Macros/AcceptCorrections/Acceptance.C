#include "helperFunctions.h"


void ReadFile(TString bin_path, Double_t *b, Double_t *x, TString type,
	      Int_t nBins){
  ifstream f(bin_path);
  
  if (type=="xN"||type=="xPi"||type=="pT"||type=="mass"||type=="rad")
    b[0] = 0.0;
  else if (type=="xF") b[0] = -1.0;
  else if (type=="vxZ_upstream") b[0] = -294.5;
  else if (type=="vxZ_downstream") b[0] = -219.5;
  else {
    std::cout << "Invalid type: " << type << " in leftright::leftright"
	      << std::endl;
    exit(EXIT_FAILURE);
  }
  
  TString boundaries = type; boundaries += " bin boundaries";
  TString averages = type; averages += " bin averages";
  
  string line;
  bool first = false, found = false;
  Int_t iter = 1;
  
  while (!f.eof()) {
    getline(f, line);
    TString tline (line);

    if (tline == boundaries) {
      found = true; first = true;
      continue;
    }
    else if (tline == averages) {
      first = false; iter = 0;
      continue;
    }
    else if (!found) continue;

    if (first) {
      if (iter == nBins) cout << "\n Bin size problem\n" << endl;
      b[iter] = atof(line.c_str() );
      iter++;
    }
    else {
      if (iter >= nBins) break;
      x[iter] = atof(line.c_str() );
      iter++;
    }
  }//while loop
  f.close();
  
  if (type=="xN"||type=="xPi"||type=="xF") b[nBins] = 1.0;
  else if (type=="pT") b[nBins] = 5.0;
  else if (type=="mass") b[nBins] = 12.0;
  else if (type=="rad") b[nBins] = 2.0;
  else if (type=="vxZ_upstream") b[nBins] = -239.3;
  else if (type=="vxZ_downstream") b[nBins] = -164.3;
}


Bool_t BinDataCounts(unsigned long long *counts, Double_t binVal,
		     std::vector<Double_t> &binValBounds){

  Int_t iter=-1;
  for (std::vector<Double_t>::iterator it=binValBounds.begin();
       it!=binValBounds.end(); it++, iter++){
    if(binVal <= *it ) {
      if(iter==-1){
	std::cout << "bin value too low!!!!" << std::endl;
	std::cout << *it << " " << binVal << std::endl;
	std::cout << " " << std::endl;
	return false;
      }
      
      counts[iter]++;
      return true;
    }
  }

  std::cout << "bin value too high!!!!" << std::endl;
  std::cout << binValBounds.back() << " " << binVal << std::endl;
  std::cout << " " << std::endl;
  return false;
}


Bool_t BinAcc(unsigned long long *Top, unsigned long long *Bottom,
		       Double_t* Acc, Double_t *e_Acc, Int_t nBins){
  
  for (Int_t i=0; i<nBins; i++) {
    Acc[i] = 1.0*Top[i]/( 1.0*Bottom[i] );
    e_Acc[i] = RatioError(Top[i], Bottom[i]);
  }

  return true;
}//BinAcc


void Acceptance(){
  ////////////////
  //Changes Here//
  TString period = "W07"; //WAll || W07 for test
  TString subper = "sp1";
  if (period=="W09"||period=="W10") period += Form("_%s", subper.Data() );
  TString mcType ="HMDY", massRange ="HM";
  const Int_t nBins=3;
  //const Int_t nBins=1;
  Bool_t toWrite =false;
  
  TString path="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant";
  
  TString MC_path = "Presents/DATA/MC_Data/YuShiangMC";
  TFile* f_phast =
    TFile::Open(Form("%s/%s/%s/Yu_%s_%s_3bins.root", path.Data(),MC_path.Data(),
		     mcType.Data(), mcType.Data(), period.Data() ) );
  
  TString Gen_MC_path = "Presents/DATA/Gen_MC_Data/Yu_BW/";
  TFile* f_gen =
    TFile::Open(Form("%s/%s/%s/%s_%s_run0.root", path.Data(),Gen_MC_path.Data(),
		     mcType.Data(), mcType.Data(), period.Data() ) );

  TString bin_path = path+"/Presents/DATA/RealData/";
  (massRange=="HM") ? bin_path += "HMDY/" : bin_path += "JPsi/";
  bin_path += "BinValues/WAll_"+massRange+"_"+Form("%i", nBins)+"bins.txt";
  ifstream f_bins_phast(bin_path);
  cout << "\nNote: Binning used from  REAL DATA  ";
  cout << " from   WAll   for any input periods\n" << endl;
  
  if ( !(f_phast) || !(f_gen) || !(f_bins_phast) ){//Basic file checks
    cout << "Error:" << endl;
    cout << "One of the files did not open" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  else f_bins_phast.close();

  
  TTree *T_phast = (TTree*)f_phast->Get("pT_Weighted");
  Double_t ph_x_beam, ph_x_target, ph_x_feynman, ph_pT, ph_Mmumu;
  Double_t ph_vx_z, ph_PhiS_simple;
  Int_t ph_trigMask;
  T_phast->SetBranchAddress("x_beam", &ph_x_beam);
  T_phast->SetBranchAddress("x_target", &ph_x_target);
  T_phast->SetBranchAddress("x_feynman", &ph_x_feynman);
  T_phast->SetBranchAddress("q_transverse", &ph_pT);
  T_phast->SetBranchAddress("Mmumu", &ph_Mmumu);
  T_phast->SetBranchAddress("vx_z", &ph_vx_z);
  T_phast->SetBranchAddress("PhiS_simple", &ph_PhiS_simple);
  T_phast->SetBranchAddress("trigMask", &ph_trigMask);

    
  TTree *T_gen = (TTree*)f_gen->Get("Event");
  Double_t gen_x_beam, gen_x_target, gen_x_feynman, gen_pT, gen_Mmumu;
  Double_t gen_vx_z, gen_vx_x, gen_vx_y;
  Double_t gen_PhiS_simple;
  T_gen->SetBranchAddress("x_beam", &gen_x_beam);
  T_gen->SetBranchAddress("x_target", &gen_x_target);
  T_gen->SetBranchAddress("x_feynman", &gen_x_feynman);
  T_gen->SetBranchAddress("q_transverse", &gen_pT);
  T_gen->SetBranchAddress("Mmumu", &gen_Mmumu);
  T_gen->SetBranchAddress("vx_z", &gen_vx_z);
  T_gen->SetBranchAddress("vx_x", &gen_vx_x);
  T_gen->SetBranchAddress("vx_y", &gen_vx_y);
  T_gen->SetBranchAddress("PhiS_simple", &gen_PhiS_simple);


  Double_t b_xN[nBins+1], x_xN[nBins];
  Double_t b_xPi[nBins+1], x_xPi[nBins];
  Double_t b_xF[nBins+1], x_xF[nBins];
  Double_t b_pT[nBins+1], x_pT[nBins];
  Double_t b_M[nBins+1], x_M[nBins];
  ReadFile(bin_path, b_xN, x_xN, "xN", nBins);
  ReadFile(bin_path, b_xPi, x_xPi, "xPi", nBins);
  ReadFile(bin_path, b_xF, x_xF, "xF", nBins);
  ReadFile(bin_path, b_pT, x_pT, "pT", nBins);
  ReadFile(bin_path, b_M, x_M, "mass", nBins);
  vector<Double_t> Bounds;
  
  ////////////////
  //Changes Here//
  Double_t *ph_phys=&ph_x_feynman, *gen_phys=&gen_x_feynman, *physB=b_xF, *xval=x_xF;
  Double_t xMin =-0.2, xMax =0.8; TString physType="xF";//x_feynman
  Int_t nHistBins = 100;
  

  for (Int_t i=0; i<nBins+1; i++, physB++) Bounds.push_back(*physB);
  const Int_t nHist = 18;
  TString accType[nHist] = {"", "UpS", "DownS", "Left", "Right",
			    "UpS_Left", "UpS_Right",
			    "DownS_Left", "DownS_Right",
			    "LL", "LO", "LL_LO",
			    "LL_Left", "LL_Right",
			    "LO_Left", "LO_Right",
			    "LL_LO_Left", "LL_LO_Right"};
  TH1D* h_phast[nHist]; TH1D* h_gen[nHist];
  for (Int_t i=0; i<nHist; i++) {
    h_phast[i] = new TH1D(Form("h_phast%s",accType[i].Data() ),
			  Form("h_phast%s",accType[i].Data() ),
			  nHistBins, xMin, xMax);
    h_gen[i] = new TH1D(Form("h_gen%s",accType[i].Data() ),
			Form("h_gen%s",accType[i].Data() ),
			nHistBins, xMin, xMax);
  }

  unsigned long long ph_counts[nHist][nBins] = {0};
  Int_t ph_entries = T_phast->GetEntries();
  for (Int_t ev=0; ev<ph_entries; ev++) {
    T_phast->GetEntry(ev);
    h_phast[0]->Fill(*ph_phys);//""
    BinDataCounts(ph_counts[0], *ph_phys, Bounds);

    Bool_t upStream=false, downStream=false;
    if (ph_vx_z>-294.5 && ph_vx_z<-239.3)                         upStream=true;
    else if (ph_vx_z>-219.5 && ph_vx_z<-164.3)                  downStream=true;

    Bool_t left=false, right=false;
    if (ph_PhiS_simple>0 && ph_PhiS_simple<=TMath::Pi() )             left=true;
    else if (ph_PhiS_simple>-1.0*TMath::Pi() &&ph_PhiS_simple<=0.0 ) right=true;

    Bool_t LL=false, LO=false;
    if ( ((ph_trigMask >> 8) & 1) )                                   LL = true;
    if ( ((ph_trigMask >> 2) & 1) )                                   LO = true;
    

    if (upStream){ h_phast[1]->Fill(*ph_phys);//UpS
      BinDataCounts(ph_counts[1], *ph_phys, Bounds);
      
      if (left) {
	h_phast[5]->Fill(*ph_phys);//UpS_Left
	BinDataCounts(ph_counts[5], *ph_phys, Bounds);
      }
      else if (right) {
	h_phast[6]->Fill(*ph_phys);//UpS_Right
	BinDataCounts(ph_counts[6], *ph_phys, Bounds);
      }
    }
    else if (downStream){ h_phast[2]->Fill(*ph_phys);//DownS
      BinDataCounts(ph_counts[2], *ph_phys, Bounds);
      
      if (left) {
	h_phast[7]->Fill(*ph_phys);//DownS_Left
	BinDataCounts(ph_counts[7], *ph_phys, Bounds);
      }
      else if (right) {
	h_phast[8]->Fill(*ph_phys);//DownS_Right
	BinDataCounts(ph_counts[8], *ph_phys, Bounds);
      }
    }

    if (left) {
      h_phast[3]->Fill(*ph_phys);//Left
      BinDataCounts(ph_counts[3], *ph_phys, Bounds);
    }
    else if (right) {
      h_phast[4]->Fill(*ph_phys);//Right
      BinDataCounts(ph_counts[4], *ph_phys, Bounds);
    }

    if (LL && !LO) {//LAST_LAST only
      h_phast[9]->Fill(*ph_phys);
      BinDataCounts(ph_counts[9], *ph_phys, Bounds);
      
      if (left){
	h_phast[12]->Fill(*ph_phys);//Left
	BinDataCounts(ph_counts[12], *ph_phys, Bounds);
      }
      else if (right){
	h_phast[13]->Fill(*ph_phys);//Right
	BinDataCounts(ph_counts[13], *ph_phys, Bounds);
      }
    }
    else if (LO && !LL){//LAST_Outer only
      h_phast[10]->Fill(*ph_phys);
      BinDataCounts(ph_counts[10], *ph_phys, Bounds);

      if (left){
	h_phast[14]->Fill(*ph_phys);//Left
	BinDataCounts(ph_counts[14], *ph_phys, Bounds);
      }
      else if (right){
	h_phast[15]->Fill(*ph_phys);//Right
	BinDataCounts(ph_counts[15], *ph_phys, Bounds);
      }
    }

    if (LL && LO){//LAST_LAST and LAST_Outer
      h_phast[11]->Fill(*ph_phys);
      BinDataCounts(ph_counts[11], *ph_phys, Bounds);

      if (left){
	h_phast[16]->Fill(*ph_phys);//Left
	BinDataCounts(ph_counts[16], *ph_phys, Bounds);
      }
      else if (right){
	h_phast[17]->Fill(*ph_phys);//Right
	BinDataCounts(ph_counts[17], *ph_phys, Bounds);
      }
    }
  }
  

  unsigned long long gen_counts[nHist][nBins] = {0};
  Int_t gen_entries = T_gen->GetEntries();
  for (Int_t ev=0; ev<gen_entries; ev++) {
    T_gen->GetEntry(ev);

    if ((gen_vx_z<-294.5||gen_vx_z>-239.3) &&(gen_vx_z<-219.5||gen_vx_z>-164.3))
      continue;
    if (gen_vx_x*gen_vx_x + gen_vx_y*gen_vx_y > 1.9*1.9) continue;
    if (gen_pT<0.4 || gen_pT>5.0) continue;
    if (gen_Mmumu<4.3 || gen_Mmumu>8.5) continue;
    
    h_gen[0]->Fill(*gen_phys);//""
    BinDataCounts(gen_counts[0], *gen_phys, Bounds);

    Bool_t upStream=false, downStream=false;
    if (gen_vx_z>-294.5 && gen_vx_z<-239.3)                       upStream=true;
    else if (gen_vx_z>-219.5 && gen_vx_z<-164.3)                downStream=true;

    Bool_t left=false, right=false;
    if (gen_PhiS_simple>0 && gen_PhiS_simple<=TMath::Pi() )           left=true;
    else if (gen_PhiS_simple>-1.0*TMath::Pi()&&gen_PhiS_simple<=0.0) right=true;


    if (upStream){ h_gen[1]->Fill(*gen_phys);//UpS
      BinDataCounts(gen_counts[1], *gen_phys, Bounds);
      
      if (left) {
	h_gen[5]->Fill(*gen_phys);//UpS_Left
	BinDataCounts(gen_counts[5], *gen_phys, Bounds);
      }
      else if (right) {
	h_gen[6]->Fill(*gen_phys);//UpS_Right
	BinDataCounts(gen_counts[6], *gen_phys, Bounds);
      }
    }
    else if (downStream){ h_gen[2]->Fill(*gen_phys);//DownS
      BinDataCounts(gen_counts[2], *gen_phys, Bounds);
      
      if (left) {
	h_gen[7]->Fill(*gen_phys);//DownS_Left
	BinDataCounts(gen_counts[7], *gen_phys, Bounds);
      }
      else if (right) {
	h_gen[8]->Fill(*gen_phys);//DownS_Right
	BinDataCounts(gen_counts[8], *gen_phys, Bounds);
      }
    }

    if (left) {
      h_gen[3]->Fill(*gen_phys);//Left
      BinDataCounts(gen_counts[3], *gen_phys, Bounds);

      h_gen[12]->Fill(*gen_phys);//LAST_LAST only Left 
      BinDataCounts(gen_counts[12], *gen_phys, Bounds);
      h_gen[14]->Fill(*gen_phys);//LAST_Outer only Left
      BinDataCounts(gen_counts[14], *gen_phys, Bounds);
      h_gen[16]->Fill(*gen_phys);//LAST_LAST and LAST_Outer Left
      BinDataCounts(gen_counts[16], *gen_phys, Bounds);
    }
    else if (right) {
      h_gen[4]->Fill(*gen_phys);//Right
      BinDataCounts(gen_counts[4], *gen_phys, Bounds);

      h_gen[13]->Fill(*gen_phys);//LAST_LAST only Right 
      BinDataCounts(gen_counts[13], *gen_phys, Bounds);
      h_gen[15]->Fill(*gen_phys);//LAST_Outer only Right
      BinDataCounts(gen_counts[15], *gen_phys, Bounds);
      h_gen[17]->Fill(*gen_phys);//LAST_LAST and LAST_Outer Right
      BinDataCounts(gen_counts[17], *gen_phys, Bounds);
    }

    h_gen[9]->Fill(*gen_phys);//LAST_LAST only
    BinDataCounts(gen_counts[9], *gen_phys, Bounds);
    h_gen[10]->Fill(*gen_phys);//LAST_Outer only
    BinDataCounts(gen_counts[10], *gen_phys, Bounds);
    h_gen[11]->Fill(*gen_phys);//LAST_LAST and LAST_Outer
    BinDataCounts(gen_counts[11], *gen_phys, Bounds);
  }

  
  TH1D *h_acc[nHist];
  for (Int_t i=0; i<nHist; i++) {
    h_phast[i]->Sumw2();
    h_acc[i] = (TH1D*)h_phast[i]->Clone();

    h_acc[i]->Divide(h_gen[i]);
    SetUpTH1(h_acc[i]); SetUpTH1(h_gen[i]); SetUpTH1(h_phast[i]);
    h_acc[i]->GetYaxis()->SetRangeUser(0, 1);
  }
  
  
  gStyle->SetOptStat(1111); gStyle->SetStatFont(22);
  TCanvas* c1 = new TCanvas();
  c1->Divide(2);

  c1->cd(1);
  h_gen[0]->SetLineColor(kBlack); h_gen[0]->Draw();
  h_phast[0]->Draw("sames");
  h_phast[9]->SetLineColor(kRed); h_phast[9]->Draw("same");
  h_phast[10]->SetLineColor(kGreen); h_phast[10]->Draw("same");
  h_phast[11]->SetLineColor(6); h_phast[11]->Draw("same");
  gPad->SetLogy();
  
  c1->cd(2); h_acc[0]->Draw();
  h_acc[9]->SetLineColor(kRed); h_acc[9]->Draw("same");
  h_acc[10]->SetLineColor(kGreen); h_acc[10]->Draw("same");
  h_acc[11]->SetLineColor(6); h_acc[11]->Draw("same");


  TCanvas* cR = new TCanvas("updown_lr"); cR->Divide(2);
  TRatioPlot *r_updown = new TRatioPlot(h_acc[1], h_acc[2] );
  cR->cd(1);  r_updown->Draw();

  TRatioPlot *r_lr = new TRatioPlot(h_acc[3], h_acc[4] );
  cR->cd(2);  r_lr->Draw();


  TCanvas* cS = new TCanvas(); cS->Divide(2);
  TRatioPlot *r_ups_lr = new TRatioPlot(h_acc[5], h_acc[6] );
  cS->cd(1);  r_ups_lr->Draw();

  TRatioPlot *r_ds_lr = new TRatioPlot(h_acc[7], h_acc[8] );
  cS->cd(2);  r_ds_lr->Draw();

  
  TCanvas* c4 = new TCanvas(); c4->Divide(2,2);
  TRatioPlot *r_ud_l = new TRatioPlot(h_acc[5], h_acc[7] );
  c4->cd(1); r_ud_l->Draw();

  TRatioPlot *r_ud_r = new TRatioPlot(h_acc[6], h_acc[8] );
  c4->cd(2); r_ud_r->Draw();


  TRatioPlot *r_ud_lr = new TRatioPlot(h_acc[5], h_acc[8] );
  c4->cd(3); r_ud_lr->Draw();

  TRatioPlot *r_ud_rl = new TRatioPlot(h_acc[6], h_acc[7] );
  c4->cd(4);  r_ud_rl->Draw();


  Double_t b_acc[nHist][nBins], eb_acc[nHist][nBins];
  Double_t ex[nBins]={0.};
  Double_t dx = 0.02;
  Double_t offset[nHist] = {0., 0., 0., 0., dx,
			    0., dx,
			    0., dx,
			    0., 0., 0.,
			    0., dx,
			    0., dx,
			    0., dx};
  TGraphErrors *g_acc[nHist], *g_acc_ratio[6];
  Double_t acc_ratio[6][nBins], e_acc_ratio[6][nBins];
  for (Int_t i=0, r=0; i<nHist; i++) {
    BinAcc(ph_counts[i], gen_counts[i], b_acc[i], eb_acc[i], nBins);
    g_acc[i] = new TGraphErrors(nBins, xval, b_acc[i], ex, eb_acc[i]);
    SetUpTGraph(g_acc[i], offset[i], nBins);

    if ( accType[i].Contains("Left") ){
      Double_t *yvals = g_acc[i]->GetY();
      Double_t *e_yvals = g_acc[i]->GetEY();
      for (Int_t j=0; j<nBins; j++) {
	acc_ratio[r][j] = yvals[j];
	e_acc_ratio[r][j] = e_yvals[j];
      }
    }
    else if ( accType[i].Contains("Right") ){
      Double_t *acc_right = g_acc[i]->GetY();
      Double_t *e_acc_right = g_acc[i]->GetEY();
      
      for (Int_t j=0; j<nBins; j++) {
	Double_t acc_left = acc_ratio[r][j];
	Double_t e_acc_left = e_acc_ratio[r][j];
	
	e_acc_ratio[r][j] = RatioError(acc_left, acc_right[j],
				       e_acc_left, e_acc_right[j] );
	acc_ratio[r][j] = acc_ratio[r][j]/acc_right[j];
      }
      
      g_acc_ratio[r] = new TGraphErrors(nBins, xval, acc_ratio[r], ex,
					e_acc_ratio[r]);
      SetUpTGraph(g_acc_ratio[r], r*dx/5.0, nBins);
      g_acc_ratio[r]->GetYaxis()->SetRangeUser(0.85, 1.15);
      r++;
    }
  }

  TCanvas* cG = new TCanvas();
  g_acc[0]->Draw("AP");

  TCanvas* cAlr = new TCanvas(); cAlr->Divide(1,2);
  cAlr->cd(1);
  g_acc[3]->Draw("AP");
  g_acc[3]->SetMarkerColor(kBlue); g_acc[3]->GetYaxis()->SetRangeUser(0.0, 0.4);
  g_acc[4]->Draw("Psame"); g_acc[4]->SetMarkerColor(kBlue);
  g_acc[12]->Draw("Psame"); g_acc[12]->SetMarkerColor(kRed);
  g_acc[13]->Draw("Psame"); g_acc[13]->SetMarkerColor(kRed);
  g_acc[14]->Draw("Psame"); g_acc[14]->SetMarkerColor(kGreen);
  g_acc[15]->Draw("Psame"); g_acc[15]->SetMarkerColor(kGreen);
  g_acc[16]->Draw("Psame"); g_acc[16]->SetMarkerColor(6);
  g_acc[17]->Draw("Psame"); g_acc[17]->SetMarkerColor(6);
  cAlr->cd(2);
  g_acc_ratio[0]->Draw("AP"); g_acc_ratio[0]->SetMarkerColor(kBlue);
  g_acc_ratio[3]->Draw("Psame"); g_acc_ratio[3]->SetMarkerColor(kRed);
  g_acc_ratio[4]->Draw("Psame"); g_acc_ratio[4]->SetMarkerColor(kGreen);
  g_acc_ratio[5]->Draw("Psame"); g_acc_ratio[5]->SetMarkerColor(6);

  Double_t min_x = g_acc_ratio[0]->GetXaxis()->GetXmin();	
  Double_t max_x = g_acc_ratio[0]->GetXaxis()->GetXmax();	
  TLine* li_acc_ratio = new TLine(min_x, 1.0, max_x, 1.0);
  li_acc_ratio->SetLineColor(1);
  li_acc_ratio->SetLineStyle(8);
  li_acc_ratio->Draw("same");


  TCanvas* cTrLR = new TCanvas(); cTrLR->Divide(1,2);
  cTrLR->cd(1);
  g_acc[5]->Draw("AP"); g_acc[5]->SetMarkerColor(kRed);
  g_acc[6]->Draw("Psame"); g_acc[6]->SetMarkerColor(kRed);
  g_acc[7]->Draw("Psame"); g_acc[7]->SetMarkerColor(kBlue);
  g_acc[8]->Draw("Psame"); g_acc[8]->SetMarkerColor(kBlue);

  cTrLR->cd(2);
  g_acc_ratio[1]->Draw("AP"); g_acc_ratio[1]->SetMarkerColor(kRed);
  g_acc_ratio[2]->Draw("Psame"); g_acc_ratio[2]->SetMarkerColor(kBlue);
  min_x = g_acc_ratio[1]->GetXaxis()->GetXmin();	
  max_x = g_acc_ratio[1]->GetXaxis()->GetXmax();	
  TLine* li_a_r_updown = new TLine(min_x, 1.0, max_x, 1.0);
  li_a_r_updown->SetLineColor(1);
  li_a_r_updown->SetLineStyle(8);
  li_a_r_updown->Draw("same");
  
  
  TString fname = "Accept/Acceptance_"; fname += physType; fname += "_";
  fname+= period; fname += "_"; fname += Form("%i", nBins); fname+= ".root";
  if (toWrite){
    TFile *fOutput = new TFile(fname, "RECREATE");
    for (Int_t i=0; i<nHist; i++) {
      g_acc[i]->Write(Form("acc_%s", accType[i].Data() ) );
    }
    fOutput->Close();
    cout <<" " << endl;
    cout << fname << "    was written" << endl;
  }
  else{
    cout <<" " << endl;
    cout << fname << "    was NOT written" << endl;
  }
  
  cout << " " << endl;
  cout << "Period is:    "  << period << endl;
  cout << "Acceptance from:  " << physType << endl;
  cout  << "xMin: " << xMin << "     xMax: " << xMax << endl;
  cout << "Number of output bins is:  " << nBins << endl;
  cout << " " << endl;
}
