#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/helperFunctions.h"
#include "/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros/include/finalSetup.h"


void DrawLegend(TGraphErrors *g){
  TString name = g->GetName();
  Double_t sigma;
  Double_t avg = WeightedAvgAndError(g, &sigma);

  TLegend *leg = new TLegend(0.25,0.9,0.7,0.99);

  //leg->AddEntry(name, Form("#bar{A}_{N} = %0.2f #pm %0.2f", avg, sigma), "p");

  //fa2TargJuraSaleve
  leg->AddEntry(name, Form("#bar{A}_{N,false} = %0.2f #pm %0.2f", avg, sigma), "p");
  
  leg->SetBorderSize(0);
  leg->SetTextFont(133); leg->SetTextSize(20);
  leg->Draw("same");
}


void OffSetLocal(TGraphErrors *g, Double_t *xvals, Int_t per){
  TAxis *xaxis = g->GetXaxis();
  Double_t range = xaxis->GetXmax() - xaxis->GetXmin();
  Double_t *change_xval = g->GetX();
  for (Int_t i=0; i<g->GetN(); i++) {
    change_xval[i] = xvals[i] + range*per/50.0;
  }
}


void FinalSetupLocal(TGraphErrors *g, TString xName, Double_t yMax,
		     Bool_t same=false){
  g->GetYaxis()->SetRangeUser(-yMax, yMax);
  FinalSetup(g);

  g->SetMarkerColor(kRed);
  g->SetMarkerStyle(20);
  
  g->SetMarkerSize(1.8);
  g->SetFillStyle(0);
  g->SetLineWidth(2);
  g->SetTitle("");

  SetTitleName(g, xName, "x");
}


void FinalLocalIntegrated(TGraphErrors *g){
  //SetTitleName(g, "A_{N}", "y");//allPhysBinned4Targ
  SetTitleName(g, "A_{N,Fa}", "y");//fa2TargJuraSaleve
  
  g->GetXaxis()->SetLabelSize(0.0);
  g->GetXaxis()->SetTickSize(0.0);
}


void KinematicsCanvas(TString inputFile=""){
  //Basic Setup
  TString localPath
    ="/Users/robertheitz/Documents/Research/DrellYan/Analysis/TGeant/\
Local_LeftRight_Analysis/Macros";
  const Int_t nPhysBinned =6;
  Double_t yMax =1.1;
  TString xNames[nPhysBinned] =
    {"", "x_{N}", "x_{#pi}", "x_{F}", "q_{T} (GeV/c)","M_{#mu#mu} (GeV/c^{2})"};
  TString periods[9] =
    {"W07", "W08", "W09", "W10", "W11", "W12", "W13", "W14", "W15"};

  Bool_t toWrite =true;

  //Aesthetics setup
  TCanvas* c1 = new TCanvas(); c1->SetLeftMargin(0.2);
  c1->Divide(nPhysBinned, 1, 0, 0);
  
  //TLegend *legend = new TLegend(0.27, 0.8, 0.99, 0.99);
  //legend->SetNColumns(3); //allPhysBinned4Targ

  fstream fIn(inputFile);
  TString fname, gname;
  Int_t iPad;
  Int_t lastPad =0, first =0;
  Int_t integrated =1, per=0;
  Double_t *xvals;
  while (fIn >> fname >> gname >> iPad){
    TFile *f1 = OpenFile(localPath+fname);
    TGraphErrors *g =(TGraphErrors*)f1->Get(gname);
    if (!g) {
      cout << "TGraph does not exist" << endl;
      exit(EXIT_FAILURE);
    }

    //Conditions for draw first
    if (lastPad ==0){
      lastPad =1;
      first =1;
    }
    else if (iPad==lastPad){ first =0; }
    else{
      integrated =0;
      first =1;
      lastPad =iPad;
      per =0;
    }
    c1->cd(iPad); gPad->SetFrameLineWidth(2);
    FinalSetupLocal(g, xNames[iPad-1], yMax);
    
    g->GetYaxis()->SetRangeUser(-0.3, 0.05);//fa2TargJuraSaleve

    //X-axis limits (must be before DrawLine)
    if (integrated){
      //allPhysBinned4Targ
      //legend->AddEntry(g, Form("%s", periods[per].Data()), "p");   
      
      g->GetXaxis()->SetLimits(0.16, 0.19);
    }
    else if (xNames[iPad-1] == "x_{N}")
      g->GetXaxis()->SetLimits(0.085, 0.315);
    else if (xNames[iPad-1] == "x_{#pi}")
      g->GetXaxis()->SetLimits(0.25, 0.9);
    else if (xNames[iPad-1] == "x_{F}")
      g->GetXaxis()->SetLimits(0.05, 0.81);
    else if (xNames[iPad-1] == "q_{T} (GeV/c)")
      g->GetXaxis()->SetLimits(0.5, 2.4);
    else if (xNames[iPad-1] == "M_{#mu#mu} (GeV/c^{2})")
      g->GetXaxis()->SetLimits(4.3, 7.2);
    
    if (first){
      g->Draw("AP");
      DrawLine(g, 0.0);
      first =0;
      xvals = g->GetX();

      if (integrated){
	FinalLocalIntegrated(g);
      }
    }
    else {
      g->Draw("Psame");
    }

    //Shift x-axis by periods
    if (integrated)
      OffSetLocal(g, xvals, 2*per);      
    else
      OffSetLocal(g, xvals, per);

    DrawLegend(g);//fa2TargJuraSaleve

    per++;
  }

  //Legend //allPhysBinned4Targ
  //c1->cd(1);
  //legend->SetBorderSize(0); legend->SetTextFont(132); legend->SetTextSize(0.08);
  //legend->Draw("same");

  //Write output
  if (toWrite){
    TFile *fOut = new TFile("output.root", "RECREATE");
    c1->Write();
  }

  cout << "\nIf there are glitches in pdf output" << endl;
  cout << "Save file as .jpg" << endl;
  cout << " " << endl;
}
