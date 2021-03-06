#pragma once
#ifndef SETUP_H
#define SETUP_H

#include "common.h"

void SetupTGraph(TGraphErrors* gr, TString title, TString xTitle,
			  Double_t ymax=1.0);


void SetupTLine(TLine*);

void SetupHist(TH1D* h);

#endif
