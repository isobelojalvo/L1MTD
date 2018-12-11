#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TF1.h"
#include <math.h>
#include <TEfficiency.h>
#include <TMath.h>
#include <iostream>
#include <string>

#include <iostream>
#include <cmath>
#include "TLegend.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TF1.h"
#include "THStack.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TTree.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TGaxis.h"


void plotRatesIsoVersions(TString type, TString variable, TString header, TString folder, TString inFile, TString outFile, TString name = "");

void tdrGrid(bool gridOn);
void setTDRStyle() ;
void setMaxErrorTo1(TGraphAsymmErrors * hist);
void applyPadStyle(TPad* pad1);
void applyLegStyle(TLegend *leg);
TH1F* plotRateReturnTH1F(TString variable, TString fileName, TString outFileName, TString l1Cut, TString xaxis, TString text, TString treePath, int bins, float low, float high, Color_t color, bool variableBin=false, int linestyle=1, TString yaxis = "Rate");
int getNEntries(TString fileName, TString histName);
