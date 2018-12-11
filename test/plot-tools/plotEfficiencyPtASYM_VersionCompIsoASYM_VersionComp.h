
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
#include <TFormula.h>
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
#include "TGraphAsymmErrors.h"


void applyPadStyle(TPad* pad1);
void applyLegStyle(TLegend *leg);
TGraphAsymmErrors* plotEfficiencyReturnASYM(TString variable, TString fileName, TString outFileName, TString l1Cut, TString recoCut, TString xaxis, TString text, TString treePath, int bins, float low, float high, Color_t color, bool variableBin=false, int markerstyle = 20, int inputfitparam = 88);
void plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp(TString inputVariable, TString lepton, TString DM, TString folder, TString inFile, TString outFile, TString isoTime, TString isoNoTime, TString label);

void setTDRStyle();
void tdrGrid(bool gridOn);
