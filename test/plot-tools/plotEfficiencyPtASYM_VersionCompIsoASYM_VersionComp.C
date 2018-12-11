
#include "plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp.h"

void plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp(TString inputVariable, TString lepton, TString DM, TString folder, TString inFile, TString outFile, TString isoTime, TString isoNoTime, TString label){
  
  
  //gROOT->LoadMacro("CMS_lumi.C");
  setTDRStyle();

  TFile *tauFile = new TFile("dummy");
  TCanvas *Tcan2 = new TCanvas("Tcan2","",100,20,800,600); Tcan2->cd();  Tcan2->SetFillColor(0);
  TPad* pad1     = new TPad("pad1","The pad",0,0,0.98,1);
  applyPadStyle(pad1);

  TLegend *leg = new TLegend(0.60,0.21,0.87,0.44);
  applyLegStyle(leg);

  TString recoCut;
  TString l1Cut;
  TString tree ;

  Int_t color0 = TColor::GetColor("#C65D00"); //dark Orange
  Int_t color1 = TColor::GetColor("#283593"); //dark blue
  Int_t color2 = TColor::GetColor("#0288D1"); //medium blue

  Int_t color3 = TColor::GetColor("#00695C"); //green blue
  Int_t color4 = TColor::GetColor("#660066"); //Purple

  if(inputVariable!="pt"&&inputVariable!="eta"){
    std::cout<<"InputVariable "<< inputVariable<<" not found!! Exiting...."<<std::endl;
    exit(0);
  }

  TString variable;
  int bin;
  float xmin;
  float xmax;
  float ymax = 1.2;
  TString xlabel;
  bool variableBin = false;


  if(inputVariable=="pt"){
    variableBin=true;
    // DM="allDM";
    variable = "recoPt";
    bin = 20;
    xmin = 0;
    xmax = 80; 
    xlabel = lepton+" p_{T} [GeV]";
    recoCut = "recoPt>20&&abs(recoEta)<2&&l1Pt>20";
    l1Cut   = "recoPt>20&&abs(recoEta)<2&&l1Pt>20";
  //xlabel = "Reco #tau_{h} #eta";
  //xlabel = "Reco #tau_{h} #phi";
  }

  if(inputVariable=="eta"){
    variableBin=false;
    variable = "recoEta";
    bin = 20;
    xmin = -2;
    xmax = 2; 
    xlabel = lepton+" #eta [GeV]";
    recoCut = "recoPt>30&&abs(recoEta)<2&&l1Pt>30";
    l1Cut   = "recoPt>30&&abs(recoEta)<2&&l1Pt>30";
  }

  tree = folder+"/"+lepton+"Tree";
  TGraphAsymmErrors* eff_pt30_NoIso    = plotEfficiencyReturnASYM(variable, inFile, outFile+"/tmp/l1Pt", l1Cut+"&&recoPt>0", recoCut+"&&recoPt>0", xlabel, "NoIso", tree, bin,xmin,xmax,  color0  , variableBin);
  TGraphAsymmErrors* eff_pt30_VLoose   = plotEfficiencyReturnASYM(variable, inFile, outFile+"/tmp/l1PtVLoose_garbage_", l1Cut+"&&((l1Iso)/l1Pt)<"+isoNoTime, recoCut, xlabel, "p_{T} > 30GeV", tree, bin,xmin,xmax, 0  , variableBin);
  TGraphAsymmErrors* eff_pt30_VLoose_2 = plotEfficiencyReturnASYM(variable, inFile, outFile+"/tmp/l1PtNoTiming", l1Cut+"&&((l1Iso)/l1Pt)<"+isoNoTime, recoCut, xlabel, "No Timing", tree, bin,xmin,xmax, color1   , variableBin);
  TGraphAsymmErrors* eff_pt30_Loose    = plotEfficiencyReturnASYM(variable, inFile, outFile+"/tmp/l1PtTiming", l1Cut+"&&recoPt>0&&((l1Iso_time)/l1Pt)<"+isoTime, recoCut+"&&recoPt>0", xlabel, "Timing", tree, bin,xmin,xmax,  color2  , variableBin);

  eff_pt30_VLoose_2->Print();
  eff_pt30_Loose->Print();
  Tcan2->cd();
  pad1->cd();
  pad1->SetGrid(5,5); 

  eff_pt30_VLoose->Draw();
  eff_pt30_NoIso->Draw("p");
  //eff_pt30_Tight->Draw("p");
  //eff_pt30_Medium->Draw("p");
  eff_pt30_Loose->Draw("p");
  eff_pt30_VLoose_2->Draw("p");

  eff_pt30_VLoose->SetMaximum(ymax);

  leg->AddEntry(eff_pt30_NoIso, "#bf{No Iso}","pe");
  leg->AddEntry(eff_pt30_VLoose_2, "#bf{No Timing WP "+isoNoTime+"}","pe");
  leg->AddEntry(eff_pt30_Loose,  "#bf{Timing WP "+isoTime+"}","pe");
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);  
  leg->SetMargin(0.1);
  leg->Draw();

  TPaveText *pt = new TPaveText(4,.95,40,1.05);
  pt->AddText("#splitline{CMS}{#bf{#it{Preliminary}}}");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetShadowColor(0);
  //pt->Draw("E");

  eff_pt30_VLoose->GetXaxis()->SetRangeUser(0,eff_pt30_VLoose->GetXaxis()->GetXmax()-10);


  Tcan2->SaveAs(outFile+"/reco"+lepton+"-"+variable+label+".png");
  Tcan2->SaveAs(outFile+"/reco"+lepton+"-"+variable+label+".pdf");


}


void setMaxErrorTo1(TGraphAsymmErrors * hist){

  for(int i = 1; i<hist->GetN(); i++){
    Double_t errorY = hist->GetErrorY(i);
    Double_t pointx, pointy;
    if(hist->GetPoint(i,pointx,pointy)<0)
      std::cout<<"error getting point "<<std::endl;
    Double_t errorUp = pointy+errorY;
    Double_t errorLow = pointy-errorY/2;
    if(errorUp>1){
      hist->SetPointEYhigh(i, 1-pointy);
    }
  }

}


void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.2);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.1);
  //pad1->SetGrid(); 
  pad1->SetGrid(10,10); 
}
void applyLegStyle(TLegend *leg){
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetHeader("");
  leg->SetShadowColor(0);
}

using namespace std;

//added xaxis
TGraphAsymmErrors* plotEfficiencyReturnASYM(TString variable, TString fileName, TString outFileName, TString l1Cut, TString recoCut, TString xaxis, TString text, TString treePath, int bins, float low, float high, Color_t color, bool variableBin, int markerstyle, int inputfitparam){

  std::cout<<"variable "<<variable<<" filename "<<fileName<< " outFileName "<< outFileName<<" l1Cut "<< l1Cut<<" recoCut " << recoCut <<" xaxis "<< xaxis <<" text "<< text << " treePath "<< treePath<< " bins "<< bins<<" low "<< low << " high "<< high <<" color "<< color <<std::endl;

  setTDRStyle();

  TFile *tauFile    = new TFile(fileName);

  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }

  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600); Tcan->cd();  Tcan->SetFillColor(0);

  TPad* pad1 = new TPad("pad1","The pad",0,0,0.98,1);

  applyPadStyle(pad1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  TLegend *leg = new TLegend(0.25,0.75,0.55,0.99);//  Leg = new TLegend(0.74,0.65,0.99,0.95);                                                                                                                                                 

  applyLegStyle(leg);

  TTree* tauTree = (TTree*)tauFile->Get(treePath);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }

  TH1F* Denom;
  Float_t xbins[11]= {0,5,10,15,20,25,30,35,40,55,80};
  float xmax = 140;

  if(variableBin){
    Denom = new TH1F("Denom","Denom",10,xbins);
  }
  else
    Denom = new TH1F("Denom","Denom",bins,low,high);

  Denom->Sumw2();
  tauTree->Draw(variable+">>+Denom",recoCut);
  
  TH1F* Num;
  if(variableBin){
     Num = new TH1F("Num","Num",10,xbins);

  }
  else
    Num = new TH1F("Num","Num",bins,low,high);


  tauTree->Draw(variable+">>+Num",l1Cut);

  Num->Divide(Denom);
  gStyle->SetErrorX(0.5);

  TGraphAsymmErrors* asym = new TGraphAsymmErrors(Num);

  asym->GetXaxis()->SetTitle(xaxis);
  asym->GetYaxis()->SetTitle("L1 Efficiency ");
  asym->GetYaxis()->SetTitleOffset(0.8);

  TPaveText *pt = new TPaveText(150,.05,250,0.15);
  pt->AddText(text);
  //pt->AddText(text1.c_str());
  //pt->AddText(text2.c_str());
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetShadowColor(0);
  asym->SetMarkerStyle(markerstyle);
  asym->SetMarkerColor(color);
  asym->SetMarkerSize(1.5);
  asym->Draw("P");

  pt->Draw("E");

  if(color==0){
    asym->SetLineColorAlpha(color,1);  
    asym->SetMarkerColorAlpha(color,1);
    asym->SetFillColorAlpha(color,1);
    asym->SetMarkerStyle(1);
  }
  else{
    asym->SetLineColor(color);  
    asym->SetFillColor(color);
  }
  asym->SetFillStyle(1001);
  asym->SetLineWidth(2);
  asym->SetMaximum(1.1);
  asym->SetMinimum(0);
  setMaxErrorTo1(asym);


  gStyle->SetOptStat(0);

  //Tcan->cd();
  //Tcan->SaveAs(outFileName+".pdf");

  return asym;
}



void tdrGrid(bool gridOn) {
  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");
  tdrStyle->SetPadGridX(gridOn);
  tdrStyle->SetPadGridY(gridOn);
}

// fixOverlay: Redraws the axis

//void fixOverlay() {
//  gPad->RedrawAxis();
//}

void setTDRStyle() {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.1);
  
  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.16);
  tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  // tdrStyle->SetTitleX(0); // Set the position of the title box
  // tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  //tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(0.9);
  tdrStyle->SetTitleYOffset(1.25);
  // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(kTRUE);
  tdrStyle->SetTickLength(0.03, "XYZ");
  tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

// Postscript options:
  tdrStyle->SetPaperSize(20.,20.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();

}
