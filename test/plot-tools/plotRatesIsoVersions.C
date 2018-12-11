
#include "plotRatesIsoVersions.h"

void plotRatesIsoVersions(TString type, TString variable, TString header, TString folder, TString inFile, TString outFile, TString name){
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/Level1/tauSeedStudy/plotRateReturnTH1F.C");
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/Level1/tauSeedStudy/plotResolution.C");
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/Level1/tauSeedStudy/plotTProfile.C");
  //gROOT->ProcessLine(".L ~/Documents/work/Analysis/Level1/tauSeedStudy/getNEntries.C");
  
  setTDRStyle();
  TFile *tauFile    = new TFile("dummy");
  TCanvas *Tcan2= new TCanvas("Tcan2","",100,20,600,600); Tcan2->cd();  Tcan2->SetFillColor(0);
  TPad* pad1 = new TPad("pad1","The pad",0,0.25,0.98,1);
  applyPadStyle(pad1);
  pad1->SetLeftMargin(0.15);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.025);

  TLegend *leg = new TLegend(0.64,0.512,0.941,0.870);//  Leg = new TLegend(0.74,0.65,0.99,0.95); 
  applyLegStyle(leg);

  Int_t color0 = TColor::GetColor("#c12600"); //dark Orange
  Int_t color1 = TColor::GetColor("#609dff"); //light blue
  Int_t color2 = TColor::GetColor("#ff9d00"); //dark blue
  Int_t color3 = TColor::GetColor("#0042ad"); //light orange
  Int_t color4 = TColor::GetColor("#c17700"); //dark orange

  TString file, Date, OptA, L1TypeA, L1TypeB, key;
  TString xaxis_label = "L1 p_{T} [GeV]";
  TString l1Cut = "";
  TString y_label = "Rate [kHz]";
  

//l1Ele_pt
  std::cout<<"variable to plot: "<<type+"_"+variable<<std::endl;
  TH1F* rate = plotRateReturnTH1F(type+"_"+variable,inFile,outFile+"/Rates/tmp/"+type+"_"+variable, "", xaxis_label, "NoIso", folder+"/", 300, -1 , 300, color0, true, 1, y_label);
  
  std::cout<<"variable to plot: "<<type+"_"+variable<<std::endl;  

  TH1F* IsoRel90 = plotRateReturnTH1F(type+"_IsoRel90_"+variable,inFile,outFile+"/Rates/tmp/"+type+"_IsoRel90_"+variable, "", xaxis_label, "IsoRel90", folder+"/", 300, -1 , 300, color1, true, 1, y_label);

  TH1F* IsoRel90_time = plotRateReturnTH1F(type+"_IsoRel90_time_"+variable,inFile,outFile+"/Rates/tmp/"+type+"_IsoRel90_time"+variable, "IsoRel90_time", xaxis_label, "IsoRel90_time", folder+"/", 300, -1 , 300, color2, true, 1, y_label);

  TH1F* IsoRel95 = plotRateReturnTH1F(type+"_IsoRel95_"+variable,inFile,outFile+"/Rates/tmp/"+type+"_IsoRel_95_"+variable, "", xaxis_label, "IsoRel95", folder+"/", 300, -1 , 300, color3, true, 1, y_label);

  TH1F* IsoRel95_time = plotRateReturnTH1F(type+"_IsoRel95_time_"+variable,inFile,outFile+"/Rates/tmp/"+type+"_IsoRel95_time_"+variable, "", xaxis_label, "IsoRel95_time", folder+"/", 300, -1 , 300, color4, true, 1, y_label);
  
  int nEvents  =34400; //= getNEntries(inFile, "L1MTDAnalyzer/nEvents");

  Tcan2->cd();
  pad1->cd();
  gPad->SetLogy();
  pad1->SetGrid(5,5); 

  float scaledRate = 0;
  float scaleOverall = 40000/5; // scale by 40 MHz
  float firstBin = 0;

  //Level 1
  float onlineRate = 1;

  firstBin   = rate->GetBinContent(1);
  scaledRate = firstBin/nEvents;
  std::cout << "bin content 1 " << firstBin << " scaled rate " << scaledRate<<std::endl;
  rate->Scale( 1/firstBin );
  rate->Scale( 1*onlineRate*scaledRate*scaleOverall );
  std::cout<<"first bin final "<< rate->GetBinContent(1)<<std::endl;

  firstBin   = IsoRel90->GetBinContent(1);
  scaledRate = firstBin/nEvents;
  std::cout << "bin content 1 " << firstBin << " scaled rate " << scaledRate<<std::endl;
  IsoRel90->Scale( 1/firstBin );
  IsoRel90->Scale( 1*onlineRate*scaledRate*scaleOverall );
  std::cout<<"first bin final "<< IsoRel90->GetBinContent(1)<<std::endl;

  firstBin   = IsoRel90_time->GetBinContent(1);
  scaledRate = firstBin/nEvents;
  std::cout << "bin content 1 " << firstBin << " scaled rate " << scaledRate<<std::endl;
  IsoRel90_time->Scale( 1/firstBin );
  IsoRel90_time->Scale( 1*onlineRate*scaledRate*scaleOverall );
  std::cout<<"first bin final "<< IsoRel90_time->GetBinContent(1)<<std::endl;

  firstBin   = IsoRel95->GetBinContent(1);
  scaledRate = firstBin/nEvents;
  std::cout << "bin content 1 " << firstBin << " scaled rate " << scaledRate<<std::endl;
  IsoRel95->Scale( 1/firstBin );
  IsoRel95->Scale( 1*onlineRate*scaledRate*scaleOverall );
  std::cout<<"first bin final "<< IsoRel95->GetBinContent(1)<<std::endl;

  firstBin   = IsoRel95_time->GetBinContent(1);
  scaledRate = firstBin/nEvents;
  std::cout << "bin content 1 " << firstBin << " scaled rate " << scaledRate<<std::endl;
  IsoRel95_time->Scale( 1/firstBin );
  IsoRel95_time->Scale( 1*onlineRate*scaledRate*scaleOverall );
  std::cout<<"first bin final "<< IsoRel95_time->GetBinContent(1)<<std::endl;

  IsoRel90->Draw();
  IsoRel90->SetMinimum(1);


  IsoRel90->GetYaxis()->SetTitleOffset(0.85);

  rate->Draw("same");
  IsoRel90->Draw("same");
  IsoRel90_time->Draw("same");
  IsoRel95->Draw("same");
  IsoRel95_time->Draw("same");

  rate->Rebin(2);
  rate->Scale(0.5);
  IsoRel90->Rebin(2);
  IsoRel90->Scale(0.5);
  IsoRel95->Rebin(2);
  IsoRel95->Scale(0.5);
  IsoRel90_time->Rebin(2);
  IsoRel90_time->Scale(0.5);
  IsoRel95_time->Rebin(2);
  IsoRel95_time->Scale(0.5);
  IsoRel90->GetXaxis()->SetRange(2,20);
  IsoRel90->SetMaximum(scaleOverall*100);
  IsoRel90->SetLineWidth(3);
  IsoRel95->SetLineWidth(3);
  IsoRel90_time->SetLineWidth(3);
  IsoRel95_time->SetLineWidth(3);


  TH1F* ratio = (TH1F*) IsoRel95->Clone();
  ratio->Sumw2();
  ratio->Add(IsoRel95_time,-1);
  ratio->GetXaxis()->SetRange(2,20);
  ratio->Divide(IsoRel95_time);
  ratio->SetLineStyle(7);
  ratio->SetLineColor(kGray+2);
  ratio->SetLineWidth(2);
  
  leg->SetHeader(header);
  leg->AddEntry(rate, "Non Isolated", "l");
  leg->AddEntry(IsoRel90, "WP90", "l");
  leg->AddEntry(IsoRel95, "WP95", "l");
  leg->AddEntry(IsoRel90_time,  "WP90 + timing ",  "l");
  leg->AddEntry(IsoRel95_time,  "WP95 + timing ",  "l");

  leg->Draw();


  Tcan2->cd();
  TPad* pad2 = new TPad("pad2","The lower pad",0,0,0.98,0.25);
  applyPadStyle(pad2);
  pad2->SetLeftMargin(0.15);  pad2->SetBottomMargin(0.13); pad2->SetRightMargin(0.025);
  pad2->cd();
  pad2->SetGrid(0,0); 

  ratio->Draw("p");

  ratio->GetYaxis()->SetTitleOffset(0.95);
  ratio->GetXaxis()->SetTitle("");
  ratio->GetYaxis()->SetTitle("#splitline{(Rate - Time Rate)}{          Rate}");
  //ratio->GetYaxis()->SetTitleSize(0.75);
  ratio->GetXaxis()->SetLabelSize(0.1);
  ratio->GetYaxis()->SetLabelSize(0.1);
  ratio->GetXaxis()->SetNdivisions(10);
  ratio->GetYaxis()->SetNdivisions(502);

  float low = 1;
  float high = 39;
  TLine *line0 = new TLine(low,0,high,0);
  line0->SetLineColor(kBlue);
  line0->Draw();

  TLine *line1 = new TLine(low,0.5,high,0.5);
  line1->SetLineColor(kBlue-7);
  line1->SetLineStyle(8);
  line1->Draw();

  TLine *line2 = new TLine(low,-0.5,high,-0.5);
  line2->SetLineColor(kBlue-7);
  line2->SetLineStyle(8);
  line2->Draw();

  ratio->Draw("psame");

  ratio->SetMaximum(0.7);
  ratio->SetMinimum(-0.7);

  Tcan2->SaveAs(outFile+"/Rates/"+type+variable+name+".png");
  Tcan2->SaveAs(outFile+"/Rates/"+type+variable+name+".pdf");


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


void setMaxErrorTo1(TGraphAsymmErrors * hist){

  for(int i = 1; i<hist->GetN(); i++){
    Double_t errorY = hist->GetErrorY(i);
    Double_t pointx, pointy;
    if(hist->GetPoint(i,pointx,pointy)<0)
      std::cout<<"error getting point "<<std::endl;
    Double_t errorUp = pointy+errorY/2;
    Double_t errorLow = pointy-errorY/2;
    std::cout<<"error Low "<< errorLow<<" errorUp "<<errorUp<< " geterrory "<< hist->GetErrorY(i)<<" geterrorx "<< hist->GetErrorX(i)<<std::endl;
    if(errorUp>1){
      hist->SetPointEYlow(i, errorLow);
      hist->SetPointEYhigh(i, 1-pointy);
    }
  }

}

void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.2);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.1);
  pad1->SetGrid(); 
}
void applyLegStyle(TLegend *leg){
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetHeader("");
  leg->SetShadowColor(0);
}

using namespace std;


Float_t xbins[16]= {0,10,16,20,24,28,32,36,44,52,60,76,100,140,200,280};
//added xaxis
TH1F* plotRateReturnTH1F(TString variable, TString fileName, TString outFileName, TString l1Cut, TString xaxis, TString text, TString treePath, int bins, float low, float high, Color_t color, bool variableBin=false, int linestyle=1, TString yaxis = "Rate"){

  setTDRStyle();
  //tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");

  TFile *tauFile    = new TFile(fileName);

  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }

  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600); Tcan->cd();  Tcan->SetFillColor(0);

  TPad* pad1 = new TPad("pad1","The pad",0,0,0.98,1);
  applyPadStyle(pad1);
  gPad->SetLogy();
  gStyle->SetOptStat(0);

  TLegend *leg = new TLegend(0.25,0.75,0.55,0.99);// 
  applyLegStyle(leg);
  TH1F* histo = (TH1F*)tauFile->Get(treePath + variable);
  if(histo == 0){
    std::cout<<"ERROR Tau histo is "<< histo <<" NOT FOUND; EXITING"<<std::endl;
    return 0;
  }
  //histo->Rebin(bins);
  
  TH1F* rateHisto = new TH1F("rateHisto","rateHisto",bins,low,high);
  Double_t Sum = 0;

  for(int i = bins; i > 0; i-- ){
    Sum += histo->GetBinContent(i);
    rateHisto->SetBinContent(i, Sum);
  }

  pad1->cd();
  rateHisto->Draw("L");

  gStyle->SetErrorX(0.25);
  rateHisto->Draw();
  rateHisto->GetXaxis()->SetTitle(xaxis);
  rateHisto->GetYaxis()->SetTitle(yaxis);

  
  TPaveText *pt = new TPaveText(150,.05,250,0.15);
  pt->AddText(text);
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetShadowColor(0);
  rateHisto->SetMarkerStyle(20);
  rateHisto->SetMarkerColor(color+2);
  pt->Draw("E");

  rateHisto->SetLineColor(color);  rateHisto->SetFillColor(0);
  rateHisto->SetFillStyle(0);
  rateHisto->SetLineWidth(2);
  rateHisto->SetLineStyle(linestyle);
  //rateHisto->Fit(g1,"R");
  //rateHisto->Fit(g2,"R+");
  //rateHisto->Fit(g0,"R+");
  //rateHisto->Fit(g,"R+");
  //rateHisto->Sumw2();
  //rateHisto->SetStats(0);
  rateHisto->SetMaximum(10000000);
  rateHisto->SetMinimum(10);
  //setMaxErrorTo1(Num);


  gStyle->SetOptStat(0);
  //
  Tcan->cd();
  Tcan->SaveAs(outFileName+".png");
  std::cout<<outFileName+".png"<<std::endl;
  return rateHisto;
}

int getNEntries(TString fileName, TString histName){
  TFile *file = new TFile(fileName);
  TH1F* nEventsHist =  (TH1F*)file->Get(histName);
  int nEvents = nEventsHist->GetEntries();
  std::cout<<"file: "<<fileName<<" nEvents: "<<nEvents<<std::endl;
  return nEvents;
};
