{
  TString inFile = "/scratch/ojalvo/phase2-trigger-Dec3/SingleNeutrino-200PU.root";
  TString outFile = "/afs/hep.wisc.edu/home/ojalvo/www/11_28_18/";
  gROOT->ProcessLine(".L plotRatesIsoVersions.C");
  TString IsoType = "Iso";
  //plotRatesIsoVersions("l1Ele", "pt", "Electrons", inFile, outFile);
  TString folder="L1MTDAnalyzer100ps1sigma";
  //plotRatesIsoVersions("l1Ele", "pt", "Electrons", folder, inFile, outFile);
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer100ps1sigma", inFile, outFile,"100ps1sigma");
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer100ps2sigma", inFile, outFile,"100ps2sigma");
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer100ps", inFile, outFile,"100ps3sigma");

  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer30ps1sigma", inFile, outFile,"30ps1sigma");
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer30ps2sigma", inFile, outFile,"30ps2sigma");
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer30ps", inFile, outFile,"30ps3sigma");

  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer50ps1sigma", inFile, outFile,"50ps1sigma");
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer50ps2sigma", inFile, outFile,"50ps2sigma");
  plotRatesIsoVersions("l1Mu", "pt", "Muons", "L1MTDAnalyzer50ps", inFile, outFile,"50ps3sigma");

}
