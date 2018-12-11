{
  gROOT->ProcessLine(".L plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp.C");
  gROOT->ProcessLine(".L tdrstyle.C");
  // variable,Rel,DM
  //     variable can be 'pt'
  //     Rel can be 'Rel',''
  //     DM can be 'allDM', '1prong', '3prong', '1prongPi0',
  //     Last two options are inFile and outFile

  TString inFile = "/scratch/ojalvo/phase2-trigger-MTD-Nov29/ZMM.root";
  TString outFile = "~/www/11_28_18/efficiency-mu";

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer50ps",inFile, outFile, "0.06", "0.1","WP_0.90_50ps3sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer50ps",inFile, outFile, "0.06", "0.1","WP_0.90_50ps3sigma");

  /*
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer50ps",inFile, outFile, "0.08", "0.12","WP_0.95_50ps3sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer50ps",inFile, outFile, "0.08", "0.12","WP_0.95_50ps3sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer50ps1sigma",inFile, outFile, "0.05", "0.125","WP_0.90_50ps1sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer50ps1sigma",inFile, outFile, "0.05", "0.125","WP_0.90_50ps1sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer50ps1sigma",inFile, outFile, "0.05", "0.125","WP_0.95_50ps1sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer50ps1sigma",inFile, outFile, "0.05", "0.125","WP_0.95_50ps1sigma");

  //2 sigma //90
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer50ps2sigma",inFile, outFile, "0.01", "0.08","WP_0.90_50ps2sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer50ps2sigma",inFile, outFile, "0.01", "0.08","WP_0.90_50ps2sigma");
  //95
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer50ps2sigma",inFile, outFile, "0.06", "0.12","WP_0.95_50ps2sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer50ps2sigma",inFile, outFile, "0.06", "0.12","WP_0.95_50ps2sigma");


  //////////////////// 30ps
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer30ps",inFile, outFile, "0.02", "0.08","WP_0.90_30ps3sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer30ps",inFile, outFile, "0.02", "0.08","WP_0.90_30ps3sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer30ps",inFile, outFile, "0.065", "0.12","WP_0.95_30ps3sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer30ps",inFile, outFile, "0.065", "0.12","WP_0.95_30ps3sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer30ps1sigma",inFile, outFile, "0.01", "0.08","WP_0.90_30ps1sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer30ps1sigma",inFile, outFile, "0.01", "0.08","WP_0.90_30ps1sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer30ps1sigma",inFile, outFile, "0.01", "0.12","WP_0.95_30ps1sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer30ps1sigma",inFile, outFile, "0.01", "0.12","WP_0.95_30ps1sigma");

  //2 sigma //90
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer30ps2sigma",inFile, outFile, "0.01", "0.08","WP_0.90_30ps2sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer30ps2sigma",inFile, outFile, "0.01", "0.08","WP_0.90_30ps2sigma");
  //95
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer30ps2sigma",inFile, outFile, "0.06", "0.12","WP_0.95_30ps2sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer30ps2sigma",inFile, outFile, "0.06", "0.12","WP_0.95_30ps2sigma");


  //////////////////// 100ps
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer100ps",inFile, outFile, "0.065", "0.08","WP_0.90_100ps3sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer100ps",inFile, outFile, "0.065", "0.08","WP_0.90_100ps3sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer100ps",inFile, outFile, "0.09", "0.12","WP_0.95_100ps3sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer100ps",inFile, outFile, "0.09", "0.12","WP_0.95_100ps3sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer100ps1sigma",inFile, outFile, "0.01", "0.08","WP_0.90_100ps1sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer100ps1sigma",inFile, outFile, "0.01", "0.08","WP_0.90_100ps1sigma");

  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer100ps1sigma",inFile, outFile, "0.055", "0.12","WP_0.95_100ps1sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer100ps1sigma",inFile, outFile, "0.055", "0.12","WP_0.95_100ps1sigma");

  //2 sigma //90
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer100ps2sigma",inFile, outFile, "0.01", "0.08","WP_0.90_100ps2sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer100ps2sigma",inFile, outFile, "0.01", "0.08","WP_0.90_100ps2sigma");
  //95
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("pt", "muon","","L1MTDAnalyzer100ps2sigma",inFile, outFile, "0.055", "0.12","WP_0.95_100ps2sigma");
  plotEfficiencyPtASYM_VersionCompIsoASYM_VersionComp("eta","muon","","L1MTDAnalyzer100ps2sigma",inFile, outFile, "0.055", "0.12","WP_0.95_100ps2sigma");
  */  
}
