// system include files
#include <memory>

// Math Include
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <memory>
#include <math.h>
#include <vector>
#include <list>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//Vertex and gen particle
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1Trigger/interface/L1PFTau.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/Phase2L1ParticleFlow/interface/PFCandidate.h"

#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticle.h"
#include "DataFormats/L1TrackTrigger/interface/L1TkEtMissParticleFwd.h"

#include "L1Trigger/L1MTD/plugins/helpers.h"
#include "helpers.h"

using namespace l1t;
using namespace edm;
using namespace std;

struct object{float pt; float eta; float phi; float time; float isoSum; float isoSumTime;};
struct jet{float pt; float eta; float phi; float time;};

struct genVisTau{ reco::Candidate::LorentzVector p4; int decayMode;};

struct l1PFCandIso{l1t::PFCandidate object;
                   l1t::PFCandidateCollection isoCands;
};

class L1MTDPFAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

  explicit L1MTDPFAnalyzer(const edm::ParameterSet&);

  ~L1MTDPFAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //custom class functions
  void calculate_charged_iso_sum(l1t::PFCandidate object, int object_index, Handle< l1t::PFCandidateCollection > pfCands,  double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime);
  void calculate_charged_iso_sum_tau(l1t::L1PFTau object, Handle<l1t::PFCandidateCollection>  pfCandidates, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime);
  //void calculate_charged_iso_sum_tau(l1t::L1PFTau, Handle< l1t::PFCandidateCollection > pfCands,  double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime);
  void fillGenParticleVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<reco::GenParticle> &selectedGenParticles,int pdgID);
  void fillGenTauVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<genVisTau> &selectedGenParticles);
  // void fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates,PFCandidateCollection &pfElectrons,PFCandidateCollection &pfPhotons,PFCandidateCollection &pfMuons);
  void fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates,PFCandidateCollection &pfElectrons, std::vector<int> & l1PFElectrons_indicies, PFCandidateCollection &pfPhotons, std::vector<int> & l1PFPhotons_indicies, PFCandidateCollection &pfMuons, std::vector<int> & l1PFMuons_indicies);
  void findJets(Handle<l1t::PFCandidateCollection>  pfCandidates, float time_cut, std::vector<jet> &foundjets);
  void findLLJets(Handle<l1t::PFCandidateCollection>  pfCandidates, float time_cut, float minimum_offset, std::vector<jet> &foundjets);
  void zeroElectronTreeVariables();
  void zeroPhotonTreeVariables();
  void zeroMuonTreeVariables();
  void zeroTauTreeVariables();
  void fillLLRates(l1t::PFCandidate l1PFCand, float time_eta_calibrated);
  void fillRates( std::vector<object> l1Electrons, std::vector<object> l1Photons, std::vector<object> l1Muons, std::vector<object> l1Taus, std::vector<jet> l1Jets, std::vector<jet> l1Jets_time);

  //void fillL1TauTreeVariables(Handle<L1PFTauCollection> l1PFTaus, std::vector<genVisTau> genTaus, Handle<L1TkPrimaryVertexCollection> L1VertexHandle);
  void fillL1TauTreeVariables(Handle<L1PFTauCollection> l1PFTaus, Handle< l1t::PFCandidateCollection > l1PFCandidates, std::vector<genVisTau> genTaus, Handle<L1TkPrimaryVertexCollection> L1VertexHandle);

  EDGetTokenT< PFCandidateCollection >          L1PFCandsToken_;
  EDGetTokenT< L1PFTauCollection >              L1PFTausToken_;
  EDGetTokenT< L1TkEtMissParticleCollection >   L1PFMETToken_;
  EDGetTokenT< L1TkEtMissParticleCollection >   L1PFMETTimeToken_;
  EDGetTokenT< std::vector<pat::Electron> >     recoElectronsToken_; 
  EDGetTokenT< std::vector<pat::Photon> >       recoPhotonsToken_; 
  EDGetTokenT< std::vector<pat::Muon> >         recoMuonsToken_;     
  EDGetTokenT< std::vector<pat::Tau> >          recoTausToken_;
  EDGetTokenT< std::vector<pat::MET> >          recoMetToken_;
  EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
  EDGetTokenT< std::vector<reco::GenJet> >      genJetsToken_;
  EDGetTokenT< L1TkPrimaryVertexCollection >    pvToken_;

  double time_cut_;
  double iso_cone_deltaR_;
  double iso_cone_deltaZ_;
  // Declare all the variables for the trees
  //electron tree
  TTree* electronTree;
  double eleRecoPt, eleRecoEta, eleRecoPhi, eleGenPt, eleGenEta, eleGenPhi, eleL1Pt, eleL1Eta, eleL1Phi, eleL1Time, eleL1Iso, eleL1Iso_time;
  //photon tree
  TTree* photonTree;
  double gammaRecoPt, gammaRecoEta, gammaRecoPhi, gammaGenPt, gammaGenEta, gammaGenPhi, gammaL1Pt, gammaL1Eta, gammaL1Phi, gammaL1Time, gammaL1Iso, gammaL1Iso_time;
  //muon tree
  TTree* muonTree;
  double muRecoPt,  muRecoEta,  muRecoPhi,  muGenPt,  muGenEta,  muGenPhi,  muL1Pt,  muL1Eta,  muL1Phi, muL1Time, muL1Iso,  muL1Iso_time;
  //tau tree
  TTree* tauTree;
  double tauRecoPt, tauRecoEta, tauRecoPhi, tauGenPt, tauGenEta, tauGenPhi, tauL1Pt, tauL1Eta, tauL1Phi, tauL1Time, tauL1Iso, tauL1Iso_time;

  //tau tree
  TTree* L1TauTree;
  double track12DZ, track13DZ, track1PVDZ, track2PVDZ, track3PVDZ;
  double track1nStubs, track2nStubs, track3nStubs;
  double track1Time, track2Time,track3Time;
  double l1DM;
  double track1ChiSquared, track2ChiSquared, track3ChiSquared;
  double zVTX;
  double track1Z, track2Z, track3Z;
  double tauL1StripPt, tauL1StripEta, tauL1StripPhi, tauL1StripDR;
  double pfCand1HoE, pfCand2HoE, pfCand3HoE, tauL1nEG, tauL1EGPt, l1TauEGTime;
  

  //met tree
  TTree* metTree;
  double metL1Et, metL1Phi, metL1Et_time, metRecoEt, metRecoPhi, metGenEt, metGenPhi;
  
  TTree* jetTree;
  TTree* jetLLTree;
  double jetRecoPt, jetRecoEta, jetRecoPhi, jetGenPt, jetGenEta, jetGenPhi, jetL1Pt, jetL1Eta, jetL1Phi, jetL1Pt_time, jetL1Eta_time, jetL1Phi_time, jetL1Time;

  int run, lumi, event;

  TH1F* l1PFCandidates_TH1F;     
  TH1F* l1PFNeutral_TH1F;        
  TH1F* l1PFMuon_TH1F;        
			  
  TH1F* l1PFCandidates_1ns_TH1F; 
  TH1F* l1PFNeutral_1ns_TH1F; 
  TH1F* l1PFMuon_1ns_TH1F;           
			  
  TH1F* l1PFCandidates_2ns_TH1F; 
  TH1F* l1PFNeutral_2ns_TH1F;    
  TH1F* l1PFMuon_2ns_TH1F;           
			  
  TH1F* l1PFCandidates_3ns_TH1F; 
  TH1F* l1PFNeutral_3ns_TH1F;
  TH1F* l1PFMuon_3ns_TH1F;               
			  
  TH1F* l1PFCandidates_4ns_TH1F; 
  TH1F* l1PFNeutral_4ns_TH1F;
  TH1F* l1PFMuon_4ns_TH1F;               
			  
  TH1F* l1PFCandidates_5ns_TH1F; 
  TH1F* l1PFNeutral_5ns_TH1F;    
  TH1F* l1PFMuon_5ns_TH1F;           
			  
  TH1F* l1PFCandidates_6ns_TH1F; 
  TH1F* l1PFNeutral_6ns_TH1F;    
  TH1F* l1PFMuon_6ns_TH1F;           
			  
  TH1F* l1PFCandidates_7ns_TH1F; 
  TH1F* l1PFNeutral_7ns_TH1F;    
  TH1F* l1PFMuon_7ns_TH1F;           
			  
  TH1F* l1PFCandidates_8ns_TH1F; 
  TH1F* l1PFNeutral_8ns_TH1F;    
  TH1F* l1PFMuon_8ns_TH1F;           
			  
  TH1F* l1PFCandidates_9ns_TH1F; 
  TH1F* l1PFNeutral_9ns_TH1F;    
			  
  TH1F* l1PFCandidates_10ns_TH1F;
  TH1F* l1PFNeutral_10ns_TH1F;   

  TH1F* nEvents;
  /*
  TH2F* l1PFCandidates_TH2F;    
  TH2F* l1PFNeutral_TH2F;       
  TH2F* l1PFCharged_TH2F;       
  TH2F* l1PFNeutralHadrons_TH2F;
  TH2F* l1PFChargedHadrons_TH2F;
  TH2F* l1PFMuons_TH2F;         
  TH2F* l1PFElectrons_TH2F;     
  TH2F* l1PFPhotons_TH2F;       
  TH2F* l1PFTaus_TH2F;          
  TH2F* l1PFJets_TH2F;          

  TH2F* l1PFCandidatesFlat_TH2F;    
  TH2F* l1PFNeutralFlat_TH2F;       
  TH2F* l1PFChargedFlat_TH2F;       
  TH2F* l1PFNeutralHadronsFlat_TH2F;
  TH2F* l1PFChargedHadronsFlat_TH2F;
  TH2F* l1PFMuonsFlat_TH2F;         
  TH2F* l1PFElectronsFlat_TH2F;     
  TH2F* l1PFPhotonsFlat_TH2F;       
  TH2F* l1PFTausFlat_TH2F;          
  TH2F* l1PFJetsFlat_TH2F;          

  TH2F* l1PFCandidates10_TH2F;    
  TH2F* l1PFNeutral10_TH2F;       
  TH2F* l1PFCharged10_TH2F;       
  TH2F* l1PFNeutralHadrons10_TH2F;
  TH2F* l1PFChargedHadrons10_TH2F;
  TH2F* l1PFMuons10_TH2F;         
  TH2F* l1PFElectrons10_TH2F;     
  TH2F* l1PFPhotons10_TH2F;       
  TH2F* l1PFTaus10_TH2F;          
  TH2F* l1PFJets10_TH2F;          

  TH2F* l1PFCandidates20_TH2F;    
  TH2F* l1PFNeutral20_TH2F;       
  TH2F* l1PFCharged20_TH2F;       
  TH2F* l1PFNeutralHadrons20_TH2F;
  TH2F* l1PFChargedHadrons20_TH2F;
  TH2F* l1PFMuons20_TH2F;         
  TH2F* l1PFElectrons20_TH2F;     
  TH2F* l1PFPhotons20_TH2F;       
  TH2F* l1PFTaus20_TH2F;          
  TH2F* l1PFJets20_TH2F;          
  */

  TH1F* l1PFTauRate_Barrel      ;
  TH1F* l1DiPFTauRate_Barrel      ;
  TH1F* l1PFMuonRate_Barrel     ;
  TH1F* l1PFElectronRate_Barrel ;
  TH1F* l1PFPhotonRate_Barrel   ;

  TH1F* l1PFTauRate_Iso90_Barrel      ;
  TH1F* l1DiPFTauRate_Iso90_Barrel      ;
  TH1F* l1PFMuonRate_Iso90_Barrel     ;
  TH1F* l1PFElectronRate_Iso90_Barrel ;
  TH1F* l1PFPhotonRate_Iso90_Barrel   ;

  TH1F* l1PFTauRate_Iso90TIME_Barrel      ;
  //    l1DiPFTauRate_Iso90Time_Barrel
  TH1F* l1DiPFTauRate_Iso90TIME_Barrel      ;
  TH1F* l1PFMuonRate_Iso90TIME_Barrel     ;
  TH1F* l1PFElectronRate_Iso90TIME_Barrel ;
  TH1F* l1PFPhotonRate_Iso90TIME_Barrel   ;

  TH1F* l1PFTauRate_Iso95_Barrel      ;
  TH1F* l1DiPFTauRate_Iso95_Barrel      ;
  TH1F* l1PFMuonRate_Iso95_Barrel     ;
  TH1F* l1PFElectronRate_Iso95_Barrel ;
  TH1F* l1PFPhotonRate_Iso95_Barrel   ;

  TH1F* l1PFTauRate_Iso95TIME_Barrel      ;
  TH1F* l1DiPFTauRate_Iso95TIME_Barrel      ;
  TH1F* l1PFMuonRate_Iso95TIME_Barrel     ;
  TH1F* l1PFElectronRate_Iso95TIME_Barrel ;
  TH1F* l1PFPhotonRate_Iso95TIME_Barrel   ;

  TH1F* l1PFJetRate_Barrel      ;
  TH1F* l1PFJetRate_time_Barrel      ;
  TH1F* l1PFJet2nsRate_Barrel   ;
  TH1F* l1PFJet3nsRate_Barrel   ;
  TH1F* l1PFJet4nsRate_Barrel   ;
  TH1F* l1PFJet5nsRate_Barrel   ;
  TH1F* l1PFJet6nsRate_Barrel   ;
  TH1F* l1PFJet7nsRate_Barrel   ;


  TH1F* l1PFJetRate_calibrated_Barrel      ;
  TH1F* l1PFJetRate_time_calibrated_Barrel      ;
  TH1F* l1PFJet2nsRate_calibrated_Barrel   ;
  TH1F* l1PFJet3nsRate_calibrated_Barrel   ;
  TH1F* l1PFJet4nsRate_calibrated_Barrel   ;
  TH1F* l1PFJet5nsRate_calibrated_Barrel   ;
  TH1F* l1PFJet6nsRate_calibrated_Barrel   ;
  TH1F* l1PFJet7nsRate_calibrated_Barrel   ;

  TH1F* l1PFTauRate_Iso90TIME      ;
  TH1F* l1DiPFTauRate_Iso90TIME      ;
  TH1F* l1PFMuonRate_Iso90TIME     ;
  TH1F* l1PFElectronRate_Iso90TIME ;
  TH1F* l1PFPhotonRate_Iso90TIME   ;

  TH1F* l1PFTauRate_Iso90      ;
  TH1F* l1DiPFTauRate_Iso90      ;
  TH1F* l1PFMuonRate_Iso90     ;
  TH1F* l1PFElectronRate_Iso90 ;
  TH1F* l1PFPhotonRate_Iso90   ;

  TH1F* l1PFTauRate_Iso95TIME      ;
  TH1F* l1DiPFTauRate_Iso95TIME      ;
  TH1F* l1PFMuonRate_Iso95TIME     ;
  TH1F* l1PFElectronRate_Iso95TIME ;
  TH1F* l1PFPhotonRate_Iso95TIME   ;

  TH1F* l1PFTauRate_Iso95      ;
  TH1F* l1DiPFTauRate_Iso95      ;
  TH1F* l1PFMuonRate_Iso95     ;
  TH1F* l1PFElectronRate_Iso95 ;
  TH1F* l1PFPhotonRate_Iso95   ;

  TH1F* l1PFTauRate      ;
  TH1F* l1DiPFTauRate      ;
  TH1F* l1PFMuonRate     ;
  TH1F* l1PFElectronRate ;
  TH1F* l1PFPhotonRate   ;

  TH1F* l1PFJetRate      ;
  TH1F* l1PFJetRate_time ;
  TH1F* l1PFJet2nsRate   ;
  TH1F* l1PFJet3nsRate   ;
  TH1F* l1PFJet4nsRate   ;
  TH1F* l1PFJet5nsRate   ;
  TH1F* l1PFJet6nsRate   ;
  TH1F* l1PFJet7nsRate   ;

};

//Constructor
L1MTDPFAnalyzer::L1MTDPFAnalyzer(const edm::ParameterSet &cfg) :  //inputs L1PFCands, reco electrons, reco muons, reco taus, gen particles
  //timingValuesToken_(  consumes<edm::ValueMap<float> >(            cfg.getParameter<edm::InputTag>("timingValuesNominal"))),
  pvToken_(consumes<L1TkPrimaryVertexCollection>(                  cfg.getParameter<InputTag>("L1VertexInputTag"))),
  L1PFCandsToken_(     consumes< std::vector<l1t::PFCandidate>  >( cfg.getParameter<InputTag>("l1PFCands")     )),
  L1PFTausToken_(      consumes< L1PFTauCollection >(              cfg.getParameter<InputTag>("l1PFTaus")      )),
  L1PFMETToken_(       consumes< L1TkEtMissParticleCollection >(   cfg.getParameter<InputTag>("l1PFMET")       )),
  L1PFMETTimeToken_(   consumes< L1TkEtMissParticleCollection >(   cfg.getParameter<InputTag>("l1PFMETTime")   )),
  recoElectronsToken_( consumes< std::vector<pat::Electron> >(     cfg.getParameter<InputTag>("recoElectrons") )),
  recoPhotonsToken_(    consumes< std::vector<pat::Photon> >(      cfg.getParameter<InputTag>("recoPhotons")   )),
  recoMuonsToken_(     consumes< std::vector<pat::Muon> >(         cfg.getParameter<InputTag>("recoMuons")     )),
  recoTausToken_(      consumes< std::vector<pat::Tau> >(          cfg.getParameter<InputTag>("recoTaus")      )),
  recoMetToken_(       consumes< std::vector<pat::MET> >(          cfg.getParameter<InputTag>("recoMet")       )),
  genParticlesToken_(  consumes< std::vector<reco::GenParticle> >( cfg.getParameter<InputTag>("genParticles")  )),
  genJetsToken_(       consumes< vector<reco::GenJet> >(           cfg.getParameter<InputTag>("genJets")       )),
  time_cut_(           cfg.getParameter<double>("time_cut")),
  iso_cone_deltaR_(    cfg.getParameter<double>("isoConeDeltaR")),
  iso_cone_deltaZ_(    cfg.getParameter<double>("isoConeDeltaZ"))
{
  //services
  usesResource("TFileService");
  Service<TFileService> fs;

  
  //Trees
  electronTree = fs->make<TTree>("electronTree","Electron Efficiency Tree" );
  electronTree->Branch("run",          &run,               "run/I"         );
  electronTree->Branch("lumi",         &lumi,              "lumi/I"        );
  electronTree->Branch("event",        &event,             "event/I"       );
  electronTree->Branch("recoPt",       &eleRecoPt,         "recoPt/D"      );
  electronTree->Branch("recoEta",      &eleRecoEta,        "recoEta/D"     );
  electronTree->Branch("recoPhi",      &eleRecoPhi,        "recoPhi/D"     );
  electronTree->Branch("genPt",        &eleGenPt,          "genPt/D"       );
  electronTree->Branch("genEta",       &eleGenEta,         "genEta/D"      );
  electronTree->Branch("genPhi",       &eleGenPhi,         "genPhi/D"      );
  electronTree->Branch("l1Pt",         &eleL1Pt,           "l1Pt/D"        );
  electronTree->Branch("l1Eta",        &eleL1Eta,          "l1Eta/D"       );
  electronTree->Branch("l1Phi",        &eleL1Phi,          "l1Phi/D"       );
  electronTree->Branch("l1Time",       &eleL1Time,         "l1Time/D"      );
  electronTree->Branch("l1Iso",        &eleL1Iso,          "l1Iso/D"       );
  electronTree->Branch("l1Iso_time",   &eleL1Iso_time,     "l1Iso_time/D"  );

  photonTree = fs->make<TTree>("photonTree","Photon Efficiency Tree" );
  photonTree->Branch("run",          &run,                 "run/I"         );
  photonTree->Branch("lumi",         &lumi,                "lumi/I"        );
  photonTree->Branch("event",        &event,               "event/I"       );
  photonTree->Branch("recoPt",       &gammaRecoPt,         "recoPt/D"      );
  photonTree->Branch("recoEta",      &gammaRecoEta,        "recoEta/D"     );
  photonTree->Branch("recoPhi",      &gammaRecoPhi,        "recoPhi/D"     );
  photonTree->Branch("genPt",        &gammaGenPt,          "genPt/D"       );
  photonTree->Branch("genEta",       &gammaGenEta,         "genEta/D"      );
  photonTree->Branch("genPhi",       &gammaGenPhi,         "genPhi/D"      );
  photonTree->Branch("l1Pt",         &gammaL1Pt,           "l1Pt/D"        );
  photonTree->Branch("l1Eta",        &gammaL1Eta,          "l1Eta/D"       );
  photonTree->Branch("l1Phi",        &gammaL1Phi,          "l1Phi/D"       );
  photonTree->Branch("l1Time",       &gammaL1Time,         "l1Time/D"      );
  photonTree->Branch("l1Iso",        &gammaL1Iso,          "l1Iso/D"       );
  photonTree->Branch("l1Iso_time",   &gammaL1Iso_time,     "l1Iso_time/D"  );
  
  muonTree     = fs->make<TTree>("muonTree"    ,"Muon Efficiency Tree"     );
  muonTree->Branch("run",              &run,               "run/D"         );
  muonTree->Branch("lumi",             &lumi,              "lumi/D"        );
  muonTree->Branch("event",            &event,             "event/D"       );
  muonTree->Branch("recoPt",           &muRecoPt,          "recoPt/D"      );
  muonTree->Branch("recoEta",          &muRecoEta,         "recoEta/D"     );
  muonTree->Branch("recoPhi",          &muRecoPhi,         "recoPhi/D"     );
  muonTree->Branch("genPt",            &muGenPt,           "genPt/D"       );
  muonTree->Branch("genEta",           &muGenEta,          "genEta/D"      );
  muonTree->Branch("genPhi",           &muGenPhi,          "genPhi/D"      );
  muonTree->Branch("l1Pt",             &muL1Pt,            "l1Pt/D"        );
  muonTree->Branch("l1Eta",            &muL1Eta,           "l1Eta/D"       );
  muonTree->Branch("l1Phi",            &muL1Phi,           "l1Phi/D"       );
  muonTree->Branch("l1Time",           &muL1Time,          "l1Time/D"      );
  muonTree->Branch("l1Iso",            &muL1Iso,           "l1Iso/D"       );
  muonTree->Branch("l1Iso_time",       &muL1Iso_time,      "l1Iso_time/D"  );
  
  tauTree      = fs->make<TTree>("tauTree"     ,"Tau Efficiency Tree"      );
  tauTree->Branch("run",               &run,               "run/D"         );
  tauTree->Branch("lumi",              &lumi,              "lumi/D"        );
  tauTree->Branch("event",             &event,             "event/D"       );
  tauTree->Branch("recoPt",            &tauRecoPt,         "recoPt/D"      );
  tauTree->Branch("recoEta",           &tauRecoEta,        "recoEta/D"     );
  tauTree->Branch("recoPhi",           &tauRecoPhi,        "recoPhi/D"     );
  tauTree->Branch("genPt",             &tauGenPt,          "genPt/D"       );
  tauTree->Branch("genEta",            &tauGenEta,         "genEta/D"      );
  tauTree->Branch("genPhi",            &tauGenPhi,         "genPhi/D"      );
  tauTree->Branch("l1Pt",              &tauL1Pt,           "l1Pt/D"        );
  tauTree->Branch("l1Eta",             &tauL1Eta,          "l1Eta/D"       );
  tauTree->Branch("l1Phi",             &tauL1Phi,          "l1Phi/D"       );
  tauTree->Branch("l1Time",            &tauL1Time,         "l1Time/D"      );
  tauTree->Branch("l1Iso",             &tauL1Iso,          "l1Iso/D"       );
  tauTree->Branch("l1Iso_time",        &tauL1Iso_time,     "l1Iso_time/D"  );

  L1TauTree      = fs->make<TTree>("L1TauTree"     ,"L1 Tau Tree"      );
  L1TauTree->Branch("run",               &run,               "run/D"         );
  L1TauTree->Branch("lumi",              &lumi,              "lumi/D"        );
  L1TauTree->Branch("event",             &event,             "event/D"       );
  L1TauTree->Branch("genPt",             &tauGenPt,          "genPt/D"       );
  L1TauTree->Branch("genEta",            &tauGenEta,         "genEta/D"      );
  L1TauTree->Branch("genPhi",            &tauGenPhi,         "genPhi/D"      );
  L1TauTree->Branch("l1Pt",              &tauL1Pt,           "l1Pt/D"        );
  L1TauTree->Branch("l1Eta",             &tauL1Eta,          "l1Eta/D"       );
  L1TauTree->Branch("l1Phi",             &tauL1Phi,          "l1Phi/D"       );
  L1TauTree->Branch("l1Time",            &tauL1Time,         "l1Time/D"      );
  L1TauTree->Branch("l1Iso",             &tauL1Iso,          "l1Iso/D"       );
  L1TauTree->Branch("l1Iso_time",        &tauL1Iso_time,     "l1Iso_time/D"  );
  L1TauTree->Branch("track12DZ",         &track12DZ,         "track12DZ/D"   );
  L1TauTree->Branch("track13DZ",         &track13DZ,         "track13DZ/D"   );
  L1TauTree->Branch("track1PVDZ",             &track1PVDZ,          "track1PVDZ/D"       );
  L1TauTree->Branch("track2PVDZ",             &track2PVDZ,          "track2PVDZ/D"       );
  L1TauTree->Branch("track3PVDZ",             &track3PVDZ,          "track3PVDZ/D"       );
  L1TauTree->Branch("track1nStubs",             &track1nStubs,          "track1nStubs/D"       );
  L1TauTree->Branch("track2nStubs",             &track2nStubs,          "track2nStubs/D"       );
  L1TauTree->Branch("track3nStubs",             &track3nStubs,          "track3nStubs/D"       );
  L1TauTree->Branch("track1Time",             &track1Time,          "track1Time/D"       );
  L1TauTree->Branch("track2Time",             &track2Time,          "track2Time/D"       );
  L1TauTree->Branch("track3Time",             &track3Time,          "track3Time/D"       );
  L1TauTree->Branch("l1DecayMode",             &l1DM,          "l1DecayMode/D"       );
  L1TauTree->Branch("track1ChiSquared",             &track1ChiSquared,          "track1ChiSquared/D"       );
  L1TauTree->Branch("track2ChiSquared",             &track2ChiSquared,          "track2ChiSquared/D"       );
  L1TauTree->Branch("track3ChiSquared",             &track3ChiSquared,          "track3ChiSquared/D"       );
  L1TauTree->Branch("zVTX", &zVTX, "zVTX/D");
  L1TauTree->Branch("track1Z", &track1Z, "track1Z/D");
  L1TauTree->Branch("track2Z", &track2Z, "track2Z/D");
  L1TauTree->Branch("track3Z", &track3Z, "track3Z/D");
  L1TauTree->Branch("tauL1StripPt", &tauL1StripPt, "tauL1StripPt/D");
  L1TauTree->Branch("tauL1StripEta", &tauL1StripEta, "tauL1StripEta/D");
  L1TauTree->Branch("tauL1StripPhi", &tauL1StripPhi, "tauL1StripPhi/D");
  L1TauTree->Branch("tauL1StripDR", &tauL1StripDR, "tauL1StripDR/D");
  L1TauTree->Branch("pfCand1HoE", &pfCand1HoE, "pfCand1HoE/D");
  L1TauTree->Branch("pfCand2HoE", &pfCand2HoE, "pfCand2HoE/D");
  L1TauTree->Branch("pfCand3HoE", &pfCand3HoE, "pfCand3HoE/D");
  L1TauTree->Branch("tauL1nEG", &tauL1nEG, "tauL1nEG/D");
  L1TauTree->Branch("tauL1EGPt", &tauL1EGPt, "tauL1EGPt/D");
  L1TauTree->Branch("l1TauEGTime", &l1TauEGTime, "l1TauEGTime/D");
  //L1TauTree->Branch("",             &,          "/D"       );
  //

  jetTree      = fs->make<TTree>("jetTree"     ,"Jet Efficiency Tree"      );
  jetTree->Branch("run",               &run,               "run/D"         );
  jetTree->Branch("lumi",              &lumi,              "lumi/D"        );
  jetTree->Branch("event",             &event,             "event/D"       );
  jetTree->Branch("recoPt",            &jetRecoPt,         "recoPt/D"      );
  jetTree->Branch("recoEta",           &jetRecoEta,        "recoEta/D"     );
  jetTree->Branch("recoPhi",           &jetRecoPhi,        "recoPhi/D"     );
  jetTree->Branch("genPt",             &jetGenPt,          "genPt/D"       );
  jetTree->Branch("genEta",            &jetGenEta,         "genEta/D"      );
  jetTree->Branch("genPhi",            &jetGenPhi,         "genPhi/D"      );
  jetTree->Branch("l1Pt",              &jetL1Pt,           "l1Pt/D"        );
  jetTree->Branch("l1Eta",             &jetL1Eta,          "l1Eta/D"       );
  jetTree->Branch("l1Phi",             &jetL1Phi,          "l1Phi/D"       );
  jetTree->Branch("l1Pt_time",         &jetL1Pt_time,      "l1Pt_time/D"   );
  jetTree->Branch("l1Eta_time",        &jetL1Eta_time,     "l1Eta_time/D"  );
  jetTree->Branch("l1Phi_time",        &jetL1Phi_time,     "l1Phi_time/D"  );
  jetTree->Branch("l1Time",            &jetL1Time,         "l1Time/D"      );

  jetLLTree      = fs->make<TTree>("jetLLTree"     ,"Jet LL Efficiency Tree" );
  jetLLTree->Branch("run",               &run,               "run/D"         );
  jetLLTree->Branch("lumi",              &lumi,              "lumi/D"        );
  jetLLTree->Branch("event",             &event,             "event/D"       );
  jetLLTree->Branch("recoPt",            &jetRecoPt,         "recoPt/D"      );
  jetLLTree->Branch("recoEta",           &jetRecoEta,        "recoEta/D"     );
  jetLLTree->Branch("recoPhi",           &jetRecoPhi,        "recoPhi/D"     );
  jetLLTree->Branch("genPt",             &jetGenPt,          "genPt/D"       );
  jetLLTree->Branch("genEta",            &jetGenEta,         "genEta/D"      );
  jetLLTree->Branch("genPhi",            &jetGenPhi,         "genPhi/D"      );
  jetLLTree->Branch("l1Pt",              &jetL1Pt,           "l1Pt/D"        );
  jetLLTree->Branch("l1Eta",             &jetL1Eta,          "l1Eta/D"       );
  jetLLTree->Branch("l1Phi",             &jetL1Phi,          "l1Phi/D"       );
  jetLLTree->Branch("l1Time",            &jetL1Time,         "l1Time/D"      );

  metTree      = fs->make<TTree>("metTree"     ,"Met Efficiency Tree"      );
  metTree->Branch("run",               &run,               "run/D"         );
  metTree->Branch("lumi",              &lumi,              "lumi/D"        );
  metTree->Branch("event",             &event,             "event/D"       );
  metTree->Branch("recoEt",            &metRecoEt,         "recoEt/D"      );
  metTree->Branch("recoPhi",           &metRecoPhi,        "recoPhi/D"     );
  metTree->Branch("genEt",             &metGenEt,          "genEt/D"       );
  metTree->Branch("genPhi",            &metGenPhi,         "genPhi/D"      );
  metTree->Branch("l1Et",              &metL1Et,           "l1Et/D"        );
  metTree->Branch("l1Phi",             &metL1Phi,          "l1Phi/D"       );
  metTree->Branch("l1Et_time",         &metL1Et_time,      "metL1Et_time/D");

  nEvents                  = fs->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  l1PFCandidates_TH1F      = fs->make<TH1F>( "l1PFCandidatesTH1F"      ,"l1PFCandidatesTH1F",     100,  0., 100. );
  l1PFNeutral_TH1F         = fs->make<TH1F>( "l1PFNeutralTH1F"         ,"l1PFNeutralTH1F",        100,  0., 100. );
  l1PFMuon_TH1F            = fs->make<TH1F>( "l1PFMuonTH1F"         ,"l1PFMuonTH1F",        100,  0., 100. );

  l1PFCandidates_1ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_1nsTH1F"  ,"l1PFCandidates_1nsTH1F", 100,  0., 100. );
  l1PFNeutral_1ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_1nsTH1F"     ,"l1PFNeutral_1nsTH1F",    100,  0., 100. );
  l1PFMuon_1ns_TH1F     = fs->make<TH1F>( "l1PFMuon_1nsTH1F"     ,"l1PFMuon_1nsTH1F",    100,  0., 100. );

  l1PFCandidates_2ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_2nsTH1F"  ,"l1PFCandidates_2nsTH1F", 100,  0., 100. );
  l1PFNeutral_2ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_2nsTH1F"     ,"l1PFNeutral_2nsTH1F",    100,  0., 100. );
  l1PFMuon_2ns_TH1F     = fs->make<TH1F>( "l1PFMuon_2nsTH1F"     ,"l1PFMuon_2nsTH1F",    100,  0., 100. );

  l1PFCandidates_3ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_3nsTH1F"  ,"l1PFCandidates_3nsTH1F", 100,  0., 100. );
  l1PFNeutral_3ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_3nsTH1F"     ,"l1PFNeutral_3nsTH1F",    100,  0., 100. );
  l1PFMuon_3ns_TH1F     = fs->make<TH1F>( "l1PFMuon_3nsTH1F"     ,"l1PFMuon_3nsTH1F",    100,  0., 100. );

  l1PFCandidates_4ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_4nsTH1F"  ,"l1PFCandidates_4nsTH1F", 100,  0., 100. );
  l1PFNeutral_4ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_4nsTH1F"     ,"l1PFNeutral_4nsTH1F",    100,  0., 100. );
  l1PFMuon_4ns_TH1F     = fs->make<TH1F>( "l1PFMuon_4nsTH1F"     ,"l1PFMuon_4nsTH1F",    100,  0., 100. );

  l1PFCandidates_5ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_5nsTH1F"  ,"l1PFCandidates_5nsTH1F", 100,  0., 100. );
  l1PFNeutral_5ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_5nsTH1F"     ,"l1PFNeutral_5nsTH1F",    100,  0., 100. );
  l1PFMuon_5ns_TH1F     = fs->make<TH1F>( "l1PFMuon_5nsTH1F"     ,"l1PFMuon_5nsTH1F",    100,  0., 100. );

  l1PFCandidates_6ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_6nsTH1F"  ,"l1PFCandidates_6nsTH1F", 100,  0., 100. );
  l1PFNeutral_6ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_6nsTH1F"     ,"l1PFNeutral_6nsTH1F",    100,  0., 100. );
  l1PFMuon_6ns_TH1F     = fs->make<TH1F>( "l1PFMuon_6nsTH1F"     ,"l1PFMuon_6nsTH1F",    100,  0., 100. );

  l1PFCandidates_7ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_7nsTH1F"  ,"l1PFCandidates_7nsTH1F", 100,  0., 100. );
  l1PFNeutral_7ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_7nsTH1F"     ,"l1PFNeutral_7nsTH1F",    100,  0., 100. );
  l1PFMuon_7ns_TH1F     = fs->make<TH1F>( "l1PFMuon_7nsTH1F"     ,"l1PFMuon_7nsTH1F",    100,  0., 100. );

  l1PFCandidates_8ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_8nsTH1F"  ,"l1PFCandidates_8nsTH1F", 100,  0., 100. );
  l1PFNeutral_8ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_8nsTH1F"     ,"l1PFNeutral_8nsTH1F",    100,  0., 100. );
  l1PFMuon_8ns_TH1F     = fs->make<TH1F>( "l1PFMuon_8nsTH1F"     ,"l1PFMuon_8nsTH1F",    100,  0., 100. );

  l1PFCandidates_9ns_TH1F  = fs->make<TH1F>( "l1PFCandidates_9nsTH1F"  ,"l1PFCandidates_9nsTH1F", 100,  0., 100. );
  l1PFNeutral_9ns_TH1F     = fs->make<TH1F>( "l1PFNeutral_9nsTH1F"     ,"l1PFNeutral_9nsTH1F",    100,  0., 100. );

  l1PFCandidates_10ns_TH1F = fs->make<TH1F>( "l1PFCandidates_10nsTH1F" ,"l1PFCandidates_10nsTH1F",100,  0., 100. );
  l1PFNeutral_10ns_TH1F    = fs->make<TH1F>( "l1PFNeutral_10nsTH1F"    ,"l1PFNeutral_10nsTH1F",   100,  0., 100. );


  /*
  l1PFCandidates_TH2F     = fs->make<TH2F>( "l1PFCandidates",     "l1PFCandidates",     50, -4, 4, 50, 0, 22);
  l1PFNeutral_TH2F        = fs->make<TH2F>( "l1PFNeutrals",       "l1PFNeutrals",       50, -4, 4, 50, 0, 22);
  l1PFCharged_TH2F        = fs->make<TH2F>( "l1PFCharged",        "l1PFCharged",        50, -4, 4, 50, 0, 22);
  l1PFNeutralHadrons_TH2F = fs->make<TH2F>( "l1PFNeutralHadrons", "l1PFNeutralHadrons", 50, -4, 4, 50, 0, 22);
  l1PFChargedHadrons_TH2F = fs->make<TH2F>( "l1PFChargedHadrons", "l1PFChargedHadrons", 50, -4, 4, 50, 0, 22);
  l1PFMuons_TH2F          = fs->make<TH2F>( "l1PFMuons",          "l1PFMuons",          50, -4, 4, 50, 0, 22);
  l1PFElectrons_TH2F      = fs->make<TH2F>( "l1PFElectrons",      "l1PFElectrons",      50, -4, 4, 50, 0, 22);
  l1PFPhotons_TH2F        = fs->make<TH2F>( "l1PFPhotons",        "l1PFPhotons",        50, -4, 4, 50, 0, 22); 
  l1PFTaus_TH2F           = fs->make<TH2F>( "l1PFTaus",           "l1PFTaus",           50, -4, 4, 50, 0, 22);
  l1PFJets_TH2F           = fs->make<TH2F>( "l1PFJets",           "l1PFJets",           50, -4, 4, 50, 0, 22);  

  l1PFCandidatesFlat_TH2F     = fs->make<TH2F>( "l1PFCandidatesFlat",     "l1PFCandidatesFlat",     50, -4, 4, 50, 0, 22);
  l1PFNeutralFlat_TH2F        = fs->make<TH2F>( "l1PFNeutralsFlat",       "l1PFNeutralsFlat",       50, -4, 4, 50, 0, 22);
  l1PFChargedFlat_TH2F        = fs->make<TH2F>( "l1PFChargedFlat",        "l1PFChargedFlat",        50, -4, 4, 50, 0, 22);
  l1PFNeutralHadronsFlat_TH2F = fs->make<TH2F>( "l1PFNeutralHadronsFlat", "l1PFNeutralHadronsFlat", 50, -4, 4, 50, 0, 22);
  l1PFChargedHadronsFlat_TH2F = fs->make<TH2F>( "l1PFChargedHadronsFlat", "l1PFChargedHadronsFlat", 50, -4, 4, 50, 0, 22);
  l1PFMuonsFlat_TH2F          = fs->make<TH2F>( "l1PFMuonsFlat",          "l1PFMuonsFlat",          50, -4, 4, 50, 0, 22);
  l1PFElectronsFlat_TH2F      = fs->make<TH2F>( "l1PFElectronsFlat",      "l1PFElectronsFlat",      50, -4, 4, 50, 0, 22);
  l1PFPhotonsFlat_TH2F        = fs->make<TH2F>( "l1PFPhotonsFlat",        "l1PFPhotonsFlat",        50, -4, 4, 50, 0, 22); 
  l1PFTausFlat_TH2F           = fs->make<TH2F>( "l1PFTausFlat",           "l1PFTausFlat",           50, -4, 4, 50, 0, 22);
  l1PFJetsFlat_TH2F           = fs->make<TH2F>( "l1PFJetsFlat",           "l1PFJetsFlat",           50, -4, 4, 50, 0, 22);  

  l1PFCandidates10_TH2F     = fs->make<TH2F>( "l1PFCandidates10",     "l1PFCandidates10",     50, -4, 4, 50, 0, 22);
  l1PFNeutral10_TH2F        = fs->make<TH2F>( "l1PFNeutrals10",       "l1PFNeutrals10",       50, -4, 4, 50, 0, 22);
  l1PFCharged10_TH2F        = fs->make<TH2F>( "l1PFCharged10",        "l1PFCharged10",        50, -4, 4, 50, 0, 22);
  l1PFNeutralHadrons10_TH2F = fs->make<TH2F>( "l1PFNeutralHadrons10", "l1PFNeutralHadrons10", 50, -4, 4, 50, 0, 22);
  l1PFChargedHadrons10_TH2F = fs->make<TH2F>( "l1PFChargedHadrons10", "l1PFChargedHadrons10", 50, -4, 4, 50, 0, 22);
  l1PFMuons10_TH2F          = fs->make<TH2F>( "l1PFMuons10",          "l1PFMuons10",          50, -4, 4, 50, 0, 22);
  l1PFElectrons10_TH2F      = fs->make<TH2F>( "l1PFElectrons10",      "l1PFElectrons10",      50, -4, 4, 50, 0, 22);
  l1PFPhotons10_TH2F        = fs->make<TH2F>( "l1PFPhotons10",        "l1PFPhotons10",        50, -4, 4, 50, 0, 22); 
  l1PFTaus10_TH2F           = fs->make<TH2F>( "l1PFTaus10",           "l1PFTaus10",           50, -4, 4, 50, 0, 22);
  l1PFJets10_TH2F           = fs->make<TH2F>( "l1PFJets10",           "l1PFJets10",           50, -4, 4, 50, 0, 22);  

  l1PFCandidates20_TH2F     = fs->make<TH2F>( "l1PFCandidates20",     "l1PFCandidates20",     50, -4, 4, 50, 0, 22);
  l1PFNeutral20_TH2F        = fs->make<TH2F>( "l1PFNeutrals20",       "l1PFNeutrals20",       50, -4, 4, 50, 0, 22);
  l1PFCharged20_TH2F        = fs->make<TH2F>( "l1PFCharged20",        "l1PFCharged20",        50, -4, 4, 50, 0, 22);
  l1PFNeutralHadrons20_TH2F = fs->make<TH2F>( "l1PFNeutralHadrons20", "l1PFNeutralHadrons20", 50, -4, 4, 50, 0, 22);
  l1PFChargedHadrons20_TH2F = fs->make<TH2F>( "l1PFChargedHadrons20", "l1PFChargedHadrons20", 50, -4, 4, 50, 0, 22);
  l1PFMuons20_TH2F          = fs->make<TH2F>( "l1PFMuons20",          "l1PFMuons20",          50, -4, 4, 50, 0, 22);
  l1PFElectrons20_TH2F      = fs->make<TH2F>( "l1PFElectrons20",      "l1PFElectrons20",      50, -4, 4, 50, 0, 22);
  l1PFPhotons20_TH2F        = fs->make<TH2F>( "l1PFPhotons20",        "l1PFPhotons20",        50, -4, 4, 50, 0, 22); 
  l1PFTaus20_TH2F           = fs->make<TH2F>( "l1PFTaus20",           "l1PFTaus20",           50, -4, 4, 50, 0, 22);
  l1PFJets20_TH2F           = fs->make<TH2F>( "l1PFJets20",           "l1PFJets20",           50, -4, 4, 50, 0, 22);  
  */
  // Barrel rate trees: Jets, Muons, Electrons, Taus, Photons
  l1PFTauRate_Barrel      = fs->make<TH1F>( "l1PFTauRate_Barrel",      "l1PFTauRate_Barrel",      300,  0., 300. );
  l1DiPFTauRate_Barrel     = fs->make<TH1F>( "l1DiPFTauRate_Barrel",      "l1DiPFTauRate_Barrel",      300,  0., 300. );
  l1PFMuonRate_Barrel     = fs->make<TH1F>( "l1PFMuonRate_Barrel",     "l1PFMuonRate_Barrel",     300,  0., 300. );
  l1PFElectronRate_Barrel = fs->make<TH1F>( "l1PFElectronRate_Barrel", "l1PFElectronRate_Barrel", 300,  0., 300. );
  l1PFPhotonRate_Barrel   = fs->make<TH1F>( "l1PFPhotonRate_Barrel",   "l1PFPhotonRate_Barrel",   300,  0., 300. );

  l1PFTauRate_Iso90_Barrel      = fs->make<TH1F>( "l1PFTauRate_Iso90_Barrel",      "l1PFTauRate_Iso90_Barrel",      300,  0., 300. );
  l1DiPFTauRate_Iso90_Barrel      = fs->make<TH1F>( "l1DiPFTauRate_Iso90_Barrel",      "l1DiPFTauRate_Iso90_Barrel",      300,  0., 300. );
  l1PFMuonRate_Iso90_Barrel     = fs->make<TH1F>( "l1PFMuonRate_Iso90_Barrel",     "l1PFMuonRate_Iso90_Barrel",     300,  0., 300. );
  l1PFElectronRate_Iso90_Barrel = fs->make<TH1F>( "l1PFElectronRate_Iso90_Barrel", "l1PFElectronRate_Iso90_Barrel", 300,  0., 300. );
  l1PFPhotonRate_Iso90_Barrel   = fs->make<TH1F>( "l1PFPhotonRate_Iso90_Barrel",   "l1PFPhotonRate_Iso90_Barrel",   300,  0., 300. );

  l1PFTauRate_Iso90TIME_Barrel      = fs->make<TH1F>( "l1PFTauRate_Iso90TIME_Barrel",      "l1PFTauRate_Iso90TIME_Barrel",      300,  0., 300. );
  l1DiPFTauRate_Iso90TIME_Barrel      = fs->make<TH1F>( "l1DiPFTauRate_Iso90TIME_Barrel",      "l1DiPFTauRate_Iso90TIME_Barrel",      300,  0., 300. );
  l1PFMuonRate_Iso90TIME_Barrel     = fs->make<TH1F>( "l1PFMuonRate_Iso90TIME_Barrel",     "l1PFMuonRate_Iso90TIME_Barrel",     300,  0., 300. );
  l1PFElectronRate_Iso90TIME_Barrel = fs->make<TH1F>( "l1PFElectronRate_Iso90TIME_Barrel", "l1PFElectronRate_Iso90TIME_Barrel", 300,  0., 300. );
  l1PFPhotonRate_Iso90TIME_Barrel   = fs->make<TH1F>( "l1PFPhotonRate_Iso90TIME_Barrel",   "l1PFPhotonRate_Iso90TIME_Barrel",   300,  0., 300. );

  l1PFTauRate_Iso95_Barrel      = fs->make<TH1F>( "l1PFTauRate_Iso95_Barrel",      "l1PFTauRate_Iso95_Barrel",      300,  0., 300. );
  l1DiPFTauRate_Iso95_Barrel      = fs->make<TH1F>( "l1DiPFTauRate_Iso95_Barrel",      "l1DiPFTauRate_Iso95_Barrel",      300,  0., 300. );
  l1PFMuonRate_Iso95_Barrel     = fs->make<TH1F>( "l1PFMuonRate_Iso95_Barrel",     "l1PFMuonRate_Iso95_Barrel",     300,  0., 300. );
  l1PFElectronRate_Iso95_Barrel = fs->make<TH1F>( "l1PFElectronRate_Iso95_Barrel", "l1PFElectronRate_Iso95_Barrel", 300,  0., 300. );
  l1PFPhotonRate_Iso95_Barrel   = fs->make<TH1F>( "l1PFPhotonRate_Iso95_Barrel",   "l1PFPhotonRate_Iso95_Barrel",   300,  0., 300. );

  l1PFTauRate_Iso95TIME_Barrel      = fs->make<TH1F>( "l1PFTauRate_Iso95TIME_Barrel",      "l1PFTauRate_Iso95TIME_Barrel",      300,  0., 300. );
  l1DiPFTauRate_Iso95TIME_Barrel      = fs->make<TH1F>( "l1DiPFTauRate_Iso95TIME_Barrel",      "l1PDiFTauRate_Iso95TIME_Barrel",      300,  0., 300. );
  l1PFMuonRate_Iso95TIME_Barrel     = fs->make<TH1F>( "l1PFMuonRate_Iso95TIME_Barrel",     "l1PFMuonRate_Iso95TIME_Barrel",     300,  0., 300. );
  l1PFElectronRate_Iso95TIME_Barrel = fs->make<TH1F>( "l1PFElectronRate_Iso95TIME_Barrel", "l1PFElectronRate_Iso95TIME_Barrel", 300,  0., 300. );
  l1PFPhotonRate_Iso95TIME_Barrel   = fs->make<TH1F>( "l1PFPhotonRate_Iso95TIME_Barrel",   "l1PFPhotonRate_Iso95TIME_Barrel",   300,  0., 300. );

  l1PFJetRate_Barrel      = fs->make<TH1F>( "l1PFJetRate_Barrel",      "l1PFJetRate_Barrel",      2000,  0., 2000. );
  l1PFJetRate_time_Barrel      = fs->make<TH1F>( "l1PFJetRate_time_Barrel",      "l1PFJetRate_time_Barrel",      2000,  0., 2000. );
  l1PFJet2nsRate_Barrel   = fs->make<TH1F>( "l1PFJet2nsRate_Barrel",   "l1PFJet2nsRate_Barrel",   2000,  0., 2000. );
  l1PFJet3nsRate_Barrel   = fs->make<TH1F>( "l1PFJet3nsRate_Barrel",   "l1PFJet3nsRate_Barrel",   2000,  0., 2000. );
  l1PFJet4nsRate_Barrel   = fs->make<TH1F>( "l1PFJet4nsRate_Barrel",   "l1PFJet4nsRate_Barrel",   2000,  0., 2000. );
  l1PFJet5nsRate_Barrel   = fs->make<TH1F>( "l1PFJet5nsRate_Barrel",   "l1PFJet5nsRate_Barrel",   2000,  0., 2000. );
  l1PFJet6nsRate_Barrel   = fs->make<TH1F>( "l1PFJet6nsRate_Barrel",   "l1PFJet6nsRate_Barrel",   2000,  0., 2000. );
  l1PFJet7nsRate_Barrel   = fs->make<TH1F>( "l1PFJet7nsRate_Barrel",   "l1PFJet7nsRate_Barrel",   2000,  0., 2000. );

  l1PFJetRate_calibrated_Barrel      = fs->make<TH1F>( "l1PFJetRate_calibrated_Barrel",      "l1PFJetRate_calibrated_Barrel",      2000,  0., 2000. );
  l1PFJetRate_time_calibrated_Barrel      = fs->make<TH1F>( "l1PFJetRate_time_calibrated_Barrel",      "l1PFJetRate_time_calibrated_Barrel",      2000,  0., 2000. );
  l1PFJet2nsRate_calibrated_Barrel   = fs->make<TH1F>( "l1PFJet2nsRate_calibrated_Barrel",   "l1PFJet2nsRate_calibrated_Barrel",   2000,  0., 2000. );
  l1PFJet3nsRate_calibrated_Barrel   = fs->make<TH1F>( "l1PFJet3nsRate_calibrated_Barrel",   "l1PFJet3nsRate_calibrated_Barrel",   2000,  0., 2000. );
  l1PFJet4nsRate_calibrated_Barrel   = fs->make<TH1F>( "l1PFJet4nsRate_calibrated_Barrel",   "l1PFJet4nsRate_calibrated_Barrel",   2000,  0., 2000. );
  l1PFJet5nsRate_calibrated_Barrel   = fs->make<TH1F>( "l1PFJet5nsRate_calibrated_Barrel",   "l1PFJet5nsRate_calibrated_Barrel",   2000,  0., 2000. );
  l1PFJet6nsRate_calibrated_Barrel   = fs->make<TH1F>( "l1PFJet6nsRate_calibrated_Barrel",   "l1PFJet6nsRate_calibrated_Barrel",   2000,  0., 2000. );
  l1PFJet7nsRate_calibrated_Barrel   = fs->make<TH1F>( "l1PFJet7nsRate_calibrated_Barrel",   "l1PFJet7nsRate_calibrated_Barrel",   2000,  0., 2000. );

  // Total rate trees: Jets, Muons, Electrons, Taus, Photons
  l1PFTauRate      = fs->make<TH1F>( "l1PFTauRate",      "l1PFTauRate",      300,  0., 300. );
  l1DiPFTauRate      = fs->make<TH1F>( "l1DiPFTauRate",      "l1DiPFTauRate",      300,  0., 300. );
  l1PFMuonRate     = fs->make<TH1F>( "l1PFMuonRate",     "l1PFMuonRate",     300,  0., 300. );
  l1PFElectronRate = fs->make<TH1F>( "l1PFElectronRate", "l1PFElectronRate", 300,  0., 300. );
  l1PFPhotonRate   = fs->make<TH1F>( "l1PFPhotonRate",   "l1PFPhotonRate",   300,  0., 300. );

  l1PFTauRate_Iso90      = fs->make<TH1F>( "l1PFTauRate_Iso90",      "l1PFTauRate_Iso90",      300,  0., 300. );
  l1DiPFTauRate_Iso90      = fs->make<TH1F>( "l1DiPFTauRate_Iso90",      "l1DiPFTauRate_Iso90",      300,  0., 300. );
  l1PFMuonRate_Iso90     = fs->make<TH1F>( "l1PFMuonRate_Iso90",     "l1PFMuonRate_Iso90",     300,  0., 300. );
  l1PFElectronRate_Iso90 = fs->make<TH1F>( "l1PFElectronRate_Iso90", "l1PFElectronRate_Iso90", 300,  0., 300. );
  l1PFPhotonRate_Iso90   = fs->make<TH1F>( "l1PFPhotonRate_Iso90",   "l1PFPhotonRate_Iso90",   300,  0., 300. );

  l1PFTauRate_Iso90TIME      = fs->make<TH1F>( "l1PFTauRate_Iso90TIME",      "l1PFTauRate_Iso90TIME",      300,  0., 300. );
  l1DiPFTauRate_Iso90TIME      = fs->make<TH1F>( "l1DiPFTauRate_Iso90TIME",      "l1DiPFTauRate_Iso90TIME",      300,  0., 300. );
  l1PFMuonRate_Iso90TIME     = fs->make<TH1F>( "l1PFMuonRate_Iso90TIME",     "l1PFMuonRate_Iso90TIME",     300,  0., 300. );
  l1PFElectronRate_Iso90TIME = fs->make<TH1F>( "l1PFElectronRate_Iso90TIME", "l1PFElectronRate_Iso90TIME", 300,  0., 300. );
  l1PFPhotonRate_Iso90TIME   = fs->make<TH1F>( "l1PFPhotonRate_Iso90TIME",   "l1PFPhotonRate_Iso90TIME",   300,  0., 300. );

  l1PFTauRate_Iso95      = fs->make<TH1F>( "l1PFTauRate_Iso95",      "l1PFTauRate_Iso95",      300,  0., 300. );
  l1DiPFTauRate_Iso95      = fs->make<TH1F>( "l1DiPFTauRate_Iso95",      "l1DiPFTauRate_Iso95",      300,  0., 300. );
  l1PFMuonRate_Iso95     = fs->make<TH1F>( "l1PFMuonRate_Iso95",     "l1PFMuonRate_Iso95",     300,  0., 300. );
  l1PFElectronRate_Iso95 = fs->make<TH1F>( "l1PFElectronRate_Iso95", "l1PFElectronRate_Iso95", 300,  0., 300. );
  l1PFPhotonRate_Iso95   = fs->make<TH1F>( "l1PFPhotonRate_Iso95",   "l1PFPhotonRate_Iso95",   300,  0., 300. );

  l1PFTauRate_Iso95TIME      = fs->make<TH1F>( "l1PFTauRate_Iso95TIME",      "l1PFTauRate_Iso95TIME",      300,  0., 300. );
  l1DiPFTauRate_Iso95TIME      = fs->make<TH1F>( "l1DiPFTauRate_Iso95TIME",      "l1DiPFTauRate_Iso95TIME",      300,  0., 300. );
  l1PFMuonRate_Iso95TIME     = fs->make<TH1F>( "l1PFMuonRate_Iso95TIME",     "l1PFMuonRate_Iso95TIME",     300,  0., 300. );
  l1PFElectronRate_Iso95TIME = fs->make<TH1F>( "l1PFElectronRate_Iso95TIME", "l1PFElectronRate_Iso95TIME", 300,  0., 300. );
  l1PFPhotonRate_Iso95TIME   = fs->make<TH1F>( "l1PFPhotonRate_Iso95TIME",   "l1PFPhotonRate_Iso95TIME",   300,  0., 300. );

  l1PFJetRate      = fs->make<TH1F>( "l1PFJetRate",      "l1PFJetRate",      2000,  0., 2000. );
  l1PFJetRate_time      = fs->make<TH1F>( "l1PFJetRate_time",      "l1PFJetRate_time",      2000,  0., 2000. );
  l1PFJet2nsRate   = fs->make<TH1F>( "l1PFJet2nsRate",   "l1PFJet2nsRate",   2000,  0., 2000. );
  l1PFJet3nsRate   = fs->make<TH1F>( "l1PFJet3nsRate",   "l1PFJet3nsRate",   2000,  0., 2000. );
  l1PFJet4nsRate   = fs->make<TH1F>( "l1PFJet4nsRate",   "l1PFJet4nsRate",   2000,  0., 2000. );
  l1PFJet5nsRate   = fs->make<TH1F>( "l1PFJet5nsRate",   "l1PFJet5nsRate",   2000,  0., 2000. );
  l1PFJet6nsRate   = fs->make<TH1F>( "l1PFJet6nsRate",   "l1PFJet6nsRate",   2000,  0., 2000. );
  l1PFJet7nsRate   = fs->make<TH1F>( "l1PFJet7nsRate",   "l1PFJet7nsRate",   2000,  0., 2000. );


}

//destructor
L1MTDPFAnalyzer::~L1MTDPFAnalyzer()
{
}

//deltaR, deltaZ, time_Cut 
void 
L1MTDPFAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  run   = iEvent.id().run();
  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();
  nEvents->Fill(1);
  //edm::Handle<edm::ValueMap<float> > timingValues;
  //iEvent.getByToken(timingValuesToken_,timingValues);

  Handle<L1TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken_,L1VertexHandle);

  //make gen vectors
  std::vector<reco::GenParticle> genElectrons, genPhotons, genMuons, genMet;
  std::vector<genVisTau> genTaus;
  // make reco handle vectors
  Handle< std::vector<pat::Electron> > recoElectrons;
  Handle< std::vector<pat::Photon> >   recoPhotons;
  Handle< std::vector<pat::Muon> >     recoMuons;
  Handle< std::vector<pat::Tau> >      recoTaus;
  Handle< std::vector<pat::MET> >      recoMet;

  // make l1 handle vectors
  Handle< l1t::PFCandidateCollection > l1PFCandidates;

  // make the L1 vectors electrons, photons, muons, taus, met 
  PFCandidateCollection l1PFElectrons, l1PFPhotons, l1PFMuons, l1PFTaus_Iso; 
  std::vector<int> l1PFElectrons_indicies, l1PFPhotons_indicies, l1PFMuons_indicies;

  Handle<L1PFTauCollection> l1PFTaus;

  //put l1 met collection here
  Handle<L1TkEtMissParticleCollection> l1Met;
  Handle<L1TkEtMissParticleCollection> l1Met_time;

  if(!iEvent.getByToken( L1PFCandsToken_, l1PFCandidates))    std::cout<<"No L1PF Cands Found!!"     <<std::endl;
  if(!iEvent.getByToken( L1PFTausToken_,  l1PFTaus))          std::cout<<"No L1PF Taus Found!!"      <<std::endl;
  if(!iEvent.getByToken( L1PFMETToken_,   l1Met))             std::cout<<"No L1PF Mets Found!!"      <<std::endl;

  //std::vector<l1t::PFCandidate*> l1PFCandidates_sort;
  //for(auto &cand : *l1PFCandidates){
  //l1PFCandidates_sort.push_back(*cand);
  //}

  //Get the reco objects
  if(!iEvent.getByToken( recoElectronsToken_, recoElectrons)) std::cout<<"No Reco ELECTRONS Found!!!"<<std::endl;
  if(!iEvent.getByToken( recoPhotonsToken_,   recoPhotons))   std::cout<<"No Reco PHOTONS Found!!!"  <<std::endl;
  if(!iEvent.getByToken( recoMuonsToken_,     recoMuons))     std::cout<<"No Reco MUONS Found!!!"    <<std::endl;
  if(!iEvent.getByToken( recoTausToken_,      recoTaus))      std::cout<<"No Reco TAUS Found!!!"     <<std::endl;
  if(!iEvent.getByToken( recoMetToken_,       recoMet))       std::cout<<"No Reco MET Found!!!"      <<std::endl;

  // sort by pt would be useful for rates, but doing this object by object for now
  //std::sort(l1PFCandidates->begin(), l1PFCandidates->end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   

  //get gen objects
  Handle<std::vector<reco::GenParticle> > genParticles;
  if(!iEvent.getByToken(genParticlesToken_,genParticles)) std::cout<<"No Gen Particles Found!!!"<<std::endl;
  
  Handle<std::vector<reco::GenJet> > genJets;
  if(!iEvent.getByToken(genJetsToken_,genJets)) std::cout<<"No Gen Jets Found!!!"<<std::endl;
  else
    std::cout<<"NGen Jets found: "<<genJets->size()<<std::endl;

  //Fill each gen particle vector with the corresponding gen particles
  //electrons pdgid == 11
  fillGenParticleVector(genParticles, genElectrons, 11);

  //photons   pdgid == 22
  fillGenParticleVector(genParticles, genPhotons,   22);

  //muons     pdgid == 13 
  fillGenParticleVector(genParticles, genMuons,     13);

  //taus      pdgid == 15 // but this is a special lepton, using visible decay products
  fillGenTauVector(genParticles, genTaus);

  //Fill Each L1 PF Cand vector with the L1 PFCands Electrons, Photons and Muons
  fillL1PFCandCollections(l1PFCandidates, l1PFElectrons, l1PFElectrons_indicies, l1PFPhotons, l1PFPhotons_indicies, l1PFMuons, l1PFMuons_indicies);

  //Create the L1 jets
  std::vector<jet> L1Jets_time;
  findJets(l1PFCandidates, 0.5, L1Jets_time);
  std::vector<jet> L1Jets;
  findJets(l1PFCandidates, -1, L1Jets);

  //Create the L1 LL jets
  std::vector<jet> L1LLJets;
  findLLJets(l1PFCandidates, 1, 2, L1LLJets); //1ns window, 2ns offset from central bunch crossing
 
  std::vector<object> tL1Electrons, tL1Photons, tL1Muons, tL1Taus;

  fillL1TauTreeVariables( l1PFTaus, l1PFCandidates, genTaus, L1VertexHandle);
  //std::cout<<"electron tree"<<std::endl;
  //Fill the Electron Tree
  //for each reco electron 
  //              find a matched gen
  //              find a matched L1
  size_t idx = 0;
  for( std::vector<pat::Electron>::const_iterator cand  = recoElectrons->begin(); 
                                                  cand != recoElectrons->end(); 
                                                  ++ cand, ++ idx) {
    double isoSum = 0;
    double isoSumTime = 0;    

    //Zero the tree
    zeroElectronTreeVariables();
    eleRecoPt  = cand->pt();
    eleRecoEta = cand->eta();
    eleRecoPhi = cand->phi();

    int idx = 0;
    //Get matched L1 Electron
    for(auto l1Cand : l1PFElectrons){
      int object_idx = l1PFElectrons_indicies.at(idx);
      idx++;
      if(reco::deltaR(cand->eta(),cand->phi(),l1Cand.eta(),l1Cand.phi())<0.1){
	//calculate iso Sum
	calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	eleL1Iso      = isoSum;
	eleL1Iso_time = isoSumTime;
	//std::cout<<"matched electron found!! eleL1Iso: "<<eleL1Iso<<" eleL1Iso_time: "<<eleL1Iso_time<<std::endl;
	eleL1Pt       = l1Cand.pt();
	eleL1Eta      = l1Cand.eta();
	eleL1Phi      = l1Cand.phi();
	eleL1Time     = l1Cand.time();
	object temp;
	temp.pt = l1Cand.pt(); temp.eta = l1Cand.eta(); temp.phi = l1Cand.phi(); temp.time = l1Cand.time(); temp.isoSum = isoSum; temp.isoSumTime = isoSumTime;
	break;
      }
    }
    //Get matched gen Electron
    for(auto genCand : genElectrons){
      //matchGen
      if(reco::deltaR(cand->eta(), cand->phi(), genCand.eta(), genCand.phi())<0.1){
	eleGenPt      = genCand.pt();
	eleGenEta     = genCand.eta();
	eleGenPhi     = genCand.phi();
	break;
      }
    }
    electronTree->Fill();
  }

  /// fill l1 elecrtron collection
  idx = 0;
  for(auto l1Cand : l1PFElectrons){
    int object_idx = l1PFElectrons_indicies.at(idx);
    idx++;
    //calculate iso Sum
    double isoSum = 0;
    double isoSumTime = 0;    

    calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
    object temp;
    temp.pt = l1Cand.pt(); temp.eta = l1Cand.eta(); temp.phi = l1Cand.phi(); temp.time = l1Cand.time(); temp.isoSum = isoSum; temp.isoSumTime = isoSumTime;
    tL1Electrons.push_back(temp);
  }
  std::sort(tL1Electrons.begin(), tL1Electrons.end(), [](object i,object j){return(i.pt > j.pt);});   

  //Fill the Photons
  //std::cout<<"Checking Photons"<<std::endl;
  idx = 0;
  for( std::vector<pat::Photon>::const_iterator cand  = recoPhotons->begin(); 
                                                  cand != recoPhotons->end(); 
                                                  ++ cand) {
    double isoSum = 0;
    double isoSumTime = 0;    

    if(cand->photonID("cutBasedPhotonID-Spring16-V2p2-loose")==0)
      continue;

    //std::cout<<"passed photonID"<<std::endl;
    //Zero the tree
    zeroPhotonTreeVariables();
    gammaRecoPt  = cand->pt();
    gammaRecoEta = cand->eta();
    gammaRecoPhi = cand->phi();

    int idx = 0;
    //Get matched L1 Photon
    for(auto l1Cand : l1PFPhotons){
      int object_idx = l1PFPhotons_indicies.at(idx);
      idx++;
      if(reco::deltaR(cand->eta(),cand->phi(),l1Cand.eta(),l1Cand.phi())<0.2){
	//calculate iso Sum
	calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	gammaL1Iso      = isoSum;
	gammaL1Iso_time = isoSumTime;
	//std::cout<<"matched photon found!! gammaL1Iso: "<<gammaL1Iso<<" gammaL1Iso_time: "<<gammaL1Iso_time<<std::endl;
	gammaL1Pt       = l1Cand.pt();
	gammaL1Eta      = l1Cand.eta();
	gammaL1Phi      = l1Cand.phi();
	gammaL1Time     = l1Cand.time();
	break;
      }
    }
    //Get matched gen Photon
    for(auto genCand : genPhotons){
      //matchGen
      if(reco::deltaR(cand->eta(), cand->phi(), genCand.eta(), genCand.phi())<0.1){
	gammaGenPt      = genCand.pt();
	gammaGenEta     = genCand.eta();
	gammaGenPhi     = genCand.phi();
	break;
      }
    }
    photonTree->Fill();
  }

  idx = 0;
  /// fill l1 photon collection
  for(auto l1Cand : l1PFPhotons){
    int object_idx = l1PFPhotons_indicies.at(idx);
    idx++;
    //calculate iso Sum
    double isoSum = 0;
    double isoSumTime = 0;    
    calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
    object temp;
    temp.pt = l1Cand.pt(); temp.eta = l1Cand.eta(); temp.phi = l1Cand.phi(); temp.time = l1Cand.time(); temp.isoSum = isoSum; temp.isoSumTime = isoSumTime;
  }

  std::sort(tL1Photons.begin(), tL1Photons.end(), [](object i,object j){return(i.pt > j.pt);});   

  //std::cout<<"Checking Muons"<<std::endl;
    //Fill the Muons
  for( std::vector<pat::Muon>::const_iterator cand  = recoMuons->begin(); 
                                                  cand != recoMuons->end(); 
                                                  ++ cand) {
    double isoSum = 0;
    double isoSumTime = 0;    

    //Zero the tree
    zeroMuonTreeVariables();
    muRecoPt  = cand->pt();
    muRecoEta = cand->eta();
    muRecoPhi = cand->phi();

    int idx = 0;
    //Get matched L1 Muon
    for(auto l1Cand : l1PFMuons){
      int object_idx = l1PFMuons_indicies.at(idx);
      idx++;
      if(reco::deltaR(cand->eta(),cand->phi(),l1Cand.eta(),l1Cand.phi())<0.1){
	//calculate iso Sum
	calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	muL1Iso      = isoSum;
	muL1Iso_time = isoSumTime;
	//std::cout<<"matched muon found!! muL1Iso: "<<muL1Iso<<" muL1Iso_time: "<<muL1Iso_time<<std::endl;
	muL1Pt       = l1Cand.pt();
	muL1Eta      = l1Cand.eta();
	muL1Phi      = l1Cand.phi();
	muL1Time     = l1Cand.time();
	break;
      }
    }
    //Get matched gen Muon
    for(auto genCand : genMuons){
      //matchGen
      if(reco::deltaR(cand->eta(), cand->phi(), genCand.eta(), genCand.phi())<0.1){
	muGenPt      = genCand.pt();
	muGenEta     = genCand.eta();
	muGenPhi     = genCand.phi();
	break;
      }
    }
    muonTree->Fill();
  }
  /// fill l1 muon n
  idx = 0;
  for(auto l1Cand : l1PFMuons){
    int object_idx = l1PFMuons_indicies.at(idx);
    idx++;
    //calculate iso Sum
    double isoSum = 0;
    double isoSumTime = 0;    

    calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
    object temp;
    temp.pt = l1Cand.pt(); temp.eta = l1Cand.eta(); temp.phi = l1Cand.phi(); temp.time = l1Cand.time(); temp.isoSum = isoSum; temp.isoSumTime = isoSumTime;
    tL1Muons.push_back(temp);
  }

  std::sort(tL1Muons.begin(), tL1Muons.end(), [](object i,object j){return(i.pt > j.pt);});   
  //std::cout<<"Checking Taus"<<std::endl;


  //for(auto &l1Cand : *l1PFTaus) 
  //l1PFTaus_TH2F->Fill( l1Cand.eta() , l1Cand.time() );

  //Fill the Taus
  for( std::vector<pat::Tau>::const_iterator cand  = recoTaus->begin(); 
                                                  cand != recoTaus->end(); 
                                                  ++ cand) {
    double isoSum = 0;
    double isoSumTime = 0;    
    
    //Zero the tree
    zeroTauTreeVariables();
    tauRecoPt  = cand->pt();
    tauRecoEta = cand->eta();
    tauRecoPhi = cand->phi();

    int idx = 0;
    //Get matched L1 Tau
    for(auto &l1Cand : *l1PFTaus){
      idx++;
      if(l1Cand.pt()<1)
	continue;

      if(reco::deltaR(cand->eta(),cand->phi(),l1Cand.eta(),l1Cand.phi())<0.1){
	//calculate iso Sum
	calculate_charged_iso_sum_tau(l1Cand, l1PFCandidates,  time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	tauL1Iso      = isoSum;
	tauL1Iso_time = isoSumTime;
	//std::cout<<"matched tau found!! tauL1Iso: "<<tauL1Iso<<" tauL1Iso_time: "<<tauL1Iso_time<<std::endl;
	tauL1Pt       = l1Cand.pt();
	tauL1Eta      = l1Cand.eta();
	tauL1Phi      = l1Cand.phi();
	tauL1Time     = l1Cand.time();
	break;
      }
    }
    //Get matched gen Tau
    for(auto genCand : genTaus){
      //matchGen
      if(reco::deltaR(cand->eta(), cand->phi(), genCand.p4.eta(), genCand.p4.phi())<0.1){
	tauGenPt      = genCand.p4.pt();
	tauGenEta     = genCand.p4.eta();
	tauGenPhi     = genCand.p4.phi();
	break;
      }
    }
    tauTree->Fill();
  }
  //std::cout<<"Here5"<<std::endl;
  /// fill l1 photon n
  for(auto &l1Cand : *l1PFTaus){
    //calculate iso Sum
    double isoSum = 0;
    double isoSumTime = 0;    

    calculate_charged_iso_sum_tau(l1Cand, l1PFCandidates,  time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
    object temp;
    temp.pt = l1Cand.pt(); temp.eta = l1Cand.eta(); temp.phi = l1Cand.phi(); temp.time = l1Cand.time(); temp.isoSum = isoSum; temp.isoSumTime = isoSumTime;
    tL1Taus.push_back(temp);
  }

  std::sort(tL1Taus.begin(), tL1Taus.end(), [](object i,object j){return(i.pt > j.pt);});   

  //fill the jets
  for(auto &genCand : *genJets){
    //matchGen
    jetGenPt  = genCand.pt();
    jetGenEta = genCand.eta();
    jetGenPhi = genCand.phi();
    jetL1Pt   = 0;
    jetL1Eta  = 0;
    jetL1Phi  = 0;
    jetL1Time = 0;

    for(auto cand : L1Jets_time){
      if(reco::deltaR(cand.eta, cand.phi, genCand.eta(), genCand.phi())<0.2){
	jetL1Pt_time   = cand.pt;
	jetL1Eta_time  = cand.eta;
	jetL1Phi_time  = cand.phi;
	jetL1Time = cand.time;
	break;
      }
    }

    for(auto cand : L1Jets){
      if(reco::deltaR(cand.eta, cand.phi, genCand.eta(), genCand.phi())<0.2){
	jetL1Pt   = cand.pt;
	jetL1Eta  = cand.eta;
	jetL1Phi  = cand.phi;
	break;
      }
    }

    jetTree->Fill();
  }

  //fill the jets
  for(auto &genCand : *genJets){
    //matchGen
    jetGenPt  = genCand.pt();
    jetGenEta = genCand.eta();
    jetGenPhi = genCand.phi();
    jetL1Pt   = 0;
    jetL1Eta  = 0;
    jetL1Phi  = 0;
    jetL1Time = 0;
    for(auto cand : L1LLJets){
      if(reco::deltaR(cand.eta, cand.phi, genCand.eta(), genCand.phi())<0.2){
	jetL1Pt   = cand.pt;
	jetL1Eta  = cand.eta;
	jetL1Phi  = cand.phi;
	jetL1Time = cand.time;
	break;
      }
    }
    jetLLTree->Fill();
  }

  fillRates(tL1Electrons, tL1Photons, tL1Muons, tL1Taus, L1Jets, L1Jets_time);
    //Fill the MET
  /*
  if(pfMet->size()>0){
    metRecoEt    = pfMet->at(0).et();
    metRecoPhi   = pfMet->at(0).phi();
    metGenEt     = ;
    metGenPhi    = ;
    metL1Et      = ;
    metL1Et_time = ;
    }*/
}

//void L1MTDPFAnalyzer::fillValidationHistograms(Handle< l1t::PFCandidateCollection > l1PFCandidates){
//  for( const auto & l1pf: *l1PFCandidates){
//    
//  }
//
//};

void L1MTDPFAnalyzer::calculate_charged_iso_sum(l1t::PFCandidate object, int object_index, Handle<l1t::PFCandidateCollection> pfCandidates, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime){
    double iso_sum = 0;
    double iso_sum_time = 0;
    int idx = 0;
    int kdx = 0;
    for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
                                               pfCand != pfCandidates->end(); 
	                                       ++pfCand, ++idx){    

      //if(idx == object_index)
      //continue;
      
      //first check if it is a charged candidate
      if(pfCand->id() == l1t::PFCandidate::Muon || pfCand->id() == l1t::PFCandidate::Electron || pfCand->id() == l1t::PFCandidate::ChargedHadron){
      
	// check deltaR match
	if(reco::deltaR(object.eta(), object.phi(), pfCand->eta(), pfCand->phi()) < deltaR){

	  if(object.pfTrack().isNull()||pfCand->pfTrack().isNull())
	    continue;
	  
	  if(fabs(object.pfTrack()->track()->getPOCA().z()-pfCand->pfTrack()->track()->getPOCA().z()) < deltaZ){
	    
	    iso_sum += pfCand->pt();
	    float objectTime = 0;
	    float pfCandTime = 0;
	    
	    if(kdx<20){
	      //std::cout<<"objectTime: "<<objectTime<<"pfCandTime: "<<pfCandTime<< " Delta Time: "<<fabs(objectTime-pfCandTime) <<std::endl;
	      kdx++;
	    } 
	    objectTime = object.time();
	    pfCandTime = pfCand->time();
	    
	    //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > objectTime(pfCand->pfTrack(), track_index);
	  //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > pfCandTime(object.pfTrack(), track_index);
	    
	  //std::cout<<"objectTime: "<<objectTime<<" pfCandTime:"<<pfCandTime<<std::endl;
	  // then check time match

	    //0 represents an invalid value
	    //if object time is 0 then always add it to the isolation cone. 

	    if(objectTime == 0 || pfCandTime == 0 || fabs(objectTime - pfCandTime) < time_cut){
	      iso_sum_time += pfCand->pt();
	    }
	  }
	}
      }
    }
    isoSum = iso_sum;
    isoSumTime = iso_sum_time;
}

void L1MTDPFAnalyzer::calculate_charged_iso_sum_tau(l1t::L1PFTau object, Handle<l1t::PFCandidateCollection>  pfCandidates, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime){
    double iso_sum = 0;
    double iso_sum_time = 0;
    int idx = 0;
    int kdx = 0;
    for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
                                               pfCand != pfCandidates->end(); 
	                                       ++pfCand, ++idx){    
      
      //if(reco::deltaR(object.eta(), object.phi(), pfCand->eta(), pfCand->phi()) > 0.1){
	
      //first check if it is a charged candidate
      if(pfCand->id() == l1t::PFCandidate::Muon || pfCand->id() == l1t::PFCandidate::Electron || pfCand->id() == l1t::PFCandidate::ChargedHadron){
	
	// check deltaR match
	if(reco::deltaR(object.eta(), object.phi(), pfCand->eta(), pfCand->phi()) < deltaR){
	  
	  if(object.pfRef().at(0).pfTrack().isNull()||pfCand->pfTrack().isNull())
	    continue;
	  
	  if(fabs(object.pfRef().at(0).pfTrack()->track()->getPOCA().z()-pfCand->pfTrack()->track()->getPOCA().z()) < deltaZ){
	    
	    iso_sum += pfCand->pt();
	    float objectTime = 0;
	    float pfCandTime = 0;
	    
	    objectTime = object.time();
	    pfCandTime = pfCand->time();

	    if(kdx<20){
	      //std::cout<<"objectTime: "<<objectTime<<"pfCandTime: "<<pfCandTime<< " Delta Time: "<<fabs(objectTime-pfCandTime) <<std::endl;
	      kdx++;
	    } 
	    //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > objectTime(pfCand->pfTrack(), track_index);
	    //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > pfCandTime(object.pfTrack(), track_index);
	    
	    //std::cout<<"objectTime: "<<objectTime<<" pfCandTime:"<<pfCandTime<<std::endl;
	    // then check time match
	    
	    //0 represents an invalid value
	    //if object time is 0 then always add it to the isolation cone. 
	    
	    if(objectTime == 0 || pfCandTime == 0 || fabs(objectTime - pfCandTime) < time_cut){
	      iso_sum_time += pfCand->pt();
	    }
	  }
	}
      }
    }
    isoSum = iso_sum;
    isoSumTime = iso_sum_time;
}

void L1MTDPFAnalyzer::fillGenParticleVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<reco::GenParticle> &selectedGenParticles,int pdgID){
  for( std::vector<reco::GenParticle>::const_iterator cand  = genParticles->begin(); 
                                                      cand != genParticles->end(); 
                                                      ++ cand) {
    if(abs(cand->pdgId()) == pdgID){
      selectedGenParticles.push_back(*cand);
    }
  }
}

void L1MTDPFAnalyzer::fillGenTauVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<genVisTau> &selectedGenParticles){
  const std::vector<reco::GenParticle> genCands = *genParticles;
  for( std::vector<reco::GenParticle>::const_iterator cand  = genParticles->begin(); 
                                                      cand != genParticles->end(); 
                                                      ++ cand) {

    if(abs(cand->pdgId())==15){      
      const reco::GenParticle tau= *cand;
      reco::Candidate::LorentzVector visGenTau = getVisMomentum(&tau, &genCands);
      genVisTau Temp;
      int decayMode  = GetDecayMode(&tau);
      Temp.p4        = visGenTau;
      Temp.decayMode = decayMode;    
      selectedGenParticles.push_back(Temp);
    }
  }
}


void L1MTDPFAnalyzer::fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates, PFCandidateCollection &pfElectrons, std::vector<int> & l1PFElectrons_indicies, PFCandidateCollection &pfPhotons, std::vector<int> & l1PFPhotons_indicies, PFCandidateCollection &pfMuons, std::vector<int> & l1PFMuons_indicies){
  int index = 0;
  //std::cout<<"fill l1 pf cand collection"<<std::endl;
  for( PFCandidateCollection::const_iterator l1PFCand  = pfCandidates->begin();
                                             l1PFCand != pfCandidates->end(); 
       ++l1PFCand, ++index){    
    //fill the histograms
    //make pt cut histos as well
    float time_eta_calibrated = 0;

    if(fabs(l1PFCand->eta())<1.67)
      time_eta_calibrated = l1PFCand->time() - 15/6*l1PFCand->eta()*l1PFCand->eta() - 3.5;
    else
      time_eta_calibrated = l1PFCand->time() - 10;

    /*Fill the Long Lived Rates*/
    fillLLRates(*l1PFCand,time_eta_calibrated);

    //l1PFCandidates_TH2F->Fill(l1PFCand->eta() , l1PFCand->time() );
    
    //if(l1PFCand->pt()>10)
    //l1PFCandidates10_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );

    //if(l1PFCand->pt()>20)
    //l1PFCandidates20_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
    /*
    if(l1PFCand->id() == l1t::PFCandidate::ChargedHadron){
      l1PFChargedHadrons_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      l1PFCharged_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>10)
	l1PFCharged10_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>20)
	l1PFCharged20_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
    }

    if(l1PFCand->id() == l1t::PFCandidate::NeutralHadron){
      l1PFNeutralHadrons_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      l1PFNeutral_TH2F->Fill(          l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>10)
	l1PFNeutral10_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>20)
	l1PFNeutral20_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
	}
    */

    //Electrons
    if(l1PFCand->id() == l1t::PFCandidate::Electron){
      if(l1PFCand->pfTrack().isNull())
	continue;
      pfElectrons.push_back(*l1PFCand);
      l1PFElectrons_indicies.push_back(index);
      /*
      //Fill the histograms
      l1PFElectrons_TH2F->Fill( l1PFCand->eta() , l1PFCand->time() );
      l1PFCharged_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );

      if(l1PFCand->pt()>10)
	l1PFElectrons10_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>20)
	l1PFElectrons20_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      */
    }

    //Photons
    if(l1PFCand->id() == l1t::PFCandidate::Photon){
      pfPhotons.push_back(*l1PFCand);
      l1PFPhotons_indicies.push_back(index);
      /*
      //Fill the histograms
      l1PFPhotons_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      l1PFNeutral_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );

      if(l1PFCand->pt()>10)
	l1PFPhotons10_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>20)
	l1PFPhotons20_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      */
    }

    //Muons
    if(l1PFCand->id() == l1t::PFCandidate::Muon){

      if(l1PFCand->pfTrack().isNull())
	continue;

      pfMuons.push_back(*l1PFCand);
      l1PFMuons_indicies.push_back(index);
      /*
      //Fill the histograms
      l1PFMuons_TH2F->Fill(     l1PFCand->eta() , l1PFCand->time() );
      l1PFCharged_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );

      if(l1PFCand->pt()>10)
	l1PFMuons10_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      if(l1PFCand->pt()>20)
	l1PFMuons20_TH2F->Fill(   l1PFCand->eta() , l1PFCand->time() );
      */
    }
  }


}

void L1MTDPFAnalyzer::fillLLRates(l1t::PFCandidate l1PFCand, float time_eta_calibrated){

  //l1PFCandidatesFlat_TH2F->Fill(l1PFCand.eta() , time_eta_calibrated)
  l1PFCandidates_TH1F->Fill(l1PFCand.pt());
  
  if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
    l1PFNeutral_TH1F->Fill(l1PFCand.pt());
  
  if(l1PFCand.id() == l1t::PFCandidate::Muon)
    l1PFMuon_TH1F->Fill(l1PFCand.pt());
  
  
  if(time_eta_calibrated < 1){
    l1PFCandidates_1ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_1ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_1ns_TH1F->Fill(l1PFCand.pt());
    
  }
  
  if(time_eta_calibrated < 2 && time_eta_calibrated > 1){
    l1PFCandidates_2ns_TH1F->Fill(l1PFCand.pt());

    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_2ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_2ns_TH1F->Fill(l1PFCand.pt());
    
  }
  
  if(time_eta_calibrated < 3 && time_eta_calibrated > 2){
    l1PFCandidates_3ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_3ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_3ns_TH1F->Fill(l1PFCand.pt());
  }
  
  if(time_eta_calibrated < 4 && time_eta_calibrated > 3){
    l1PFCandidates_4ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_4ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_4ns_TH1F->Fill(l1PFCand.pt());
  }
  
  if(time_eta_calibrated < 5 && time_eta_calibrated > 4){
    l1PFCandidates_5ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_5ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_5ns_TH1F->Fill(l1PFCand.pt());
  }
  
  if(time_eta_calibrated < 6 && time_eta_calibrated > 5){
    l1PFCandidates_6ns_TH1F->Fill(l1PFCand.pt());

    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_6ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_6ns_TH1F->Fill(l1PFCand.pt());
  }
  
  if(time_eta_calibrated < 7 && time_eta_calibrated > 6){
    l1PFCandidates_7ns_TH1F->Fill(l1PFCand.pt());

    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_7ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_7ns_TH1F->Fill(l1PFCand.pt());
  }
  
  if(time_eta_calibrated < 8 && time_eta_calibrated > 7){
    l1PFCandidates_8ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_8ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::Muon)
      l1PFMuon_8ns_TH1F->Fill(l1PFCand.pt());
    
  }
  
  if(time_eta_calibrated < 9 && time_eta_calibrated > 8){
    l1PFCandidates_9ns_TH1F->Fill(l1PFCand.pt());
    
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_9ns_TH1F->Fill(l1PFCand.pt());
  }
  
  if( time_eta_calibrated > 9){
    l1PFCandidates_10ns_TH1F->Fill(l1PFCand.pt());
    if(l1PFCand.id() == l1t::PFCandidate::NeutralHadron)
      l1PFNeutral_10ns_TH1F->Fill(l1PFCand.pt());
  }
  //std::cout<<"Finished filling LL Rates"<<std::endl;
}

void L1MTDPFAnalyzer::findJets(Handle<l1t::PFCandidateCollection>  pfCandidates, float time_cut, std::vector<jet> &foundjets){

  PFCandidateCollection seeds;
  //find all seeds, 10 GeV 
  int nSeeds = 0;
  std::vector<jet> tempjets;
  for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
       pfCand != pfCandidates->end(); 
       ++pfCand){

    if(pfCand->pt()>10){
      seeds.push_back(*pfCand);
      nSeeds++;
    }
  }

  if(nSeeds > 10){
    std::cout<<"number of seeds greater than 10!! Maybe increase the seed threshold"<<std::endl;
  }

  for(auto seed: seeds){
    jet temp;
    temp.pt   = seed.pt();
    temp.eta  = seed.eta();
    temp.phi  = seed.phi();
    temp.time = seed.time();    

    for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
	 pfCand != pfCandidates->end(); 
	 ++pfCand){
      if(reco::deltaR( pfCand->eta(), pfCand->phi(),temp.eta, temp.phi)<0.4){
	if(time_cut!=-1 && fabs(pfCand->time()-temp.time) < time_cut )
	temp.eta = (temp.pt*temp.eta + pfCand->pt()*pfCand->eta())/(temp.pt + pfCand->pt());
	temp.phi = (temp.pt*temp.phi + pfCand->pt()*pfCand->phi())/(temp.pt + pfCand->pt());
	temp.pt += pfCand->pt();
      }
    }
    tempjets.push_back(temp);
  }  
  //cluster PF Cands near the jet - deltaR, time cut configurable
  //sort jets by pt, return top 6
  std::sort(tempjets.begin(), tempjets.end(), [](jet i,jet j){return(i.pt > j.pt);});   

  // set the first jet, highest pt jet
  foundjets.push_back(tempjets.at(0));

  //cross clean the jets
  for(auto jet : tempjets){

    bool passDeltaRCleaning = true;

    // check if temp jet is a new jet and if it is far enough away from the sorted pt jets
    for(auto nextJet : foundjets){
      if(reco::deltaR(jet.eta,jet.phi,nextJet.eta,nextJet.phi)<0.4){
	passDeltaRCleaning = false;
      }	
    }

    if(passDeltaRCleaning)
      foundjets.push_back(jet);
  }
}

void L1MTDPFAnalyzer::findLLJets(Handle<l1t::PFCandidateCollection>  pfCandidates, float time_cut, float minimum_offset, std::vector<jet> &foundjets){

  PFCandidateCollection seeds;
  std::vector<jet> tempjets;
  //find all seeds, 5 GeV 
  int nSeeds = 0;
  for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
       pfCand != pfCandidates->end(); 
       ++pfCand){

    float time_eta_calibrated = 0;

    // calculate flattened time
    if(fabs(pfCand->eta())<1.67)
      time_eta_calibrated = pfCand->time() - 15/6*pfCand->eta()*pfCand->eta() - 3.5;
    else
      time_eta_calibrated = pfCand->time() - 10;

    if((pfCand->pt() > 10) && (pfCand->time() > minimum_offset)){
      seeds.push_back(*pfCand);
      nSeeds++;
    }
  }

  if(nSeeds > 10){
    std::cout<<"number of seeds greater than 10!! Maybe increase the seed threshold"<<std::endl;
  }

  for(auto seed: seeds){
    jet temp;
    temp.pt   = seed.pt();
    temp.eta  = seed.eta();
    temp.phi  = seed.phi();
    temp.time = seed.time();    

    for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
	 pfCand != pfCandidates->end(); 
	 ++pfCand){
      if(reco::deltaR( pfCand->eta(), pfCand->phi(),temp.eta, temp.phi)<0.4){
	if(time_cut!=-1 && fabs(pfCand->time()-temp.time) < time_cut )
	temp.eta = (temp.pt*temp.eta + pfCand->pt()*pfCand->eta())/(temp.pt + pfCand->pt());
	temp.phi = (temp.pt*temp.phi + pfCand->pt()*pfCand->phi())/(temp.pt + pfCand->pt());
	temp.pt += pfCand->pt();
      }
    }
    tempjets.push_back(temp);
  }  
  //cluster PF Cands near the jet - deltaR, time cut configurable
  //sort jets by pt, return top 6
  std::sort(foundjets.begin(), foundjets.end(), [](jet i,jet j){return(i.pt > j.pt);});   

  // set the first jet, highest pt jet
  foundjets.push_back(tempjets.at(0));

  //cross clean the jets
  for(auto jet : tempjets){

    bool passDeltaRCleaning = true;

    // check if temp jet is a new jet and if it is far enough away from the sorted pt jets
    for(auto nextJet : foundjets){
      if(reco::deltaR(jet.eta,jet.phi,nextJet.eta,nextJet.phi)<0.4){
	passDeltaRCleaning = false;
      }	
    }

    if(passDeltaRCleaning)
      foundjets.push_back(jet);
  }

}

//take in electrons, photons, muons, taus, jets
void L1MTDPFAnalyzer::fillRates( std::vector<object> l1Electrons, std::vector<object> l1Photons, std::vector<object> l1Muons, std::vector<object> l1Taus, std::vector<jet> l1Jets, std::vector<jet> l1Jets_time){


  if(l1Electrons.size()>0){
    l1PFElectronRate->Fill(     l1Electrons.at(0).pt);
  }    

  bool fill        = false;
  bool fillIso90     = false;
  bool fillIso90Time = false;
  bool fillIso95     = false;
  bool fillIso95Time = false;
  
  bool fillB        = false;
  bool fillBIso90     = false;
  bool fillBIso90Time = false;
  bool fillBIso95     = false;
  bool fillBIso95Time = false;

  /* 1cm Z  
  float WP90ele      = 0.11;
  float WP90ele_time = 0.065;

  float WP95ele      = 0.16;
  float WP95ele_time = 0.1;

  float WP90mu       = 0.11;
  float WP90mu_time  = 0.065;

  float WP95mu       = 0.15;
  float WP95mu_time  = 0.09;

  float WP90tau      = 0.26;
  float WP90tau_time = 0.064;

  float WP95tau      = 0.49;
  float WP95tau_time = 0.24;
  */

  // 0.5cm Z  
  /*
  float WP90ele      = 0.075;
  float WP95ele      = 0.125;

  float WP90ele_time = 0.042;
  float WP95ele_time = 0.075;

  float WP90mu       = 0.067;
  float WP95mu       = 0.095;

  float WP90mu_time  = 0.01;
  float WP95mu_time  = 0.065;

  float WP90tau      = 0.21;
  float WP95tau      = 0.45;

  float WP90tau_time = 0.05;
  float WP95tau_time = 0.20;
*/

  // 0.3cm Z  

  float WP90ele      = 0.015;
  float WP95ele      = 0.07;

  float WP90ele_time = 0.065;
  float WP95ele_time = 0.12;

  float WP90mu       = 0.01;
  float WP95mu       = 0.053;

  float WP90mu_time  = 0.055;
  float WP95mu_time  = 0.081;

  float WP90tau      = 0.005;
  float WP95tau      = 0.16;

  float WP90tau_time = 0.2;
  float WP95tau_time = 0.37;



  for(auto object : l1Electrons){

    //electron iso time barrel
    if(fabs(object.eta)<1.6){
      if(!fillB){
	//if(object.isoSum/object.pt < 0.2)
	l1PFElectronRate_Barrel->Fill(    object.pt);
	fillB = true;
      }

      if((object.isoSumTime-object.pt)/object.pt < WP90ele_time && !fillBIso90Time){
	l1PFElectronRate_Iso90TIME_Barrel->Fill(object.pt);
	fillBIso90Time = true;
      }
      if((object.isoSumTime-object.pt)/object.pt < WP95ele_time && !fillBIso95Time){
	l1PFElectronRate_Iso95TIME_Barrel->Fill(object.pt);
	fillBIso95Time = true;
      }

      if((object.isoSum-object.pt)/object.pt < WP90ele && !fillBIso90){
	l1PFElectronRate_Iso90_Barrel->Fill(object.pt);
	fillBIso90 = true;
      }

      if((object.isoSum-object.pt)/object.pt < WP95ele && !fillBIso95){
	l1PFElectronRate_Iso95_Barrel->Fill(object.pt);
	fillBIso95 = true;
      }
    }
    
    //full detector
    if((object.isoSumTime-object.pt)/object.pt < WP90ele_time && !fillIso90Time){
      l1PFElectronRate_Iso90TIME->Fill(object.pt);
      fillIso90Time = true;
    }
    if((object.isoSumTime-object.pt)/object.pt < WP95ele_time && !fillIso95Time){
      l1PFElectronRate_Iso95TIME->Fill(object.pt);
      fillIso95Time = true;
    }
    
    if((object.isoSum-object.pt)/object.pt < WP90ele && !fillIso90){
      l1PFElectronRate_Iso90->Fill(object.pt);
      fillIso90 = true;
    }
    
    if((object.isoSum-object.pt)/object.pt < WP95ele && !fillIso95){
      l1PFElectronRate_Iso95->Fill(object.pt);
      fillIso95 = true;
    }
    
  }
  
  fill        = false;
  fillIso90     = false;
  fillIso90Time = false;
  fillIso95     = false;
  fillIso95Time = false;
  
  fillB        = false;
  fillBIso90     = false;
  fillBIso90Time = false;
  fillBIso95     = false;
  fillBIso95Time = false;
  
  if(l1Muons.size()>0){
    //if(object.isoSum/object.pt < 0.2)
    l1PFMuonRate->Fill(     l1Muons.at(0).pt);
  }    

  
  for(auto object : l1Muons){

    //electron iso time barrel
    if(fabs(object.eta)<1.6){
      if(!fillB){
	//if(object.isoSum/object.pt < 0.2)
	l1PFMuonRate_Barrel->Fill(    object.pt);
	fillB = true;
      }

      if((object.isoSumTime-object.pt)/object.pt < WP90mu_time && !fillBIso90Time){
	l1PFMuonRate_Iso90TIME_Barrel->Fill(object.pt);
	fillBIso90Time = true;
      }
      if((object.isoSumTime-object.pt)/object.pt < WP95mu_time && !fillBIso95Time){
	l1PFMuonRate_Iso95TIME_Barrel->Fill(object.pt);
	fillBIso95Time = true;
      }

      if((object.isoSum-object.pt)/object.pt < WP90mu && !fillBIso90){
	l1PFMuonRate_Iso90_Barrel->Fill(object.pt);
	fillBIso90 = true;
      }

      if((object.isoSum-object.pt)/object.pt < WP95mu && !fillBIso95){
	l1PFMuonRate_Iso95_Barrel->Fill(object.pt);
	fillBIso95 = true;
      }
    }
    
    //full detector
    if((object.isoSumTime-object.pt)/object.pt < WP90mu_time && !fillIso90Time){
      l1PFMuonRate_Iso90TIME->Fill(object.pt);
      fillIso90Time = true;
    }
    if((object.isoSumTime-object.pt)/object.pt < WP95mu_time && !fillIso95Time){
      l1PFMuonRate_Iso95TIME->Fill(object.pt);
      fillIso95Time = true;
    }
    
    if((object.isoSum-object.pt)/object.pt < WP90mu && !fillIso90){
      l1PFMuonRate_Iso90->Fill(object.pt);
      fillIso90 = true;
    }
    
    if((object.isoSum-object.pt)/object.pt < WP95mu && !fillIso95){
      l1PFMuonRate_Iso95->Fill(object.pt);
      fillIso95 = true;
    }

  }

  
  if(l1Photons.size()>0){
    l1PFPhotonRate->Fill(     l1Photons.at(0).pt);
  }

  if(l1Taus.size()>0){
    l1PFTauRate->Fill(     l1Taus.at(0).pt);
  }  

  if(l1Taus.size()>1){
    l1DiPFTauRate->Fill(     l1Taus.at(1).pt);
  }  

  fillIso90     = false;
  fillIso90Time = false;
  fillIso95     = false;
  fillIso95Time = false;

  fillB          = false;
  fillBIso90     = false;
  fillBIso90Time = false;
  fillBIso95     = false;
  fillBIso95Time = false;
  
  bool DfillIso90     = false;
  bool DfillIso90Time = false;
  bool DfillIso95     = false;
  bool DfillIso95Time = false;

  bool DfillB         = false;
  bool DfillBIso90     = false;
  bool DfillBIso90Time = false;
  bool DfillBIso95     = false;
  bool DfillBIso95Time = false;
  
  for(auto object : l1Taus){
    
    //electron iso time barrel
    if(fabs(object.eta)<1.6){
      // if barrel already filled once, fill double tau
      if(fillB){
	if(!DfillB){
	  l1DiPFTauRate_Barrel->Fill(    object.pt);
	  DfillB = true;
	}
      }
      
      if(!fillB){
	//if(object.isoSum/object.pt < 0.2)
	l1PFTauRate_Barrel->Fill(    object.pt);
	fillB = true;
      }

      // if barrel already filled once, fill double tau
      if(fillBIso90Time){
	if((object.isoSumTime-object.pt)/object.pt < WP90tau_time && !DfillBIso90Time){
	  l1DiPFTauRate_Iso90TIME_Barrel->Fill(    object.pt);
	  DfillBIso90Time = true;
	}
      }

      if((object.isoSumTime-object.pt)/object.pt < WP90tau_time && !fillBIso90Time){
	l1PFTauRate_Iso90TIME_Barrel->Fill(object.pt);
	fillBIso90Time = true;
      }
      
      // if barrel already filled once, fill double tau
      if(fillBIso95Time){
	if((object.isoSumTime-object.pt)/object.pt < WP95tau_time && !DfillBIso95Time){
	  l1DiPFTauRate_Iso95TIME_Barrel->Fill(    object.pt);
	  DfillBIso95Time = true;
	}
      }
      
      if((object.isoSumTime-object.pt)/object.pt < WP95tau_time && !fillBIso95Time){
	l1PFTauRate_Iso95TIME_Barrel->Fill(object.pt);
	fillBIso95Time = true;
      }
      
      // if barrel already filled once, fill double tau
      if(fillBIso90){
	if((object.isoSum-object.pt)/object.pt < WP90tau && !DfillBIso90){
	  l1DiPFTauRate_Iso90_Barrel->Fill(object.pt);
	  DfillBIso90 = true;
	}
      }
      
      if((object.isoSum-object.pt)/object.pt < WP90tau && !fillBIso90){
	l1PFTauRate_Iso90_Barrel->Fill(object.pt);
	fillBIso90 = true;
      }

      
      // if barrel already filled once, fill double tau
      if(fillBIso95){
	if((object.isoSum-object.pt)/object.pt < WP95tau && !DfillBIso95){
	  l1DiPFTauRate_Iso95_Barrel->Fill(object.pt);
	  DfillBIso95 = true;
	}
      }
      
      if((object.isoSum-object.pt)/object.pt < WP95tau && !fillBIso95){
	l1PFTauRate_Iso95_Barrel->Fill(object.pt);
	fillBIso95 = true;
      }
    }
    
    //full detector
    // if already filled once, fill double tau
    if(fillIso90Time){
      if((object.isoSumTime-object.pt)/object.pt < WP90tau_time && !DfillIso90Time){
	l1DiPFTauRate_Iso90TIME->Fill(    object.pt);
	DfillIso90Time = true;
      }
    }
    
    if((object.isoSumTime-object.pt)/object.pt < WP90tau_time && !fillIso90Time){
      l1PFTauRate_Iso90TIME->Fill(object.pt);
      fillIso90Time = true;
    }
    
    // if barrel already filled once, fill double tau
    if(fillIso95Time){
      if((object.isoSumTime-object.pt)/object.pt < WP95tau_time && !DfillIso95Time){
	l1DiPFTauRate_Iso95TIME->Fill(    object.pt);
	DfillIso95Time = true;
      }
    }
    
    if((object.isoSumTime-object.pt)/object.pt < WP95tau_time && !fillIso95Time){
      l1PFTauRate_Iso95TIME->Fill(object.pt);
      fillIso95Time = true;
    }

    // if barrel already filled once, fill double tau
    if(fillIso90){
      if((object.isoSum-object.pt)/object.pt < WP90tau && !DfillIso90){
	l1DiPFTauRate_Iso90->Fill(object.pt);
	DfillIso90 = true;
      }
    }
    
    if((object.isoSum-object.pt)/object.pt < WP90tau && !fillIso90){
      l1PFTauRate_Iso90->Fill(object.pt);
      fillIso90 = true;
    }

    // if barrel already filled once, fill double tau
    if(fillIso95){
      if((object.isoSum-object.pt)/object.pt < WP95tau && !DfillIso95){
	l1DiPFTauRate_Iso95->Fill(object.pt);
	DfillIso95 = true;
      }
    }
    
    
    if((object.isoSum-object.pt)/object.pt < WP95tau && !fillIso95){
      l1PFTauRate_Iso95->Fill(object.pt);
      fillIso95 = true;
    }
  }
  
  bool fill0 = false;
  bool fill2 = false;
  bool fill3 = false;
  bool fill4 = false;
  bool fill5 = false;
  bool fill6 = false;
  bool fill7 = false;
  
  bool fillB0 = false;
  bool fillB0_time = false;
  bool fillB2 = false;
  bool fillB3 = false;
  bool fillB4 = false;
  bool fillB5 = false;
  bool fillB6 = false;
  bool fillB7 = false;
  
  
  if(l1Jets.size()>0){
    l1PFJetRate->Fill(l1Jets.at(0).pt);
    //fill0=true;
  }

  if(l1Jets_time.size()>0){
    l1PFJetRate_time->Fill(l1Jets_time.at(0).pt);
    //fill0=true;
  }


  for(auto object : l1Jets){
    if(fabs(object.eta)<1.6 && fillB0_time == false){
      l1PFJetRate_Barrel->Fill(object.pt);
      l1PFJetRate_calibrated_Barrel->Fill(object.pt*0.65);
      fillB0_time = true;
    }
  }

  for(auto object : l1Jets_time){  

    float time_eta_calibrated = 0;

    // calculate flattened time
    if(fabs(object.eta)<1.67)
      time_eta_calibrated = object.time - 15/6*object.eta*object.eta - 3.5;
    else
      time_eta_calibrated = object.time - 10;
    
    
    if(time_eta_calibrated > 2 && fill2 == false){
      l1PFJet2nsRate->Fill(object.pt);
      fill2 = true;
    }

    if(time_eta_calibrated > 3 && fill3 ==false){      
      l1PFJet3nsRate->Fill(object.pt);
      fill3 = true;
    }

    if(time_eta_calibrated > 4 && fill4 ==false){
      l1PFJet4nsRate->Fill(object.pt);
      fill4 = true;
    }

    if(time_eta_calibrated > 5 && fill5 == false){
      l1PFJet5nsRate->Fill(object.pt);
      fill5 = true;
    }

    if(time_eta_calibrated > 6 && fill6 == false){
      l1PFJet6nsRate->Fill(object.pt);
      fill6 = true;
    }

    if(time_eta_calibrated > 7 && fill7 ==false ){
      l1PFJet7nsRate->Fill(object.pt);
      fill7 = true;
    }



    if(fabs(object.eta)<1.6 ){

      if(fillB0==false){
	l1PFJetRate_time_Barrel->Fill(object.pt);
	l1PFJetRate_time_calibrated_Barrel->Fill(object.pt*0.65);
	fillB0 = true;
      }

      if(time_eta_calibrated > 2 && fillB2 ==false){
	l1PFJet2nsRate_Barrel->Fill(object.pt);
	l1PFJet2nsRate_calibrated_Barrel->Fill(object.pt*0.65);
	fillB2 = true;
      }

      if(time_eta_calibrated > 3 && fillB3 ==false){
	l1PFJet3nsRate_Barrel->Fill(object.pt);
	l1PFJet3nsRate_calibrated_Barrel->Fill(object.pt*0.65);
	fillB3 = true;
      }
      
      if(time_eta_calibrated > 4 && fillB4 == false){
	l1PFJet4nsRate_Barrel->Fill(object.pt);
	l1PFJet4nsRate_calibrated_Barrel->Fill(object.pt*0.65);
	fillB4 = true;
      }

      if(time_eta_calibrated > 5 && fillB5 ==false){
	l1PFJet5nsRate_Barrel->Fill(object.pt);
	l1PFJet5nsRate_calibrated_Barrel->Fill(object.pt*0.65);
	fillB5 = true;
      }

      if(time_eta_calibrated > 6 && fillB6 ==false){
	l1PFJet6nsRate_Barrel->Fill(object.pt);
	l1PFJet6nsRate_calibrated_Barrel->Fill(object.pt*0.65);
	fillB6 = true;
      }

      if(time_eta_calibrated > 7 && fillB7 ==false){
	l1PFJet7nsRate_Barrel->Fill(object.pt);
	l1PFJet7nsRate_calibrated_Barrel->Fill(object.pt*0.65);
	fillB7 = true;
      }
    }
  }
}
// Finish me!!
void L1MTDPFAnalyzer::fillL1TauTreeVariables(Handle<L1PFTauCollection> l1PFTaus, Handle< l1t::PFCandidateCollection > l1PFCandidates, std::vector<genVisTau> genTaus, Handle<L1TkPrimaryVertexCollection> L1VertexHandle){
  for(auto &l1Cand : *l1PFTaus){

    zeroTauTreeVariables();

    l1DM = l1Cand.tauType();
    tauL1Pt = l1Cand.pt();
    tauL1Eta = l1Cand.eta();
    tauL1Phi = l1Cand.phi();
    tauL1Time = l1Cand.time();
    tauL1StripPt  = l1Cand.strip_p4().pt();
    tauL1StripEta = l1Cand.strip_p4().eta();
    tauL1StripPhi = l1Cand.strip_p4().phi();
    tauL1StripDR  = reco::deltaR(l1Cand.pfRef().at(0).eta(), l1Cand.pfRef().at(0).phi(), l1Cand.strip_p4().eta(), l1Cand.strip_p4().phi());   
    track1Time = l1Cand.time();

    int nEG = 0;
    l1t::PFCandidateCollection egammas;
    if(l1Cand.egPFRef().size()>0){
      for(auto eg : l1Cand.egPFRef()){
	egammas.push_back(eg);
	nEG++;
      }
      l1TauEGTime = l1Cand.egPFRef().at(0).time();
    }
    tauL1nEG = nEG;

    std::sort(egammas.begin(), egammas.end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   
    if(l1Cand.egPFRef().size()>0){
      tauL1EGPt = egammas.at(0).pt();
    }

    //calculate iso Sum
    double isoSum = 0;
    double isoSumTime = 0;    
    calculate_charged_iso_sum_tau(l1Cand, l1PFCandidates,  time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);

    tauL1Iso = isoSum;
    tauL1Iso_time = isoSumTime;
    //for (vtxIter = L1VertexHandle->begin(); vtxIter != L1VertexHandle->end(); ++vtxIter) {
    
    zVTX = 0;
    if(L1VertexHandle->size()>0){
      zVTX = L1VertexHandle->at(0).getZvertex();
    }


    //get tracks for tau
    //for(auto pf : l1Cand.pfRef()){
    if(l1Cand.pfRef().size()>0){
      track1Z      =  l1Cand.pfRef().at(0).pfTrack()->track()->getPOCA().z();
      track1nStubs =  l1Cand.pfRef().at(0).pfTrack()->track()->getStubRefs().size();
      track1ChiSquared = l1Cand.pfRef().at(0).pfTrack()->track()->getChi2();

      if(! l1Cand.pfRef().at(0).pfCluster().isNull())
	pfCand1HoE   =  l1Cand.pfRef().at(0).pfCluster()->hOverE();

      track1PVDZ = (l1Cand.pfRef().at(0).pfTrack()->track()->getPOCA().z() - zVTX);
    }

    if(l1Cand.pfRef().size()>2){
      track2Z      =  l1Cand.pfRef().at(1).pfTrack()->track()->getPOCA().z();
      track2nStubs =  l1Cand.pfRef().at(1).pfTrack()->track()->getStubRefs().size();
      track2ChiSquared = l1Cand.pfRef().at(1).pfTrack()->track()->getChi2();
      track2PVDZ = (l1Cand.pfRef().at(1).pfTrack()->track()->getPOCA().z() - zVTX);

      if(! l1Cand.pfRef().at(1).pfCluster().isNull())
	pfCand2HoE =  l1Cand.pfRef().at(1).pfCluster()->hOverE();

      track2Time = l1Cand.pfRef().at(1).time();

      track3Z      =  l1Cand.pfRef().at(2).pfTrack()->track()->getPOCA().z();
      track3nStubs =  l1Cand.pfRef().at(2).pfTrack()->track()->getStubRefs().size();
      track3ChiSquared = l1Cand.pfRef().at(2).pfTrack()->track()->getChi2();
      track3PVDZ = (l1Cand.pfRef().at(2).pfTrack()->track()->getPOCA().z() - zVTX);
      if(! l1Cand.pfRef().at(2).pfCluster().isNull())
	pfCand3HoE =  l1Cand.pfRef().at(2).pfCluster()->hOverE();
      track3Time = l1Cand.pfRef().at(2).time();
      track12DZ = track1Z - track2Z;
      track13DZ = track1Z - track3Z;

    }

      //edm::Ref< L1TkPrimaryVertexCollection > vtxRef( L1VertexHandle, ivtx );
      //ivtx ++;
      //float tVTX = vtxIter->getTime();
      //}


    
    for(auto genTau : genTaus){
      if(reco::deltaR(l1Cand.eta(), l1Cand.phi(), genTau.p4.eta(), genTau.p4.phi())<0.1){
	tauGenPt      = genTau.p4.pt();
	tauGenEta     = genTau.p4.eta();
	tauGenPhi     = genTau.p4.phi();
	break;
      }
    }


    L1TauTree->Fill();
  }

}

void L1MTDPFAnalyzer::zeroElectronTreeVariables(){
  eleRecoPt = -10, eleRecoEta = -10, eleRecoPhi = -10;
  eleGenPt  = -10, eleGenEta  = -10, eleGenPhi  = -10; 
  eleL1Pt   = -10, eleL1Eta   = -10, eleL1Phi   = -10; 
  eleL1Iso = -10, eleL1Iso_time = -10;
}

void L1MTDPFAnalyzer::zeroPhotonTreeVariables(){
  gammaRecoPt = -10, gammaRecoEta = -10, gammaRecoPhi = -10;
  gammaGenPt  = -10, gammaGenEta  = -10, gammaGenPhi  = -10; 
  gammaL1Pt   = -10, gammaL1Eta   = -10, gammaL1Phi   = -10; 
  gammaL1Iso = -10, gammaL1Iso_time = -10;
}

void L1MTDPFAnalyzer::zeroMuonTreeVariables(){
  muRecoPt = -10, muRecoEta = -10, muRecoPhi = -10;
  muGenPt  = -10, muGenEta  = -10, muGenPhi  = -10; 
  muL1Pt   = -10, muL1Eta   = -10, muL1Phi   = -10; 
  muL1Iso = -10, muL1Iso_time = -10;
}

void L1MTDPFAnalyzer::zeroTauTreeVariables(){
  tauRecoPt = -10, tauRecoEta = -10, tauRecoPhi = -10;
  tauGenPt  = -10, tauGenEta  = -10, tauGenPhi  = -10; 
  tauL1Pt   = -10, tauL1Eta   = -10, tauL1Phi   = -10; 
  tauL1Iso = -10, tauL1Iso_time = -10; 
  track12DZ = -10, track13DZ = -10; 
  track1PVDZ = -10, track2PVDZ = -10, track3PVDZ = -10;
  track1nStubs = -10, track2nStubs = -10, track3nStubs = -10;
  track1Time = -10, track2Time = -10, track3Time = -10;
  l1DM = -10; 
  track1ChiSquared = -10, track2ChiSquared = -10, track3ChiSquared = -10;
  zVTX = -10;
  track1Z = -10, track2Z = -10, track3Z = -10;
  tauL1StripPt = -10, tauL1StripEta = -10, tauL1StripPhi = -10, tauL1StripDR = -10;
  pfCand1HoE = -10, tauL1nEG = -10, tauL1EGPt = -10;
}


// ------------ method called once each job just before starting event loop  ------------
void
L1MTDPFAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1MTDPFAnalyzer::endJob()
{
}

void
L1MTDPFAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(L1MTDPFAnalyzer);
