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

using namespace l1t;
using namespace edm;
using namespace std;
struct genVisTau{ reco::Candidate::LorentzVector p4; int decayMode;};

class mtdIsoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

  explicit mtdIsoAnalyzer(const edm::ParameterSet&);

  ~mtdIsoAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //custom class functions
  void calculate_charged_iso_sum(l1t::PFCandidate object, int object_index, Handle< l1t::PFCandidateCollection > pfCands, edm::Handle<edm::ValueMap<float> > timingValues, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime);
  void fillGenParticleVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<reco::GenParticle> &selectedGenParticles,int pdgID);
  //void fillGenTauVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<genVisTau> &selectedGenParticles);
  // void fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates,PFCandidateCollection &pfElectrons,PFCandidateCollection &pfPhotons,PFCandidateCollection &pfMuons);
  void fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates,PFCandidateCollection &pfElectrons, std::vector<int> & l1PFElectrons_indicies, PFCandidateCollection &pfPhotons, std::vector<int> & l1PFPhotons_indicies, PFCandidateCollection &pfMuons, std::vector<int> &l1PFMuons_indicies);
  void zeroElectronTreeVariables();
  void zeroPhotonTreeVariables();
  void zeroMuonTreeVariables();
  void zeroTauTreeVariables();

  EDGetTokenT<edm::ValueMap<float> >            timingValuesToken_;
  EDGetTokenT< PFCandidateCollection >          L1PFCandsToken_;
  EDGetTokenT< L1PFTauCollection >              L1PFTausToken_;
  EDGetTokenT< L1TkEtMissParticleCollection >   L1PFMETToken_;
  EDGetTokenT< L1TkEtMissParticleCollection >   L1PFMETTimeToken_;
  EDGetTokenT< std::vector<pat::Electron> >     recoElectronsToken_; 
  EDGetTokenT< std::vector<pat::Muon> >         recoMuonsToken_;     
  EDGetTokenT< std::vector<pat::Tau> >          recoTausToken_;
  EDGetTokenT< std::vector<pat::MET> >          recoMetToken_;
  EDGetTokenT< std::vector<reco::GenParticle> > genParticlesToken_;
  double time_cut_;
  double iso_cone_deltaR_;
  double iso_cone_deltaZ_;
  // Declare all the variables for the trees
  //electron tree
  TTree* electronTree;
  double eleRecoPt, eleRecoEta, eleRecoPhi, eleGenPt, eleGenEta, eleGenPhi, eleL1Pt, eleL1Eta, eleL1Phi, eleL1Iso, eleL1Iso_time;
  //photon tree
  TTree* photonTree;
  double gammaRecoPt, gammaRecoEta, gammaRecoPhi, gammaGenPt, gammaGenEta, gammaGenPhi, gammaL1Pt, gammaL1Eta, gammaL1Phi, gammaL1Iso, gammaL1Iso_time;
  //muon tree
  TTree* muonTree;
  double muRecoPt,  muRecoEta,  muRecoPhi,  muGenPt,  muGenEta,  muGenPhi,  muL1Pt,  muL1Eta,  muL1Phi,  muL1Iso,  muL1Iso_time;
  //tau tree
  TTree* tauTree;
  double tauRecoPt, tauRecoEta, tauRecoPhi, tauGenPt, tauGenEta, tauGenPhi, tauL1Pt, tauL1Eta, tauL1Phi,tauL1Iso, tauL1Iso_time;
  //met tree
  TTree* metTree;
  double metL1Et, metL1Phi, metL1Et_time, metRecoEt, metRecoPhi, metGenEt, metGenPhi;
  int run, lumi, event;
};

//Constructor
mtdIsoAnalyzer::mtdIsoAnalyzer(const edm::ParameterSet &cfg) :  //inputs L1PFCands, reco electrons, reco muons, reco taus, gen particles
  timingValuesToken_(  consumes<edm::ValueMap<float> >(                       cfg.getParameter<edm::InputTag>("timingValuesNominal"))),
  L1PFCandsToken_(     consumes< std::vector<l1t::PFCandidate>  >(            cfg.getParameter<InputTag>("l1PFCands")     )),
  L1PFTausToken_(      consumes< L1PFTauCollection >(                         cfg.getParameter<InputTag>("l1PFTaus")      )),
  L1PFMETToken_(       consumes< L1TkEtMissParticleCollection >(              cfg.getParameter<InputTag>("l1PFMET")       )),
  L1PFMETTimeToken_(   consumes< L1TkEtMissParticleCollection >(              cfg.getParameter<InputTag>("l1PFMETTime")   )),
  recoElectronsToken_( consumes< std::vector<pat::Electron> >(                cfg.getParameter<InputTag>("recoElectrons") )),
  recoMuonsToken_(     consumes< std::vector<pat::Muon> >(                    cfg.getParameter<InputTag>("recoMuons")     )),
  recoTausToken_(      consumes< std::vector<pat::Tau> >(                     cfg.getParameter<InputTag>("recoTaus")      )),
  recoMetToken_(       consumes< std::vector<pat::MET> >(                     cfg.getParameter<InputTag>("recoMet")       )),
  genParticlesToken_(  consumes< std::vector<reco::GenParticle> >(            cfg.getParameter<InputTag>("genParticles")  )),
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
  electronTree->Branch("l1Iso",        &eleL1Iso,          "l1Iso/D"       );
  electronTree->Branch("l1Iso_time",   &eleL1Iso_time,     "l1Iso_time/D"  );
  
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
  tauTree->Branch("l1Iso",             &tauL1Iso,          "l1Iso/D"       );
  tauTree->Branch("l1Iso_time",        &tauL1Iso_time,     "l1Iso_time/D"  );

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
  
}

//destructor
mtdIsoAnalyzer::~mtdIsoAnalyzer()
{
}

//deltaR, deltaZ, time_Cut 
void 
mtdIsoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  run   = iEvent.id().run();
  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();

  edm::Handle<edm::ValueMap<float> > timingValues;
  iEvent.getByToken(timingValuesToken_,timingValues);

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
  l1t::PFCandidateCollection l1PFElectrons, l1PFPhotons, l1PFMuons;
  std::vector<int> l1PFElectrons_indicies, l1PFPhotons_indicies, l1PFMuons_indicies;
  Handle<L1PFTauCollection> l1PFTaus;

  //put l1 met collection here
  Handle<L1TkEtMissParticleCollection> l1Met;
  Handle<L1TkEtMissParticleCollection> l1Met_time;

  if(!iEvent.getByToken( L1PFCandsToken_, l1PFCandidates))    std::cout<<"No L1PF Cands Found!!"     <<std::endl;
  if(!iEvent.getByToken( L1PFTausToken_,  l1PFTaus))          std::cout<<"No L1PF Taus Found!!"      <<std::endl;
  if(!iEvent.getByToken( L1PFMETToken_,   l1Met))             std::cout<<"No L1PF Mets Found!!"      <<std::endl;

  //Get the reco objects
  if(!iEvent.getByToken( recoElectronsToken_, recoElectrons)) std::cout<<"No Reco ELECTRONS Found!!!"<<std::endl;
  //if(!iEvent.getByToken( recoPhotonsToken_,   recoPhotons))   std::cout<<"No Reco PHOTONS Found!!!"  <<std::endl;
  if(!iEvent.getByToken( recoMuonsToken_,     recoMuons))     std::cout<<"No Reco MUONS Found!!!"    <<std::endl;
  if(!iEvent.getByToken( recoTausToken_,      recoTaus))      std::cout<<"No Reco TAUS Found!!!"     <<std::endl;
  if(!iEvent.getByToken( recoMetToken_,       recoMet))       std::cout<<"No Reco MET Found!!!"      <<std::endl;

  // sort by pt not necessary for now... perhaps needed for an emulator
  //std::sort(l1PFCandidates_sort.begin(), l1PFCandidates_sort.end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   

  //get gen objects
  Handle<std::vector<reco::GenParticle> > genParticles;
  if(!iEvent.getByToken(genParticlesToken_,genParticles)) std::cout<<"No Gen Particles Found!!!"<<std::endl;

  //Fill each gen particle vector with the corresponding gen particles
  //electrons pdgid == 11
  fillGenParticleVector(genParticles, genElectrons, 11);

  //photons   pdgid == 22
  fillGenParticleVector(genParticles, genPhotons,   22);

  //muons     pdgid == 13 
  fillGenParticleVector(genParticles, genMuons,     13);

  //taus      pdgid == 15 // but this is a special lepton, using visible decay products
  //fillGenTauVector(genParticles, genTaus);

  //met       // composite object

  //Fill Each L1 PF Cand vector with the L1 PFCands Electrons, Photons and Muons
  fillL1PFCandCollections(l1PFCandidates, l1PFElectrons, l1PFElectrons_indicies, l1PFPhotons, l1PFPhotons_indicies, l1PFMuons, l1PFMuons_indicies);
  
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
	calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, timingValues, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	eleL1Iso      = isoSum;
	eleL1Iso_time = isoSumTime;
	std::cout<<"matched electron found!! eleL1Iso: "<<eleL1Iso<<" eleL1Iso_time: "<<eleL1Iso_time<<std::endl;
	eleL1Pt       = l1Cand.pt();
	eleL1Eta      = l1Cand.eta();
	eleL1Phi      = l1Cand.phi();
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

  //Fill the Photons
  /*
  idx = 0;
  for( std::vector<pat::Photon>::const_iterator cand  = recoPhotons->begin(); 
                                                  cand != recoPhotons->end(); 
                                                  ++ cand, ++ idx) {
    double isoSum = 0;
    double isoSumTime = 0;    

    //Zero the tree
    zeroPhotonTreeVariables();
    eleRecoPt  = cand->pt();
    eleRecoEta = cand->eta();
    eleRecoPhi = cand->phi();

    //Get matched L1 Photon
    for(auto l1Cand : l1PFPhotons){
      if(reco::deltaR(cand->eta(),cand->phi(),l1Cand.eta(),l1Cand.phi())<0.2){
	//calculate iso Sum
	calculate_charged_iso_sum(l1Cand, l1PFCandidates, timingValues, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	eleL1Iso      = isoSum;
	eleL1Iso_time = isoSumTime;
	std::cout<<"matched photon found!! eleL1Iso: "<<eleL1Iso<<" eleL1Iso_time: "<<eleL1Iso_time<<std::endl;
	eleL1Pt       = l1Cand.pt();
	eleL1Eta      = l1Cand.eta();
	eleL1Phi      = l1Cand.phi();
	break;
      }
    }
    //Get matched gen Photon
    for(auto genCand : genPhotons){
      //matchGen
      if(reco::deltaR(cand->eta(), cand->phi(), genCand.eta(), genCand.phi())<0.1){
	eleGenPt      = genCand.pt();
	eleGenEta     = genCand.eta();
	eleGenPhi     = genCand.phi();
	break;
      }
    }
    photonTree->Fill();
  }
  */

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
	calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, timingValues, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);
	muL1Iso      = isoSum;
	muL1Iso_time = isoSumTime;
	std::cout<<"matched muon found!! muL1Iso: "<<muL1Iso<<" muL1Iso_time: "<<muL1Iso_time<<std::endl;
	muL1Pt       = l1Cand.pt();
	muL1Eta      = l1Cand.eta();
	muL1Phi      = l1Cand.phi();
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

    //Fill the Taus

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

void mtdIsoAnalyzer::calculate_charged_iso_sum(l1t::PFCandidate object, int object_index, Handle< l1t::PFCandidateCollection > pfCandidates, edm::Handle<edm::ValueMap<float> > timingValues, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime){
    double iso_sum = 0;
    double iso_sum_time = 0;
    int idx = 0;
    for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
                                               pfCand != pfCandidates->end(); 
	                                       ++pfCand, ++idx){    

      if(idx == object_index)
	continue;
      
      //first check if it is a charged candidate
      if(pfCand->id() == l1t::PFCandidate::Muon || pfCand->id() == l1t::PFCandidate::Electron || pfCand->id() == l1t::PFCandidate::ChargedHadron){
      
	// first check deltaR match
	if(reco::deltaR(object.eta(), object.phi(), pfCand->eta(), pfCand->phi()) < deltaR){

	  if(object.pfTrack().isNull()||pfCand->pfTrack().isNull())
	    continue;
	  
	  if(fabs(object.pfTrack()->track()->getPOCA().z()-pfCand->pfTrack()->track()->getPOCA().z()) < deltaZ){
	    
	    iso_sum += pfCand->pt();
	    float objectTime = 0;
	    float pfCandTime = 0;
	    
	    
	    if(!object.pfTrack().isNull()){
	      objectTime = (*timingValues)[object.pfTrack()->track()];
	    }
	    if(!pfCand->pfTrack().isNull()){
	      pfCandTime = (*timingValues)[pfCand->pfTrack()->track()];
	    }
	    //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > objectTime(pfCand->pfTrack(), track_index);
	  //edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > pfCandTime(object.pfTrack(), track_index);
	    
	  //std::cout<<"objectTime: "<<objectTime<<" pfCandTime:"<<pfCandTime<<std::endl;
	  // then check time match
	    if(fabs(objectTime - pfCandTime) < time_cut){
	      iso_sum_time += pfCand->pt();
	    }
	  }
	}
      }
    }
    isoSum = iso_sum;
    isoSumTime = iso_sum_time;
}

void mtdIsoAnalyzer::fillGenParticleVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<reco::GenParticle> &selectedGenParticles,int pdgID){
  for( std::vector<reco::GenParticle>::const_iterator cand  = genParticles->begin(); 
                                                      cand != genParticles->end(); 
                                                      ++ cand) {
    if(abs(cand->pdgId()) == pdgID){
      selectedGenParticles.push_back(*cand);
    }
  }
}
/*
void mtdIsoAnalyzer::fillGenTauVector(Handle<std::vector<reco::GenParticle> > genParticles, std::vector<genVisTau> &selectedGenParticles){
  for( std::vector<reco::GenParticle>::const_iterator cand  = genParticles->begin(); 
                                                      cand != genParticles->end(); 
                                                      ++ cand) {

    if(abs(cand->pdgId())==15){      
      reco::Candidate::LorentzVector visGenTau = getVisMomentum(cand, genParticles);
      genVisTau Temp;
      int decayMode  = GetDecayMode(*cand);
      Temp.p4        = visGenTau;
      Temp.decayMode = decayMode;    
      selectedGenParticles.push_back(Temp);
    }
  }
}
*/
void mtdIsoAnalyzer::fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates,PFCandidateCollection &pfElectrons, std::vector<int> & l1PFElectrons_indicies, PFCandidateCollection &pfPhotons, std::vector<int> & l1PFPhotons_indicies, PFCandidateCollection &pfMuons, std::vector<int> & l1PFMuons_indicies){
  int index = 0;
  for( PFCandidateCollection::const_iterator l1PFCand  = pfCandidates->begin();
                                             l1PFCand != pfCandidates->end(); 
       ++l1PFCand, ++index){    
    if(l1PFCand->id() == l1t::PFCandidate::Electron){
      pfElectrons.push_back(*l1PFCand);
      l1PFElectrons_indicies.push_back(index);
    }
    if(l1PFCand->id() == l1t::PFCandidate::Photon){
      pfPhotons.push_back(*l1PFCand);
      l1PFPhotons_indicies.push_back(index);
    }
    if(l1PFCand->id() == l1t::PFCandidate::Muon){
      pfMuons.push_back(*l1PFCand);
      l1PFMuons_indicies.push_back(index);
    }
  }
}

void mtdIsoAnalyzer::zeroElectronTreeVariables(){
  eleRecoPt = -10, eleRecoEta = -10, eleRecoPhi = -10;
  eleGenPt  = -10, eleGenEta  = -10, eleGenPhi  = -10; 
  eleL1Pt   = -10, eleL1Eta   = -10, eleL1Phi   = -10; 
  eleL1Iso = -10, eleL1Iso_time = -10;
}

void mtdIsoAnalyzer::zeroPhotonTreeVariables(){
  gammaRecoPt = -10, gammaRecoEta = -10, gammaRecoPhi = -10;
  gammaGenPt  = -10, gammaGenEta  = -10, gammaGenPhi  = -10; 
  gammaL1Pt   = -10, gammaL1Eta   = -10, gammaL1Phi   = -10; 
  gammaL1Iso = -10, gammaL1Iso_time = -10;
}

void mtdIsoAnalyzer::zeroMuonTreeVariables(){
  muRecoPt = -10, muRecoEta = -10, muRecoPhi = -10;
  muGenPt  = -10, muGenEta  = -10, muGenPhi  = -10; 
  muL1Pt   = -10, muL1Eta   = -10, muL1Phi   = -10; 
  muL1Iso = -10, muL1Iso_time = -10;
}

void mtdIsoAnalyzer::zeroTauTreeVariables(){
  tauRecoPt = -10, tauRecoEta = -10, tauRecoPhi = -10;
  tauGenPt  = -10, tauGenEta  = -10, tauGenPhi  = -10; 
  tauL1Pt   = -10, tauL1Eta   = -10, tauL1Phi   = -10; 
  tauL1Iso = -10, tauL1Iso_time = -10;
}


// ------------ method called once each job just before starting event loop  ------------
void
mtdIsoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
mtdIsoAnalyzer::endJob()
{
}

void
mtdIsoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(mtdIsoAnalyzer);
