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
struct simple_object{
  double pt;
  double eta;
  double relativeIso;
  double relativeIso_time;
};

class mtdIsoAnalyzerRates : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

  explicit mtdIsoAnalyzerRates(const edm::ParameterSet&);

  ~mtdIsoAnalyzerRates();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  //custom class functions
  void calculate_charged_iso_sum(l1t::PFCandidate object, int object_index, Handle< l1t::PFCandidateCollection > pfCands, edm::Handle<edm::ValueMap<float> > timingValues, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime);

  void fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates,PFCandidateCollection &pfElectrons, std::vector<int> & l1PFElectrons_indicies, PFCandidateCollection &pfPhotons, std::vector<int> & l1PFPhotons_indicies, PFCandidateCollection &pfMuons, std::vector<int> &l1PFMuons_indicies);

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
  double WP_90_ele_;
  double WP_time_90_ele_;
  double WP_95_ele_;
  double WP_time_95_ele_;
  double WP_90_mu_;
  double WP_time_90_mu_;
  double WP_95_mu_;
  double WP_time_95_mu_;
  // Declare all the variables for the rates
  //electron rate
  TH1F* l1Ele_pt;
  TH1F* l1Ele_IsoRel90_pt;
  TH1F* l1Ele_IsoRel90_time_pt;
  TH1F* l1Ele_IsoRel95_pt;
  TH1F* l1Ele_IsoRel95_time_pt;

  //photon rate
  TH1F* l1Gamma_pt;
  TH1F* l1Gamma_IsoRel90_pt;
  TH1F* l1Gamma_IsoRel90_time_pt;
  TH1F* l1Gamma_IsoRel95_pt;
  TH1F* l1Gamma_IsoRel95_time_pt;

  //muon rate
  TH1F* l1Mu_pt;
  TH1F* l1Mu_IsoRel90_pt;
  TH1F* l1Mu_IsoRel90_time_pt;
  TH1F* l1Mu_IsoRel95_pt;
  TH1F* l1Mu_IsoRel95_time_pt;

  //tau rate
  TH1F* l1Tau_pt;
  TH1F* l1Tau_IsoRel90_pt;
  TH1F* l1Tau_IsoRel90_time_pt;
  TH1F* l1Tau_IsoRel95_pt;
  TH1F* l1Tau_IsoRel95_time_pt;    

  //met rate

};

//Constructor
mtdIsoAnalyzerRates::mtdIsoAnalyzerRates(const edm::ParameterSet &cfg) :  //inputs L1PFCands, reco electrons, reco muons, reco taus, gen particles
  timingValuesToken_( consumes<edm::ValueMap<float> >(             cfg.getParameter<edm::InputTag>("timingValuesNominal"))),
  L1PFCandsToken_(    consumes< std::vector<l1t::PFCandidate> >(   cfg.getParameter<InputTag>("l1PFCands")     )),
  L1PFTausToken_(     consumes< L1PFTauCollection >(               cfg.getParameter<InputTag>("l1PFTaus")      )),
  L1PFMETToken_(      consumes< L1TkEtMissParticleCollection >(    cfg.getParameter<InputTag>("l1PFMET")       )),
  L1PFMETTimeToken_(  consumes< L1TkEtMissParticleCollection >(    cfg.getParameter<InputTag>("l1PFMETTime")   )),
  time_cut_(          cfg.getParameter<double>("time_cut"     )),
  iso_cone_deltaR_(   cfg.getParameter<double>("isoConeDeltaR")),
  iso_cone_deltaZ_(   cfg.getParameter<double>("isoConeDeltaZ")),
  WP_90_ele_(         cfg.getParameter<double>("WP90ele")),
  WP_time_90_ele_(    cfg.getParameter<double>("WPtime90ele")),
  WP_95_ele_(         cfg.getParameter<double>("WP95ele")),
  WP_time_95_ele_(    cfg.getParameter<double>("WPtime95ele")),
  WP_90_mu_(          cfg.getParameter<double>("WP90mu")),
  WP_time_90_mu_(     cfg.getParameter<double>("WPtime90mu")),
  WP_95_mu_(          cfg.getParameter<double>("WP95mu")),
  WP_time_95_mu_(     cfg.getParameter<double>("WPtime95mu"))   
{
  //services
  usesResource("TFileService");
  Service<TFileService> fs;

  l1Ele_pt		   = fs->make<TH1F>( "l1Ele_pt"                 , "p_{t}", 300,  0., 300. );
  l1Ele_IsoRel90_pt	   = fs->make<TH1F>( "l1Ele_IsoRel90_pt"        , "p_{t}", 300,  0., 300. );
  l1Ele_IsoRel90_time_pt   = fs->make<TH1F>( "l1Ele_IsoRel90_time_pt"   , "p_{t}", 300,  0., 300. );
  l1Ele_IsoRel95_pt	   = fs->make<TH1F>( "l1Ele_IsoRel95_pt"        , "p_{t}", 300,  0., 300. );
  l1Ele_IsoRel95_time_pt   = fs->make<TH1F>( "l1Ele_IsoRel95_time_pt"   , "p_{t}", 300,  0., 300. );
  
  //photon rate		  
  l1Gamma_pt		   = fs->make<TH1F>( "l1Gamma_pt"               , "p_{t}", 300,  0., 300. );
  l1Gamma_IsoRel90_pt	   = fs->make<TH1F>( "l1Gamma_IsoRel90_pt"      , "p_{t}", 300,  0., 300. );
  l1Gamma_IsoRel90_time_pt = fs->make<TH1F>( "l1Gamma_IsoRel90_time_pt" , "p_{t}", 300,  0., 300. );
  l1Gamma_IsoRel95_pt	   = fs->make<TH1F>( "l1Gamma_IsoRel95_pt"      , "p_{t}", 300,  0., 300. );
  l1Gamma_IsoRel95_time_pt = fs->make<TH1F>( "l1Gamma_IsoRel95_time_pt" , "p_{t}", 300,  0., 300. );
  
   //muon rate			  
  l1Mu_pt		   = fs->make<TH1F>( "l1Mu_pt"                  , "p_{t}", 300,  0., 300. );
  l1Mu_IsoRel90_pt	   = fs->make<TH1F>( "l1Mu_IsoRel90_pt"         , "p_{t}", 300,  0., 300. );
  l1Mu_IsoRel90_time_pt    = fs->make<TH1F>( "l1Mu_IsoRel90_time_pt"    , "p_{t}", 300,  0., 300. );
  l1Mu_IsoRel95_pt	   = fs->make<TH1F>( "l1Mu_IsoRel95_pt"         , "p_{t}", 300,  0., 300. );
  l1Mu_IsoRel95_time_pt    = fs->make<TH1F>( "l1Mu_IsoRel95_time_pt"    , "p_{t}", 300,  0., 300. );
   
   //tau rate			  
   l1Tau_pt		   = fs->make<TH1F>( "l1Tau_pt"                 , "p_{t}", 300,  0., 300. );
   l1Tau_IsoRel90_pt	   = fs->make<TH1F>( "l1Tau_IsoRel90_pt"        , "p_{t}", 300,  0., 300. );
   l1Tau_IsoRel90_time_pt  = fs->make<TH1F>( "l1Tau_IsoRel90_time_pt"   , "p_{t}", 300,  0., 300. );
   l1Tau_IsoRel95_pt	   = fs->make<TH1F>( "l1Tau_IsoRel95_pt"        , "p_{t}", 300,  0., 300. );
   l1Tau_IsoRel95_time_pt  = fs->make<TH1F>( "l1Tau_IsoRel95_time_pt"   , "p_{t}", 300,  0., 300. );

  
}

//destructor
mtdIsoAnalyzerRates::~mtdIsoAnalyzerRates()
{
}

//deltaR, deltaZ, time_Cut 
void 
mtdIsoAnalyzerRates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::Handle<edm::ValueMap<float> > timingValues;
  iEvent.getByToken(timingValuesToken_,timingValues);

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

  // sort by pt not necessary for now... perhaps needed for an emulator
  //std::sort(l1PFCandidates.begin(), l1PFCandidates.end(), [](l1t::PFCandidate i,l1t::PFCandidate j){return(i.pt() > j.pt());});   

  //Fill Each L1 PF Cand vector with the L1 PFCands Electrons, Photons and Muons
  fillL1PFCandCollections(l1PFCandidates, l1PFElectrons, l1PFElectrons_indicies, l1PFPhotons, l1PFPhotons_indicies, l1PFMuons, l1PFMuons_indicies);
  
  //Generate electron isolation
  std::vector<simple_object> l1Electrons;
  std::vector<simple_object> l1Electrons90;
  std::vector<simple_object> l1Electrons90_time;
  std::vector<simple_object> l1Electrons95;
  std::vector<simple_object> l1Electrons95_time;

  int idx = 0;
  for(auto l1Cand : l1PFElectrons){

    int object_idx = l1PFElectrons_indicies.at(idx);
    idx++;
    simple_object temp;
    temp.pt       = l1Cand.pt();
    temp.eta      = l1Cand.eta();
    double isoSum = 0;
    double isoSumTime = 0;    

    //calculate iso Sum
    calculate_charged_iso_sum(l1Cand, object_idx, l1PFCandidates, timingValues, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);

    double relativeIso      = isoSum/l1Cand.pt();
    double relativeIso_time = isoSumTime/l1Cand.pt();
    
    temp.relativeIso      = relativeIso;
    temp.relativeIso_time = relativeIso_time;

    l1Electrons.push_back(temp);

    //90% efficiency WP
    if(relativeIso<WP_90_ele_)
      l1Electrons90.push_back(temp);
    if(relativeIso_time<WP_time_90_ele_)
      l1Electrons90_time.push_back(temp);

    //95% efficiency WP
    if(relativeIso<WP_95_ele_)
      l1Electrons95.push_back(temp);
    if(relativeIso_time<WP_time_95_ele_)
      l1Electrons95_time.push_back(temp);
  }

  //sort by pt as we trigger on the highest pt candidate
  std::sort(l1Electrons90.begin(),      l1Electrons90.end(),      [](simple_object i,simple_object j){return(i.pt > j.pt);});   
  std::sort(l1Electrons95.begin(),      l1Electrons95.end(),      [](simple_object i,simple_object j){return(i.pt > j.pt);});   
  std::sort(l1Electrons90_time.begin(), l1Electrons90_time.end(), [](simple_object i,simple_object j){return(i.pt > j.pt);});   
  std::sort(l1Electrons95_time.begin(), l1Electrons95_time.end(), [](simple_object i,simple_object j){return(i.pt > j.pt);});   

  //fill the rates
  
  //std::cout<<"l1Electrons.size() "<<l1Electrons.size()<<std::endl;
  //std::cout<<"l1Electrons.at(0) pt "<<l1Electrons.at(0).pt<<std::endl;
  if(l1Electrons.size()>0)
    l1Ele_pt->Fill( l1Electrons.at(0).pt );

  if(l1Electrons90.size()>0)
    l1Ele_IsoRel90_pt->Fill( l1Electrons90.at(0).pt) ;

  if(l1Electrons95.size()>0)
    l1Ele_IsoRel95_pt->Fill( l1Electrons95.at(0).pt );

  if(l1Electrons90_time.size()>0)
    l1Ele_IsoRel90_time_pt->Fill( l1Electrons90_time.at(0).pt );

  if(l1Electrons95_time.size()>0)
    l1Ele_IsoRel95_time_pt->Fill( l1Electrons95_time.at(0).pt );


  //// Muons    
  //Generate muon isolation
  std::vector<simple_object> l1Muons;
  std::vector<simple_object> l1Muons90;
  std::vector<simple_object> l1Muons90_time;
  std::vector<simple_object> l1Muons95;
  std::vector<simple_object> l1Muons95_time;

  idx = 0;
  for(auto l1Cand : l1PFMuons){

    double isoSum     = 0;
    double isoSumTime = 0;    
    int object_idx    = l1PFMuons_indicies.at(idx);
    idx++;
    simple_object temp;

    temp.pt           = l1Cand.pt();
    temp.eta          = l1Cand.eta();

    //calculate iso Sum
    calculate_charged_iso_sum( l1Cand, object_idx, l1PFCandidates, timingValues, time_cut_, iso_cone_deltaR_, iso_cone_deltaZ_, isoSum, isoSumTime);

    double relativeIso      = isoSum/l1Cand.pt();
    double relativeIso_time = isoSumTime/l1Cand.pt();
    
    temp.relativeIso        = relativeIso;
    temp.relativeIso_time   = relativeIso_time;

    l1Muons.push_back(temp);

    //90% efficiency WP
    if(relativeIso<WP_90_mu_){
      l1Muons90.push_back(temp);
    }
    if((relativeIso_time<WP_time_90_mu_ && l1Cand.pt()>15)||(relativeIso<WP_90_mu_ && l1Cand.pt()<15 )){
      l1Muons90_time.push_back(temp);
    }

    //95% efficiency WP
    if(relativeIso<WP_95_mu_)
      l1Muons95.push_back(temp);

    if((relativeIso_time<WP_time_95_mu_ && l1Cand.pt()>15)||(relativeIso<WP_95_mu_ && l1Cand.pt()<15 )){
      l1Muons95_time.push_back(temp);
    }
  }

  std::sort(l1Muons90.begin(),      l1Muons90.end(),      [](simple_object i,simple_object j){return(i.pt > j.pt);});   
  std::sort(l1Muons95.begin(),      l1Muons95.end(),      [](simple_object i,simple_object j){return(i.pt > j.pt);});   
  std::sort(l1Muons90_time.begin(), l1Muons90_time.end(), [](simple_object i,simple_object j){return(i.pt > j.pt);});   
  std::sort(l1Muons95_time.begin(), l1Muons95_time.end(), [](simple_object i,simple_object j){return(i.pt > j.pt);});   

  //fill the rates
  if(l1Muons.size()>0)
    l1Mu_pt->Fill( l1Muons.at(0).pt );

  if(l1Muons90.size()>0)
    l1Mu_IsoRel90_pt->Fill( l1Muons90.at(0).pt );

  if(l1Muons95.size()>0)
    l1Mu_IsoRel95_pt->Fill( l1Muons95.at(0).pt );

  if(l1Muons90_time.size()>0)
    l1Mu_IsoRel90_time_pt->Fill( l1Muons90_time.at(0).pt );

  if(l1Muons95_time.size()>0)
    l1Mu_IsoRel95_time_pt->Fill( l1Muons95_time.at(0).pt );

}


void mtdIsoAnalyzerRates::calculate_charged_iso_sum(l1t::PFCandidate object, int object_index, Handle< l1t::PFCandidateCollection > pfCandidates, edm::Handle<edm::ValueMap<float> > timingValues, double time_cut, double deltaR, double deltaZ, double &isoSum, double &isoSumTime){

    double iso_sum = 0;
    double iso_sum_time = 0;
    int idx = 0;

    for( PFCandidateCollection::const_iterator pfCand  = pfCandidates->begin();
                                               pfCand != pfCandidates->end(); 
	                                       ++pfCand, ++idx){    
      //don't double count the object
      if(idx == object_index)
	continue;
      
      //first check if it is a charged candidate
      if(pfCand->id() == l1t::PFCandidate::Muon || pfCand->id() == l1t::PFCandidate::Electron || pfCand->id() == l1t::PFCandidate::ChargedHadron){
	
	// first check deltaR match
	if(reco::deltaR(object.eta(), object.phi(), pfCand->eta(), pfCand->phi()) < deltaR){

	  if(object.pfTrack().isNull()||pfCand->pfTrack().isNull())
	    continue;
	  
	  if(fabs( object.pfTrack()->track()->getPOCA().z() - pfCand->pfTrack()->track()->getPOCA().z() ) < deltaZ){
	    
	    iso_sum += pfCand->pt();
	    float objectTime = 0;
	    float pfCandTime = 0;
	    
	    
	    if(!object.pfTrack().isNull()){
	      objectTime = (*timingValues)[object.pfTrack()->track()];
	    }
	    if(!pfCand->pfTrack().isNull()){
	      pfCandTime = (*timingValues)[pfCand->pfTrack()->track()];
	    }
	    
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


void mtdIsoAnalyzerRates::fillL1PFCandCollections(Handle<PFCandidateCollection> pfCandidates, PFCandidateCollection &pfElectrons, std::vector<int> & l1PFElectrons_indicies, PFCandidateCollection &pfPhotons, std::vector<int> & l1PFPhotons_indicies, PFCandidateCollection &pfMuons, std::vector<int> & l1PFMuons_indicies){

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


// ------------ method called once each job just before starting event loop  ------------
void
mtdIsoAnalyzerRates::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
mtdIsoAnalyzerRates::endJob()
{
}

void
mtdIsoAnalyzerRates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(mtdIsoAnalyzerRates);
