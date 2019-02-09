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
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "DataFormats/ForwardDetId/interface/MTDDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"
#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CLHEP/Units/GlobalPhysicalConstants.h"

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
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"

#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"

using namespace edm;
using namespace std;


struct MTDinfo {

  float sim_energy;
  float sim_time;
  float sim_x;
  float sim_y;
  float sim_z;

  uint32_t digi_row[2];
  uint32_t digi_col[2];
  uint32_t digi_charge[2];
  uint32_t digi_time1[2];
  uint32_t digi_time2[2];

  float ureco_charge[2];
  float ureco_time[2];

  float reco_energy;
  float reco_time;

};

struct {
  l1extra::L1JetParticle l1object;
  FTLCluster cluster;
} l1_time_match; 

class L1MTDLLPAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

  explicit L1MTDLLPAnalyzer(const edm::ParameterSet&);

  ~L1MTDLLPAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  void zeroTreeValues();
  void calculateRingValues(float eta, float time);
  void averageRings();

  const float btlMinEnergy_;

  typedef std::vector<reco::GenParticle> GenParticleCollectionType;

  const MTDGeometry* geom_; 
 
  EDGetTokenT<std::vector<reco::GenParticle> > genToken_;

  EDGetTokenT< BTLDigiCollection >    btlDigisToken_;
  EDGetTokenT< ETLDigiCollection >    etlDigisToken_;

  EDGetTokenT< FTLRecHitCollection >  btlRecHitToken_;
  EDGetTokenT< FTLRecHitCollection >  etlRecHitToken_;
  
  EDGetTokenT< FTLClusterCollection > btlClusterToken_;
  EDGetTokenT< FTLClusterCollection > etlClusterToken_;

  EDGetTokenT< vector<PSimHit> >      btlSimHitsToken_;
  EDGetTokenT<edm::HepMCProduct>      HepMCProductToken_;
  EDGetTokenT< vector<l1extra::L1JetParticle> > l1TausToken_;

  InputTag genSrc_;
  double time_cut_;


  const HepMC::GenParticle* GetFinal(const HepMC::GenParticle* p){ // includes mixing (assuming mixing is not occurring more than 5 times back and forth)
    const HepMC::GenParticle* aPart = p;
    for (unsigned int iMix = 0; iMix < 10; iMix++) {
      bool foundSimilar = false;
      if(aPart->end_vertex()){ 
	if(aPart->end_vertex()->particles_out_size()!=0){ 
	  for(HepMC::GenVertex::particles_out_const_iterator d=aPart->end_vertex()->particles_out_const_begin(); d!=aPart->end_vertex()->particles_out_const_end();d++){ 
	    if(abs((*d)->pdg_id())==abs(aPart->pdg_id())){ 
	      aPart = *d;         
	      foundSimilar = true;
	      break;
	    } 
	  } 
	}
	if (!foundSimilar) break;
      } 
    } 
    return aPart;
  }




  double genPt, genEta, genPhi, genEnergy, genPdgId, genStatus;
  double matchedClusterTime;
  double sumMatchedClusterTime;
  double etaRingAvgTime;
  double eventAvgTime;
  double genLLPLTau, genLLPLifetime, genMuLifetime, genPVTime, genSVTime, genSVLTau, genMuLTau, genMuLTau_noLLP, genMuLifetime_noLLP;
  double sumTimeRing0p4;
  double sumTimeRing0p7;
  double sumTimeRing0p8;
  double sumTimeRing0p9;
  double sumTimeRing1p0;
  double sumTimeRing1p1;
  double sumTimeRing1p2;
  double sumTimeRing1p3;
  double sumTimeRing1p4;
  double sumTimeRing1p5;

  double avgTimeRing0p4;
  double avgTimeRing0p7;
  double avgTimeRing0p8;
  double avgTimeRing0p9;
  double avgTimeRing1p0;
  double avgTimeRing1p1;
  double avgTimeRing1p2;
  double avgTimeRing1p3;
  double avgTimeRing1p4;
  double avgTimeRing1p5;

  double nHitsCluster;
  double nHitsRing0p4;
  double nHitsRing0p7;
  double nHitsRing0p8;
  double nHitsRing0p9;
  double nHitsRing1p0;
  double nHitsRing1p1;
  double nHitsRing1p2;
  double nHitsRing1p3;
  double nHitsRing1p4;
  double nHitsRing1p5;

  bool genIsPromptFinalState;
  double eta, phi;
  TTree* recHitTree;
  TTree* clusterTree;
  TTree* digiTree;
  TTree* genParticleTree;
  double time, isTimeValid, timeError, energy;
  double x;           
  double y;
  double z;
  double genLLP_m;       
  double genLLP_p;       	      
  double genLLP_rho;     
  double genLLP_beta;    
  double genLLP_SV_x;    
  double genLLP_SV_y;    
  double genLLP_SV_z;    
  double genLLP_PV_x;    
  double genLLP_PV_y;    
  double genLLP_PV_z;    
 
  double gen_mu_m;       
  double gen_mu_eta;     
  double gen_mu_phi;     
  double gen_mu_p;       
  double gen_mu_rho;     
  double gen_mu_beta;    
  double gen_mu_tof;
  double gen_mu_x;       
  double gen_mu_y;       
  double gen_mu_z;       

  //double time;        
  //double timeError;   
  double size;        
  double sizeX;       
  double sizeY;       
  double hitOffset;   
  double hitEnergy;   
  double hitTime;     
  double hitTimeError;

  std::unordered_map<uint32_t, MTDinfo> btl_hits;
  std::unordered_map<uint32_t, MTDinfo> etl_hits[2];
  
};

L1MTDLLPAnalyzer::L1MTDLLPAnalyzer(const edm::ParameterSet &cfg) :  
  btlMinEnergy_( cfg.getParameter<double>("BTLMinimumEnergy") ),
  geom_(nullptr),
  btlDigisToken_(      consumes< BTLDigiCollection >   ( cfg.getParameter<InputTag>("FTLBarrel"))),
  etlDigisToken_(      consumes< ETLDigiCollection >   ( cfg.getParameter<InputTag>("FTLEndcap"))),
  btlRecHitToken_(     consumes< FTLRecHitCollection > ( cfg.getParameter<InputTag>("recHitBarrel"))), // finish me
  etlRecHitToken_(     consumes< FTLRecHitCollection > ( cfg.getParameter<InputTag>("recHitEndcap"))), // finish me
  btlClusterToken_(    consumes< FTLClusterCollection > ( cfg.getParameter<InputTag>("mtdClusterBarrel"))),
  etlClusterToken_(    consumes< FTLClusterCollection > ( cfg.getParameter<InputTag>("mtdClusterEndcap"))),
  btlSimHitsToken_(    consumes< vector<PSimHit> >      ( cfg.getParameter<InputTag>("BTLSimHits"))),
  HepMCProductToken_(  consumes< edm::HepMCProduct >    ( cfg.getParameter<InputTag>("HepMCProduct"))),
  l1TausToken_(        consumes< vector<l1extra::L1JetParticle> >( cfg.getParameter<InputTag>("l1Taus"))),
  genSrc_ ((           cfg.getParameter<edm::InputTag>( "genParticles"))),
  //vector<l1extra::L1JetParticle>        "l1extraParticles"          "Tau"             "RECO"
  //genSrc_ ((           cfg.getParameter<edm::InputTag>( "genParticles"))),
  time_cut_(           cfg.getParameter<double>("time_cut"))
{
    //services
  usesResource("TFileService");
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
  Service<TFileService> fs;
  
  recHitTree = fs->make<TTree>("recHitTree","recHit Tree" );
  recHitTree->Branch("time",        &time,         "time/D");
  recHitTree->Branch("isTimeValid", &isTimeValid, "isTimeValid/D");
  recHitTree->Branch("timeError",   &timeError,   "timeError/D");
  recHitTree->Branch("energy",      &energy,       "energy/D");

  clusterTree = fs->make<TTree>("clusterTree","cluster Tree" );
  clusterTree->Branch("eta",           &eta,            "eta/D");           
  clusterTree->Branch("phi",           &phi,            "phi/D");           
  clusterTree->Branch("x",             &x,              "x/D");           
  clusterTree->Branch("y",             &y,           	"y/D");           
  clusterTree->Branch("time",          &time,        	"time/D");        
  clusterTree->Branch("timeError",     &timeError,   	"timeError/D");   
  clusterTree->Branch("size",          &size,        	"size/D");        
  clusterTree->Branch("sizeX",         &sizeX,       	"sizeX/D");       
  clusterTree->Branch("sizeY",         &sizeY,       	"sizeY/D");       
  clusterTree->Branch("energy",        &energy,      	"energy/D");      
  clusterTree->Branch("hitOffset",     &hitOffset,   	"hitOffset/D");   
  clusterTree->Branch("hitEnergy",     &hitEnergy,   	"hitEnergy/D");   
  clusterTree->Branch("hitTime",       &hitTime,     	"hitTime/D");     
  clusterTree->Branch("hitTimeError",  &hitTimeError,	"hitTimeError/D");

  digiTree = fs->make<TTree>("digiTree","digi Tree" );
  digiTree->Branch("eta",           &eta,            "eta/D");           
  digiTree->Branch("phi",           &phi,            "phi/D");           
  digiTree->Branch("x",             &x,              "x/D");           
  digiTree->Branch("y",             &y,           	"y/D");           
  digiTree->Branch("time",          &time,        	"time/D");        
  digiTree->Branch("timeError",     &timeError,   	"timeError/D");   
  digiTree->Branch("size",          &size,        	"size/D");        
  digiTree->Branch("sizeX",         &sizeX,       	"sizeX/D");       
  digiTree->Branch("sizeY",         &sizeY,       	"sizeY/D");       
  digiTree->Branch("energy",        &energy,      	"energy/D");      
  digiTree->Branch("hitOffset",     &hitOffset,   	"hitOffset/D");   
  digiTree->Branch("hitEnergy",     &hitEnergy,   	"hitEnergy/D");   
  digiTree->Branch("hitTime",       &hitTime,     	"hitTime/D");     
  digiTree->Branch("hitTimeError",  &hitTimeError,	"hitTimeError/D");

  genParticleTree = fs->make<TTree>("genParticleTree","generated particle Tree" );
  genParticleTree->Branch("gen_pt",                 &genPt,               "gen_pt/D");           
  genParticleTree->Branch("gen_eta",                &genEta,              "gen_eta/D");           
  genParticleTree->Branch("gen_phi",                &genPhi,              "gen_phi/D");           
  genParticleTree->Branch("gen_Energy",             &genEnergy,           "gen_Energy/D");           
  genParticleTree->Branch("gen_pdgId",              &genPdgId,            "gen_pdgId/D"); 
          
  genParticleTree->Branch("gen_LLP_Lifetime",       &genLLPLifetime,      "gen_LLP_Lifetime/D");
  genParticleTree->Branch("gen_LLP_LTau",           &genLLPLTau,          "gen_LLP_LTau/D");

  genParticleTree->Branch("gen_LLP_m",              &genLLP_m,            "gen_LLP_m/D");
  genParticleTree->Branch("gen_LLP_p",              &genLLP_p,            "gen_LLP_p/D");
  genParticleTree->Branch("gen_LLP_rho",            &genLLP_rho,          "gen_LLP_rho/D");
  genParticleTree->Branch("gen_LLP_beta",           &genLLP_beta,         "gen_LLP_beta/D");

  genParticleTree->Branch("gen_LLP_SV_x",           &genLLP_SV_x,         "gen_LLP_SV_x/D");
  genParticleTree->Branch("gen_LLP_SV_y",           &genLLP_SV_y,         "gen_LLP_SV_y/D");
  genParticleTree->Branch("gen_LLP_SV_z",           &genLLP_SV_z,         "gen_LLP_SV_z/D");

  genParticleTree->Branch("gen_LLP_PV_x",           &genLLP_PV_x,         "gen_LLP_PV_x/D");
  genParticleTree->Branch("gen_LLP_PV_y",           &genLLP_PV_y,         "gen_LLP_PV_y/D");
  genParticleTree->Branch("gen_LLP_PV_z",           &genLLP_PV_z,         "gen_LLP_PV_z/D");

  genParticleTree->Branch("gen_mu_Lifetime",        &genMuLifetime,       "gen_mu_Lifetime/D");
  genParticleTree->Branch("gen_mu_Ltau",            &genMuLTau,           "gen_mu_LTau/D");

  genParticleTree->Branch("gen_mu_m",               &gen_mu_m,            "gen_mu_m/D");
  genParticleTree->Branch("gen_mu_eta",             &gen_mu_eta,          "gen_mu_eta/D");
  genParticleTree->Branch("gen_mu_phi",             &gen_mu_phi,          "gen_mu_phi/D");
  genParticleTree->Branch("gen_mu_p",               &gen_mu_p,            "gen_mu_p/D");
  genParticleTree->Branch("gen_mu_rho",             &gen_mu_rho,          "gen_mu_rho/D");
  genParticleTree->Branch("gen_mu_beta",            &gen_mu_beta,         "gen_mu_beta/D");
  genParticleTree->Branch("gen_mu_tof",             &gen_mu_tof,          "gen_mu_tof/D");  

  genParticleTree->Branch("gen_mu_x",               &gen_mu_x,            "gen_mu_x/D");
  genParticleTree->Branch("gen_mu_y",               &gen_mu_y,            "gen_mu_y/D");
  genParticleTree->Branch("gen_mu_z",               &gen_mu_z,            "gen_mu_z/D");

  genParticleTree->Branch("gen_mu_Lifetime_noLLP",  &genMuLifetime_noLLP, "gen_mu_Lifetime_noLLP/D");
  genParticleTree->Branch("gen_mu_Ltau_noLLP",      &genMuLTau_noLLP,     "gen_mu_Ltau_noLLP/D");

  //genParticleTree->Branch("gen_SV_time",            &genSVTime,        "gen_SV_time/D");
  //genParticleTree->Branch("gen_SV_Ltau",            &genSVLTau,        "gen_SV_LTau/D");
  //genParticleTree->Branch("gen_status",             &genStatus,         "gen_status/D");           
  //genParticleTree->Branch("gen_isPromptFinalState", &genIsPromptFinalState,      "gen_isPromptFinalState/b");           
  genParticleTree->Branch("matchedClusterTime",     &matchedClusterTime, "matchedClusterTime/D");
  genParticleTree->Branch("sumMatchedClusterTime",     &sumMatchedClusterTime, "sumMatchedClusterTime/D");
  genParticleTree->Branch("eventAvgTime",     &eventAvgTime,     "eventAvgTime/D" );

  genParticleTree->Branch("sumTimeRing0p4", &sumTimeRing0p4 ,"sumTimeRing0p4/D" );
  genParticleTree->Branch("sumTimeRing0p7", &sumTimeRing0p7 ,"sumTimeRing0p7/D" );
  genParticleTree->Branch("sumTimeRing0p8", &sumTimeRing0p8 ,"sumTimeRing0p8/D" );
  genParticleTree->Branch("sumTimeRing0p9", &sumTimeRing0p9 ,"sumTimeRing0p9/D" );
  genParticleTree->Branch("sumTimeRing1p0", &sumTimeRing1p0 ,"sumTimeRing1p0/D" );
  genParticleTree->Branch("sumTimeRing1p1", &sumTimeRing1p1 ,"sumTimeRing1p1/D" );
  genParticleTree->Branch("sumTimeRing1p2", &sumTimeRing1p2 ,"sumTimeRing1p2/D" );
  genParticleTree->Branch("sumTimeRing1p3", &sumTimeRing1p3 ,"sumTimeRing1p3/D" );
  genParticleTree->Branch("sumTimeRing1p4", &sumTimeRing1p4 ,"sumTimeRing1p4/D" );
  genParticleTree->Branch("sumTimeRing1p5", &sumTimeRing1p5 ,"sumTimeRing1p5/D" );

  genParticleTree->Branch("avgTimeRing0p4", &avgTimeRing0p4 ,"avgTimeRing0p4/D" );
  genParticleTree->Branch("avgTimeRing0p7", &avgTimeRing0p7 ,"avgTimeRing0p7/D" );
  genParticleTree->Branch("avgTimeRing0p8", &avgTimeRing0p8 ,"avgTimeRing0p8/D" );
  genParticleTree->Branch("avgTimeRing0p9", &avgTimeRing0p9 ,"avgTimeRing0p9/D" );
  genParticleTree->Branch("avgTimeRing1p0", &avgTimeRing1p0 ,"avgTimeRing1p0/D" );
  genParticleTree->Branch("avgTimeRing1p1", &avgTimeRing1p1 ,"avgTimeRing1p1/D" );
  genParticleTree->Branch("avgTimeRing1p2", &avgTimeRing1p2 ,"avgTimeRing1p2/D" );
  genParticleTree->Branch("avgTimeRing1p3", &avgTimeRing1p3 ,"avgTimeRing1p3/D" );
  genParticleTree->Branch("avgTimeRing1p4", &avgTimeRing1p4 ,"avgTimeRing1p4/D" );
  genParticleTree->Branch("avgTimeRing1p5", &avgTimeRing1p5 ,"avgTimeRing1p5/D" );

  genParticleTree->Branch("nHitsCluster", &nHitsCluster ,"nHitsCluster/D" );
  genParticleTree->Branch("nHitsRing0p4", &nHitsRing0p4 ,"nHitsRing0p4/D" );
  genParticleTree->Branch("nHitsRing0p7", &nHitsRing0p7 ,"nHitsRing0p7/D" );
  genParticleTree->Branch("nHitsRing0p8", &nHitsRing0p8 ,"nHitsRing0p8/D" );
  genParticleTree->Branch("nHitsRing0p9", &nHitsRing0p9 ,"nHitsRing0p9/D" );
  genParticleTree->Branch("nHitsRing1p0", &nHitsRing1p0 ,"nHitsRing1p0/D" );
  genParticleTree->Branch("nHitsRing1p1", &nHitsRing1p1 ,"nHitsRing1p1/D" );
  genParticleTree->Branch("nHitsRing1p2", &nHitsRing1p2 ,"nHitsRing1p2/D" );
  genParticleTree->Branch("nHitsRing1p3", &nHitsRing1p3 ,"nHitsRing1p3/D" );
  genParticleTree->Branch("nHitsRing1p4", &nHitsRing1p4 ,"nHitsRing1p4/D" );
  genParticleTree->Branch("nHitsRing1p5", &nHitsRing1p5 ,"nHitsRing1p5/D" );
    
}
//destructor
L1MTDLLPAnalyzer::~L1MTDLLPAnalyzer()
{
}


void 
L1MTDLLPAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::ESHandle<MTDGeometry> geom;
  iSetup.get<MTDDigiGeometryRecord>().get(geom);
  geom_ = geom.product();
  //std::cout<<"got geometry"<<std::endl;
   

  edm::Handle<BTLDigiCollection> h_BTL_digi;
  iEvent.getByToken(btlDigisToken_,h_BTL_digi);

  edm::Handle<ETLDigiCollection> etlDigis;
  iEvent.getByToken(etlDigisToken_,etlDigis);

  edm::Handle< FTLRecHitCollection > btlRecHits;
  iEvent.getByToken(btlRecHitToken_,btlRecHits);

  edm::Handle< FTLRecHitCollection > etlRecHits;
  iEvent.getByToken(etlRecHitToken_,etlRecHits);

  edm::Handle< FTLClusterCollection > btlClusters;
  iEvent.getByToken(btlClusterToken_,btlClusters);

  edm::Handle< FTLClusterCollection > etlClusters;
  iEvent.getByToken(etlClusterToken_,etlClusters);

  edm::Handle< std::vector<PSimHit> > btlSimHits;
  iEvent.getByToken(btlSimHitsToken_,btlSimHits);
  
  edm::Handle< vector<l1extra::L1JetParticle> > l1Taus;
  iEvent.getByToken(l1TausToken_, l1Taus);


  if( btlRecHits->size() > 0 ){
    for(const auto& recHit : *btlRecHits){

      time        = (double)recHit.time();
      isTimeValid = (double)recHit.isTimeValid();
      timeError   = (double)recHit.timeError();;
      energy      = (double)recHit.energy();
      recHitTree->Fill();
    }
  }

  if( btlClusters->size() > 0 ){
    for(const auto& detIds : *btlClusters){
      for(const auto& cluster : detIds){
	//eta          = (double)cluster.eta();
	//phi          = (double)cluster.phi();
	x            = (double)cluster.x();
	y            = (double)cluster.y();
	time         = (double)cluster.time();
	timeError    = (double)cluster.timeError();
	size         = (double)cluster.size();
	sizeX        = (double)cluster.sizeX();
	sizeY        = (double)cluster.sizeY();
	energy       = (double)cluster.energy();
	hitOffset    = (double)cluster.hitOffset().at(0);
	hitEnergy    = (double)cluster.hitENERGY().at(0);
	hitTime      = (double)cluster.hitTIME().at(0);
	hitTimeError = (double)cluster.hitTIME_ERROR().at(0);
	clusterTree->Fill();
      }
    }
  }

  for(const auto& l1object : *l1Taus){
    std::cout<<"x: "<<l1object.p4().x()<<std::endl;
    std::cout<<"y: "<<l1object.p4().y()<<std::endl;
  }

  unsigned int n_digi_btl[2] = {0,0};

  btl_hits.clear();
  
  if (h_BTL_digi->size() > 0 ) {

    //std::cout << " ----------------------------------------" << std::endl;
    //std::cout << " BTL DIGI collection:" << std::endl;
    
    for (const auto& dataFrame: *h_BTL_digi) {
      // in case print outs are needed for debugging this can be uncommented
      /*
      // --- detector element ID:
      std::cout << "   det ID:  det = " << dataFrame.id().det() 
		<< "  subdet = "  << dataFrame.id().mtdSubDetector() 
		<< "  side = "    << dataFrame.id().mtdSide() 
		<< "  rod = "     << dataFrame.id().mtdRR() 
		<< "  mod = "     << dataFrame.id().module() 
		<< "  type = "    << dataFrame.id().modType() 
		<< "  crystal = " << dataFrame.id().crystal() 
		<< std::endl;


      // --- loop over the dataFrame samples
      for (int isample = 0; isample<dataFrame.size(); ++isample){

	const auto& sample = dataFrame.sample(isample);

	std::cout << "       sample " << isample << ":"; 
	if ( sample.data()==0 && sample.toa()==0 ) {
	  std::cout << std::endl;
	  continue;
	}
	std::cout << "  amplitude = " << sample.data() 
		  << "  time1 = " <<  sample.toa() 
		  << "  time2 = " <<  sample.toa2() << std::endl;
      */


      DetId id =  dataFrame.id();
      
      const auto& sample_L = dataFrame.sample(0);
      const auto& sample_R = dataFrame.sample(1);

      btl_hits[id.rawId()].digi_row[0] = sample_L.row();
      btl_hits[id.rawId()].digi_row[1] = sample_R.row();
      btl_hits[id.rawId()].digi_col[0] = sample_L.column();
      btl_hits[id.rawId()].digi_col[1] = sample_R.column();

      btl_hits[id.rawId()].digi_charge[0] = sample_L.data();
      btl_hits[id.rawId()].digi_charge[1] = sample_R.data();
      btl_hits[id.rawId()].digi_time1[0]  = sample_L.toa();
      btl_hits[id.rawId()].digi_time1[1]  = sample_R.toa();
      btl_hits[id.rawId()].digi_time2[0]  = sample_L.toa2();
      btl_hits[id.rawId()].digi_time2[1]  = sample_R.toa2();

      if ( sample_L.data() > 0 )
	n_digi_btl[0]++;

      if ( sample_R.data() > 0 )
	n_digi_btl[1]++;
      
      
    } // digi loop

  } // if (h_BTL_digi->size() > 0 )

  unsigned int n_reco_btl = 0;

  if (btlRecHits->size() > 0 ) {

    for (const auto& recHit: *btlRecHits) {

      DetId id = recHit.id();

      btl_hits[id.rawId()].reco_energy = recHit.energy();
      btl_hits[id.rawId()].reco_time   = recHit.time();

      if ( recHit.energy() > 0. )
	n_reco_btl++;


    } // recHit loop

  } // if ( h_BTL_reco->size() > 0 )

  Local3DPoint muonEntry;
  PSimHit simHit_muon;
  for(const auto& simHit : *btlSimHits){
    if(simHit.particleType()==13){
      //std::cout<<"Sim Muon Hit time of flight: "<<simHit.timeOfFlight()<<std::endl;
      muonEntry   = simHit.entryPoint();
      gen_mu_tof  = simHit.timeOfFlight();
      simHit_muon = simHit;
    }
  }
 
  //// hep mc info
   edm::Handle<edm::HepMCProduct> HepMCHandle;
   iEvent.getByToken( HepMCProductToken_, HepMCHandle) ;
   const HepMC::GenEvent* Evt = HepMCHandle->GetEvent() ;

   int i = 0;
   //float timeGenLLP = 0;
   TVector3 MuonVertex;

   for ( HepMC::GenEvent::vertex_const_iterator itVtx=Evt->vertices_begin(); itVtx!=Evt->vertices_end(); ++itVtx ){
     i++;

     for ( HepMC::GenVertex::particles_out_const_iterator itPartOut=(*itVtx)->particles_out_const_begin(); itPartOut!=(*itVtx)->particles_out_const_end(); ++itPartOut ) {
       if((*itPartOut)->pdg_id()==1000013){

	 //std::cout<< "pdgID: "<<(*itPartOut)->pdg_id() <<" mass: "<< (*itPartOut)->momentum().m()<< " px: "<<(*itPartOut)->momentum.px() <<" time: "<<time<<std::endl;
	 const HepMC::GenParticle* pf = GetFinal(*itPartOut); // inlcude mixing

	 if((*itPartOut)->production_vertex() && (*itPartOut)->end_vertex()){

	   TVector3 PV((*itPartOut)->production_vertex()->point3d().x(), (*itPartOut)->production_vertex()->point3d().y(), (*itPartOut)->production_vertex()->point3d().z()); 
	   TVector3 SV((*itPartOut)->end_vertex()->point3d().x(),(*itPartOut)->end_vertex()->point3d().y(),(*itPartOut)->end_vertex()->point3d().z()); 
	   
	   TVector3 DL=SV-PV; 
	   
	   double c(2.99792458E8);
	   double Ltau(DL.Mag()/100)/*cm->m*/;
	   double beta((*itPartOut)->momentum().rho()/(*itPartOut)->momentum().m()); 
	   double lt=Ltau/(c*beta);

	   if(pf->end_vertex()){
	     TVector3 SVf(pf->end_vertex()->point3d().x(),pf->end_vertex()->point3d().y(),pf->end_vertex()->point3d().z());
	     DL=SVf-PV;
	     Ltau=DL.Mag()/100;
	     lt=Ltau/(c*beta);
	     
	     //std::cout<<"LTau "<<Ltau<<" lifetime "<<lt<< " lifetime, no beta: "<<Ltau/(c) <<std::endl;
	     //timeGenLLP = lt;
	     //if(lt>1E-16)lifetime_final->Fill(log10(lt),weight);
	     
	     if(!(muonEntry.x()==0 && muonEntry.y()==0 && muonEntry.z()==0 )){
	       TVector3 muon(muonEntry.x(), muonEntry.y(), muonEntry.z());
	       TVector3 dl        = muon-SV;
	       double   Ltau_muon = dl.Mag()/100;
	       double   lt_muon   = Ltau_muon/c;
	       
	       TVector3 dl_noLLP      = muon-PV;
	       double Ltau_muon_noLLP = dl_noLLP.Mag()/100;
	       double lt_muon_noLLP   = Ltau_muon_noLLP/c;
	       //std::cout<<"Ltau_muon: "<<Ltau_muon<<" lt_muon: "<<lt_muon<< " LTau_muon_noLLP: "<< Ltau_muon_noLLP<< " lt_muon_noLLP: "<< lt_muon_noLLP<< " Time Offset = "<< (lt+lt_muon)-lt_muon_noLLP <<std::endl;
	       
	       genLLPLifetime  = lt;
	       genLLPLTau      = Ltau;
	       genLLP_m        = (*itPartOut)->momentum().m();
	       genLLP_p        = sqrt((*itPartOut)->momentum().px()*(*itPartOut)->momentum().px() + 
				      (*itPartOut)->momentum().py()*(*itPartOut)->momentum().py() + 
				      (*itPartOut)->momentum().pz()*(*itPartOut)->momentum().pz() );
	       genLLP_rho      = (*itPartOut)->momentum().rho();
	       genLLP_beta     = (*itPartOut)->momentum().rho()/(*itPartOut)->momentum().m();
	       genLLP_SV_x     = (*itPartOut)->end_vertex()->point3d().x();
	       genLLP_SV_y     = (*itPartOut)->end_vertex()->point3d().y();
	       genLLP_SV_z     = (*itPartOut)->end_vertex()->point3d().z();
	       genLLP_PV_x     = (*itPartOut)->production_vertex()->point3d().x();
	       genLLP_PV_y     = (*itPartOut)->production_vertex()->point3d().y();
	       genLLP_PV_z     = (*itPartOut)->production_vertex()->point3d().z();
	       
	       genMuLifetime   = lt_muon; 
	       genMuLTau       = Ltau_muon;
	       
	       gen_mu_x        = muonEntry.x();
	       gen_mu_y        = muonEntry.y();
	       gen_mu_z        = muonEntry.z();
	       
	       genMuLifetime_noLLP = lt_muon_noLLP;
	       genMuLTau_noLLP = Ltau_muon_noLLP;
	     }
	   }
	 }
       }
     }
   }

   
   
   //// end hep mc info
   
   edm::Handle<GenParticleCollectionType> genParticles_;

   if(!iEvent.getByToken(genToken_,genParticles_))
     std::cout<<"No gen Particles Found "<<std::endl;
    
   GenParticleCollectionType motherParticles;
   int pdgID_mother = 1000013;
   std::vector<const reco::Candidate*> theMuons;

   // std::cout<<"Number of Daughters for PDGID "<<pdgID_mother<<" : "<<genParticles_->at(0).daughter(0)->pdgId()<<std::endl;   

   for(unsigned int i = 0; i < genParticles_->size(); i++){
     reco::GenParticle genParticle = genParticles_->at(i);

     if(genParticle.pdgId()==pdgID_mother){
       motherParticles.push_back(genParticle);

       if(genParticle.numberOfDaughters()>1){
	  if(genParticles_->at(i).daughter(1)->pdgId()==13){
	    theMuons.push_back(genParticles_->at(i).daughter(1));
	  }

	 //std::cout<<"Second Daughter PDGID: "<<genParticles_->at(i).daughter(1)->pdgId()<< " pt: "<< genParticles_->at(i).daughter(1)->pt() <<" eta: "<< genParticles_->at(i).daughter(1)->eta() <<" phi: "<< genParticles_->at(i).daughter(1)->phi() <<std::endl;
       }
       //std::cout<<"Number of Daughters for PDGID "<<pdgID_mother<<" : "<<genParticle.numberOfDaughters()<<" First Daughter PDGID: "<<genParticles_->at(i).daughter(0)->pdgId()<<std::endl;
     }
   }
   
   zeroTreeValues();

   
   for (auto const& hit: btl_hits) {

      BTLDetId detId(hit.first); 
      DetId geoId = BTLDetId(detId.mtdSide(),detId.mtdRR(),detId.module()+14*(detId.modType()-1),0,1);
      // this seg faults from time to time for specific events! 
      // it appears to originate in the external MTD geometry function...
      // this should be solved with Lindsey and other developers.
      // this portion is needed to convert the detID to eta/phi

      const MTDGeomDet* thedet           = geom_->idToDet(geoId);
      const ProxyMTDTopology& topoproxy  = static_cast<const ProxyMTDTopology&>(thedet->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

      if ( (hit.second).reco_energy < btlMinEnergy_ ) continue;

      Local3DPoint simscaled( 0.1*(hit.second).sim_x, 
			      0.1*(hit.second).sim_y, 
			      0.1*(hit.second).sim_z );

      simscaled              = topo.pixelToModuleLocalPoint( simscaled, detId.row( topo.nrows()), detId.column(topo.nrows()) );

      // Finally get an object in terms of eta, phi
      const auto& global_pos = thedet->toGlobal(simscaled);

      eta          = (double)global_pos.eta();
      phi          = (double)global_pos.phi();
      x            = (double)global_pos.x();
      y            = (double)global_pos.y();
      z            = (double)global_pos.z();
      time         = (double)(hit.second).reco_time;
      energy       = (double)(hit.second).reco_energy;

      for(auto theMuon : theMuons){

	// 0.02 seems to be a reasonable matching parameter, make this configurable
	// it appears that many 
	if(fabs(theMuon->eta()-eta)<0.02&&fabs(theMuon->phi()-phi)<0.02 && time<19){
	  sumMatchedClusterTime = time;
	  nHitsCluster          = 1;
	}
      }

      calculateRingValues(eta, time);

      digiTree->Fill();
    }

    for(auto mu : theMuons){
      genPt      = mu->pt();    
      genEta     = mu->eta();
      genPhi     = mu->phi();   
      genEnergy  = mu->energy();
      genPdgId   = 13; 
     
      gen_mu_m          = mu->mass();
      gen_mu_eta         = mu->eta();
      gen_mu_phi         = mu->phi();
      gen_mu_p           = mu->p();
      gen_mu_rho         = mu->momentum().rho();
      gen_mu_beta        = sqrt(mu->vx()*mu->vx() + mu->vy()*mu->vy() + mu->vz()*mu->vz())/2.99792458E8;
      gen_mu_tof         = simHit_muon.timeOfFlight();
      matchedClusterTime = sumMatchedClusterTime/nHitsCluster;
          
    }

    averageRings();

    genParticleTree->Fill();

}



//
void 
L1MTDLLPAnalyzer::zeroTreeValues(){
  genPt = -10;    
  genEta = -10;
  genPhi = -10;   
  genEnergy = -10;
  genPdgId = -10; 
  matchedClusterTime = -10;
  eventAvgTime = -10;
  
  sumMatchedClusterTime = 0;
  nHitsCluster = 0;
  
  sumTimeRing0p4 = -10;
  sumTimeRing0p7 = -10;
  sumTimeRing0p8 = -10;
  sumTimeRing0p9 = -10;
  sumTimeRing1p0 = -10;
  sumTimeRing1p1 = -10;
  sumTimeRing1p2 = -10;
  sumTimeRing1p3 = -10;
  sumTimeRing1p4 = -10;
  sumTimeRing1p5 = -10;
  
  nHitsRing0p4 = 0;
  nHitsRing0p7 = 0;
  nHitsRing0p8 = 0;
  nHitsRing0p9 = 0;
  nHitsRing1p0 = 0;
  nHitsRing1p1 = 0;
  nHitsRing1p2 = 0;
  nHitsRing1p3 = 0;
  nHitsRing1p4 = 0;
  nHitsRing1p5 = 0;
  
}

//
void 
L1MTDLLPAnalyzer::calculateRingValues(float eta, float time){

      //ring 0
      if(fabs(eta)<0.4){
	sumTimeRing0p4 += time;
	nHitsRing0p4++;
      }

      //ring 1
      if(fabs(eta)>0.4&&fabs(eta)<0.7){
	sumTimeRing0p7 += time;
	nHitsRing0p7++;
      }

      //ring 2
      if(fabs(eta)>0.7&&fabs(eta)<0.8){
	sumTimeRing0p8 += time;
	nHitsRing0p8++;
      }

      if(fabs(eta)>0.8&&fabs(eta)<0.9){
	sumTimeRing0p9 += time;
	nHitsRing0p9++;
      }

      if(fabs(eta)>0.9&&fabs(eta)<1.0){
	sumTimeRing1p0 += time;
	nHitsRing1p0++;
      }

      if(fabs(eta)>1.0&&fabs(eta)<1.1){
	sumTimeRing1p1 += time;
	nHitsRing1p1++;
      }

      if(fabs(eta)>1.1&&fabs(eta)<1.2){
	sumTimeRing1p2 += time;
	nHitsRing1p2++;
      }

      if(fabs(eta)>1.2&&fabs(eta)<1.3){
	sumTimeRing1p3 += time;
	nHitsRing1p3++;
      }

      if(fabs(eta)>1.3&&fabs(eta)<1.4){
	sumTimeRing1p4 += time;
	nHitsRing1p4++;
      }

      if(fabs(eta)>1.4&&fabs(eta)<1.5){
	sumTimeRing1p5 += time;
	nHitsRing1p5++;
      }
}

void
L1MTDLLPAnalyzer::averageRings(){

  if(nHitsRing0p4>0)
    avgTimeRing0p4 = sumTimeRing0p4/nHitsRing0p4;
  if(nHitsRing0p7>0)
    avgTimeRing0p7 = sumTimeRing0p7/nHitsRing0p7;
  if(nHitsRing0p8>0)
    avgTimeRing0p8 = sumTimeRing0p8/nHitsRing0p8;
  if(nHitsRing0p9>0)
    avgTimeRing0p9 = sumTimeRing0p9/nHitsRing0p9;
  
  if(nHitsRing1p0>0)
    avgTimeRing1p0 = sumTimeRing1p0/nHitsRing1p0;
  if(nHitsRing1p1>0)
    avgTimeRing1p1 = sumTimeRing1p1/nHitsRing1p1;
  if(nHitsRing1p2>0)
    avgTimeRing1p2 = sumTimeRing1p2/nHitsRing1p2;
  if(nHitsRing1p3>0)
    avgTimeRing1p3 = sumTimeRing1p3/nHitsRing1p3;
  if(nHitsRing1p4>0)
    avgTimeRing1p4 = sumTimeRing1p4/nHitsRing1p4;
  if(nHitsRing1p5>0)
    avgTimeRing1p5 = sumTimeRing1p5/nHitsRing1p5;
  
}

// ------------ method called once each job just before starting event loop  ------------
void
L1MTDLLPAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1MTDLLPAnalyzer::endJob()
{
}

void
L1MTDLLPAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(L1MTDLLPAnalyzer);
