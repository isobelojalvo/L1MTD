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

  const float btlMinEnergy_;

  typedef std::vector<reco::GenParticle> GenParticleCollectionType;

  const MTDGeometry* geom_; 
 
  EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  
  EDGetTokenT< BTLDigiCollection >  btlDigisToken_;
  EDGetTokenT< ETLDigiCollection >  etlDigisToken_;

  EDGetTokenT< FTLRecHitCollection > btlRecHitToken_;
  EDGetTokenT< FTLRecHitCollection > etlRecHitToken_;
  
  EDGetTokenT< FTLClusterCollection > btlClusterToken_;
  EDGetTokenT< FTLClusterCollection > etlClusterToken_;
  EDGetTokenT< vector<l1extra::L1JetParticle> > l1TausToken_;
  InputTag genSrc_;
  double time_cut_;

  double genPt, genEta, genPhi, genEnergy, genPdgId, genStatus;
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
  //double time;        
  //double timeError;   
  double size;        
  double sizeX;       
  double sizeY;       
  //double energy;      
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
  l1TausToken_(        consumes< vector<l1extra::L1JetParticle> >( cfg.getParameter<InputTag>("l1Taus"))),
  //vector<l1extra::L1JetParticle>        "l1extraParticles"          "Tau"             "RECO"
  genSrc_ ((           cfg.getParameter<edm::InputTag>( "genParticles"))),
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
  genParticleTree->Branch("gen_pt",            &genPt,            "gen_pt/D");           
  genParticleTree->Branch("gen_eta",           &genEta,           "gen_eta/D");           
  genParticleTree->Branch("gen_phi",           &genPhi,           "gen_phi/D");           
  genParticleTree->Branch("gen_Energy",        &genEnergy,        "gen_Energy/D");           
  genParticleTree->Branch("gen_pdgId",         &genPdgId,         "gen_pdgId/D");           
  genParticleTree->Branch("gen_status",        &genStatus,         "gen_status/D");           
  genParticleTree->Branch("gen_isPromptFinalState",      &genIsPromptFinalState,      "gen_isPromptFinalState/b");           


}
//destructor
L1MTDLLPAnalyzer::~L1MTDLLPAnalyzer()
{
}


void 
L1MTDLLPAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  edm::ESHandle<MTDGeometry> geom;
  if( geom_ == nullptr ) {
    iSetup.get<MTDDigiGeometryRecord>().get(geom);
    geom_ = geom.product();
  }

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
  
  edm::Handle< vector<l1extra::L1JetParticle> > l1Taus;
  iEvent.getByToken(l1TausToken_, l1Taus);

  edm::Handle< std::vector< reco::GenParticle> > genParticles;
  iEvent.getByToken(genToken_, genParticles);


  // *** Generator Particles
  if( genParticles->size() > 0 ){
    for(const auto& genParticle : *genParticles){
      genPt        = (double)genParticle.p4().pt();
      genEta       = (double)genParticle.p4().eta();
      genPhi       = (double)genParticle.p4().phi();
      genEnergy    = (double)genParticle.p4().e();
      genPdgId     = (double)genParticle.pdgId();
      genStatus    = (double)genParticle.status();
      genIsPromptFinalState  = (bool)genParticle.isPromptFinalState();

      genParticleTree->Fill();

    }
  }


  // *** BTL RecHits
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

  std::cout<<"l1Taus size "<<l1Taus->size()<<std::endl;
  for(const auto& l1object : *l1Taus){
    std::cout<<"x: "<<l1object.p4().x()<<std::endl;
    std::cout<<"y: "<<l1object.p4().y()<<std::endl;
  }

  unsigned int n_digi_btl[2] = {0,0};
  
  if (h_BTL_digi->size() > 0 ) {
    
    std::cout << " ----------------------------------------" << std::endl;
    std::cout << " BTL DIGI collection:" << std::endl;
    
    for (const auto& dataFrame: *h_BTL_digi) {
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

  std::cout<<"number of btl_hits "<<btl_hits.size()<<std::endl;
    for (auto const& hit: btl_hits) {
      
      BTLDetId detId(hit.first); 
      
      DetId geoId = BTLDetId(detId.mtdSide(),detId.mtdRR(),detId.module()+14*(detId.modType()-1),0,1);
      const MTDGeomDet* thedet = geom_->idToDet(geoId);
      const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
      const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());
      
      if ( (hit.second).reco_energy < btlMinEnergy_ ) continue;

      Local3DPoint simscaled(0.1*(hit.second).sim_x,0.1*(hit.second).sim_y,0.1*(hit.second).sim_z);
      simscaled = topo.pixelToModuleLocalPoint(simscaled,detId.row(topo.nrows()),detId.column(topo.nrows()));
      const auto& global_pos = thedet->toGlobal(simscaled);
      
      eta          = (double)global_pos.eta();
      phi          = (double)global_pos.phi();
      x            = (double)global_pos.x();
      y            = (double)global_pos.y();
      z            = (double)global_pos.z();
      time         = (double)(hit.second).reco_time;
      //timeError    = (double);
      energy       = (double)(hit.second).reco_energy;
      
      digiTree->Fill();
    }
  

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
