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

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

//#include "DataFormats/FTLDigi/interface/FTLDigiCollections.h"

//#include "DataFormats/FTLRecHit/interface/FTLRecHit.h"
//#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"

//#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"

using namespace edm;
using namespace std;

class L1MTDGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
  public:

  explicit L1MTDGenAnalyzer(const edm::ParameterSet&);

  ~L1MTDGenAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  typedef std::vector<reco::GenParticle> GenParticleCollectionType;
  
  EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  EDGetTokenT<edm::HepMCProduct> HepMCProductToken_;

  //EDGetTokenT< BTLDigiCollection >  btlDigisToken_;
  //EDGetTokenT< ETLDigiCollection >  etlDigisToken_;

  //EDGetTokenT< FTLRecHitCollection > btlRecHitToken_;
  //EDGetTokenT< FTLRecHitCollection > etlRecHitToken_;
  
  //EDGetTokenT< FTLClusterCollection > btlClusterToken_;
  //EDGetTokenT< FTLClusterCollection > etlClusterToken_;

  double time_cut_;

  InputTag genSrc_;

  double genPt, genEta, genPhi;

  TTree* recHitTree;
  TTree* clusterTree;
  double time, isTimeValid, timeError, energy;
  double x;           
  double y;           
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
  
};

L1MTDGenAnalyzer::L1MTDGenAnalyzer(const edm::ParameterSet &cfg) :  
  //btlDigisToken_(      consumes< BTLDigiCollection >   ( cfg.getParameter<InputTag>("FTLBarrel"))),
  //etlDigisToken_(      consumes< ETLDigiCollection >   ( cfg.getParameter<InputTag>("FTLEndcap"))),
  //btlRecHitToken_(     consumes< FTLRecHitCollection > ( cfg.getParameter<InputTag>("recHitBarrel"))), // finish me
  //etlRecHitToken_(     consumes< FTLRecHitCollection > ( cfg.getParameter<InputTag>("recHitEndcap"))), // finish me
  //btlClusterToken_(    consumes< FTLClusterCollection > ( cfg.getParameter<InputTag>("mtdClusterBarrel"))),
  //etlClusterToken_(    consumes< FTLClusterCollection > ( cfg.getParameter<InputTag>("mtdClusterEndcap"))),
  HepMCProductToken_(    consumes< edm::HepMCProduct >    ( cfg.getParameter<InputTag>("HepMCProduct"))),
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
}
//destructor
L1MTDGenAnalyzer::~L1MTDGenAnalyzer()
{
}


void 
L1MTDGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //edm::Handle<BTLDigiCollection> btlDigis;
  //iEvent.getByToken(btlDigisToken_,btlDigis);

  //edm::Handle<ETLDigiCollection> etlDigis;
  //iEvent.getByToken(etlDigisToken_,etlDigis);

  //edm::Handle< FTLRecHitCollection > btlRecHits;
  //iEvent.getByToken(btlRecHitToken_,btlRecHits);

  //edm::Handle< FTLRecHitCollection > etlRecHits;
  //iEvent.getByToken(etlRecHitToken_,etlRecHits);

  //edm::Handle< FTLClusterCollection > btlClusters;
  //iEvent.getByToken(btlClusterToken_,btlClusters);

  //edm::Handle< FTLClusterCollection > etlClusters;
  //iEvent.getByToken(etlClusterToken_,etlClusters);

  //edm::Handle<std::vector<reco::GenParticle> > genParticleHandle;
  //for(unsigned int i = 0; i< genParticleHandle->size(); i++){
  //edm::Ptr<reco::GenParticle> ptr(genParticleHandle, i);
  //}

   edm::Handle<edm::HepMCProduct> HepMCHandle;
   iEvent.getByToken( HepMCProductToken_, HepMCHandle) ;
   const HepMC::GenEvent* Evt = HepMCHandle->GetEvent() ;
   //
   // this is an example loop over the hierarchy of vertices
   //
   std::cout<<"Here 1"<<std::endl;
   int i = 0;
   for ( HepMC::GenEvent::vertex_const_iterator
	   itVtx=Evt->vertices_begin(); itVtx!=Evt->vertices_end(); ++itVtx )
     {
       std::cout<<"Here 2"<<std::endl;
       i++;
       //
       // this is an example loop over particles coming out of each vertex in the loop
       //
       int j = 0; 
       for ( HepMC::GenVertex::particles_out_const_iterator
	       itPartOut=(*itVtx)->particles_out_const_begin();
	     itPartOut!=(*itVtx)->particles_out_const_end(); ++itPartOut )
	 {
	   std::cout<<"Here 3 "<<std::endl;
	   double time = (*itPartOut)->momentum().t();
	   if(j<20)
	     std::cout<< "pdgID: "<<(*itPartOut)->pdg_id() <<" mass: "<< (*itPartOut)->momentum().m()<< " px: "<<(*itPartOut)->momentum.px() <<" time: "<<time<<std::endl;
	   j++;
	   std::cout<<"Here 4"<<std::endl;
	 }
     }
   std::cout<<"n HepMC Vertices: "<<i<<std::endl;

  //if( btlRecHits->size() > 0 ){
  //for(const auto& recHit : *btlRecHits){

  //  time        = (double)recHit.time();
  //  isTimeValid = (double)recHit.isTimeValid();
  //  timeError   = (double)recHit.timeError();;
  //  energy      = (double)recHit.energy();
  //  recHitTree->Fill();
  //}
  // }

  //if( btlClusters->size() > 0 ){
  //for(const auto& detIds : *btlClusters){
  //  for(const auto& cluster : detIds){
  //	x            = (double)cluster.x();
  //	y            = (double)cluster.y();
  //	time         = (double)cluster.time();
  //	timeError    = (double)cluster.timeError();
  //	size         = (double)cluster.size();
  //	sizeX        = (double)cluster.sizeX();
  //	sizeY        = (double)cluster.sizeY();
  //	energy       = (double)cluster.energy();
  //	hitOffset    = (double)cluster.hitOffset().at(0);
  //	hitEnergy    = (double)cluster.hitENERGY().at(0);
  //	hitTime      = (double)cluster.hitTIME().at(0);
  //	hitTimeError = (double)cluster.hitTIME_ERROR().at(0);
  //	clusterTree->Fill();
  //  }
  //}
  //}

  /*
  if (btlDigis->size() > 0 ) {

    std::cout << " ----------------------------------------" << std::endl;
    std::cout << " BTL DIGI collection:" << std::endl;
  
    for (const auto& dataFrame: *btlDigis) {

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

      } // isaple loop

    } // digi loop

  } // if (btlDigis->size() > 0 )
  */

}

// ------------ method called once each job just before starting event loop  ------------
void
L1MTDGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1MTDGenAnalyzer::endJob()
{
}

void
L1MTDGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(L1MTDGenAnalyzer);
