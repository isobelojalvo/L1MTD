// -*- C++ -*-
//
// Package:    L1Trigger/L1MTDAnalyzer
// Class:      L1MTDAnalyzer
//
/**\class L1MTDAnalyzer L1MTDAnalyzer.cc L1Trigger/L1MTDAnalyzer/plugins/L1MTDAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Wed, 10 Oct 2018 19:40:40 GMT
//
//


// system include files
#include <memory>
#include <utility>

// Math Include
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"

// General tool includes for TrackStub
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"

// Geometry includes for TrackStub
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class L1MTDAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit L1MTDAnalyzer(const edm::ParameterSet&);
      ~L1MTDAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_; 
      edm::EDGetTokenT<edm::ValueMap<float> > timingValuesToken_;
      edm::EDGetTokenT<edm::ValueMap<float> > timingValuesSmearedToken_;
      bool saveStubs; // option to save also stubs in the ntuples (makes them large...) // BBT 10-18-18
      edm::EDGetTokenT< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > ttStubToken_;   // BBT 10-18-18
      edm::EDGetTokenT< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > ttStubMCTruthToken_;   // BBT 10-18-18
      edm::InputTag L1StubInputTag; // BBT 10-18-18
      edm::InputTag MCTruthStubInputTag; // BBT 10-18-18
      edm::InputTag primarySimVertexPosTag; // BBT 10-26-18
      edm::InputTag primarySimVertexTimeTag; // BBT 10-26-18
      edm::EDGetTokenT< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > primarySimVertexPosToken_;
      edm::EDGetTokenT< float > primarySimVertexTimeToken_;
  
  // --------------- output tree and members ---------------
  TTree* trackerTree;
  
      std::vector<double> v_trackPt; 
      std::vector<double> v_trackEta; 
      std::vector<double> v_trackTime; 
      std::vector<double> v_trackTime_smeared; 
      std::vector<double> v_trackTime_PVcorrected; 

      // Vertex 10-26-18 BBT
      double d_vertexX; 
      double d_vertexY; 
      double d_vertexZ; 
      double d_vertexTime; 

      // SimHit 10-26-18 BBT
      std::vector<double> v_simHitX; 
      std::vector<double> v_simHitY; 
      std::vector<double> v_simHitZ; 
      std::vector<double> v_simHitTime; 

      // ALL stubs
      double nStubsPerEvent;
      double nStubsPerEvent_2p5;
      double nStubsPerEvent_5;
      std::vector<float>* m_allstub_x;
      std::vector<float>* m_allstub_y;
      std::vector<float>* m_allstub_z;

      std::vector<int>* m_allstub_isBarrel; // stub is in barrel (1) or in disk (0)
      std::vector<int>* m_allstub_layer;
      std::vector<int>* m_allstub_isPSmodule;

      std::vector<float>* m_allstub_trigDisplace;
      std::vector<float>* m_allstub_trigOffset;
      std::vector<float>* m_allstub_trigPos;
      std::vector<float>* m_allstub_trigBend;
     
      // stub associated with tracking particle ?
      std::vector<int>*   m_allstub_matchTP_pdgid; // -999 if not matched
      std::vector<float>* m_allstub_matchTP_pt;    // -999 if not matched
      std::vector<float>* m_allstub_matchTP_eta;   // -999 if not matched
      std::vector<float>* m_allstub_matchTP_phi; // -999 if not matched

      std::vector<int>* m_allstub_genuine;
    
      // stubs matched to simHits, BBT 10-24-18
      std::vector<float>* m_allstub_matchSimHit_pt; 
      std::vector<float>* m_allstub_matchSimHit_eta;
      std::vector<float>* m_allstub_matchSimHit_phi;


      // --------------- histograms ---------------
      TH1F* track_pt;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1MTDAnalyzer::L1MTDAnalyzer(const edm::ParameterSet& cfg):
  ttTrackToken_(  consumes< std::vector< TTTrack < Ref_Phase2TrackerDigi_ > > >(cfg.getParameter<edm::InputTag>("L1TrackInputTag"))),
  timingValuesToken_( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("timingValuesNominal"))),
  timingValuesSmearedToken_( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("timingValuesSmeared")))
{
  // track stubs, BBT 10-18-18
  saveStubs = cfg.getParameter< bool >("saveStubs");
  L1StubInputTag = cfg.getParameter<edm::InputTag>("L1StubInputTag");
  MCTruthStubInputTag = cfg.getParameter<edm::InputTag>("MCTruthStubInputTag");
  ttStubToken_ = consumes< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > >(L1StubInputTag);
  ttStubMCTruthToken_ = consumes< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > >(MCTruthStubInputTag);
  primarySimVertexPosTag = cfg.getParameter<edm::InputTag>("genParticles_xyz");
  primarySimVertexTimeTag = cfg.getParameter<edm::InputTag>("genParticles_t");
  primarySimVertexPosToken_ = consumes< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> >(primarySimVertexPosTag);
  primarySimVertexTimeToken_ = consumes< float >(primarySimVertexTimeTag);

  // output tree and members
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  trackerTree = fs->make<TTree>("trackerTree", "track information");
  trackerTree->Branch("track_pt",    &v_trackPt);
  trackerTree->Branch("track_eta",    &v_trackEta);
  trackerTree->Branch("track_time",    &v_trackTime);
  trackerTree->Branch("track_time_smeared",    &v_trackTime_smeared);
  trackerTree->Branch("track_time_PVcorrected",    &v_trackTime_PVcorrected);

  trackerTree->Branch("vertex_X",    &d_vertexX);
  trackerTree->Branch("vertex_Y",    &d_vertexY);
  trackerTree->Branch("vertex_Z",    &d_vertexZ);
  trackerTree->Branch("vertex_time", &d_vertexTime);

  trackerTree->Branch("nStubsPerEvent",    &nStubsPerEvent);
  trackerTree->Branch("nStubsPerEvent_2p5",    &nStubsPerEvent_2p5);
  trackerTree->Branch("nStubsPerEvent_5",    &nStubsPerEvent_5);

  // histograms
  track_pt             = fs->make<TH1F>( "track_pt"  , "p_{t}", 200,  0., 200. );

}
		  

L1MTDAnalyzer::~L1MTDAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1MTDAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //std::cout << " !!! Running L1MTDAnalyzer::analyze() !!!" << std::endl;

   using namespace edm;
   std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks;
   //std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks_time;
   // L1 tracks
   edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > l1trackHandle;
   iEvent.getByToken(ttTrackToken_, l1trackHandle);

   // L1 stubs, BBT 10-18-18
   edm::Handle< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_ > > > TTStubHandle;
   if (saveStubs) iEvent.getByToken(ttStubToken_, TTStubHandle);

   // more for TTStubs, BBT 10-18-18
   edm::ESHandle<TrackerGeometry> geometryHandle;
   iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

   edm::ESHandle<TrackerTopology> tTopoHandle;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

   edm::ESHandle<TrackerGeometry> tGeomHandle;
   iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);

   const TrackerTopology* const tTopo = tTopoHandle.product();
   const TrackerGeometry* const theTrackerGeom = tGeomHandle.product();

   // MC truth association maps, BBT 10-18-18
   edm::Handle< TTStubAssociationMap< Ref_Phase2TrackerDigi_ > > MCTruthTTStubHandle;
   iEvent.getByToken(ttStubMCTruthToken_, MCTruthTTStubHandle);

   // Timing
   edm::Handle<edm::ValueMap<float> > timingValues;
   edm::Handle<edm::ValueMap<float> > timingValues_smeared;
   iEvent.getByToken(timingValuesToken_,timingValues);
   iEvent.getByToken(timingValuesSmearedToken_,timingValues_smeared);

   // pair stuff
   std::vector< std::pair < TTTrack< Ref_Phase2TrackerDigi_ >, double > > l1Tracks_pair;

   // Get vertex info, BBT 10-26-18
   edm::Handle< ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> > primarySimVertexPos;
   iEvent.getByToken(primarySimVertexPosToken_, primarySimVertexPos);
   edm::Handle< float > primarySimVertexTime;
   iEvent.getByToken(primarySimVertexTimeToken_, primarySimVertexTime);

   // -----------------------------------------------------------------------------------------------

   // variables for output
   l1Tracks.clear();
   //l1Tracks_time.clear();
   v_trackPt.clear();
   v_trackEta.clear();
   v_trackTime.clear();
   v_trackTime_smeared.clear();
   v_trackTime_PVcorrected.clear();
   d_vertexTime = -999;
   d_vertexX = -999;
   d_vertexY = -999;
   d_vertexZ = -999;
   nStubsPerEvent=0;
   nStubsPerEvent_2p5=0;
   nStubsPerEvent_5=0;

   if (saveStubs) {
     m_allstub_x->clear();
     m_allstub_y->clear();
     m_allstub_z->clear();

     m_allstub_isBarrel->clear();
     m_allstub_layer->clear();
     m_allstub_isPSmodule->clear();

     m_allstub_trigDisplace->clear();
     m_allstub_trigOffset->clear();
     m_allstub_trigPos->clear();
     m_allstub_trigBend->clear();

     m_allstub_matchTP_pdgid->clear();
     m_allstub_matchTP_pt->clear();
     m_allstub_matchTP_eta->clear();
     m_allstub_matchTP_phi->clear();
     
     // match SimHit, BBT 10-24-18
     m_allstub_matchSimHit_pt->clear();
     m_allstub_matchSimHit_eta->clear();
     m_allstub_matchSimHit_phi->clear();

     m_allstub_genuine->clear();
   }

   // VERTEX stuff, BBT 10-26-18
   d_vertexTime = (*primarySimVertexTime) ;
   d_vertexX    = (*primarySimVertexPos).X() ;
   d_vertexY    = (*primarySimVertexPos).Y() ;
   d_vertexZ    = (*primarySimVertexPos).Z() ;
   
   //Find and sort the tracks
   for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
       
       edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, track_index);
       double value = (*timingValues)[ptr];
       double value_smeared = (*timingValues_smeared)[ptr];
       
       //if (track_index = 0)
       //std::cout << "unsorted " << track_index << "th l1track, value: "<< value << " , pt: " << ptr->getMomentum().perp() <<std::endl;
       double pt  = ptr->getMomentum().perp();
       double eta = ptr->getMomentum().eta();
       v_trackPt.push_back(pt);
       v_trackEta.push_back(eta);
       v_trackTime.push_back( value );     
       v_trackTime_smeared.push_back( value_smeared );     
       v_trackTime_PVcorrected.push_back( value_smeared - (*primarySimVertexTime) );     
       //only using tracks with eta less than 1.5 and pt greater than 2.5 GeV
       //if(abs(eta)<1.5 && pt > 2.5)
       l1Tracks.push_back(l1trackHandle->at(track_index));       
       //l1Tracks_time.push_back(value);       

       l1Tracks_pair.push_back( std::make_pair(l1trackHandle->at(track_index), value) );
     }

   //std::cout << "0th l1Track before sorting: " << l1Tracks.at(0).getMomentum().perp() << std::endl;
   //std::cout << "1st l1Track before sorting: " << l1Tracks.at(1).getMomentum().perp() << std::endl;


   //std::sort(l1Tracks.begin(), l1Tracks.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   


   //std::sort(l1Tracks_pair.begin(), l1Tracks_pair.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   
   //std::cout << "0th l1Track after sorting: "  << l1Tracks.at(0).getMomentum().perp() << std::endl;
   //std::cout << "1st l1Track before sorting: " << l1Tracks.at(1).getMomentum().perp() << std::endl;
   
   /*
   // store pt ordered l1Tracks
   for(unsigned int i = 0; i < l1Tracks.size(); i++) {
     v_trackPt.push_back( l1Tracks.at(i).getMomentum().perp() );
     v_trackEta.push_back( l1Tracks.at(i).getMomentum().eta() );

     edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, i);
     double value = (*timingValues)[ptr];
     v_trackTime.push_back( value );     
     std::cout<<"sorted " << i << "th l1track, value: "<< value << " , pt: " << l1Tracks.at(i).getMomentum().perp() <<std::endl;

     //std::cout<<"2sorted " << i << "th l1track, value: "<< l1Tracks_pair.at(i).second << " , pt: " << l1Tracks.at(i).getMomentum().perp() <<std::endl;
   }
   */


   // ----------------------------------------------------------------------------------------------
   // loop over L1 stubs
   // ----------------------------------------------------------------------------------------------

   if (saveStubs) {
     
     for (auto gd=theTrackerGeom->dets().begin(); gd != theTrackerGeom->dets().end(); gd++) {
       
       DetId detid = (*gd)->geographicalId();
       if(detid.subdetId()!=StripSubdetector::TOB && detid.subdetId()!=StripSubdetector::TID ) continue;
       if(!tTopo->isLower(detid) ) continue; // loop on the stacks: choose the lower arbitrarily
       DetId stackDetid = tTopo->stack(detid); // Stub module detid
       
       if (TTStubHandle->find( stackDetid ) == TTStubHandle->end() ) continue;
       
       // Get the DetSets of the Clusters
       edmNew::DetSet< TTStub< Ref_Phase2TrackerDigi_ > > stubs = (*TTStubHandle)[ stackDetid ];
       const GeomDetUnit* det0 = theTrackerGeom->idToDetUnit( detid );
       const PixelGeomDetUnit* theGeomDet = dynamic_cast< const PixelGeomDetUnit* >( det0 );
       const PixelTopology* topol = dynamic_cast< const PixelTopology* >( &(theGeomDet->specificTopology()) );
       
       // loop over stubs
       for ( auto stubIter = stubs.begin();stubIter != stubs.end();++stubIter ) {
	 nStubsPerEvent++;

	 edm::Ref< edmNew::DetSetVector< TTStub< Ref_Phase2TrackerDigi_  > >, TTStub< Ref_Phase2TrackerDigi_  > >
	   tempStubPtr = edmNew::makeRefTo( TTStubHandle, stubIter );
	 
	 int isBarrel = 0;
	 int layer=-999999;
	 if ( detid.subdetId()==StripSubdetector::TOB ) {
	   isBarrel = 1;
	   layer  = static_cast<int>(tTopo->layer(detid));
	 }
	 else if ( detid.subdetId()==StripSubdetector::TID ) {
	   isBarrel = 0;
	   layer  = static_cast<int>(tTopo->layer(detid));
	 }
	 else {
	   std::cout << "WARNING -- neither TOB or TID stub, shouldn't happen..." << std::endl;
	   layer = -1;
	 }
	 
	 int isPSmodule=0;
	 if (topol->nrows() == 960) isPSmodule=1;
	 
	 MeasurementPoint coords = tempStubPtr->getClusterRef(0)->findAverageLocalCoordinatesCentered();      
	 LocalPoint clustlp = topol->localPosition(coords);
	 GlobalPoint posStub  =  theGeomDet->surface().toGlobal(clustlp);
	 
	 double tmp_stub_x=posStub.x();
	 double tmp_stub_y=posStub.y();
	 double tmp_stub_z=posStub.z();
	 
	 float trigDisplace = tempStubPtr->getTriggerDisplacement();
	 float trigOffset = tempStubPtr->getTriggerOffset();
	 float trigPos = tempStubPtr->getTriggerPosition();
	 float trigBend = tempStubPtr->getTriggerBend();
	 
	 m_allstub_x->push_back(tmp_stub_x);
	 m_allstub_y->push_back(tmp_stub_y);
	 m_allstub_z->push_back(tmp_stub_z);
	 
	 m_allstub_isBarrel->push_back(isBarrel);
	 m_allstub_layer->push_back(layer);
	 m_allstub_isPSmodule->push_back(isPSmodule);
	 
	 m_allstub_trigDisplace->push_back(trigDisplace);
	 m_allstub_trigOffset->push_back(trigOffset);
	 m_allstub_trigPos->push_back(trigPos);
	 m_allstub_trigBend->push_back(trigBend);
	 
	 // matched to tracking particle? 
	 edm::Ptr< TrackingParticle > my_tp = MCTruthTTStubHandle->findTrackingParticlePtr(tempStubPtr);

	 int myTP_pdgid = -999;
	 float myTP_pt  = -999;
	 float myTP_eta = -999;
	 float myTP_phi = -999;
	 
	 if (my_tp.isNull() == false) {
	   int tmp_eventid = my_tp->eventId().event();
	   
	   if (tmp_eventid > 0) continue; // this means stub from pileup track
	   
	   myTP_pdgid = my_tp->pdgId();
	   myTP_pt = my_tp->p4().pt();
	   myTP_eta = my_tp->p4().eta();
	   myTP_phi = my_tp->p4().phi();
	 }

	 if (myTP_pt > 2.5) nStubsPerEvent_2p5++;
	 if (myTP_pt > 5) nStubsPerEvent_5++;

	 m_allstub_matchTP_pdgid->push_back(myTP_pdgid);
	 m_allstub_matchTP_pt->push_back(myTP_pt);
	 m_allstub_matchTP_eta->push_back(myTP_eta);
	 m_allstub_matchTP_phi->push_back(myTP_phi);
	 
	 int tmp_stub_genuine = 0;
	 if (MCTruthTTStubHandle->isGenuine(tempStubPtr)) tmp_stub_genuine = 1;
	 
	 m_allstub_genuine->push_back(tmp_stub_genuine);
	
	 // try to match stubs to SimHits, BBT 10-24-18
	 for(unsigned int track_idx = 0; track_idx < l1Tracks.size(); track_idx++) {
	   double etaTrack = l1Tracks.at(track_idx).getMomentum().eta();
	   double phiTrack = l1Tracks.at(track_idx).getMomentum().phi();
	   //double etaStub  = posStub.eta();
	   //double phiStub  = posStub.phi();
	   double etaStub  = myTP_eta;
	   double phiStub  = myTP_phi;
	   double deltaR   = sqrt( (etaTrack - etaStub)*(etaTrack - etaStub) + (phiTrack - phiStub)*(phiTrack - phiStub));
	   if (deltaR < 0.1) {
	     m_allstub_matchSimHit_pt->push_back(myTP_pt);
	     m_allstub_matchSimHit_eta->push_back(myTP_eta);
	     m_allstub_matchSimHit_phi->push_back(myTP_phi);
	   }
	 }
	 
       } // end loop over stubs
     }
   } // end saveStubs
   
   // ***. Fill tree
   trackerTree->Fill();


   
}


// ------------ method called once each job just before starting event loop  ------------
void
L1MTDAnalyzer::beginJob()
{
  std::cout << "STARTING JOB!!!" << std::endl;

  m_allstub_x = new std::vector<float>;
  m_allstub_y = new std::vector<float>;
  m_allstub_z = new std::vector<float>;

  m_allstub_isBarrel = new std::vector<int>;
  m_allstub_layer    = new std::vector<int>;
  m_allstub_isPSmodule = new std::vector<int>;

  m_allstub_trigDisplace = new std::vector<float>;
  m_allstub_trigOffset   = new std::vector<float>;
  m_allstub_trigPos      = new std::vector<float>;
  m_allstub_trigBend     = new std::vector<float>;

  m_allstub_matchTP_pdgid = new std::vector<int>;
  m_allstub_matchTP_pt    = new std::vector<float>;
  m_allstub_matchTP_eta   = new std::vector<float>;
  m_allstub_matchTP_phi = new std::vector<float>;

  m_allstub_matchSimHit_pt    = new std::vector<float>;
  m_allstub_matchSimHit_eta   = new std::vector<float>;
  m_allstub_matchSimHit_phi = new std::vector<float>;

  m_allstub_genuine = new std::vector<int>;

  if (saveStubs) {
    trackerTree->Branch("allstub_x", &m_allstub_x);
    trackerTree->Branch("allstub_y", &m_allstub_y);
    trackerTree->Branch("allstub_z", &m_allstub_z);

    trackerTree->Branch("allstub_isBarrel",   &m_allstub_isBarrel);
    trackerTree->Branch("allstub_layer",      &m_allstub_layer);
    trackerTree->Branch("allstub_isPSmodule", &m_allstub_isPSmodule);

    trackerTree->Branch("allstub_trigDisplace", &m_allstub_trigDisplace);
    trackerTree->Branch("allstub_trigOffset",   &m_allstub_trigOffset);
    trackerTree->Branch("allstub_trigPos",      &m_allstub_trigPos);
    trackerTree->Branch("allstub_trigBend", &m_allstub_trigBend);

    trackerTree->Branch("allstub_matchTP_pdgid", &m_allstub_matchTP_pdgid);
    trackerTree->Branch("allstub_matchTP_pt",    &m_allstub_matchTP_pt);
    trackerTree->Branch("allstub_matchTP_eta",   &m_allstub_matchTP_eta);
    trackerTree->Branch("allstub_matchTP_phi", &m_allstub_matchTP_phi);

    trackerTree->Branch("allstub_matchSimHit_pt",    &m_allstub_matchSimHit_pt);
    trackerTree->Branch("allstub_matchSimHit_eta",   &m_allstub_matchSimHit_eta);
    trackerTree->Branch("allstub_matchSimHit_phi", &m_allstub_matchSimHit_phi);

    trackerTree->Branch("allstub_genuine", &m_allstub_genuine);
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void
L1MTDAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1MTDAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(L1MTDAnalyzer);
