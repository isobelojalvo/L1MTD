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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
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
      TTree* trackerTree;

      std::vector<double> v_trackPt; 
      std::vector<double> v_trackEta; 
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
  timingValuesToken_( consumes<edm::ValueMap<float> >(cfg.getParameter<edm::InputTag>("timingValuesNominal")))
{
  usesResource("TFileService");
  edm::Service<TFileService> fs;
  trackerTree = fs->make<TTree>("trackerTree", "track information");
  trackerTree->Branch("track_pt",    &v_trackPt);
  trackerTree->Branch("track_eta",    &v_trackEta);

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
   using namespace edm;
   std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks;
   edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > l1trackHandle;
   iEvent.getByToken(ttTrackToken_, l1trackHandle);
   l1Tracks.clear();
   v_trackPt.clear();
   v_trackEta.clear();

   edm::Handle<edm::ValueMap<float> > timingValues;
   iEvent.getByToken(timingValuesToken_,timingValues);

   //Find and sort the tracks
   for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
     {
       
       edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, track_index);
       double value = (*timingValues)[ptr];
       std::cout<<"value: "<< value<<std::endl;
       double pt  = ptr->getMomentum().perp();
       double eta = ptr->getMomentum().eta();
       v_trackPt.push_back(pt);
       v_trackEta.push_back(eta);

       //only using tracks with eta less than 1.5 and pt greater than 2.5 GeV
       //if(abs(eta)<1.5 && pt > 2.5)
       l1Tracks.push_back(l1trackHandle->at(track_index));       
       
     }
   
   std::sort(l1Tracks.begin(), l1Tracks.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   
   
   // ***. Fill tree
   trackerTree->Fill();
   
}


// ------------ method called once each job just before starting event loop  ------------
void
L1MTDAnalyzer::beginJob()
{
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
