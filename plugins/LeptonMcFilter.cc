// -*- C++ -*-
//
// Package:    AnalyzeZphill/ZphiTupler
// Class:      LeptonMcFilter
// 
/**\class LeptonMcFilter LeptonMcFilter.cc AnalyzeZphill/ZphiTupler/plugins/LeptonMcFilter.cc

 Description: apply BF filter 1, 3, 4, & 5

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Artemio Ayala Sanchez
//         Created:  Tue, 18 Jun 2019 15:02:36 GMT
//         Updated: v7-12/09/20 not currently used 
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// specific include files
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
//Data fotmats includes
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/Common/interface/View.h>
#include <DataFormats/HepMCCandidate/interface/GenParticle.h>
#include "TMath.h"
#include "Math/VectorUtil.h"

//
// class declaration
//

class LeptonMcFilter : public edm::stream::EDProducer<> {
   public:
      explicit LeptonMcFilter(const edm::ParameterSet&);
      ~LeptonMcFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
      edm::EDGetTokenT<reco::VertexCollection> pVToken_;
      
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
LeptonMcFilter::LeptonMcFilter(const edm::ParameterSet& iConfig):
muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons")))
{
   //muonToken_   = consumes<pat::MuonCollection>(iConfig.getParameter< edm::InputTag>("muons"));
   pVToken_     = consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"));
   
   produces<pat::CompositeCandidateCollection>("dilepton");
}


LeptonMcFilter::~LeptonMcFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
LeptonMcFilter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
//   std::auto_ptr<pat::MuonCollection> selectedCollection(new pat::MuonCollection);
   std::unique_ptr<pat::CompositeCandidateCollection> dileptonColl(new pat::CompositeCandidateCollection);

//   edm::Handle<pat::MuonCollection> muons;
//   iEvent.getByToken(muonToken_,muons);
   edm::Handle< View<pat::Muon> > muons;
   iEvent.getByToken(muonToken_,muons);
   edm::ESHandle<TransientTrackBuilder> theTTBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(pVToken_, vertices);

   if (vertices->empty()) return; //if there is no vertices, return
   reco::VertexCollection::const_iterator PV = vertices->end();
   //loop over all the possible vertices in the container, starting with the one with higger PT
   for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
    if ( !(vtx->isFake())                                  // if is not fake
         && vtx->ndof()>=4. && vtx->position().Rho() < 2.0 // and number of degres of freedeom (fit) is > than 4
         && fabs(vtx->position().Z()) < 24.0) {            // and is near the beeam pipe we select it as our primary vertex PV
         PV = vtx;
         break;
    }
  }
  if ( PV==vertices->end() ) return;
  // loop over all the muons and check if it is a valid track
  //std::cout << "START LEPTON FILTER muon container lenght: "<< muons->size() << std::endl;
  for (View<pat::Muon>::const_iterator muon1 = muons->begin(); muon1 != muons->end(); ++muon1 ) {
        reco::TrackRef track_1 = muon1->track();
        if (track_1.isNull()) continue;
        reco::TransientTrack tt1 = theTTBuilder->build(*muon1->track());
        std::pair<bool,Measurement1D> tkPVdist1 = IPTools::absoluteImpactParameter3D(tt1,*PV);
        if (!tkPVdist1.first) continue;
        if (abs(tkPVdist1.second.significance())>4.) continue;
        //std::cout << "PASS MU 1, LOOPING OVER MU 2" << std::endl;
        for (View<pat::Muon>::const_iterator muon2 = muons->begin() ;muon2 !=  muons->end(); ++muon2 ) {
                if((muon1->charge() * muon2->charge()) != -1) continue;
                if(muon1==muon2) continue; //v7 add
                if(!(track_1->quality(reco::TrackBase::highPurity))) continue; //v7
                if(!(track_1->quality(reco::TrackBase::highPurity))) continue; //v7
                //std::cout << "PASS HIGH PURITY" << std::endl;
                reco::TrackRef track_2 = muon2->track();
                if (track_2.isNull()) continue;
                //std::cout << "PASS is null"<< std::endl;
                reco::TransientTrack tt2 = theTTBuilder->build(*muon2->track());
                std::pair<bool,Measurement1D> tkPVdist2 = IPTools::absoluteImpactParameter3D(tt2,*PV);
                if (!tkPVdist2.first) continue;
                if (abs(tkPVdist2.second.significance())>4.) continue;
                pat::CompositeCandidate myCand;
                //std::cout << "Pass significance saving event " << std::endl; 
                myCand.setP4( muon1->p4() + muon2->p4());
                //Leoton 1 is negative and 2 is positive
                if(muon1->charge() == -1){
                   myCand.addDaughter( *muon1      ,"lepton1");
                   myCand.addDaughter( *muon2      ,"lepton2");
                }
                else{
                   myCand.addDaughter( *muon2      ,"lepton1");
                   myCand.addDaughter( *muon1      ,"lepton2");
                }
                dileptonColl->push_back(myCand);
        }
   }// end for muon1
   iEvent.put(std::move(dileptonColl),"dilepton");
   
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
LeptonMcFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
LeptonMcFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
LeptonMcFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LeptonMcFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LeptonMcFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LeptonMcFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LeptonMcFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonMcFilter);
