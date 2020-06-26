// -*- C++ -*-
//
// Package:    AnalyzeZll/LeptonMcFilterElectron
// Class:      LeptonMcFilterElectron
// 
/**\class LeptonMcFilterElectron LeptonMcFilterElectron.cc AnalyzeZll/LeptonMcFilterElectron/plugins/LeptonMcFilterElectron.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rogelio Reyes Almanza
//         Created:  Mon, 20 Aug 2018 09:33:34 GMT
//
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

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/PatCandidates/interface/UserData.h" //new

// DataFormat includes
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"//new
#include "DataFormats/PatCandidates/interface/Electron.h"//new
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

class LeptonMcFilterElectron : public edm::stream::EDProducer<>{
   public:
      explicit LeptonMcFilterElectron(const edm::ParameterSet&);
      ~LeptonMcFilterElectron();
//      UInt_t isTriggerMatched(const pat::MuonCollection *);

//      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      //void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      //void removeDuplicates(pat::MuonCollection&);
      Float_t getIsoVar(const pat::Electron& );       
	      
      // ----------member data ---------------------------
      edm::EDGetToken electronsMiniAODToken_;
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
LeptonMcFilterElectron::LeptonMcFilterElectron(const edm::ParameterSet& iConfig)
{

   electronsMiniAODToken_    = mayConsume<edm::View<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electronsMiniAOD"));
   pVToken_ 	= consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"));



   produces<pat::CompositeCandidateCollection>("dilepton"); 

 
}


LeptonMcFilterElectron::~LeptonMcFilterElectron()
{
 

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
LeptonMcFilterElectron::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//   typedef Candidate::LorentzVector LorentzVector;
   using namespace edm; 

   //std::auto_ptr<pat::MuonCollection> selectedCollection(new pat::MuonCollection);
   std::auto_ptr<pat::ElectronCollection > selectedCollection(new pat::ElectronCollection );
   std::unique_ptr<pat::CompositeCandidateCollection> dileptonColl(new pat::CompositeCandidateCollection); 
	
   //edm::Handle<pat::MuonCollection> muons;
   //iEvent.getByToken(muonToken_,muons);
   edm::Handle<edm::View<pat::Electron> > electrons;
   if( !electrons.isValid() ) iEvent.getByToken(electronsMiniAODToken_,electrons);

   edm::ESHandle<TransientTrackBuilder> theTTBuilder;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(pVToken_, vertices);

   
 
 
   if (vertices->empty()) return;
   reco::VertexCollection::const_iterator PV = vertices->end();
   for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
    if ( !(vtx->isFake())
         && vtx->ndof()>=4. && vtx->position().Rho() < 2.0
         && fabs(vtx->position().Z()) < 24.0) {
         PV = vtx;
         break;
    }
  }
  if ( PV==vertices->end() ) return;
 


  for ( edm::View<pat::Electron>::const_iterator el1 = electrons->begin() ; el1 != electrons->end(); ++el1) { //Remove Ele1++
        //std::cout << "fist loop Electron ok " << std::endl; 
	selectedCollection->push_back(*el1); 
    } 

	 
  for (pat::ElectronCollection::const_iterator electron1 = selectedCollection->begin() ;electron1 !=  selectedCollection->end(); ++electron1) {
	for (pat::ElectronCollection::const_iterator electron2 = electron1+1 ;electron2 !=  selectedCollection->end(); ++electron2){  
	 	if((electron1->charge() * electron2->charge()) != -1) continue; 
                //if(electron1 == electron2) continue; //try to remove overcounting
                //double dR = deltaR(*muon2, *(muon1->innerTrack()));
	  	//if (dR < drCut_) continue; 
		//reco::Candidate::LorentzVector  dilept = muon1->p4()+muon2->p4(); 
		//if ( dilept.M() > 80.5)continue; 
		
        if((electron1->charge() * electron2->charge()) != -1) continue;
        pat::CompositeCandidate myCand;
        myCand.setP4( electron1->p4() + electron2->p4());
        if(electron1->charge() == -1){
           myCand.addDaughter( *electron1      ,"lepton1");
           myCand.addDaughter( *electron2      ,"lepton2");
        }
        else{
            myCand.addDaughter( *electron2      ,"lepton1");
            myCand.addDaughter( *electron1      ,"lepton2");
        }
		dileptonColl->push_back(myCand);
        //std::cout<< "pass Electron Filter" << std::endl;
	}


   }// end for electron1
   iEvent.put(std::move(dileptonColl),"dilepton");
}
/*
Float_t LeptonMcFilterElectron::getIsoVar(const pat::Muon& mu){
		Float_t coriso = 99.0; 
		reco::MuonPFIsolation pfR03 = mu.pfIsolationR03();
	        coriso = pfR03.sumChargedHadronPt + std::max(0., pfR03.sumNeutralHadronEt+pfR03.sumPhotonEt-0.5*pfR03.sumPUPt);
return coriso; 
}
*/
Float_t LeptonMcFilterElectron::getIsoVar(const pat::Electron& el){
                reco::GsfElectron::PflowIsolationVariables pfIso1 = el.pfIsolationVariables();
                float coriso1 = pfIso1.sumChargedHadronPt + std::max(0., pfIso1.sumNeutralHadronEt+pfIso1.sumPhotonEt-0.5*pfIso1.sumPUPt);
return coriso1;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
//void
//LeptonMcFilterElectron::beginStream(edm::StreamID)
//{
//}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
//void
//LeptonMcFilterElectron::endStream() {
//}

// ------------ method called when starting to processes a run  ------------
/*
void
LeptonMcFilterElectron::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
LeptonMcFilterElectron::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
LeptonMcFilterElectron::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
LeptonMcFilterElectron::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
LeptonMcFilterElectron::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}*/

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonMcFilterElectron);
