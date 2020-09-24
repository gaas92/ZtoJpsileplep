// -*- C++ -*-
//
// Package:    AnalyzeZphill/ZphiTupler
// Class:      jpsiLepLepKmcFitter
// 
/**\class jpsiLepLepKmcFitter jpsiLepLepKmcFitter.cc AnalyzeZphill/ZphiTupler/plugins/jpsiLepLepKmcFitter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Artemio Ayala Sanchez
//         Created:  Tue, 18 Jun 2019 19:12:13 GMT
//         Updated:  v7-12/09/2020
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

// specific include files
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//Data format includes
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//Pat Data formats
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
//Kinematic Fit includes
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"

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
#include "TLorentzVector.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

//
// class declaration
//

class jpsiLepLepKmcFitter : public edm::stream::EDProducer<> {
   public:
      explicit jpsiLepLepKmcFitter(const edm::ParameterSet&);
      ~jpsiLepLepKmcFitter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
      bool IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu);
      double getIso(const pat::Muon& mu);
   private:
      bool    isAncestor(const reco::Candidate*, const reco::Candidate*);
      bool    isAncestor(int, const reco::Candidate*);
      virtual void beginStream(edm::StreamID) override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      //edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trackCollection_label; //miniAPD
      edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
      edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
      //edm::EDGetTokenT<pat::CompositeCandidateCollection> dilepton_Label;
      edm::EDGetTokenT<edm::View<pat::Muon>> leptonToken_;

      //edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
      edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;
    
      //what MC
      bool isMC4l_;
      //Cuts
      double ImparSigm_ = 0;
      double ImparSigl_ = 0;
      double dxym_      = 0;
      double dxyl_      = 0;
      double dzm_       = 0;
      double dzl_       = 0;
      double tlwm_      = 0;
      double plwm_      = 0;
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
jpsiLepLepKmcFitter::jpsiLepLepKmcFitter(const edm::ParameterSet& iConfig)
{
   //trackCollection_label = consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Tracks")); //miniAOD
   dimuon_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuon"));
   primaryVertices_Label = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"));
   //dilepton_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dilepton"));
   leptonToken_ = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("leptons"));

   //BSLabel_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotLabel"));
   genCands_ = consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"));
   packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
   
   isMC4l_ = iConfig.getParameter<bool>("isMC4l" );
    
   ImparSigm_       = iConfig.getParameter<double>("ImparSigm");
   ImparSigl_       = iConfig.getParameter<double>("ImparSigl");
   dxym_            = iConfig.getParameter<double>("dxym");
   dxyl_            = iConfig.getParameter<double>("dxyl");
   dzm_             = iConfig.getParameter<double>("dzm");
   dzl_             = iConfig.getParameter<double>("dzl");
   tlwm_            = iConfig.getParameter<double>("trackerLayersWithMeasurement");
   plwm_            = iConfig.getParameter<double>("pixelLayersWithMeasurement");

   produces<pat::CompositeCandidateCollection>("ZCandidates");
}


jpsiLepLepKmcFitter::~jpsiLepLepKmcFitter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//recursively check is a given particle is ancestor
bool jpsiLepLepKmcFitter::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool jpsiLepLepKmcFitter::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}
// ------------ method called to produce the data  ------------
void
jpsiLepLepKmcFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<pat::CompositeCandidateCollection> ZCandColl(new pat::CompositeCandidateCollection); //produce ZCandidates

   edm::ESHandle<TransientTrackBuilder> theTTB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB); 

   //edm::Handle<pat::CompositeCandidateCollection> dileptons;
   //iEvent.getByToken(dilepton_Label,dileptons); //dilepton
   edm::Handle< View<pat::Muon> > leptons;
   iEvent.getByToken(leptonToken_,leptons);
    
   edm::Handle<pat::CompositeCandidateCollection> dimuons;
   iEvent.getByToken(dimuon_Label,dimuons); //dimuon

   //edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;  //miniAOD
   //iEvent.getByToken(trackCollection_label,thePATTrackHandle);  //Tracks

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(primaryVertices_Label, vertices);    //primaryVertices
   /*is beam SpotLabel necessary ?? */
   //MC
   edm::Handle<reco::GenParticleCollection> pruned;
   iEvent.getByToken(genCands_, pruned);
   // Packed particles are all the status 1, so usable to remake jets
   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
   edm::Handle<pat::PackedGenParticleCollection> packed;
   iEvent.getByToken(packedGenToken_,packed);
    
   if (vertices->empty()) return; //if there is no vertices, return
   reco::VertexCollection::const_iterator PV = vertices->end();
   //reco::Vertex bestPtVtx;
   //loop over all the possible vertices in the container, starting with the one with higger PT
   int vertexRef_i = -1;
   double pvChi2 = -1;
   for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
      vertexRef_i++;
      if ( !(vtx->isFake())                                  // if is not fake
           && vtx->ndof()>=4. && vtx->position().Rho() < 2.0 // and number of degres of freedeom (fit) is > than 4
           && fabs(vtx->position().Z()) < 24.0) {            // and is near the beeam pipe we select it as our primary vertex PV
           PV = vtx;
           pvChi2 = vtx->normalizedChi2();
           break;
           }
   }

    
    TLorentzVector gen_z_p4,gen_jpsi_p4,gen_muon1_p4,gen_muon2_p4,gen_lepton1_p4,gen_lepton2_p4;
    TVector3       gen_z_vtx,gen_jpsi_vtx;

    gen_z_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_z_vtx.SetXYZ(0.,0.,0.);
    gen_jpsi_vtx.SetXYZ(0.,0.,0.);
    int n_Z_dau = 0;
    int Event_Cand = 1;
    if (isMC4l_){
        TLorentzVector temp_mu1P, temp_mu2P, temp_mu1N, temp_mu2N, temp_di1, temp_di2, temp_di3, temp_di4;
        temp_mu1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        temp_mu1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        temp_mu2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        temp_mu2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
        if (packed.isValid()) {
            for(size_t k=0; k<packed->size(); k++){
                const reco::Candidate * dau = &(*packed)[k];
                int stable_id = (*packed)[k].pdgId();
                if(!isAncestor(23, dau)) continue;
                if(stable_id ==  13 && temp_mu1N.M() == 0) temp_mu1N.SetPtEtaPhiM(dau->pt(), dau->eta(), dau->phi(), dau->mass());
                else if(stable_id == -13 && temp_mu1P.M() == 0) temp_mu1P.SetPtEtaPhiM(dau->pt(), dau->eta(), dau->phi(), dau->mass());
                else if(stable_id ==  13 && temp_mu2N.M() == 0) temp_mu2N.SetPtEtaPhiM(dau->pt(), dau->eta(), dau->phi(), dau->mass());
                else if(stable_id == -13 && temp_mu2P.M() == 0) temp_mu2P.SetPtEtaPhiM(dau->pt(), dau->eta(), dau->phi(), dau->mass());
            }//end for packed
            std::cout << temp_mu1N.Pt() << " | " << temp_mu2N.Pt() << " | " << temp_mu1P.Pt() << " | " << temp_mu2P.Pt() << std::endl;
            temp_di1 = temp_mu1P + temp_mu1N;
            temp_di2 = temp_mu2P + temp_mu2N;
            
            temp_di3 = temp_mu1P + temp_mu2N;
            temp_di4 = temp_mu2P + temp_mu1N;
            if (abs(temp_di1.M() - temp_di2.M()) > abs(temp_di3.M() - temp_di4.M())) {
                if (temp_di1.M() > temp_di2.M()){
                    gen_lepton1_p4 = temp_mu1N;
                    gen_lepton2_p4 = temp_mu1P;
                    gen_muon1_p4 = temp_mu2N;
                    gen_muon2_p4 = temp_mu2P;
                }
                else{
                    gen_lepton1_p4 = temp_mu2N;
                    gen_lepton2_p4 = temp_mu2P;
                    gen_muon1_p4 = temp_mu1N;
                    gen_muon2_p4 = temp_mu1P;
                }
            }
            else{
                if (temp_di3.M() > temp_di4.M()){
                    gen_lepton1_p4 = temp_mu2N;
                    gen_lepton2_p4 = temp_mu1P;
                    gen_muon1_p4 = temp_mu1N;
                    gen_muon2_p4 = temp_mu2P;
                }
                else{
                    gen_lepton1_p4 = temp_mu1N;
                    gen_lepton2_p4 = temp_mu2P;
                    gen_muon1_p4 = temp_mu1N;
                    gen_muon2_p4 = temp_mu2P;
                }
            }
            gen_z_p4 = gen_lepton1_p4 + gen_lepton2_p4 + gen_muon1_p4 + gen_muon2_p4;
            gen_jpsi_p4 = gen_muon1_p4 + gen_muon2_p4;
        }// end if packed valid
        std::cout << "found Z: "<< gen_z_p4.M() << std::endl;
    }
    else{
        //std::cout << "test" << std::endl;
        if ( pruned.isValid() ) {
          int foundit = 0;
          //std::cout << "MC ok " << std::endl;

          for (size_t i=0; i<pruned->size(); i++) {
             foundit = 0;
             const reco::Candidate *dau = &(*pruned)[i];
             ///ndau = dau->numberOfDaughters();

             if ( (abs(dau->pdgId()) == 23) ) { //&& (dau->status() == 2) ) { //found Z
                foundit++;
                //const reco::Candidate * Zboson = dau;
                gen_z_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
                gen_z_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
                n_Z_dau = dau->numberOfDaughters();
                if (n_Z_dau!=3) continue;
                //std::cout << " Z daugh: " << dau->numberOfDaughters() << std::endl;
                for (size_t k=0; k<dau->numberOfDaughters(); k++) {
                  const reco::Candidate *gdau = dau->daughter(k);
                  //std::cout << "MC Z daughter pdgID: " << gdau->pdgId() << std::endl;
                  if (gdau->pdgId()==443 ) { //&& gdau->status()==2) {   //found jpsi
                    foundit++;
                    gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
                    gen_jpsi_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                    int nm = 0;
                    for (size_t k=0; k<packed->size(); k++) {
                       //const reco::Candidate * dauInPrunedColl = (*packed)[k].mother(0);
                       const reco::Candidate * dauInPrunedColl = &(*packed)[k];
                       int stable_id = (*packed)[k].pdgId();
                       if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) {
                          //if (ndau<1) std::cout << (*packed)[k].pdgId() << " ";
                          //std::cout<<" phi = "<< gdau->pdgId()<< " daughter ID " << stable_id << std::endl;
                          if(stable_id == 13) { //found muon-
                                  gen_muon1_p4.SetPtEtaPhiM(dauInPrunedColl->pt(),dauInPrunedColl->eta(),dauInPrunedColl->phi(),dauInPrunedColl->mass());
                                  nm++;
                                 // std::cout<< "works K+ "<< dauInPrunedColl->mass() <<std::endl;
                          }
                           if(stable_id == -13){ //found muon+
                               gen_muon2_p4.SetPtEtaPhiM(dauInPrunedColl->pt(),dauInPrunedColl->eta(),dauInPrunedColl->phi(),dauInPrunedColl->mass());
                                  nm++;
                                 // std::cout<< "works K- "<< dauInPrunedColl->mass() << std::endl;
                          }
                       }
                    }
        
                  }//end found pai
                  if (gdau->pdgId()==13 ) {// pdgid for muon=13
                     foundit++;
                 gen_lepton1_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                 }
                  if (gdau->pdgId()==-13 ) {// pdgid for muon+=13
                     foundit++;
                     gen_lepton2_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                 }
                }// end number of daughters
             } //endif found Z
          }//end for

        } //end pruned
    }//end mc jpsi
    //NEW for muons and psi pairs
    int nonia = dimuons->size();
    int nmuons = leptons->size();
    int nPV    = vertices->size();
    
   //double chiVtxSqdProb = ChiSquaredProbability((double)(PV.chi2()),(double)(PV.ndof())); 
   //int breaker = 0;
   // We Cycle over dileptons for the Kfit
   //for (pat::CompositeCandidateCollection::const_iterator dilepton = dileptons->begin(); dilepton != dileptons->end() /*test && breaker < 10*/; ++dilepton){
   for (View<pat::Muon>::const_iterator lepton1 = leptons->begin(); lepton1 != leptons->end(); ++lepton1 ) {
     reco::TrackRef track_1 = lepton1->track();
     if (track_1.isNull()) continue;
     //std::cout << "PASS MU 1, LOOPING OVER MU 2" << std::endl;
     for (View<pat::Muon>::const_iterator lepton2 = leptons->begin() ; lepton2 !=  leptons->end(); ++lepton2 ) {
         if((lepton1->charge() * lepton2->charge()) != -1) continue;
         if(lepton1==lepton2) continue; //v7 add
         reco::TrackRef track_2 = lepton2->track();
         if (track_2.isNull()) continue;
         if(!(track_1->quality(reco::TrackBase::highPurity))) continue; //v7
         if(!(track_2->quality(reco::TrackBase::highPurity))) continue; //v7
         const pat::Muon *lept1 = 0;
         const pat::Muon *lept2 = 0;
         if(lepton1->charge() == -1 && lepton2->charge() == 1){ lept1 = &(*lepton1); lept2 = &(*lepton2);}
         else if (lepton1->charge() == 1 && lepton2->charge() == -1) {lept1 = &(*lepton2); lept2 = &(*lepton1);}
         else continue;
        //const pat::Muon* lept1 = dynamic_cast<const pat::Muon*>(dilepton->daughter("lepton1"));
        //const pat::Muon* lept2 = dynamic_cast<const pat::Muon*>(dilepton->daughter("lepton2"));
        //check if the muons came from PV
        //int pass1 = 0;
        //int pass2 = 0;
        // std::cout << "Reading muons" <<std::endl;
        //Check if track came from PV NOT WORKING FOR MUONS
        /*
        int l1_pass = 0;
        int l2_pass = 0;
        
        reco::VertexRef pvRef1 = lept1->vertexRef();
        reco::VertexRef pvRef2 = lept2->vertexRef();
        int key1 = (int) pvRef1.key();
        int key2 = (int) pvRef2.key();
        if (key1 == vertexRef_i ) l1_pass++;
        if (key2 == vertexRef_i ) l2_pass++;
        */
        /*
        for ( std::vector<reco::TrackBaseRef >::const_iterator tempTrack = bestPtVtx.tracks_begin();
           tempTrack != bestPtVtx.tracks_end(); ++tempTrack) {
           reco::TrackRef trackRef = tempTrack->castTo<reco::TrackRef>();
           if(lept1->track().key() == trackRef.key() ){
              pass1++;
           }
           if(lept2->track().key() == trackRef.key() ){
             pass2++;
           }
        }*/
        //std::cout << "pass: "<< pass << std::endl;
        //if (pass < 1) continue; //if at least one muon comes from PV, originaly 2 of them
        //std::cout<< "PASS"<< std::endl;
        ///////////////////////////////////////
        //////////// MY MUON ID ///////////////
        ///////////////////////////////////////
        int ZLe1Qid = 0;
        int ZLe2Qid = 0;
        // Lepton 1 (from Z)
	    if ( lept1->isGlobalMuon()){
	         ZLe1Qid = 1;
	         //std::cout << "1L is Global" << std::endl;
	    }
	    if ( lept1->isLooseMuon()){
	       ZLe1Qid += 10;
	      // std::cout << "1L is Loose " << std::endl;
	    }
	    if ( lept1->isMediumMuon()){
	       ZLe1Qid += 100;
	     //  std::cout << "1L is Medium " << std::endl;
	    }
	    if ( lept1->isTightMuon(*PV)){
	       ZLe1Qid += 1000;
	     //  std::cout << "1L is Tight " << std::endl;
	    }
	    if ( lept1->isSoftMuon(*PV)){
	       ZLe1Qid += 10000;
	     //  std::cout << "1L is Soft " << std::endl;
	    }
	    if ( lept1->isHighPtMuon(*PV)){
	       ZLe1Qid += 100000;
	     //  std::cout << "1L is HighPt " << std::endl;
	    }
	    if ( lept1->isPFMuon()){
	       ZLe1Qid += 1000000;
	     //    std::cout << "1L is ParticleFlow " << std::endl;
	    }
	    if ( lept1->isTrackerMuon()){
	       ZLe1Qid += 10000000;
	     //  std::cout << "1L is HighPt " << std::endl;
	    }
	    // Lepton 2 (from Z)
	    if ( lept2->isGlobalMuon()){
	       ZLe2Qid = 1;
	       // std::cout << "2L is Global " << std::endl;
	    }
	    if ( lept2->isLooseMuon()){
	       ZLe2Qid += 10;
	       //std::cout << "2L is Loose " << std::endl;
	    }
	    if ( lept2->isMediumMuon()){
	       ZLe2Qid += 100;
	       //std::cout << "2L is Medium " << std::endl;
	    }
	    if ( lept2->isTightMuon(*PV)){
	       ZLe2Qid += 1000;
	       //std::cout << "2L is Tight " << std::endl;
	    }
	    if ( lept2->isSoftMuon(*PV)){
	       ZLe2Qid += 10000;
	       //std::cout << "2L is Soft " << std::endl;
	    }
	    if ( lept2->isHighPtMuon(*PV)){
	       ZLe2Qid += 100000;
	       //std::cout << "2L is HighPt " << std::endl;
	    }
	    if ( lept2->isPFMuon()){
	       ZLe2Qid += 1000000;
	       //std::cout << "2L is Global " << std::endl;
	    }
	    if ( lept2->isTrackerMuon()){
	       ZLe2Qid += 10000000;
	     //  std::cout << "1L is HighPt " << std::endl;
	    }

	    int ZLe1_TrackerLWM = lept1->innerTrack()->hitPattern().trackerLayersWithMeasurement();
	    int ZLe1_PixelLWM   = lept1->innerTrack()->hitPattern().pixelLayersWithMeasurement();
	    int ZLe1_ValPixHit  = lept1->innerTrack()->hitPattern().numberOfValidPixelHits();

	    int ZLe2_TrackerLWM = lept2->innerTrack()->hitPattern().trackerLayersWithMeasurement();
	    int ZLe2_PixelLWM   = lept2->innerTrack()->hitPattern().pixelLayersWithMeasurement();
	    int ZLe2_ValPixHit  = lept2->innerTrack()->hitPattern().numberOfValidPixelHits();
         
        if (ZLe1_TrackerLWM < tlwm_) continue;
        if (ZLe2_TrackerLWM < tlwm_) continue;
        if (ZLe1_PixelLWM < plwm_) continue;
        if (ZLe2_PixelLWM < plwm_) continue;

  	    float ldxy1 = lept1->muonBestTrack()->dxy(PV->position());
	    float ldz1 = lept1->muonBestTrack()->dz(PV->position());

	    float ldxy2 = lept2->muonBestTrack()->dxy(PV->position());
	    float ldz2 = lept2->muonBestTrack()->dz(PV->position());
        if (abs(ldxy1) > dxyl_) continue;
        if (abs(ldxy2) > dxyl_) continue;
        if (abs(ldz1)  > dzl_ ) continue;
        if (abs(ldz2)  > dzl_ ) continue;
        
        reco::TrackRef glbTrack1 = lept1->track();
        reco::TrackRef glbTrack2 = lept2->track();
        if( glbTrack1.isNull() || glbTrack2.isNull() ) {
           //std::cout << "continue due to no track ref" << endl;
           continue;
        }
        //std::cout << "pass glbMuonTracks" <<std::endl;
        int trackPurity1; int trackPurity2;
        if(glbTrack1->quality(reco::TrackBase::loose))      {trackPurity1  = 1;}
        if(glbTrack1->quality(reco::TrackBase::tight))      {trackPurity1 += 10;}
        if(glbTrack1->quality(reco::TrackBase::highPurity)) {trackPurity1 += 100;}

        if(glbTrack2->quality(reco::TrackBase::loose))      {trackPurity2  = 1;}
        if(glbTrack2->quality(reco::TrackBase::tight))      {trackPurity2 += 10;}
        if(glbTrack2->quality(reco::TrackBase::highPurity)) {trackPurity2 += 100;}
	
        // Check if pair of leptons (muons) is coming from our Primary Vertex
        
        //std::vector<reco::TransientTrack> LLTTks;
        //LLTTks.push_back(theTTB->build(lept1->innerTrack()));
        //LLTTks.push_back(theTTB->build(lept2->innerTrack()));
        //std::pair<bool,Measurement1D> tkPVdistel1 = IPTools::absoluteImpactParameter3D(LLTTks.at(0),*PV);
        //std::pair<bool,Measurement1D> tkPVdistel2 = IPTools::absoluteImpactParameter3D(LLTTks.at(1),*PV);
        //if (!tkPVdistel1.first|| !tkPVdistel2.first ) continue;

        for (pat::CompositeCandidateCollection::const_iterator dimuon = dimuons->begin(); dimuon != dimuons->end() /*test && breaker < 10*/; ++dimuon){
            const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon1"));
            const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon2"));
            
            //check if the muons came from PV
            //int pass1 = 0;
            //int pass2 = 0;
            if (lept1 == muon1 || lept1==muon2 || lept2 == muon1 || lept2 == muon2  ) {
                //std::cout << "blocks repeating by same muon, it works!" << std::endl;
                continue;
            }

            if (lept1->innerTrack()==muon1->innerTrack() || lept1->innerTrack()==muon2->innerTrack() || lept2->innerTrack()==muon1->innerTrack() || lept2->innerTrack()==muon2->innerTrack() ) {
                //std:: cout << "blocks repeating by same track, it works !" << std::endl;
                continue;
                
            }
            //Isolation
            double dR1 = -1, dR2 = -1, dR3 = -1, dR4 = -1, dR5 = -1, dR6 = -1;
            double dRm1m2 = -1, dRl1l2 = -1, dRm1l1 = -1, dRm1l2 = -1, dRm2l1 = -1, dRm2l2 = -1;

            dR1 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
            dR2 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
            dR3 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
            dR4 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
            dR5 = deltaR(*(lept1->innerTrack()), *(lept2->innerTrack()));
            dR6 = deltaR(*(muon1->innerTrack()), *(muon2->innerTrack()));
             
            dRm1m2 = dR6;
            dRl1l2 = dR5;
            if (lept1->pt() > lept2->pt()){
                if(muon1->pt() > muon2->pt()){
                    dRm1l1 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
                    dRm2l1 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
                    dRm1l2 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
                    dRm2l2 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
                }
                else{
                    dRm1l1 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
                    dRm2l1 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
                    dRm1l2 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
                    dRm2l2 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
                }
            }
            else{
                if(muon1->pt() > muon2->pt()){
                    dRm1l1 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
                    dRm2l1 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
                    dRm1l2 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
                    dRm2l2 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
                }
                else{
                    dRm1l1 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
                    dRm2l1 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
                    dRm1l2 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
                    dRm2l2 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
                }
            }
            if ( dR1<0.01 || dR2<0.01 || dR3<0.01 ||dR4<0.01 ) continue;

            //////////// MY MUON ID ///////////////
            ///////////////////////////////////////
            int ZMu1Qid = 0;
            int ZMu2Qid = 0;
            // Muon 1 (from Jpsi)
            if ( muon1->isGlobalMuon()){
                ZMu1Qid = 1;
                //std::cout << "1L is Global" << std::endl;
            }
            if ( muon1->isLooseMuon()){
            ZMu1Qid += 10;
            // std::cout << "1L is Loose " << std::endl;
            }
            if ( muon1->isMediumMuon()){
            ZMu1Qid += 100;
            //  std::cout << "1L is Medium " << std::endl;
            }
            if ( muon1->isTightMuon(*PV)){
            ZMu1Qid += 1000;
            //  std::cout << "1L is Tight " << std::endl;
            }
            if ( muon1->isSoftMuon(*PV)){
            ZMu1Qid += 10000;
            //  std::cout << "1L is Soft " << std::endl;
            }
            if ( muon1->isHighPtMuon(*PV)){
            ZMu1Qid += 100000;
            //  std::cout << "1L is HighPt " << std::endl;
            }
            if ( muon1->isPFMuon()){
            ZMu1Qid += 1000000;
            //    std::cout << "1L is ParticleFlow " << std::endl;
            }
            if ( muon1->isTrackerMuon()){
            ZMu1Qid += 10000000;
            //  std::cout << "1L is HighPt " << std::endl;
            }
            // Muon 2 (from Jpsi)
            if ( muon2->isGlobalMuon()){
            ZMu2Qid = 1;
            // std::cout << "2L is Global " << std::endl;
            }
            if ( muon2->isLooseMuon()){
            ZMu2Qid += 10;
            //std::cout << "2L is Loose " << std::endl;
            }
            if ( muon2->isMediumMuon()){
            ZMu2Qid += 100;
            //std::cout << "2L is Medium " << std::endl;
            }
            if ( muon2->isTightMuon(*PV)){
            ZMu2Qid += 1000;
            //std::cout << "2L is Tight " << std::endl;
            }
            if ( muon2->isSoftMuon(*PV)){
            ZMu2Qid += 10000;
            //std::cout << "2L is Soft " << std::endl;
            }
            if ( muon2->isHighPtMuon(*PV)){
            ZMu2Qid += 100000;
            //std::cout << "2L is HighPt " << std::endl;
            }
            if ( muon2->isPFMuon()){
            ZMu2Qid += 1000000;
            //std::cout << "2L is Global " << std::endl;
            }
            if ( muon2->isTrackerMuon()){
            ZMu2Qid += 10000000;
            //  std::cout << "1L is HighPt " << std::endl;
            }
              
            /////////////////////////////////////////////////////
            //    T r a n s i e n t   T r a c k s  b u i l d e r //
            //////////////////////////////////////////////////////
            reco::TrackRef muTrack1 = muon1->innerTrack();
            reco::TrackRef muTrack2 = muon2->innerTrack();
            if(!(muTrack1->quality(reco::TrackBase::highPurity))) continue; //v7
            if(!(muTrack1->quality(reco::TrackBase::highPurity))) continue; //v7
            reco::TrackRef JpsiTk[2] = {muTrack1, muTrack2};
            std::vector<reco::TransientTrack> MuMuTTks;
            MuMuTTks.push_back(theTTB->build(&JpsiTk[0]));
            MuMuTTks.push_back(theTTB->build(&JpsiTk[1]));
            std::pair<bool,Measurement1D> tkPVdist1 = IPTools::absoluteImpactParameter3D(MuMuTTks.at(0),*PV);
            std::pair<bool,Measurement1D> tkPVdist2 = IPTools::absoluteImpactParameter3D(MuMuTTks.at(1),*PV);
            if (!tkPVdist1.first || !tkPVdist2.first ) continue;
            if (fabs(tkPVdist1.second.significance()) > ImparSigm_) continue;
            if (fabs(tkPVdist2.second.significance()) > ImparSigm_) continue;
            reco::TrackRef dilepTk[2]={lept1->innerTrack(), lept2->innerTrack()};
            //build the dimuon secondary vertex
            std::vector<reco::TransientTrack> LLTTks;
            LLTTks.push_back(theTTB->build(&dilepTk[0]));
            LLTTks.push_back(theTTB->build(&dilepTk[1]));
            std::pair<bool,Measurement1D> tkPVdistel1 = IPTools::absoluteImpactParameter3D(LLTTks.at(0),*PV);
            std::pair<bool,Measurement1D> tkPVdistel2 = IPTools::absoluteImpactParameter3D(LLTTks.at(1),*PV);
            if (!tkPVdistel1.first || !tkPVdistel2.first) continue;
            if (fabs(tkPVdistel1.second.significance()) > ImparSigl_) continue;
            if (fabs(tkPVdistel2.second.significance()) > ImparSigl_) continue;
            const ParticleMass muMass    = 0.10565837;
            float muSigma                = muMass*1.e-6;
             
            //Z Kinematic Fit
            KinematicParticleFactoryFromTransientTrack pFactory;
             
            std::vector<RefCountedKinematicParticle> ZDaughters;
            ZDaughters.push_back(pFactory.particle (MuMuTTks[0], muMass, float(0), float(0), muSigma));
            ZDaughters.push_back(pFactory.particle (MuMuTTks[1], muMass, float(0), float(0), muSigma));
            ZDaughters.push_back(pFactory.particle (LLTTks[0],  muMass, float(0), float(0), muSigma));
            ZDaughters.push_back(pFactory.particle (LLTTks[1],  muMass, float(0), float(0), muSigma));

            KinematicParticleVertexFitter ZVertexFitter;
            RefCountedKinematicTree ZTree = ZVertexFitter.fit(ZDaughters);
             
            if (ZTree->isEmpty())continue;
            ZTree->movePointerToTheTop();
            RefCountedKinematicParticle fitZ = ZTree->currentParticle();
            RefCountedKinematicVertex ZDecayVertex = ZTree->currentDecayVertex();
            //if the Fit is valid Fill the Tree and provide the Zcandidates
            int passFit = 0;
            float ZM_fit    = 0;
            float ZPx_fit   = 0;
            float ZPy_fit   = 0;
            float ZPz_fit   = 0;
            float ZVtxX_fit = 0;
            float ZVtxY_fit = 0;
            float ZVtxZ_fit = 0;
            float ZVtxP_fit = 0;
             
            if (fitZ->currentState().isValid()) {
                ZM_fit  = fitZ->currentState().mass();
                ZPx_fit = fitZ->currentState().kinematicParameters().momentum().x();
                ZPy_fit = fitZ->currentState().kinematicParameters().momentum().y();
                ZPz_fit = fitZ->currentState().kinematicParameters().momentum().z();
                ZVtxX_fit = ZDecayVertex->position().x();
                ZVtxY_fit = ZDecayVertex->position().y();
                ZVtxZ_fit = ZDecayVertex->position().z();
                //in some cases we will be saving NaN due to a negative chiSquared
                //if ( (double) ZDecayVertex->chiSquared() < 0 ) continue;
		        //std::cout << "WACHAMEEEE chisqd= "<< ZDecayVertex->chiSquared() << "dof = " << ZDecayVertex->degreesOfFreedom() <<std::endl;
                ZVtxP_fit = ChiSquaredProbability((double)(ZDecayVertex->chiSquared()),
						       (double)(ZDecayVertex->degreesOfFreedom()));
                //if (ZVtxP_fit <= 0.0) continue;
                passFit = 1;
           }
           else {
                ZM_fit    = 0;
                ZPx_fit   = 0;
                ZPy_fit   = 0;
                ZPz_fit   = 0;
                ZVtxX_fit = 0;
                ZVtxY_fit = 0;
                ZVtxZ_fit = 0;
                ZVtxP_fit = ChiSquaredProbability((double)(ZDecayVertex->chiSquared()), (double)(ZDecayVertex->degreesOfFreedom()));
           }
		   float Z_px = muon1->px()+muon2->px()+lept1->px()+lept2->px();
		   float Z_py = muon1->py()+muon2->py()+lept1->py()+lept2->py();
		   float Z_pz = muon1->pz()+muon2->pz()+lept1->pz()+lept2->pz();
           TLorentzVector m1, m2, l1, l2;
		   //k1.SetPtEtaPhiM(iTrack1->pt(), iTrack1->eta(),iTrack1->phi(),iTrack1->mass()); //miniAOD
		   //k2.SetPtEtaPhiM(iTrack2->pt(), iTrack2->eta(),iTrack2->phi(),iTrack2->mass());
		   m1.SetPtEtaPhiM(muon1->pt(), muon1->eta(),muon1->phi(),muon1->mass());
		   m2.SetPtEtaPhiM(muon2->pt(), muon2->eta(),muon2->phi(),muon2->mass());
		   l1.SetPtEtaPhiM(lept1->pt(), lept1->eta(),lept1->phi(),lept1->mass());
		   l2.SetPtEtaPhiM(lept2->pt(), lept2->eta(),lept2->phi(),lept2->mass());
		   TLorentzVector theZ = m1+m2+l1+l2;
                
  		   float mdxy1 = -1, mdz1 = -1, mdxy2 = -1, mdz2 = -1;

           mdxy1 = muon1->muonBestTrack()->dxy(PV->position());
           mdz1  = muon1->muonBestTrack()->dz(PV->position());
		   mdxy2 = muon2->muonBestTrack()->dxy(PV->position());
           mdz2  = muon2->muonBestTrack()->dz(PV->position());

	       if (abs(mdxy1) > dxym_ ) continue;
           if (abs(mdxy2) > dxym_) continue;
           if (abs(mdz1)  > dzm_ ) continue;
           if (abs(mdz2)  > dzm_ ) continue;

           if (ZM_fit < 60.0) continue;
           if (ZM_fit > 220.0) continue;
		   reco::CompositeCandidate recoZ(0, math::XYZTLorentzVector(ZPx_fit, ZPy_fit, ZPz_fit,
			      sqrt(ZM_fit*ZM_fit + ZPx_fit*ZPx_fit + ZPy_fit*ZPy_fit +
			      ZPz_fit*ZPz_fit)), math::XYZPoint(ZVtxX_fit,
			       ZVtxY_fit, ZVtxZ_fit), 23);

		   reco::CompositeCandidate msrdZ(0, math::XYZTLorentzVector( Z_px, Z_py, Z_pz,
			      sqrt(theZ.M()*theZ.M() + Z_px*Z_px + Z_py*Z_py +
			      Z_pz*Z_pz)), math::XYZPoint(ZVtxX_fit,
			       ZVtxY_fit, ZVtxZ_fit), 23);
		   pat::CompositeCandidate patMZ(msrdZ);
		   pat::CompositeCandidate patZ(recoZ);
           //New
           patZ.addUserInt("passFit_", passFit);
           patZ.addUserInt("nonia_", nonia );
           patZ.addUserInt("nmuons_",nmuons);
           patZ.addUserInt("nPV_",   nPV   );
           
           patZ.addUserFloat("vProb",ZVtxP_fit);
           patZ.addUserFloat("vChi2",ZDecayVertex->chiSquared());
            
           patZ.addUserFloat("pvChi2_", pvChi2);
           patZ.addUserInt("pvIndex", vertexRef_i);
		   patZ.addUserFloat("ZvtxX",ZVtxX_fit);
		   patZ.addUserFloat("ZvtxY",ZVtxY_fit);
		   patZ.addUserFloat("ZvtxZ",ZVtxZ_fit);
		   patZ.addUserFloat("dRm1m2",dR6);
		   patZ.addUserFloat("dRl1l2",dR5);
		   patZ.addUserFloat("dRl1m1",dR1);
		   patZ.addUserFloat("dRl1m2",dR2);
		   patZ.addUserFloat("dRl2m1",dR3);
		   patZ.addUserFloat("dRl2m2",dR4);

           patZ.addUserFloat("dRm1m2_",dRm1m2);
           patZ.addUserFloat("dRl1l2_",dRl1l2);
           patZ.addUserFloat("dRl1m1_",dRm1l1);
           patZ.addUserFloat("dRl1m2_",dRm2l1);
           patZ.addUserFloat("dRl2m1_",dRm1l2);
           patZ.addUserFloat("dRl2m2_",dRm2l2);

		   bool child = ZTree->movePointerToTheFirstChild();
		   //get first muon
		   RefCountedKinematicParticle fitMu1 = ZTree->currentParticle();
           float mu1M_fit ;
           float mu1Q_fit ;
           float mu1Px_fit;
           float mu1Py_fit;
           float mu1Pz_fit;
                 
		   if (!child){
               //std::cout << "M1" << std::endl;
               mu1M_fit  = 0;
               mu1Q_fit  = 0;
               mu1Px_fit = 0;
               mu1Py_fit = 0;
               mu1Pz_fit = 0;
           }
           else{
               mu1M_fit  = fitMu1->currentState().mass();
               mu1Q_fit  = fitMu1->currentState().particleCharge();
               mu1Px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
               mu1Py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
               mu1Pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
           }
		   reco::CompositeCandidate msrdM1(muon1->charge(), math::XYZTLorentzVector(muon1->px(), muon1->py(), muon1->pz(),
                   sqrt((muMass)*(muMass) + (muon1->px())*(muon1->px()) + (muon1->py())*(muon1->py()) +
                   (muon1->pz())*(muon1->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));

		   pat::CompositeCandidate pat_msrdM1(msrdM1);
		   reco::CompositeCandidate recoM1(mu1Q_fit, math::XYZTLorentzVector(mu1Px_fit, mu1Py_fit, mu1Pz_fit,
					     sqrt(mu1M_fit*mu1M_fit + mu1Px_fit*mu1Px_fit + mu1Py_fit*mu1Py_fit +
					     mu1Pz_fit*mu1Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));

		   pat::CompositeCandidate patM1(recoM1);

		   patM1.addUserFloat("Dxy",mdxy1);
		   patM1.addUserFloat("Dz",mdz1);
           //new
           patM1.addUserFloat("dRIso"    ,getIso( *muon1 ) );
           //v7 new
           //patM1.addUserFloat("isolationR03_",          muon1->isolationR03());
           //patM1.addUserFloat("isolationR05_",          muon1->isolationR05());
           //patM1.addUserFloat("pfIsolationR03_",        muon1->pfIsolationR03());
           //patM1.addUserFloat("pfMeanDRIsoProfileR03_", muon1->pfMeanDRIsoProfileR03());
           //patM1.addUserFloat("pfSumDRIsoProfileR03_",  muon1->pfSumDRIsoProfileR03());
           //patM1.addUserFloat("pfIsolationR04_",        muon1->pfIsolationR04());
           //patM1.addUserFloat("pfMeanDRIsoProfileR04_", muon1->pfMeanDRIsoProfileR04());
           //patM1.addUserFloat("pfSumDRIsoProfileR04_",  muon1->pfSumDRIsoProfileR04());

            
            
           patM1.addUserFloat("dIP3DSig", tkPVdist1.second.significance());
        
		   patM1.addUserFloat("dIP3D",tkPVdist1.second.value());
		   patM1.addUserFloat("dIP3DErr",tkPVdist1.second.error());

		   int psiM1_TrackerLWM = muon1->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
		   int psiM1_PixelLWM   = muon1->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
		   int psiM1_ValPixHit  = muon1->muonBestTrack()->hitPattern().numberOfValidPixelHits();
		   int psiM2_TrackerLWM = muon2->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
		   int psiM2_PixelLWM   = muon2->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
		   int psiM2_ValPixHit  = muon2->muonBestTrack()->hitPattern().numberOfValidPixelHits();

           if (psiM1_TrackerLWM <= tlwm_) continue;
           if (psiM2_TrackerLWM <= tlwm_) continue;
           if (psiM1_PixelLWM <= plwm_) continue;
           if (psiM2_PixelLWM <= plwm_) continue;
           patM1.addUserInt("ZMu1Qid_", ZMu1Qid);
		   patM1.addUserFloat("psiM1_TrackerLWM_", psiM1_TrackerLWM);
		   patM1.addUserFloat("psiM1_PixelLWM_",  psiM1_PixelLWM);
		   patM1.addUserFloat("psiM1_ValPixHit_", psiM1_ValPixHit);
           ///////////////////////////////////////////////////////
           //get second Muon
		   child = ZTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle fitMu2 = ZTree->currentParticle();
           float mu2M_fit ;
           float mu2Q_fit ;
           float mu2Px_fit;
           float mu2Py_fit;
           float mu2Pz_fit;

           if (!child){
               //std::cout << "Mu2" << std::endl;
               mu2M_fit  = 0;
               mu2Q_fit  = 0;
               mu2Px_fit = 0;
               mu2Py_fit = 0;
               mu2Pz_fit = 0;
           }
           else{
               mu2M_fit  = fitMu2->currentState().mass();
               mu2Q_fit  = fitMu2->currentState().particleCharge();
               mu2Px_fit =  fitMu2->currentState().kinematicParameters().momentum().x();
               mu2Py_fit =  fitMu2->currentState().kinematicParameters().momentum().y();
               mu2Pz_fit =  fitMu2->currentState().kinematicParameters().momentum().z();
           }
                 
		   reco::CompositeCandidate recoM2(mu2Q_fit, math::XYZTLorentzVector(mu2Px_fit, mu2Py_fit, mu2Pz_fit,
					     sqrt(mu2M_fit*mu2M_fit + mu2Px_fit*mu2Px_fit + mu2Py_fit*mu2Py_fit +
					     mu2Pz_fit*mu2Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));

           reco::CompositeCandidate msrdM2(muon2->charge(), math::XYZTLorentzVector(muon2->px(), muon2->py(), muon2->pz(),
                         sqrt((muMass)*(muMass) + (muon2->px())*(muon2->px()) + (muon2->py())*(muon2->py()) +
					     (muon2->pz())*(muon2->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));

		   pat::CompositeCandidate pat_msrdM2(msrdM2);

		   pat::CompositeCandidate patM2(recoM2);
		   patM2.addUserFloat("Dxy",mdxy2);
		   patM2.addUserFloat("Dz",mdz2);
           //New
           patM2.addUserFloat("dRIso"    ,getIso( *muon2 ) );
           //v7 new
           //patM2.addUserFloat("isolationR03_",          muon2->isolationR03());
           //patM2.addUserFloat("isolationR05_",          muon2->isolationR05());
           //patM2.addUserFloat("pfIsolationR03_",        muon2->pfIsolationR03());
           //patM2.addUserFloat("pfMeanDRIsoProfileR03_", muon2->pfMeanDRIsoProfileR03());
           //patM2.addUserFloat("pfSumDRIsoProfileR03_",  muon2->pfSumDRIsoProfileR03());
           //patM2.addUserFloat("pfIsolationR04_",        muon2->pfIsolationR04());
           //patM2.addUserFloat("pfMeanDRIsoProfileR04_", muon2->pfMeanDRIsoProfileR04());
           //patM2.addUserFloat("pfSumDRIsoProfileR04_",  muon2->pfSumDRIsoProfileR04());

           patM2.addUserFloat("dIP3DSig", tkPVdist2.second.significance());
                 
           patM2.addUserFloat("dIP3d", tkPVdist2.second.significance());
		   patM2.addUserFloat("dIP3D",tkPVdist2.second.value());
		   patM2.addUserFloat("dIP3DErr",tkPVdist2.second.error());
           
           patM2.addUserInt("ZMu2Qid_", ZMu2Qid);
		   patM2.addUserFloat("psiM2_TrackerLWM_", psiM2_TrackerLWM);
		   patM2.addUserFloat("psiM2_PixelLWM_",  psiM2_PixelLWM);
		   patM2.addUserFloat("psiM2_ValPixHit_", psiM2_ValPixHit);

           //Denfine Phi from to tracks-kaons
		   pat::CompositeCandidate psi;
		   psi.addDaughter(patM1,"muon1");
		   psi.addDaughter(patM2,"muon2");
		   psi.setP4(patM1.p4()+patM2.p4());
           //check how to acces this variables
           psi.addUserFloat("vProb", dimuon->userFloat("vProb"));
           psi.addUserFloat("vChi2", dimuon->userFloat("vNChi2"));
		   //save measured values on a different daughter
		   pat::CompositeCandidate msrd_psi;
		   msrd_psi.addDaughter(pat_msrdM1,"msrd_muon1");
		   msrd_psi.addDaughter(pat_msrdM2,"msrd_muon2");
		   msrd_psi.setP4(pat_msrdM1.p4()+pat_msrdM2.p4());
		   //get Lepton 1 (Muon1)
		   child = ZTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle fitL1 = ZTree->currentParticle();
           float L1M_fit ;
           float L1Q_fit ;
           float L1Px_fit;
           float L1Py_fit;
           float L1Pz_fit;

           if (!child){
               L1M_fit = 0;
               L1Q_fit = 0;
               L1Px_fit = 0;
               L1Py_fit = 0;
               L1Pz_fit = 0;
           }
           else{
               L1M_fit  = fitL1->currentState().mass();
               L1Q_fit  = fitL1->currentState().particleCharge();
               L1Px_fit =   fitL1->currentState().kinematicParameters().momentum().x();
               L1Py_fit =   fitL1->currentState().kinematicParameters().momentum().y();
               L1Pz_fit =   fitL1->currentState().kinematicParameters().momentum().z();
           }
		   reco::CompositeCandidate recoL1(L1Q_fit, math::XYZTLorentzVector(L1Px_fit, L1Py_fit, L1Pz_fit,
					       sqrt(L1M_fit*L1M_fit + L1Px_fit*L1Px_fit + L1Py_fit*L1Py_fit +
					       L1Pz_fit*L1Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
		   reco::CompositeCandidate msrdL1(lept1->charge(), math::XYZTLorentzVector(lept1->px(), lept1->py(), lept1->pz(),
					     sqrt((lept1->mass())*(lept1->mass()) + (lept1->px())*(lept1->px()) + (lept1->py())*(lept1->py()) +
					     (lept1->pz())*(lept1->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
		   pat::CompositeCandidate pat_msrdL1(msrdL1);

		   pat::CompositeCandidate patL1(recoL1);
		   patL1.addUserFloat("Dxy"        ,ldxy1);
		   patL1.addUserFloat("Dz"         ,ldz1);
		   patL1.addUserFloat("dRIso"      ,getIso( *lept1 ) );
           //New
           //v7 new
           //patL1.addUserFloat("isolationR03_",          lept1->isolationR03());
           //patL1.addUserFloat("isolationR05_",          lept1->isolationR05());
           //patL1.addUserFloat("pfIsolationR03_",        lept1->pfIsolationR03());
           //patL1.addUserFloat("pfMeanDRIsoProfileR03_", lept1->pfMeanDRIsoProfileR03());
           //patL1.addUserFloat("pfSumDRIsoProfileR03_",  lept1->pfSumDRIsoProfileR03());
           //patL1.addUserFloat("pfIsolationR04_",        lept1->pfIsolationR04());
           //patL1.addUserFloat("pfMeanDRIsoProfileR04_", lept1->pfMeanDRIsoProfileR04());
           //patL1.addUserFloat("pfSumDRIsoProfileR04_",  lept1->pfSumDRIsoProfileR04());

            
           patL1.addUserFloat("dIP3DSig"   ,tkPVdistel1.second.significance());
           
		   patL1.addUserFloat("dIP3D"      ,tkPVdistel1.second.value());
		   patL1.addUserFloat("dIP3DErr"   ,tkPVdistel1.second.error());

		   // MY Muon ID
		   patL1.addUserInt("ZLe1Qid_", ZLe1Qid);
		   patL1.addUserFloat("ZLe1_TrackerLWM_", ZLe1_TrackerLWM);
		   patL1.addUserFloat("ZLe1_PixelLWM_", ZLe1_PixelLWM);
		   patL1.addUserFloat("ZLe1_ValPixHit_", ZLe1_ValPixHit);
        
           //patL1.addUserFloat("pass", pass1);
           ////////////////////////////////////////////////////////
		   //get Lepton2 (Muon2)
		   child = ZTree->movePointerToTheNextChild();
		   RefCountedKinematicParticle fitL2 = ZTree->currentParticle();
           float L2M_fit ;
           float L2Q_fit ;
           float L2Px_fit;
           float L2Py_fit;
           float L2Pz_fit;

           if (!child){
               L2M_fit  = 0;
               L2Q_fit  = 0;
               L2Px_fit = 0;
               L2Py_fit = 0;
               L2Pz_fit = 0;
           }
           else {
               L2M_fit  = fitL2->currentState().mass();
               L2Q_fit  = fitL2->currentState().particleCharge();
               L2Px_fit =    fitL2->currentState().kinematicParameters().momentum().x();
               L2Py_fit =    fitL2->currentState().kinematicParameters().momentum().y();
               L2Pz_fit =    fitL2->currentState().kinematicParameters().momentum().z();
           }
		   reco::CompositeCandidate recoL2(L2Q_fit, math::XYZTLorentzVector(L2Px_fit, L2Py_fit, L2Pz_fit,
					       sqrt(L2M_fit*L2M_fit + L2Px_fit*L2Px_fit + L2Py_fit*L2Py_fit +
					       L2Pz_fit*L2Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
		   reco::CompositeCandidate msrdL2(lept2->charge(), math::XYZTLorentzVector(lept2->px(), lept2->py(), lept2->pz(),
					     sqrt((lept2->mass())*(lept2->mass()) + (lept2->px())*(lept2->px()) + (lept2->py())*(lept2->py()) +
					     (lept2->pz())*(lept2->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
		   pat::CompositeCandidate pat_msrdL2(msrdL2);

		   pat::CompositeCandidate patL2(recoL2);
		   patL2.addUserFloat("Dxy"        ,ldxy2);
		   patL2.addUserFloat("Dz"         ,ldz2);
		   patL2.addUserFloat("dRIso"      ,getIso( *lept2 ) );
           //v7 new
           //patL2.addUserFloat("isolationR03_",          lept2->isolationR03());
           //patL2.addUserFloat("isolationR05_",          lept2->isolationR05());
           //patL2.addUserFloat("pfIsolationR03_",        lept2->pfIsolationR03());
           //patL2.addUserFloat("pfMeanDRIsoProfileR03_", lept2->pfMeanDRIsoProfileR03());
           //patL2.addUserFloat("pfSumDRIsoProfileR03_",  lept2->pfSumDRIsoProfileR03());
           //patL2.addUserFloat("pfIsolationR04_",        lept2->pfIsolationR04());
           //patL2.addUserFloat("pfMeanDRIsoProfileR04_", lept2->pfMeanDRIsoProfileR04());
           //patL2.addUserFloat("pfSumDRIsoProfileR04_",  lept2->pfSumDRIsoProfileR04());
            
           patL2.addUserFloat("dIP3DSig"   ,tkPVdistel2.second.significance());
		   patL2.addUserFloat("dIP3D"      ,tkPVdistel2.second.value());
		   patL2.addUserFloat("dIP3DErr"   ,tkPVdistel2.second.error());
		   // My muon Q id
		   patL2.addUserInt("ZLe2Qid_", ZLe2Qid);

		   patL2.addUserFloat("ZLe2_TrackerLWM_", ZLe2_TrackerLWM);
		   patL2.addUserFloat("ZLe2_PixelLWM_", ZLe2_PixelLWM);
		   patL2.addUserFloat("ZLe2_ValPixHit_", ZLe2_ValPixHit);

           ///MONTECARLO
           reco::CompositeCandidate mc_Z(0, math::XYZTLorentzVector(gen_z_p4.Px(), gen_z_p4.Py(), gen_z_p4.Pz(),
               sqrt((gen_z_p4.M())*(gen_z_p4.M()) + (gen_z_p4.Px())*(gen_z_p4.Px()) + (gen_z_p4.Py())*(gen_z_p4.Py()) +
               (gen_z_p4.Pz())*(gen_z_p4.Pz()))), math::XYZPoint(gen_z_vtx.x(), gen_z_vtx.y(), gen_z_vtx.z()));
           pat::CompositeCandidate pat_mcZ(mc_Z);
           
           reco::CompositeCandidate mc_Psi(0, math::XYZTLorentzVector(gen_jpsi_p4.Px(), gen_jpsi_p4.Py(), gen_jpsi_p4.Pz(),
               sqrt((gen_jpsi_p4.M())*(gen_jpsi_p4.M()) + (gen_jpsi_p4.Px())*(gen_jpsi_p4.Px()) + (gen_jpsi_p4.Py())*(gen_jpsi_p4.Py()) +
               (gen_jpsi_p4.Pz())*(gen_jpsi_p4.Pz()))), math::XYZPoint(gen_jpsi_vtx.x(), gen_jpsi_vtx.y(), gen_jpsi_vtx.z()));
           pat::CompositeCandidate pat_mcPsi(mc_Psi);

           reco::CompositeCandidate mc_M1(1, math::XYZTLorentzVector(gen_muon1_p4.Px(), gen_muon1_p4.Py(), gen_muon1_p4.Pz(),
           sqrt((gen_muon1_p4.M())*(gen_muon1_p4.M()) + (gen_muon1_p4.Px())*(gen_muon1_p4.Px()) + (gen_muon1_p4.Py())*(gen_muon1_p4.Py()) +
               (gen_muon1_p4.Pz())*(gen_muon1_p4.Pz()))), math::XYZPoint(gen_jpsi_vtx.x(), gen_jpsi_vtx.y(), gen_jpsi_vtx.z()));
           pat::CompositeCandidate pat_mcM1(mc_M1);

           reco::CompositeCandidate mc_M2(-1, math::XYZTLorentzVector(gen_muon2_p4.Px(), gen_muon2_p4.Py(), gen_muon2_p4.Pz(),
           sqrt((gen_muon2_p4.M())*(gen_muon2_p4.M()) + (gen_muon2_p4.Px())*(gen_muon2_p4.Px()) + (gen_muon2_p4.Py())*(gen_muon2_p4.Py()) +
               (gen_muon2_p4.Pz())*(gen_muon2_p4.Pz()))), math::XYZPoint(gen_jpsi_vtx.x(), gen_jpsi_vtx.y(), gen_jpsi_vtx.z()));
           pat::CompositeCandidate pat_mcM2(mc_M2);
           
           reco::CompositeCandidate mc_L1(-1, math::XYZTLorentzVector(gen_lepton1_p4.Px(), gen_lepton1_p4.Py(), gen_lepton1_p4.Pz(),
           sqrt((gen_lepton1_p4.M())*(gen_lepton1_p4.M()) + (gen_lepton1_p4.Px())*(gen_lepton1_p4.Px()) + (gen_lepton1_p4.Py())*(gen_lepton1_p4.Py()) +
               (gen_lepton1_p4.Pz())*(gen_lepton1_p4.Pz()))), math::XYZPoint(gen_z_vtx.x(), gen_z_vtx.y(), gen_z_vtx.z()));
           pat::CompositeCandidate pat_mcL1(mc_L1);

           reco::CompositeCandidate mc_L2(1, math::XYZTLorentzVector(gen_lepton2_p4.Px(), gen_lepton2_p4.Py(), gen_lepton2_p4.Pz(),
           sqrt((gen_lepton2_p4.M())*(gen_lepton2_p4.M()) + (gen_lepton2_p4.Px())*(gen_lepton2_p4.Px()) + (gen_lepton2_p4.Py())*(gen_lepton2_p4.Py()) +
               (gen_lepton2_p4.Pz())*(gen_lepton2_p4.Pz()))), math::XYZPoint(gen_z_vtx.x(), gen_z_vtx.y(), gen_z_vtx.z()));
           pat::CompositeCandidate pat_mcL2(mc_L2);
           
        
		   patZ.addDaughter(psi,"psi");
		   patZ.addDaughter(patL1,"lepton1");
		   patZ.addDaughter(patL2,"lepton2");

        
           patZ.addDaughter(msrd_psi, "msrd_psi");
		   patZ.addDaughter(pat_msrdL1, "msrd_lepton1");
		   patZ.addDaughter(pat_msrdL2, "msrd_lepton2");
		   patZ.addDaughter(patMZ, "patMZ");
           
           patZ.addDaughter(pat_mcZ, "mcZ");
           patZ.addDaughter(pat_mcPsi, "mcPsi");
           patZ.addDaughter(pat_mcM1, "mcM1");
           patZ.addDaughter(pat_mcM2, "mcM2");
           patZ.addDaughter(pat_mcL1, "mcL1");
           patZ.addDaughter(pat_mcL2, "mcL2");
           patZ.addUserInt("Event_Cand_", Event_Cand);
           //std::cout<< "Z Cand Event: "<< Event_Cand << std::endl;
           Event_Cand++;

           ZCandColl->push_back(patZ);

        }//end of loop for dimuon

   }//end for dilepton
   }
   iEvent.put(std::move(ZCandColl),"ZCandidates");
}//end produce 

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
jpsiLepLepKmcFitter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
jpsiLepLepKmcFitter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
jpsiLepLepKmcFitter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
jpsiLepLepKmcFitter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
jpsiLepLepKmcFitter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
jpsiLepLepKmcFitter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
jpsiLepLepKmcFitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
bool jpsiLepLepKmcFitter::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool jpsiLepLepKmcFitter::IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk->eta());
  double DeltaP   = fabs(mu.p()-tk->p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
double jpsiLepLepKmcFitter::getIso(const pat::Muon& mu){
                Float_t coriso = 99.0;
                reco::MuonPFIsolation pfR03 = mu.pfIsolationR03();
                coriso = pfR03.sumChargedHadronPt + std::max(0., pfR03.sumNeutralHadronEt+pfR03.sumPhotonEt-0.5*pfR03.sumPUPt);
return coriso;
}

//define this as a plug-in
DEFINE_FWK_MODULE(jpsiLepLepKmcFitter);
