// -*- C++ -*-
//
// Package:    AnalyzeZphill/ZphiTupler
// Class:      jpsi4LepLepKmcFitter
// 
/**\class jpsi4LepLepKmcFitter jpsi4LepLepKmcFitter.cc AnalyzeZphill/ZphiTupler/plugins/jpsi4LepLepKmcFitter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Artemio Ayala Sanchez
//         Created:  Tue, 18 Jun 2019 19:12:13 GMT
//
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

#include <iostream>  
#include<string>  
#include <unordered_map> 

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

//
// class declaration
//

class jpsi4LepLepKmcFitter : public edm::stream::EDProducer<> {
   public:
      explicit jpsi4LepLepKmcFitter(const edm::ParameterSet&);
      ~jpsi4LepLepKmcFitter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
      bool IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu);
      double getIso(const pat::Muon& mu);
   private:
      void printMCtree(const reco::Candidate *, int);
      //recursive function to analyze a decay and match values to any of the 2-4 muon electrons and return the kind of decay channel
      void analyzeDecay(const reco::Candidate*, TLorentzVector&, TLorentzVector&, TLorentzVector&,
                        TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, TLorentzVector&, int&, int);
      std::string printName(int);
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
      edm::EDGetTokenT<pat::CompositeCandidateCollection> dilepton_Label;
      //edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;

      //edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
      edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;

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
jpsi4LepLepKmcFitter::jpsi4LepLepKmcFitter(const edm::ParameterSet& iConfig)
{
   //trackCollection_label = consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Tracks")); //miniAOD
    dimuon_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuon"));
    primaryVertices_Label = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"));
    dilepton_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dilepton"));
    //muonToken_ = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));

   //BSLabel_ = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotLabel"));
   genCands_ = consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"));
   packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
   
   produces<pat::CompositeCandidateCollection>("ZCandidates");
}


jpsi4LepLepKmcFitter::~jpsi4LepLepKmcFitter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
//Try to print mc Tree
std::string jpsi4LepLepKmcFitter::printName(int pdgid){
    /*
    std::string retstr;
    if(pdgid == 11){
        retstr = "e-";
    }
    else if (pdgid == -11){
        retstr = "e+";
    }
    else if (pdgid == 22){
        retstr = "gamma";
    }
    else if (pdgid == 13){
        retstr = "mu-";
    }
    else if (pdgid == -13){
        retstr = "mu+";
    }
    else if (pdgid == 23){
        retstr = "Z";
    }
    else {
        retstr = std::to_string(pdgid);
    }
    */
    std::unordered_map<int, std::string> umap; 
    umap[11]  = "e-";
    umap[-11] = "e+";
    umap[13]  = "mu-";
    umap[-13] = "mu+";
    umap[22]  = "gamma";
    umap[23]  = "Z";
    umap[12]  = "Ve";
    umap[14]  = "Vmu";
    umap[-14] = "-Vmu";
    umap[15]  = "tau-";
    umap[-15] = "tau+";
    umap[16]  = "Vtau";
    umap[-16] = "-Vtau";
    umap[111] = "pi0";
    umap[211] = "pi+";
    umap[-211]= "pi-";
    std::string retstr;
    if (umap.find(pdgid) == umap.end()){
        retstr = std::to_string(pdgid);
    }
    else{
        retstr = umap.at(pdgid);
    }

    return retstr;
}

void jpsi4LepLepKmcFitter::printMCtree(const reco::Candidate* mother, int indent=0){
    if (mother == NULL){
         std::cout << "end tree" << std::endl;
         return;
    }
    if (mother->numberOfDaughters() > 1){
        if(indent){
                std::cout << std::setw(indent) << " ";
        }  
        std::cout << printName(mother->pdgId()) <<" has "<< mother->numberOfDaughters() <<" daughters " <<std::endl;
    }
    int extraIndent = 0;    
    for (size_t i=0; i< mother->numberOfDaughters(); i++){
        const reco::Candidate * daughter = mother->daughter(i);
        if (mother->numberOfDaughters() > 1){
            if(indent){
                std::cout << std::setw(indent) << " ";
            }
            std::cout << " daugter "<< i+1 <<": "<<  printName(daughter->pdgId()) << " with Pt: ";
            std::cout << daughter->pt() << " | Eta: "<< daughter->eta() << " | Phi: "<< daughter->phi() << " | mass: "<< daughter->mass() << std::endl;
            extraIndent+=4;
        }
        if (daughter->numberOfDaughters()) printMCtree(daughter, indent+extraIndent);
    }
}
void jpsi4LepLepKmcFitter::analyzeDecay(const reco::Candidate* mother, TLorentzVector& muP1, TLorentzVector& muN1, TLorentzVector& muP2, TLorentzVector& muN2,
                                        TLorentzVector& elP1, TLorentzVector& elN1, TLorentzVector& elP2, TLorentzVector& elN2, int& decay, int indent = 0){

    int momID = mother->pdgId();
    if (mother == NULL){
         std::cout << "end tree" << std::endl;
         return;
    }
    if (mother->numberOfDaughters() > 1){
        if(indent){
                std::cout << std::setw(indent) << " ";
        }
        std::cout << printName(mother->pdgId()) <<" has "<< mother->numberOfDaughters() <<" daughters " <<std::endl;
    }
    int extraIndent = 0;
    for (size_t i=0; i< mother->numberOfDaughters(); i++){
        const reco::Candidate * daughter = mother->daughter(i);
        int dauID = daughter->pdgId();
        if(indent){
            std::cout << std::setw(indent) << " ";
        }
        std::cout << " daugter "<< i+1 <<": "<<  printName(dauID) << " with Pt: ";
        std::cout << daughter->pt() << " | Eta: "<< daughter->eta() << " | Phi: "<< daughter->phi() << " | mass: "<< daughter->mass() << std::endl;
        extraIndent+=4;
        if (daughter->numberOfDaughters() == 0) {
            // decay 1  Z--> 2mu
            // decay 2  Z--> 2el
            // decay 3  Z--> 4mu via gamma
            // decay 3  Z--> 4mu via Z*
            // decay 4  Z--> 4el via gamma
            // decay 4  Z--> 4el via Z*
            // decay 5  Z--> 2mu->2el via gamma
            // decay 5  Z--> 2mu->2el via Z*
            // decay 6  Z--> 2el->2mu via gamma
            // decay 6  Z--> 2el->2mi via Z*
            // decay 11 Z--> any with tau
            //std::cout<< "final state "<<std::endl;
            if (dauID == 22) continue;
            if (dauID == 11 and elN1.M() == 0){
              //  std::cout<< "saving electron 1 ... "<<std::endl;
                elN1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            else if (dauID == -11 and elP1.M() == 0){
                //std::cout<< "saving positron 1 ... "<<std::endl;
                elP1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            else if (dauID == 13 and muN1.M() == 0){
                //std::cout<< "saving muon 1 ... "<<std::endl;
                muN1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            else if (dauID == -13 and muP1.M() == 0){
                //std::cout<< "saving anti-muon 1 ... "<<std::endl;
                muP1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            //// lepton 2
            else if (dauID == 11 and elN2.M() == 0){
                //std::cout<< "saving electron 2 ... "<< std::endl;
                elN2.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            else if (dauID == -11 and elP2.M() == 0){
                //std::cout<< "saving positron 2 ... "<< std::endl;
                elP2.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            else if (dauID == 13 and muN2.M() == 0){
                //std::cout<< "saving muon 2 ... "<< std::endl;
                muN2.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            else if (dauID == -13 and muP2.M() == 0){
                //std::cout<< "saving anti-muon 2 ... "<< std::endl;
                muP2.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
            }
            if(elN1.M() == 0 && elP1.M() == 0 && muN1.M() != 0 && muP1.M() != 0 && elN2.M() == 0 && elP2.M() == 0 && muN2.M() == 0 && muP2.M() == 0){
                decay = 1;
                //std::cout<< "is Z --> 2 mu"<< std::endl;
            }
            else if(elN1.M() != 0 && elP1.M() != 0 && muN1.M() == 0 && muP1.M() == 0 && elN2.M() == 0 && elP2.M() == 0 && muN2.M() == 0 && muP2.M() == 0){
                decay = 2;
                //std::cout<< "is Z --> 2 el"<< std::endl;
            }
            else if(elN1.M() == 0 && elP1.M() == 0 && muN1.M() != 0 && muP1.M() != 0 && elN2.M() == 0 && elP2.M() == 0 && muN2.M() != 0 && muP2.M() != 0){
                decay = 3;// or 4
                std::cout<< "is Z --> 4mu via ??? "<< std::endl;
                if (momID == 23) std::cout << "via Z" << std::endl;
                else if (momID == 22) std::cout << "via Gamma"<<std::endl;
            }
            else if(elN1.M() != 0 && elP1.M() != 0 && muN1.M() == 0 && muP1.M() == 0 && elN2.M() != 0 && elP2.M() != 0 && muN2.M() == 0 && muP2.M() == 0){
                decay = 4; //or 6
                std::cout<< "is Z --> 4el via ??? "<< std::endl;
                //if (momID == 23) std::cout << "via Z" << std::endl;
                //else if (momID == 22) std::cout << "via Gamma"<<std::endl;
            }
            else if(elN1.M() != 0 && elP1.M() != 0 && muN1.M() != 0 && muP1.M() != 0 && elN2.M() == 0 && elP2.M() == 0 && muN2.M() == 0 && muP2.M() == 0){
                //decay = 5; //or 8 9 10 2mu->2el
                std::cout<< "is Z --> 2mu 2el or 2el 2mu ??? "<< std::endl;
                //if (momID == 23) std::cout << "via Z" << std::endl;
                //else if (momID == 22) std::cout << "via Gamma"<<std::endl;
                if (isAncestor(11, daughter) || isAncestor(-11, daughter)){
                    std::cout<< "Z-> 2el -> 2mu"<<std::endl;
                    decay = 5;
                }
                else if (isAncestor(13, daughter) || isAncestor(-13, daughter)){
                    std::cout<< "Z-> 2mu -> 2el"<<std::endl;
                    decay = 6;
                }
                else std::cout<< "cant identify ..."<< std::endl;
            }
            else{
                std::cout<<"DECAY NOT IDENTIFIED !!" << std::endl;
            }
            
        }
        else{
            analyzeDecay(daughter, muP1, muN1, muP2, muN2, elP1, elN1, elP2, elN2, decay, indent+extraIndent);
        }
    }
}
//recursively check is a given particle is ancestor
bool jpsi4LepLepKmcFitter::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool jpsi4LepLepKmcFitter::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}
// ------------ method called to produce the data  ------------
void
jpsi4LepLepKmcFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   std::unique_ptr<pat::CompositeCandidateCollection> ZCandColl(new pat::CompositeCandidateCollection); //produce ZCandidates

   edm::ESHandle<TransientTrackBuilder> theTTB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTB); 

   edm::Handle<pat::CompositeCandidateCollection> dileptons;
   iEvent.getByToken(dilepton_Label,dileptons); //dilepton
    
   edm::Handle<pat::CompositeCandidateCollection> dimuons;
   iEvent.getByToken(dimuon_Label,dimuons); //dimuon

   //edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;  //miniAOD
   //iEvent.getByToken(trackCollection_label,thePATTrackHandle);  //Tracks

   //edm::Handle< View<pat::Muon> > muons;
   //iEvent.getByToken(muonToken_,muons);
    
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

   for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
   vertexRef_i++;
   if ( !(vtx->isFake())                                  // if is not fake
        && vtx->ndof()>=4. && vtx->position().Rho() < 2.0 // and number of degres of freedeom (fit) is > than 4
        && fabs(vtx->position().Z()) < 24.0) {            // and is near the beeam pipe we select it as our primary vertex PV
        PV = vtx;
        break;
        }
   }
    /*
   for (size_t i = 0; i < vertices->size(); ++i) {
   const reco::Vertex &vtx = (*vertices)[i];
   if ( !(vtx.isFake())                                  // if is not fake
        && vtx.ndof()>=4. && vtx.position().Rho() < 2.0 // and number of degres of freedeom (fit) is > than 4
        && fabs(vtx.position().Z()) < 24.0) {            // and is near the beeam pipe we select it as our primary vertex PV
        bestPtVtx = vtx;
        break;
        }
   }
*/
    
    TLorentzVector gen_mu1P, gen_mu1N, gen_mu2P, gen_mu2N, gen_el1P, gen_el1N, gen_el2P, gen_el2N;
    TVector3       gen_z_vtx,gen_jpsi_vtx;

    gen_z_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_mu1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_mu1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_mu2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_mu2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_el1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_el1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_el2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_el2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    gen_z_vtx.SetXYZ(0.,0.,0.);
    gen_jpsi_vtx.SetXYZ(0.,0.,0.);
    int n_Z_dau = 0;
    int decaychannel = 0;
    int Event_Cand = 1;
    //NEW NEW NEW MC ALV CSPM
    //std::cout << "<------------------------------------ NEW EVENT----------------------------------------->" << std::endl;

    if (pruned.isValid()){
        TLorentzVector temp_mu1P, temp_mu1N, temp_mu2P, temp_mu2N, temp_el1P, temp_el1N, temp_el2P, temp_el2N;
        for(size_t i=0; i<pruned->size(); i++){
            const reco::Candidate *mom = &(*pruned)[i];
            if (std::isnan(mom->mass())) continue;
            if(abs(mom->pdgId()) == 23){ // if generated is Z boson
                temp_mu1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                int efromZ = 0;
                int tfromZ = 0;
                int mfromZ = 0;
                for(size_t k=0; k<packed->size(); k++){
                     const reco::Candidate * stable_dau = &(*packed)[k];
                     int stable_id = (*packed)[k].pdgId();
                     if (stable_dau != nullptr && isAncestor(mom,stable_dau)) {
                         if (stable_id == 13 || stable_id == -13) { //electros 11 muons 13 as final states
                             mfromZ++;
                         }
                         else if (stable_id == 11 || stable_id == -11){
                             efromZ++;
                         }
                         else if (stable_id == 16 || stable_id == -16){
                             tfromZ++;
                         }
                     }
                 }
                if (tfromZ) break;
                //if (efromZ + mfromZ < 3) break;
                //std::cout << "<------------------------------------ PRINT DECAY----------------------------------------->" << std::endl;
                analyzeDecay(mom, temp_mu1P, temp_mu1N, temp_mu2P, temp_mu2N,temp_el1P, temp_el1N, temp_el2P, temp_el2N, decaychannel, 0);
                if(decaychannel){
                    //std::cout << "good decay" << std::endl;
                    gen_z_p4.SetPtEtaPhiM(mom->pt(), mom->eta(), mom->phi(), mom->mass());
                    break;
                }//end if good decay chain

            }//end if Z boson
        }//end loop over pruned
        //match with final states
        for(size_t k=0; k<packed->size(); k++){
            const reco::Candidate * stable_dau = &(*packed)[k];
            int stable_id = (*packed)[k].pdgId();
            // decay 1  Z--> 2mu
            // decay 2  Z--> 2el
            // decay 3  Z--> 4mu
            // decay 4  Z--> 4el
            // decay 5  Z--> 2mu->2el
            // decay 6  Z--> 2el->2mu
            if (decaychannel == 1){
                if (stable_id == 13 && deltaR(temp_mu1N.Eta(), temp_mu1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu1N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu N in Z->2mu"<<std::endl;
                }
                if (stable_id == -13 && deltaR(temp_mu1P.Eta(), temp_mu1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu1P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu P in Z->2mu"<<std::endl;
                }
            }
            if (decaychannel == 2){
                if (stable_id == 11 && deltaR(temp_el1N.Eta(), temp_el1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el1N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el N in Z->2el"<<std::endl;
                }
                if (stable_id == -11 && deltaR(temp_el1P.Eta(), temp_el1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el1P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el P in Z->2el"<<std::endl;
                }
            }
            if (decaychannel == 3){
                if (stable_id == 13 && deltaR(temp_mu1N.Eta(), temp_mu1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu1N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 1 N in Z->4mu"<<std::endl;
                }
                if (stable_id == -13 && deltaR(temp_mu1P.Eta(), temp_mu1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu1P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 1 P in Z->4mu"<<std::endl;
                }
                if (stable_id == 13 && deltaR(temp_mu2N.Eta(), temp_mu2N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu2N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 2 N in Z->4mu"<<std::endl;
                }
                if (stable_id == -13 && deltaR(temp_mu2P.Eta(), temp_mu2P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu2P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 2 P in Z->4mu"<<std::endl;
                }
            }
            if (decaychannel == 4){
                if (stable_id == 11 && deltaR(temp_el1N.Eta(), temp_el1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el1N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 1 N in Z->4el"<<std::endl;
                }
                if (stable_id == -11 && deltaR(temp_el1P.Eta(), temp_el1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el1P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 1 P in Z->4el"<<std::endl;
                }
                if (stable_id == 11 && deltaR(temp_el2N.Eta(), temp_el2N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el2N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 2 N in Z->4el"<<std::endl;
                }
                if (stable_id == -11 && deltaR(temp_el2P.Eta(), temp_el2P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el2P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 2 P in Z->4el"<<std::endl;
                }
            }
            if (decaychannel == 5){
                if (stable_id == 13 && deltaR(temp_mu1N.Eta(), temp_mu1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu1N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 1 N in Z->2mu 2el"<<std::endl;
                }
                if (stable_id == -13 && deltaR(temp_mu1P.Eta(), temp_mu1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu1P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 1 P in Z->2mu 2el"<<std::endl;
                }
                if (stable_id == 11 && deltaR(temp_el1N.Eta(), temp_el1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el2N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 2 N in Z->2mu 2el"<<std::endl;
                }
                if (stable_id == -11 && deltaR(temp_el1P.Eta(), temp_el1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el2P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 2 P in Z->2mu 2el"<<std::endl;
                }
            }
            if (decaychannel == 6){
                if (stable_id == 13 && deltaR(temp_mu1N.Eta(), temp_mu1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu2N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 2 N in Z->2el 2mu"<<std::endl;
                }
                if (stable_id == -13 && deltaR(temp_mu1P.Eta(), temp_mu1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_mu2P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched mu 2 P in Z->2el 2mu"<<std::endl;
                }
                if (stable_id == 11 && deltaR(temp_el1N.Eta(), temp_el1N.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el1N.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 1 N in Z->2el 2mu"<<std::endl;
                }
                if (stable_id == -11 && deltaR(temp_el1P.Eta(), temp_el1P.Phi(), stable_dau->eta(), stable_dau->phi()) < 0.05){
                    gen_el1P.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                    std::cout << "matched el 1 P in Z->2el 2mu"<<std::endl;
                }
            }
        }// end packed
    }//end if pruned
    //NEW MC ALV
    int tst = 0;
    int z_gen = 0;
    /*
    if (pruned.isValid()){ //if mc exist
    //std::cout<< "Valid pruned container" << std::endl;
        for(size_t i=0; i<pruned->size(); i++){ // loop over generated events        
            const reco::Candidate *mom = &(*pruned)[i];
            if (std::isnan(mom->mass())) continue;
            if(abs(mom->pdgId()) == 23){ // if generated is Z boson
                TLorentzVector temp_lep_1, temp_lep_2, temp_mu_1, temp_mu_2; //define tempotals
                temp_lep_1.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_lep_2.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu_1.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu_2.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                int mfromZ = 0;
                int efromZ = 0;
                int tfromZ = 0;
                int zfromZ = 0;
                z_gen++;
                if (mom->numberOfDaughters() == 1){
                    if (mom->daughter(0)->pdgId() == 23) zfromZ++; 
                }
                printMCtree(mom, 0);
                
                if (zfromZ) continue;
                for(size_t k=0; k<packed->size(); k++){
                    const reco::Candidate * stable_dau = &(*packed)[k];
                    int stable_id = (*packed)[k].pdgId();
                    if (stable_dau != nullptr && isAncestor(mom,stable_dau)) {
                        if (stable_id == 13 || stable_id == -13) { //electros 11 muons 13 as final states
                            mfromZ++;
                        }
                        else if (stable_id == 11 || stable_id == -11){
                            efromZ++;
                        }
                        else if (stable_id == 16 || stable_id == -16){
                            tfromZ++;
                        }
                    }
                }
                //if ((efromZ == 2 && mfromZ == 2) || mfromZ == 4) {
                if (mfromZ > 2 && tfromZ == 0){    
                    if (tst > 2) break;
                    std::cout<< "Found Z with mass: "<< mom->mass() <<", print Tree ... " << tst << std::endl;                    
                    printMCtree(mom, 0);
                    for(size_t k=0; k<packed->size(); k++){
                       const reco::Candidate * stable_dau = &(*packed)[k];
                       int stable_id = (*packed)[k].pdgId();
                       if (stable_dau != nullptr && isAncestor(mom,stable_dau)) { 
                           if (stable_id != 13 && stable_id != -13) continue;
                           if(stable_id == 13 && temp_lep_1.M() == 0){ // if muon- && not previiulsy assigned
                                std::cout << " matched 1 "  << " with Pt: ";
                                std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_lep_1.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                            }
                            else if(stable_id == -13 && temp_lep_2.M() == 0){ // if muon+ && not previusly assigned
                                std::cout << " matched 2 "  << " with Pt: ";
                                std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_lep_2.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                            }
                            else if(stable_id == 13 && temp_mu_1.M() == 0){ // if muon- && not previusly assigned
                                std::cout << " matched 3 "  << " with Pt: ";
                                std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_mu_1.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                            }
                            else if(stable_id == -13 && temp_mu_2.M() == 0){ // if muon+ && not previusly assigned
                                std::cout << " matched 4 "  << " with Pt: ";
                                std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_mu_2.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                            }
                       }
                    }
                    tst++;
                }
                if (temp_lep_1.M() != 0 && temp_lep_2.M() !=0 && temp_mu_1.M() != 0 && temp_mu_2.M() !=0){ //if 4 leptons has been found
                    gen_z_p4.SetPtEtaPhiM(mom->pt(),mom->eta(),mom->phi(),mom->mass());
                    if(!std::isnan(gen_z_p4.M())){
                       TLorentzVector dil_1 = temp_lep_1 + temp_lep_2;
                       TLorentzVector dim_1 = temp_mu_1  + temp_mu_2;

                       TLorentzVector dil_2 = temp_lep_1 + temp_mu_2;
                       TLorentzVector dim_2 = temp_lep_2 + temp_mu_1;

                       float DM1, DM2;
                       DM1 = dil_1.M() - dim_1.M();
                       DM2 = dil_2.M() - dim_2.M();
                       if(std::abs(DM1) > std::abs(DM2)){
                          if (DM1 > 0){ 
                              gen_lepton1_p4 = temp_lep_1;
                              gen_lepton2_p4 = temp_lep_2;
                              gen_muon1_p4 = temp_mu_1;
                              gen_muon2_p4 = temp_mu_2;  
                          }
                          else {
                              gen_lepton1_p4 = temp_mu_1;
                              gen_lepton2_p4 = temp_mu_2;
                              gen_muon1_p4 = temp_lep_1;
                              gen_muon2_p4 = temp_lep_2;  
                          }
                       }
                       else {
                          if (DM2 > 0){
                              gen_lepton1_p4 = temp_lep_1;
                              gen_lepton2_p4 = temp_mu_2;
                              gen_muon1_p4 = temp_mu_1;
                              gen_muon2_p4 = temp_lep_2;
                          } 
                          else {
                              gen_lepton1_p4 = temp_mu_1;
                              gen_lepton2_p4 = temp_lep_2;
                              gen_muon1_p4 = temp_lep_1;
                              gen_muon2_p4 = temp_mu_2;
                          }
                       } 
                       TLorentzVector zz = gen_lepton1_p4 + gen_lepton2_p4 + gen_muon1_p4 + gen_muon2_p4;
                       TLorentzVector gen_dilep = gen_lepton1_p4 + gen_lepton2_p4;
                       TLorentzVector gen_dimun = gen_muon1_p4 + gen_muon2_p4;
                       std::cout << "Found Z to 4l (2 mu + 2 mu), Z cand mass ~ " << gen_z_p4.M() << std::endl;
                       std::cout << "4 lep gen mass ~ " << zz.M() << std::endl;
                       std::cout << "lep1 Pt: " << gen_lepton1_p4.Pt() << std::endl;
                       std::cout << "lep2 Pt: " << gen_lepton2_p4.Pt() << std::endl;
                       std::cout << "mu1  Pt: " << gen_muon1_p4.Pt()  << std::endl;
                       std::cout << "mu2  Pt: " << gen_muon2_p4.Pt()  << std::endl;
                       std::cout << "dilep & dimuon: " << gen_dilep.M() << " & " << gen_dimun.M() << std::endl; 

                    }
                }// end if generated 4 muons
             
            }//end if generated is Z
        }//end loop over pruned 
    }//end if pruned
    */
    //if (tst) std::cout << "x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x NEW EVENT -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x" << std::endl;
    /*
    if (pruned.isValid() ){ // if mc collection exists
        for(size_t i=0; i<pruned->size(); i++){ // loop over generated events
            const reco::Candidate *mom = &(*pruned)[i];
            if (std::isnan(mom->mass())) continue;
            if(abs(mom->pdgId()) == 23){ // if generated is Z boson
                TLorentzVector temp_lep_1, temp_lep_2, temp_mu_1, temp_mu_2; //define tempotals
                temp_lep_1.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_lep_2.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu_1.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu_2.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                int fromZ = 0;
                for(size_t k=0; k<packed->size(); k++){
                    const reco::Candidate * stable_dau = &(*packed)[k];
                    int stable_id = (*packed)[k].pdgId();
                    if (stable_dau != nullptr && isAncestor(mom,stable_dau)) {
                        if (stable_id != 11 && stable_id != -11) continue;
                        fromZ++;
                    }
                }
                //std::cout << "how many daughters from Z ~  " << fromZ << std::endl;
                for(size_t k=0; k<packed->size(); k++){ //loop over stable particle collection
                    const reco::Candidate * stable_dau = &(*packed)[k];
                    int stable_id = (*packed)[k].pdgId();
                    if (stable_id != 11 && stable_id != -11) continue;
                    if (stable_dau != nullptr && isAncestor(mom,stable_dau) && fromZ > 3) { // if stable comes from Z
                        //std::cout<< "+-> Mu " << stable_dau->charge() <<" ||id & n muons " << stable_id <<" "<< fromZ <<  std::endl;
                        if(stable_id == 11 && temp_lep_1.M() == 0){ // if muon- && not previiulsy assigned
                            //std::cout << " matched 1 " << std::endl;
                            temp_lep_1.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                        }
                        else if(stable_id == -11 && temp_lep_2.M() == 0){ // if muon+ && not previusly assigned
                            //std::cout << " matched 2 " << std::endl;
                            temp_lep_2.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                        }
                        else if(stable_id == 11 && temp_mu_1.M() == 0){ // if muon- && not previusly assigned
                            //std::cout << " matched 3 " << std::endl;
                            temp_mu_1.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                        }
                        else if(stable_id == -11 && temp_mu_2.M() == 0){ // if muon+ && not previusly assigned
                            //std::cout << " matched 4 " << std::endl;
                            temp_mu_2.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                        }
                    }// end if stable particle cames from Z
                }//end loop of stable particles
                if (temp_lep_1.M() != 0 && temp_lep_2.M() !=0 && temp_mu_1.M() != 0 && temp_mu_2.M() !=0){ //if 4 leptons has been found
                    gen_z_p4.SetPtEtaPhiM(mom->pt(),mom->eta(),mom->phi(),mom->mass());
                    if(!std::isnan(gen_z_p4.M())){
                       gen_lepton1_p4 = temp_lep_1.Pt() > temp_mu_1.Pt() ? temp_lep_1 : temp_mu_1;
                       gen_lepton2_p4 = temp_lep_2.Pt() > temp_mu_2.Pt() ? temp_lep_2 : temp_mu_2;
                       gen_muon1_p4 = temp_mu_1.Pt() < temp_lep_1.Pt() ? temp_mu_1 : temp_lep_1;
                       gen_muon2_p4 = temp_mu_2.Pt() < temp_lep_2.Pt() ? temp_mu_2 : temp_lep_2;
                       gen_jpsi_p4 = temp_mu_1 + temp_mu_2;
                       gen_z_vtx.SetXYZ(mom->vx(),mom->vy(),mom->vz());
                       TLorentzVector zz = temp_lep_1 + temp_lep_2 + temp_mu_1 + temp_mu_2;
                       //std::cout << "Found Z to 4l (2 mu + 2 mu), Z cand mass ~ " << gen_z_p4.M() << std::endl;
                       //std::cout << "4 lep gen mass ~ " << zz.M() << std::endl;
                       //std::cout << "lep1 Pt: " << temp_lep_1.Pt() << std::endl;
                       //std::cout << "lep2 Pt: " << temp_lep_2.Pt() << std::endl;
                       //std::cout << "mu1  Pt: " << temp_mu_1.Pt()  << std::endl;
                       //std::cout << "mu2  Pt: " << temp_mu_2.Pt()  << std::endl;

                    }
                }
            }// end if Z
        }// end loop of generated events
    }//end pruned
    */
    //NEW for muons and psi pairs
    int nonia = dimuons->size();
    int nmuons = dileptons->size()*2;
    int nPV    = vertices->size();
    
   //double chiVtxSqdProb = ChiSquaredProbability((double)(PV.chi2()),(double)(PV.ndof())); 
   //int breaker = 0;
   // We Cycle over dileptons for the Kfit
   //if (gen_z_p4.M() != 0){   //only for mc
   if(1){
    //if ((nonia && nmuons) || gen_z_p4.M() !=0 ){
    //std::cout<< "+X+X+X+X+X+X+X+X+X+X+X+X+X+X START READING DATA X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X" << std::endl;
    //std::cout<< "dimuon length: " << nonia << std::endl;
    //std::cout<< "dilepton length:  " << nmuons << std::endl;
    for (pat::CompositeCandidateCollection::const_iterator dilepton = dileptons->begin(); dilepton != dileptons->end() /*test && breaker < 10*/; ++dilepton){
    //std::cout << " enters dilepton " << std::endl;    
    for (pat::CompositeCandidateCollection::const_iterator dimuon = dimuons->begin(); dimuon != dimuons->end() /*test && breaker < 10*/; ++dimuon){
       //std::cout << " enters dimuon " << std::endl;     
       const pat::Muon* lept1 = dynamic_cast<const pat::Muon*>(dilepton->daughter("lepton1"));
       const pat::Muon* lept2 = dynamic_cast<const pat::Muon*>(dilepton->daughter("lepton2"));
       const pat::Muon* muon1;  
       const pat::Muon* muon2; 
       
       if (dimuon->daughter("muon1")->charge() == -1 && dimuon->daughter("muon2")->charge() == 1){
            muon1 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon1"));
            muon2 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon2"));
       }
       else if (dimuon->daughter("muon1")->charge() == 1 && dimuon->daughter("muon2")->charge() == -1) {
            muon2 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon1"));
            muon1 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon2"));
       }
       else continue;
       //check muon charge 1 --> -  |||| 2 --> +

       reco::TrackRef glbTrack_l1 = lept1->track();
       reco::TrackRef glbTrack_l2 = lept2->track();
       reco::TrackRef glbTrack_m1 = muon1->track();
       reco::TrackRef glbTrack_m2 = muon2->track();
       if (glbTrack_l1.isNull() || glbTrack_l2.isNull() || glbTrack_m1.isNull() || glbTrack_m2.isNull()) continue;
       
       if (!(glbTrack_l1->quality(reco::TrackBase::highPurity))) continue;
       if (!(glbTrack_l2->quality(reco::TrackBase::highPurity))) continue;
       if (!(glbTrack_m1->quality(reco::TrackBase::highPurity))) continue;
       if (!(glbTrack_m2->quality(reco::TrackBase::highPurity))) continue;
       //std::cout << " pass track purity  " << std::endl;
       if (muon1->pt() == lept1->pt() && muon1->eta() == lept1->eta() && muon1->phi() == lept1->phi()) continue;
       if (muon2->pt() == lept2->pt() && muon1->eta() == lept1->eta() && muon1->phi() == lept1->phi()) continue;


         
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
       
  	   float ldxy1 = lept1->muonBestTrack()->dxy(PV->position());
	   float ldz1  = lept1->muonBestTrack()->dz(PV->position());
   
	   float ldxy2 = lept2->muonBestTrack()->dxy(PV->position());
	   float ldz2  = lept2->muonBestTrack()->dz(PV->position());

       // Check if pair of leptons (muons) is coming from our Primary Vertex
       /////////////////////////////////////////////////////
       //    T r a n s i e n t   T r a c k s  b u i l d e r //
       //////////////////////////////////////////////////////
       reco::TransientTrack lep1TT((*theTTB).build(glbTrack_l1));
       reco::TransientTrack lep2TT((*theTTB).build(glbTrack_l2));
       reco::TransientTrack muon1TT((*theTTB).build(glbTrack_m1));
       reco::TransientTrack muon2TT((*theTTB).build(glbTrack_m2));

       std::pair<bool,Measurement1D> tkPVdistel1 = IPTools::absoluteImpactParameter3D(lep1TT,*PV);
       std::pair<bool,Measurement1D> tkPVdistel2 = IPTools::absoluteImpactParameter3D(lep2TT,*PV);
       std::pair<bool,Measurement1D> tkPVdistem1 = IPTools::absoluteImpactParameter3D(muon1TT,*PV);
       std::pair<bool,Measurement1D> tkPVdistem2 = IPTools::absoluteImpactParameter3D(muon2TT,*PV);
         
       if(!tkPVdistel1.first || !tkPVdistel2.first || !tkPVdistem1.first || !tkPVdistem2.first) continue;
       //std::cout << " pass track valid  " << std::endl;

       //std::cout << "muon ok 520" << std::endl;

       //Isolation
       double dR1 = -1, dR2 = -1, dR3 = -1, dR4 = -1, dR5 = -1, dR6 = -1;
       dR1 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
       dR2 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
       dR3 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
       dR4 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
       dR5 = deltaR(*(lept1->innerTrack()), *(lept2->innerTrack()));
       dR6 = deltaR(*(muon1->innerTrack()), *(muon2->innerTrack()));
       if ( dR1<0.01 || dR2<0.01 || dR3<0.01 ||dR4<0.01 ) continue;
       //std::cout << " pass dR " << std::endl;

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
           //std::cout << "1L is Medium " << std::endl;
       }
       if ( muon1->isTightMuon(*PV)){
           ZMu1Qid += 1000;
           //std::cout << "1L is Tight " << std::endl;
       }
       if ( muon1->isSoftMuon(*PV)){
           ZMu1Qid += 10000;
           //std::cout << "1L is Soft " << std::endl;
       }
       if ( muon1->isHighPtMuon(*PV)){
           ZMu1Qid += 100000;
           //std::cout << "1L is HighPt " << std::endl;
       }
       if ( muon1->isPFMuon()){
           ZMu1Qid += 1000000;
           //std::cout << "1L is ParticleFlow " << std::endl;
       }
       if ( muon1->isTrackerMuon()){
           ZMu1Qid += 10000000;
           //std::cout << "1L is HighPt " << std::endl;
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

       const ParticleMass muMass    = 0.10565837;
       float muSigma                = muMass*1.e-6;
                 
       //Z Kinematic Fit
       KinematicParticleFactoryFromTransientTrack pFactory;
       float chi = 0;
       float ndf = 0;
       
       std::vector<RefCountedKinematicParticle> ZDaughters;
       try {
          ZDaughters.push_back(pFactory.particle(muon1TT,   muMass, chi, ndf, muSigma));
          ZDaughters.push_back(pFactory.particle(muon2TT,   muMass, chi, ndf, muSigma));
          ZDaughters.push_back(pFactory.particle(lep1TT,  muMass, chi, ndf, muSigma));
          ZDaughters.push_back(pFactory.particle(lep2TT,  muMass, chi, ndf, muSigma));
       }
       catch (...){
           continue;
       }
       KinematicParticleVertexFitter ZVertexFitter;
       RefCountedKinematicTree ZTree;
       try{
          ZTree = ZVertexFitter.fit(ZDaughters);
       }
       catch(...){
           continue;
       }
         
       if (ZTree->isEmpty()) continue;
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
        //std::cout << "Fiz Zpass" <<std::endl;
        //if (ZDecayVertex->chiSquared() < 0) continue;
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
            ZVtxP_fit = ChiSquaredProbability((double)(ZDecayVertex->chiSquared()),
                               (double)(ZDecayVertex->degreesOfFreedom()));
        }
       float Z_px = muon1->px()+muon2->px()+lept1->px()+lept2->px();
       float Z_py = muon1->py()+muon2->py()+lept1->py()+lept2->py();
       float Z_pz = muon1->pz()+muon2->pz()+lept1->pz()+lept2->pz();

       TLorentzVector m1, m2, l1, l2;
       m1.SetPtEtaPhiM(muon1->pt(), muon1->eta(),muon1->phi(),muon1->mass());
       m2.SetPtEtaPhiM(muon2->pt(), muon2->eta(),muon2->phi(),muon2->mass());
       l1.SetPtEtaPhiM(lept1->pt(), lept1->eta(),lept1->phi(),lept1->mass());
       l2.SetPtEtaPhiM(lept2->pt(), lept2->eta(),lept2->phi(),lept2->mass());
       TLorentzVector theZ = m1+m2+l1+l2;
                
       float mdxy1 = -1, mdz1 = -1, mdxy2 = -1, mdz2 = -1;

       mdxy1 = muon1->muonBestTrack()->dxy(PV->position());
       mdz1 =  muon1->muonBestTrack()->dz(PV->position());
       mdxy2 = muon2->muonBestTrack()->dxy(PV->position());
       mdz2 =  muon2->muonBestTrack()->dz(PV->position());

       if (ZM_fit < 60.0) continue;
       if (ZM_fit > 150.0) continue;
       //std::cout << " pass mass window " << std::endl;


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
       patZ.addUserInt("tst_", tst);
       patZ.addUserInt("z_gen_", z_gen);
                 
       patZ.addUserFloat("vProb",ZVtxP_fit);
       patZ.addUserFloat("vChi2",ZDecayVertex->chiSquared());
       patZ.addUserFloat("ZvtxX",ZVtxX_fit);
       patZ.addUserFloat("ZvtxY",ZVtxY_fit);
       patZ.addUserFloat("ZvtxZ",ZVtxZ_fit);
       patZ.addUserFloat("dRm1m2",dR6);
       patZ.addUserFloat("dRl1l2",dR5);
       patZ.addUserFloat("dRl1m1",dR1);
       patZ.addUserFloat("dRl1m2",dR2);
       patZ.addUserFloat("dRl2m1",dR3);
       patZ.addUserFloat("dRl2m2",dR4);

       bool child = ZTree->movePointerToTheFirstChild();
       /////////////////////////////////////////////////
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
       patM1.addUserFloat("dIP3DSig", tkPVdistem1.second.significance());
        
       patM1.addUserFloat("dIP3D",tkPVdistem1.second.value());
       patM1.addUserFloat("dIP3DErr",tkPVdistem1.second.error());

       int psiM1_TrackerLWM = muon1->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
       int psiM1_PixelLWM   = muon1->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
       int psiM1_ValPixHit  = muon1->muonBestTrack()->hitPattern().numberOfValidPixelHits();
       int psiM2_TrackerLWM = muon2->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
       int psiM2_PixelLWM   = muon2->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
       int psiM2_ValPixHit  = muon2->muonBestTrack()->hitPattern().numberOfValidPixelHits();

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
       patM2.addUserFloat("dIP3DSig", tkPVdistem2.second.significance());
                 
        patM2.addUserFloat("dIP3d", tkPVdistem2.second.significance());
		patM2.addUserFloat("dIP3D", tkPVdistem2.second.value());
		patM2.addUserFloat("dIP3DErr", tkPVdistem2.second.error());
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
        psi.addUserFloat("vProb", 0);
        psi.addUserFloat("vChi2", 0);
		//save measured values on a different daughter
		pat::CompositeCandidate msrd_psi;
		msrd_psi.addDaughter(pat_msrdM1,"msrd_muon1");
		msrd_psi.addDaughter(pat_msrdM2,"msrd_muon2");
		msrd_psi.setP4(pat_msrdM1.p4()+pat_msrdM2.p4());
                ////////////////////////////////////////////////
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
        patL2.addUserFloat("dIP3DSig"   ,tkPVdistel2.second.significance());
		patL2.addUserFloat("dIP3D"      ,tkPVdistel2.second.value());
		patL2.addUserFloat("dIP3DErr"   ,tkPVdistel2.second.error());
		// My muon Q id
		patL2.addUserInt("ZLe2Qid_", ZLe2Qid);

		patL2.addUserFloat("ZLe2_TrackerLWM_", ZLe2_TrackerLWM);
		patL2.addUserFloat("ZLe2_PixelLWM_", ZLe2_PixelLWM);
		patL2.addUserFloat("ZLe2_ValPixHit_", ZLe2_ValPixHit);

                //patL2.addUserFloat("pass", pass2);
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
        std::cout << "Zcand OK " << gen_z_p4.M() << " vs " << ZM_fit <<std::endl;
        std::cout << "Lep 1 Pt "  << gen_lepton1_p4.Pt()  << " vs " << patL1.pt() << std::endl;
        std::cout << "Lep 2 Pt "  << gen_lepton2_p4.Pt()  << " vs " << patL2.pt() << std::endl;
        std::cout << "Muon 1 Pt " << gen_muon1_p4.Pt()    << " vs " << patM1.pt() << std::endl;
        std::cout << "Muon 2 Pt " << gen_muon2_p4.Pt()    << " vs " << patM2.pt() << std::endl;

/*test
        ZCandColl->push_back(patZ);
        breaker +=1;
        std::cout<< "Pass on: " << breaker << std::endl;
        if(breaker > 10){std::cout<<"BREAK"<<std::endl; break;}
        continue;
*/ //endtest


             //}end of ifValidFit

    }//end loop for dimuon
    }//end loop for dilepton

   iEvent.put(std::move(ZCandColl),"ZCandidates");
   std::cout<< "+X+X+X+X+X+X+X+X+X+X+X+X+X+X END READING DATA X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X" << std::endl;
   }//end if gen only MC
}//end produce 

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
jpsi4LepLepKmcFitter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
jpsi4LepLepKmcFitter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
jpsi4LepLepKmcFitter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
jpsi4LepLepKmcFitter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
jpsi4LepLepKmcFitter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
jpsi4LepLepKmcFitter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
jpsi4LepLepKmcFitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
bool jpsi4LepLepKmcFitter::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool jpsi4LepLepKmcFitter::IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk->eta());
  double DeltaP   = fabs(mu.p()-tk->p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
double jpsi4LepLepKmcFitter::getIso(const pat::Muon& mu){
                Float_t coriso = 99.0;
                reco::MuonPFIsolation pfR03 = mu.pfIsolationR03();
                coriso = pfR03.sumChargedHadronPt + std::max(0., pfR03.sumNeutralHadronEt+pfR03.sumPhotonEt-0.5*pfR03.sumPUPt);
return coriso;
}

//define this as a plug-in
DEFINE_FWK_MODULE(jpsi4LepLepKmcFitter);
