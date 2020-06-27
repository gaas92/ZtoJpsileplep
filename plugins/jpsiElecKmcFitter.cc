// -*- C++ -*-
//
// Package:    AnalyzeZll/jpsiElecKmcFitter
// Class:      jpsiElecKmcFitter
// 
/**\class jpsiElecKmcFitter jpsiElecKmcFitter.cc AnalyzeZll/jpsiElecKmcFitter/plugins/jpsiElecKmcFitter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rogelio Reyes Almanza
//         Created:  Tue, 04 Sep 2018 09:41:00 GMT
//
//


// system include files
#include <memory>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

//ok
#include "TLorentzVector.h"

///Data formats
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h> 
#include "DataFormats/Candidate/interface/Candidate.h"
#include <DataFormats/PatCandidates/interface/PackedCandidate.h>
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include <DataFormats/MuonReco/interface/MuonFwd.h>
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
//ok
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"





//
// class declaration
//

class jpsiElecKmcFitter : public edm::stream::EDProducer<>{
   public:
      explicit jpsiElecKmcFitter(const edm::ParameterSet&);
      ~jpsiElecKmcFitter() override{};

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
   // virtual void beginStream(edm::StreamID) override;
      Float_t getIso(const pat::Muon& );
      Float_t getIsoVar(const pat::Electron&);       
      Float_t getEtaInSeed(const pat::Electron&);
    
      float ElectronRelIso(const reco::Candidate *cand, float rho);
      float MuonRelIso(const reco::Candidate *cand, float rho);
    
      bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
      bool IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu);
      bool    isAncestor(const reco::Candidate*, const reco::Candidate*);
      bool    isAncestor(int, const reco::Candidate*);
    
      int convertBinaryToDecimal(unsigned long long);

      void produce(edm::Event&, const edm::EventSetup&) override;
   //   virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
        edm::EDGetTokenT<pat::CompositeCandidateCollection> dimuon_Label;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> dielec_Label;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
    
        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;
    
        edm::EDGetTokenT<double> fixedGridRhoFastjetAll_;

};


jpsiElecKmcFitter::jpsiElecKmcFitter(const edm::ParameterSet& iConfig){
	dimuon_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuon"));
	dielec_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dilepton"));
	primaryVertices_Label = consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"));
    
    genCands_ = consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"));
    packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
    
    fixedGridRhoFastjetAll_ = consumes<double> (iConfig.getParameter <edm::InputTag>("fixedGridRhoFastjetAll"));
    
   	produces<pat::CompositeCandidateCollection>("ZCandidates"); 
 
}


//jpsiElecKmcFitter::~jpsiElecKmcFitter(){}


//
// member functions
//
//recursively check is a given particle is ancestor
bool jpsiElecKmcFitter::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool jpsiElecKmcFitter::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}


// ------------ method called to produce the data  ------------
void jpsiElecKmcFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::auto_ptr<pat::ElectronCollection > selectedCollection(new pat::ElectronCollection ); //new
  std::unique_ptr<pat::CompositeCandidateCollection> ZCandColl(new pat::CompositeCandidateCollection); 
     
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);
  
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  iEvent.getByToken(dielec_Label,dileptons);

  edm::Handle<edm::View<pat::Electron> > electrons; //new
  //if( !electrons.isValid() ) iEvent.getByToken(electronsMiniAODToken_,electrons); //new
  
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  edm::ESHandle<TransientTrackBuilder> theTTBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
 
    //MC
    edm::Handle<reco::GenParticleCollection> pruned;
    iEvent.getByToken(genCands_, pruned);
    // Packed particles are all the status 1, so usable to remake jets
    // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
    edm::Handle<pat::PackedGenParticleCollection> packed;
    iEvent.getByToken(packedGenToken_,packed);
    
  //////////////////////////////////
  //// Select  the best PV      ////
  //////////////////////////////////
  
  if (!dileptons.isValid() )return; 	

  if(primaryVertices_handle->empty())return;  
  

  reco::VertexCollection::const_iterator PV = primaryVertices_handle->end();
  for (reco::VertexCollection::const_iterator vtx = primaryVertices_handle->begin(); vtx != primaryVertices_handle->end(); ++vtx, ++PV) {
    if ( !(vtx->isFake())
         && vtx->ndof()>=4. && vtx->position().Rho() < 2.0
         && fabs(vtx->position().Z()) < 24.0) {
      PV = vtx;
      break;
    }
  }  
  if ( PV==primaryVertices_handle->end() ) return;

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
           //std::cout << " Z mass : " << dau->mass() << std::endl;
           gen_z_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
           n_Z_dau = dau->numberOfDaughters();
           if (n_Z_dau!=3) continue;
           //std::cout << " Z daugh: " << dau->numberOfDaughters() << std::endl;
           for (size_t k=0; k<dau->numberOfDaughters(); k++) {
             const reco::Candidate *gdau = dau->daughter(k);
             //std::cout << "MC Z daughter pdgID: " << gdau->pdgId() << "mass: " << gdau->mass() << std::endl;
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
                     //std::cout<<" psi = "<< gdau->pdgId()<< " daughter ID " << stable_id << std::endl;
                     if(stable_id == 13) { //found muon-
                         gen_muon1_p4.SetPtEtaPhiM(dauInPrunedColl->pt(),dauInPrunedColl->eta(),dauInPrunedColl->phi(),dauInPrunedColl->mass());
                         nm++;
                         //std::cout<< "works m- "<< dauInPrunedColl->mass() <<std::endl;
                     }
                      if(stable_id == -13){ //found muon+
                        gen_muon2_p4.SetPtEtaPhiM(dauInPrunedColl->pt(),dauInPrunedColl->eta(),dauInPrunedColl->phi(),dauInPrunedColl->mass());
                        nm++;
                        //std::cout<< "works m+ "<< dauInPrunedColl->mass() << std::endl;
                     }
                  }
               }
                 
             }//end found psi
             if (gdau->pdgId()==11 ) {// pdgid for electron=11
                foundit++;
            gen_lepton1_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                 
            }
             if (gdau->pdgId()==-11 ) {// pdgid for muon+=13
                foundit++;
                gen_lepton2_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
            }
           }// end number of daughters
        } //endif found Z
     }//end for

  } //end pruned
    
  //NEW for muons and psi pairs
  int nonia = dimuons->size();
  int nmuons = dileptons->size()*2;
  int nPV    = primaryVertices_handle->size();
    
  for (pat::CompositeCandidateCollection::const_iterator dimuon = dimuons->begin(); dimuon != dimuons->end(); ++dimuon ) {
        //Jpsi Muons
	const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon1"));
	const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon2"));

        int muon1Mu8DiEle12 = 0;
        int muon2Mu8DiEle12 = 0;
        /*
        try {
           const pat::TriggerObjectStandAloneCollection muon1HLT_Mu8DiEle12 = muon1->triggerObjectMatchesByFilter("HLTMu8DiEle12CaloIdLTrackIdLElectronlegSequence");
           const pat::TriggerObjectStandAloneCollection muon2HLT_Mu8DiEle12 = muon2->triggerObjectMatchesByFilter("HLTMu8DiEle12CaloIdLTrackIdLElectronlegSequence");
           
           //std::cout << "IT WORKS ... for now" << std::endl;
           //std::cout << "Sise of muon1: " << muon1HLT_Mu8DiEle12.size()  << std::endl;
           //std::cout << "Sise of muon2: " << muon2HLT_Mu8DiEle12.size()  << std::endl;
                     
           if (muon1HLT_Mu8DiEle12.size() > 0 ){
               //std::cout << "muon1 IsoMu HLT matched" << std::endl;
               muon1Mu8DiEle12 = 1;
               }
           if (muon2HLT_Mu8DiEle12.size() > 0 ){
               //std::cout << "muon2 IsoMu HLT matched" << std::endl;
               muon2Mu8DiEle12 = 1;
              }
        }
        catch ( ... ){
             std::cout << "Esta madre no jala" << std::endl;
             }
		*/
	  float jpsiVprob=0;
	  float jpsiChi2=0;
	  jpsiVprob = dimuon->userFloat("vProb");
	  jpsiChi2  = dimuon->userFloat("vNChi2");
    
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
      int psiM1_TrackerLWM = muon1->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
      int psiM1_PixelLWM   = muon1->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
      int psiM1_ValPixHit  = muon1->muonBestTrack()->hitPattern().numberOfValidPixelHits();
      int psiM2_TrackerLWM = muon2->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
      int psiM2_PixelLWM   = muon2->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
      int psiM2_ValPixHit  = muon2->muonBestTrack()->hitPattern().numberOfValidPixelHits();
 	
      reco::TrackRef JpsiTk[2] = {muon1->innerTrack(), muon2->innerTrack()};
	
      std::vector<reco::TransientTrack> MuMuTTks;
      MuMuTTks.push_back(theTTBuilder->build(&JpsiTk[0]));
      MuMuTTks.push_back(theTTBuilder->build(&JpsiTk[1]));
       
      std::pair<bool,Measurement1D> tkPVdist1 = IPTools::absoluteImpactParameter3D(MuMuTTks.at(0),*PV);
      std::pair<bool,Measurement1D> tkPVdist2 = IPTools::absoluteImpactParameter3D(MuMuTTks.at(1),*PV);
      if (!tkPVdist1.first|| !tkPVdist2.first ) continue;
      //if (fabs(tkPVdist1.second.significance())>4.) continue;
      //if (fabs(tkPVdist2.second.significance())>4.) continue;
      //std::cout << "works for jpsimuons" << std::endl;
      //std::cout << "diElectrons handle Size: " << dileptons->size() << " " << std::endl;
	  for (pat::CompositeCandidateCollection::const_iterator dilepton = dileptons->begin(); dilepton != dileptons->end(); ++dilepton){
                //std::cout << "Enters elec" << std::endl;
                //cambiar pues ahora son electrones	
                int lept1Ele25wpT = 0;
                int lept2Ele25wpT = 0;
                     
                int lept1Ele23_12 = 0;
                int lept2Ele23_12 = 0;
                     
                int lept1Mu8DiEle12 = 0;
                int lept2Mu8DiEle12 = 0;
          
                const pat::Electron* lept1 = dynamic_cast<const pat::Electron*>(dilepton->daughter("lepton1"));
                const pat::Electron* lept2 = dynamic_cast<const pat::Electron*>(dilepton->daughter("lepton2"));
          /*
                try {
                     const pat::TriggerObjectStandAloneCollection lept1HLT_Ele25wpT = lept1->triggerObjectMatchesByFilter("HLTEle25WPTightGsfSequence");
                     const pat::TriggerObjectStandAloneCollection lept2HLT_Ele25wpT = lept2->triggerObjectMatchesByFilter("HLTEle25WPTightGsfSequence");
                     
                     const pat::TriggerObjectStandAloneCollection lept1HLT_Ele23_12 = lept1->triggerObjectMatchesByFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter");
                     const pat::TriggerObjectStandAloneCollection lept2HLT_Ele23_12 = lept2->triggerObjectMatchesByFilter("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter");
                     
                     const pat::TriggerObjectStandAloneCollection lept1HLT_Mu8DiEle12 = lept1->triggerObjectMatchesByFilter("HLTMu8DiEle12CaloIdLTrackIdLElectronlegSequence");
                     const pat::TriggerObjectStandAloneCollection lept2HLT_Mu8DiEle12 = lept2->triggerObjectMatchesByFilter("HLTMu8DiEle12CaloIdLTrackIdLElectronlegSequence");
                     
                     //const pat::TriggerObjectStandAloneCollection muHLTMatches1_t2 = iMuon1->triggerObjectMatchesByFilter("hltJpsiTkVertexFilter");
                     //std::cout << "IT WORKS ... for now" << std::endl;
                     //std::cout << "Sise of lept1: " << lept1HLT_IsoMu24.size()  << std::endl;
                     //std::cout << "Sise of lept2: " << lept2HLT_IsoMu24.size()  << std::endl;
                     //
                    //
                     if (lept1HLT_Ele25wpT.size() > 0 ){
                         std::cout << "lept1 Ele25 HLT matched" << std::endl;
                         lept1Ele25wpT = 1;
                         }
                     if (lept2HLT_Ele25wpT.size() > 0 ){
                         std::cout << "lept2 Ele25 HLT matched" << std::endl;
                         lept2Ele25wpT = 1;
                         }
                     if (lept1HLT_Ele23_12.size() > 0 ){
                         std::cout << "lept1 Ele23_12 HLT matched" << std::endl;
                         lept1Ele23_12 = 1;
                         }
                     if (lept2HLT_Ele23_12.size() > 0 ){
                         std::cout << "lept2 Ele23_12 HLT matched" << std::endl;
                         lept2Ele23_12 = 1;
                         }
                     if (lept1HLT_Mu8DiEle12.size() > 0 ){
                         std::cout << "lept1 Mu8DiEle12 HLT matched" << std::endl;
                         lept1Mu8DiEle12 = 1;
                         }
                     if (lept2HLT_Mu8DiEle12.size() > 0 ) {
                         std::cout << "lept2 Mu8DiEle12 HLT matched" << std::endl;
                         lept2Mu8DiEle12 = 1 ;
                         } 
                     }
                catch ( ... ){
                     std::cout << "Esta madre no jala" << std::endl;
                    }
          */
		//if (lept1 == muon1 || lept1==muon2 || lept2 == muon1 || lept2== muon2  ) continue; //realmente no importa porque ahora son electrones y muones
	        ///////////////////////////////////////
                //////////// MY ELECTRON ID ///////////////
	        ///////////////////////////////////////

                int ZLe1Qid = 0;
                int ZLe2Qid = 0;
                unsigned long long ZLe1Qid_n = 0;
                unsigned long long ZLe2Qid_n = 0;

        
 	        // Lepton 1 (from Z) 2016 promt & legacy
 	        if ( lept1->electronID("cutBasedElectronID-Summer16-80X-V1-veto")==1 ){
                   ZLe1Qid_n += 1;
                   //std::cout << "1L is  80X veto Electron" << std::endl;
                }
	        if ( lept1->electronID("cutBasedElectronID-Summer16-80X-V1-loose")==1 ){
                   ZLe1Qid_n += 10;
                   //std::cout << "1L is 80X loose " << std::endl;
                }
            if ( lept1->electronID("cutBasedElectronID-Summer16-80X-V1-medium")==1 ){
                   ZLe1Qid_n += 100;
                   //std::cout << "1L is 80X Medium " << std::endl;
                }
		    if ( lept1->electronID("cutBasedElectronID-Summer16-80X-V1-tight")==1 ){
                   ZLe1Qid_n += 1000;
                   //std::cout << "1L is 80X Tight " << std::endl;
                }
          
            if ( lept1->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp80")==1 ){
                   ZLe1Qid_n += 10000;
                   //std::cout << "1L is  80X V1 wp80" << std::endl;
                }
            if ( lept1->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp90")==1 ){
                   ZLe1Qid_n += 100000;
                   //std::cout << "1L is 80X V1 wp90 " << std::endl;
                }
            // 2016 Legacy MiniAOD V2
            if ( lept1->electronID("cutBasedElectronID-Fall17-94X-V2-loose")==1 ){
                   ZLe1Qid_n += 1000000;
                   //std::cout << "1L is 94X loose " << std::endl;
                }
            if ( lept1->electronID("cutBasedElectronID-Fall17-94X-V2-medium")==1 ){
                   ZLe1Qid_n += 10000000;
                   //std::cout << "1L is 94X medium " << std::endl;
                }
            if ( lept1->electronID("cutBasedElectronID-Fall17-94X-V2-tight")==1 ){
                   ZLe1Qid_n += 100000000;
                   //std::cout << "1L is 94X tight " << std::endl;
                }
            if ( lept1->electronID("mvaEleID-Fall17-iso-V2-wp80")==1 ){
                   ZLe1Qid_n += 1000000000;
                   //std::cout << "1L is 94X iso V2 wp80 " << std::endl;
                }
            if ( lept1->electronID("mvaEleID-Fall17-iso-V2-wp90")==1 ){
                   ZLe1Qid_n += 1000000000;
                   //std::cout << "1L is 94X iso V2 wp90 " << std::endl;
                }
//////////////////////////////////////////////////////////////////////////////////////////
 	        // Lepton 2 (from Z)
//////////////////////////////////////////////////////////////////////////////////////////
            // Lepton 2 (from Z) 2016 promt & legacy
            if ( lept2->electronID("cutBasedElectronID-Summer16-80X-V1-veto")==1 ){
                    ZLe2Qid_n += 1;
                    //std::cout << "2L is  80X veto Electron" << std::endl;
                }
            if ( lept2->electronID("cutBasedElectronID-Summer16-80X-V1-loose")==1 ){
                    ZLe2Qid_n += 10;
                    //std::cout << "2L is 80X loose " << std::endl;
                }
            if ( lept2->electronID("cutBasedElectronID-Summer16-80X-V1-medium")==1 ){
                    ZLe2Qid_n += 100;
                    //std::cout << "2L is 80X Medium " << std::endl;
                }
            if ( lept2->electronID("cutBasedElectronID-Summer16-80X-V1-tight")==1 ){
                    ZLe2Qid_n += 1000;
                    //std::cout << "2L is 80X Tight " << std::endl;
                }
            
            if ( lept2->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp80")==1 ){
                    ZLe2Qid_n += 10000;
                    //std::cout << "2L is  80X V1 wp80" << std::endl;
                }
            if ( lept2->electronID("mvaEleID-Spring16-GeneralPurpose-V1-wp90")==1 ){
                    ZLe2Qid_n += 100000;
                    //std::cout << "2L is 80X V1 wp90 " << std::endl;
                }
            // 2016 Legacy MiniAOD V2
            if ( lept2->electronID("cutBasedElectronID-Fall17-94X-V2-loose")==1 ){
                    ZLe2Qid_n += 1000000;
                   //std::cout << "2L is 94X loose " << std::endl;
                }
            if ( lept2->electronID("cutBasedElectronID-Fall17-94X-V2-medium")==1 ){
                    ZLe2Qid_n += 10000000;
                    //std::cout << "2L is 94X medium " << std::endl;
                }
            if ( lept2->electronID("cutBasedElectronID-Fall17-94X-V2-tight")==1 ){
                    ZLe2Qid_n += 100000000;
                    //std::cout << "2L is 94X tight " << std::endl;
                }
            if ( lept2->electronID("mvaEleID-Fall17-iso-V2-wp80")==1 ){
                    ZLe2Qid_n += 1000000000;
                    //std::cout << "2L is 94X iso V2 wp80 " << std::endl;
                }
            if ( lept2->electronID("mvaEleID-Fall17-iso-V2-wp90")==1 ){
                     ZLe2Qid_n += 1000000000;
                     //std::cout << "2L is 94X iso V2 wp90 " << std::endl;
                }
            ZLe1Qid = convertBinaryToDecimal(ZLe1Qid_n);
            ZLe2Qid = convertBinaryToDecimal(ZLe2Qid_n);
          
            int ZLe1_TrackerLWM = lept1->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
            int ZLe1_PixelLWM   = lept1->gsfTrack()->hitPattern().pixelLayersWithMeasurement();
            int ZLe1_ValPixHit  = lept1->gsfTrack()->hitPattern().numberOfValidPixelHits();
 		
            int ZLe2_TrackerLWM = lept2->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
            int ZLe2_PixelLWM   = lept2->gsfTrack()->hitPattern().pixelLayersWithMeasurement();
            int ZLe2_ValPixHit  = lept2->gsfTrack()->hitPattern().numberOfValidPixelHits();

            //int ZLe1_ElecMissHits = 0;
            //int ZLe2_ElecMissHits = 0;
            int ZLe1_ElecMissHits = lept1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
            //int ZLe1_ElecMissHits = lept1->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
            int ZLe2_ElecMissHits = lept2->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
            //int ZLe2_ElecMissHits = lept2->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
            //std::cout<< "ZLe1_ElecMissHits: " << ZLe1_ElecMissHits << std::endl;
            //std::cout<< "ZLe2_ElecMissHits: " << ZLe2_ElecMissHits << std::endl;
            /*
            double dR1 = deltaR(*(lept1->innerTrack()), *(muon1->innerTrack()));
		    double dR2 = deltaR(*(lept1->innerTrack()), *(muon2->innerTrack()));
		    double dR3 = deltaR(*(lept2->innerTrack()), *(muon1->innerTrack()));
		    double dR4 = deltaR(*(lept2->innerTrack()), *(muon2->innerTrack()));
		    double dR5 = deltaR(*(lept1->innerTrack()), *(lept2->innerTrack()));
		    double dR6 = deltaR(*(muon1->innerTrack()), *(muon2->innerTrack()));
	        */
            double dRm1m2   = deltaR(*(muon1->innerTrack()), *(muon2->innerTrack()));
            double dRel1el2 = deltaR(*(lept1->gsfTrack()), *(lept2->gsfTrack())    );
            double dRel1mu1 = deltaR(*(lept1->gsfTrack()), *(muon1->innerTrack()));
            double dRel1mu2 = deltaR(*(lept1->gsfTrack()), *(muon2->innerTrack()));
            double dRel2mu1 = deltaR(*(lept2->gsfTrack()), *(muon1->innerTrack()));
            double dRel2mu2 = deltaR(*(lept2->gsfTrack()), *(muon2->innerTrack()));

		    //////////////////////////////////////////////////////
		    //  I m p a c t   P a ra m e t e r s                //
		    //////////////////////////////////////////////////////
		    float mdxy1 = muon1->muonBestTrack()->dxy(PV->position());
           	float mdz1 = muon1->muonBestTrack()->dz(PV->position());	
		
 		    float mdxy2 = muon2->muonBestTrack()->dxy(PV->position());
           	float mdz2 = muon2->muonBestTrack()->dz(PV->position());	
		
		    float ldxy1 = lept1->bestTrack()->dxy(PV->position());
           	float ldz1 = lept1->bestTrack()->dz(PV->position());	
		
 		    float ldxy2 = lept2->bestTrack()->dxy(PV->position());
           	float ldz2 = lept2->bestTrack()->dz(PV->position());	
		

            //std::cout<< "enters cycle for non resonant electrons" << std::endl;
       
  		    //if ( dR1<0.02 || dR2<0.02 || dR3<0.02 ||dR4<0.02 ) continue; // agregar corte en R5 y R6
		    if ( dRel1mu1 <0.01 || dRel1mu2<0.01 || dRel2mu1 <0.01 || dRel2mu2<0.01 ) continue;
            /////////////////////////////////////////////////////
		    //    T r a n s i e n t   T r a c k s  b u i l d e r //
		    //////////////////////////////////////////////////////
			
            //std::cout << "test track" << std::endl;
            //std::cout << "ElectronID: " << ZLe2Qid << std::endl;
            ///ORIGINAL CODE
            //reco::TransientTrack tt1 = theTTBuilder->build(lept1->gsfTrack());
            //reco::TransientTrack tt2 = theTTBuilder->build(lept2->gsfTrack());
            //std::pair<bool,Measurement1D> tkPVdistel1 = IPTools::absoluteImpactParameter3D(tt1,*PV);
             //std::pair<bool,Measurement1D> tkPVdistel2 = IPTools::absoluteImpactParameter3D(tt2,*PV);
            
            ////NEW CODE
            //reco::TrackRef dilepTk[2]={lept1->gsfTrack(), lept2->gsfTrack()};

            reco::TrackRef dilepTk[2]={lept1->closestCtfTrackRef(), lept2->closestCtfTrackRef()};

            ///https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGsfElectronObject#Other_tracks
          
            std::vector<reco::TransientTrack> LLTTks;
            //std::cout << "Ok " << std::endl;
            try{
                LLTTks.push_back(theTTBuilder->build(&dilepTk[0]));
                LLTTks.push_back(theTTBuilder->build(&dilepTk[1]));
                std::cout << "ok" << std::endl;
            }
            catch (...){
                std::cout<<"Bad electron, i.e. no track" << std::endl;
                continue;
            }
            std::pair<bool,Measurement1D> tkPVdistel1 = IPTools::absoluteImpactParameter3D(LLTTks.at(0),*PV);
            std::pair<bool,Measurement1D> tkPVdistel2 = IPTools::absoluteImpactParameter3D(LLTTks.at(1),*PV);
          
            if (!tkPVdistel1.first || !tkPVdistel2.first) continue;

            //std::cout << "pass Track builder " << std::endl;
		        
            /////////////////////////////////////////////////////
            //    S t a r t   t h e   K i n e m a t i c   F i t //
            //////////////////////////////////////////////////////
            const ParticleMass elMass(0.000511);
            const ParticleMass muMass(0.10565837);
            float muSigma = 1.e-6;
            float elSigma = 1.e-6;

            KinematicParticleFactoryFromTransientTrack pFactory;

			std::vector<RefCountedKinematicParticle> ZDaughters;
			ZDaughters.push_back(pFactory.particle (MuMuTTks[0], muMass, float(0), float(0), muSigma));
			ZDaughters.push_back(pFactory.particle (MuMuTTks[1], muMass, float(0), float(0), muSigma));
            //ZDaughters.push_back(pFactory.particle (tt1 , elMass, float(0), float(0), elSigma));
            //ZDaughters.push_back(pFactory.particle (tt2 , elMass, float(0), float(0), elSigma));
          
            ZDaughters.push_back(pFactory.particle (LLTTks[0] , elMass, float(0), float(0), elSigma));
            ZDaughters.push_back(pFactory.particle (LLTTks[1] , elMass, float(0), float(0), elSigma));

			KinematicParticleVertexFitter ZVertexFitter;
			RefCountedKinematicTree ZTree = ZVertexFitter.fit(ZDaughters);
					
            //std::cout << "is ZTree empty? " << ZTree->isEmpty() << std::endl;
			if (ZTree->isEmpty())continue;
			ZTree->movePointerToTheTop();
            RefCountedKinematicParticle fitZ = ZTree->currentParticle();
			RefCountedKinematicVertex ZDecayVertex = ZTree->currentDecayVertex();

            //std::cout << "is FitZ Valid ??" << fitZ->currentState().isValid() << std::endl;
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
                //if (ZDecayVertex->chiSquared() < 0) continue;
                ZM_fit  = fitZ->currentState().mass();
                ZPx_fit = fitZ->currentState().kinematicParameters().momentum().x();
                ZPy_fit = fitZ->currentState().kinematicParameters().momentum().y();
                ZPz_fit = fitZ->currentState().kinematicParameters().momentum().z();
                ZVtxX_fit = ZDecayVertex->position().x();
                ZVtxY_fit = ZDecayVertex->position().y();
                ZVtxZ_fit = ZDecayVertex->position().z();
                ZVtxP_fit = ChiSquaredProbability((double)(ZDecayVertex->chiSquared()), (double)(ZDecayVertex->degreesOfFreedom()));
                //if (ZVtxP_fit <= 0.0) continue;
                passFit = 1;
            }
            else{
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
				TLorentzVector m1;
				TLorentzVector m2;
				TLorentzVector l1;
				TLorentzVector l2;
                m1.SetPtEtaPhiM(muon1->pt(), muon1->eta(),muon1->phi(),muon1->mass());
                m2.SetPtEtaPhiM(muon2->pt(), muon2->eta(),muon2->phi(),muon2->mass());
                l1.SetPtEtaPhiM(lept1->pt(), lept1->eta(),lept1->phi(),lept1->mass());
                l2.SetPtEtaPhiM(lept2->pt(), lept2->eta(),lept2->phi(),lept2->mass());
                TLorentzVector theZ = m1+m2+l1+l2;
                // Get the Z boson
                std::cout << "Z Fit diff: " << ZM_fit - gen_z_p4.M() << std::endl;
                std::cout << "Z msrd diff: "<< theZ.M() - gen_z_p4.M() << std::endl;
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
				patZ.addUserFloat("ZvtxX",ZVtxX_fit);
				patZ.addUserFloat("ZvtxY",ZVtxY_fit);
				patZ.addUserFloat("ZvtxZ",ZVtxZ_fit);
				patZ.addUserFloat("dRm1m2",dRm1m2);	
				patZ.addUserFloat("dRl1l2",dRel1el2);	
				patZ.addUserFloat("dRl1m1",dRel1mu1);	
				patZ.addUserFloat("dRl1m2",dRel1mu2);	
				patZ.addUserFloat("dRl2m1",dRel2mu1);	
				patZ.addUserFloat("dRl2m2",dRel2mu2);	

                ////////////////////////////////
                ///////// MY MUON Q ID /////////
                ////////////////////////////////
				
                bool child = ZTree->movePointerToTheFirstChild();
				//get first muon
                RefCountedKinematicParticle fitMu1 = ZTree->currentParticle();
                float mu1M_fit;
                float mu1Q_fit ;
                float mu1Px_fit;
                float mu1Py_fit;
                float mu1Pz_fit;
          
            	if (!child){
                    std::cout << "Mu1" << std::endl;
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

			    //int muId ;
				//if (mu1Q_fit > 0 ) muId = 13;
				//else  muId = -13;	 
				reco::CompositeCandidate recoMu1(mu1Q_fit, math::XYZTLorentzVector(mu1Px_fit, mu1Py_fit, mu1Pz_fit, 
		                                             sqrt(mu1M_fit*mu1M_fit + mu1Px_fit*mu1Px_fit + mu1Py_fit*mu1Py_fit + 
		                                             mu1Pz_fit*mu1Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
                                
                //save measured values
                reco::CompositeCandidate msrdMu1(muon1->charge(), math::XYZTLorentzVector(muon1->px(), muon1->py(), muon1->pz(),
                                                 sqrt((muon1->mass())*(muon1->mass()) + (muon1->px())*(muon1->px()) + (muon1->py())*(muon1->py()) +
                                                 (muon1->pz())*(muon1->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
                pat::CompositeCandidate pat_msrdMu1(msrdMu1);

			    pat::CompositeCandidate patMu1(recoMu1);
           		//patMu1.addUserFloat("mu1Q_", muon1->charge());
                patMu1.addUserFloat("Dxy",mdxy1);
           		patMu1.addUserFloat("Dz",mdz1);
                //new
                patMu1.addUserFloat("dRIso"    ,getIso( *muon1 ) );
                patMu1.addUserFloat("dIP3DSig", tkPVdist1.second.significance());
                //patMu1.addUserFloat("dR")
          
           		patMu1.addUserFloat("dIP3D",tkPVdist1.second.value());
           		patMu1.addUserFloat("dIP3DErr",tkPVdist1.second.error());
                ////////////////////////////////
                ///////// MY MUON Q ID /////////
                ////////////////////////////////
                //patMu1.addUserFloat("JpM1Qid_", JpM1Qid);
                                
                patMu1.addUserFloat("muon1Mu8DiEle12_", muon1Mu8DiEle12);
                       
                patMu1.addUserInt("JpM1Qid_", ZMu1Qid);
                //patMu1.addUserFloat("JpM1Qid_TP_", JpM1Qid_TP);
  
                patMu1.addUserFloat("psiM1_TrackerLWM_", psiM1_TrackerLWM);
                patMu1.addUserFloat("psiM1_PixelLWM_",  psiM1_PixelLWM);
                patMu1.addUserFloat("psiM1_ValPixHit_", psiM1_ValPixHit);
                //get second muon
                child = ZTree->movePointerToTheNextChild();
			    RefCountedKinematicParticle fitMu2 = ZTree->currentParticle();
                float mu2M_fit ;
                float mu2Q_fit ;
                float mu2Px_fit;
                float mu2Py_fit;
                float mu2Pz_fit;
          
		        if (!child){
                    std::cout << "Mu2" << std::endl;
                    mu2M_fit  = 0;
                    mu2Q_fit  = 0;
                    mu2Px_fit = 0;
                    mu2Py_fit = 0;
                    mu2Pz_fit = 0;
                }
                else{
			        mu2M_fit  = fitMu2->currentState().mass();
		            mu2Q_fit  = fitMu2->currentState().particleCharge();
		            mu2Px_fit =    fitMu2->currentState().kinematicParameters().momentum().x();
			        mu2Py_fit =    fitMu2->currentState().kinematicParameters().momentum().y();
			        mu2Pz_fit =    fitMu2->currentState().kinematicParameters().momentum().z();
                }
		        reco::CompositeCandidate recoMu2(mu2Q_fit, math::XYZTLorentzVector(mu2Px_fit, mu2Py_fit, mu2Pz_fit,
                                                 sqrt(mu2M_fit*mu2M_fit + mu2Px_fit*mu2Px_fit + mu2Py_fit*mu2Py_fit +
		                                         mu2Pz_fit*mu2Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
		    		
                reco::CompositeCandidate msrdMu2(muon2->charge(), math::XYZTLorentzVector(muon2->px(), muon2->py(), muon2->pz(),
                                                 sqrt((muon2->mass())*(muon2->mass()) + (muon2->px())*(muon2->px()) + (muon2->py())*(muon2->py()) +
                                                 (muon2->pz())*(muon2->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
                pat::CompositeCandidate pat_msrdMu2(msrdMu2);

			    pat::CompositeCandidate patMu2(recoMu2);
                patMu2.addUserFloat("Dxy",mdxy2);
                patMu2.addUserFloat("Dz",mdz2);
                //New
                patMu2.addUserFloat("dRIso"    ,getIso( *muon2 ) );
                patMu2.addUserFloat("dIP3DSig", tkPVdist2.second.significance());
          
                patMu2.addUserFloat("dIP3D",tkPVdist2.second.value());
                patMu2.addUserFloat("dIP3DErr",tkPVdist2.second.error());
		
                //patMu2.addUserFloat("JpM2Qid_", JpM2Qid);

                patMu2.addUserFloat("muon2Mu8DiEle12_", muon2Mu8DiEle12);
		             
                patMu2.addUserInt("JpM2Qid_", ZMu2Qid);
                //patMu2.addUserFloat("JpM2Qid_TP_", JpM2Qid_TP);

                patMu2.addUserFloat("psiM2_TrackerLWM_", psiM2_TrackerLWM);
                patMu2.addUserFloat("psiM2_PixelLWM_",  psiM2_PixelLWM);
                patMu2.addUserFloat("psiM2_ValPixHit_", psiM2_ValPixHit);
	
                //Define Onia from two muon and the information of Onia2Mumu
                pat::CompositeCandidate jpsi;
                jpsi.addDaughter(patMu1,"muon1");
                jpsi.addDaughter(patMu2,"muon2");
                jpsi.setP4(patMu1.p4()+patMu2.p4());
				jpsi.addUserFloat("vProb",jpsiVprob);
				jpsi.addUserFloat("vChi2",jpsiChi2);
                //save measured values on a different daughter
                pat::CompositeCandidate msrd_jpsi;
                msrd_jpsi.addDaughter(pat_msrdMu1,"msrd_muon1");
                msrd_jpsi.addDaughter(pat_msrdMu2,"msrd_muon2");
                msrd_jpsi.setP4(pat_msrdMu1.p4()+pat_msrdMu2.p4());
				msrd_jpsi.addUserFloat("vProb",jpsiVprob);
				msrd_jpsi.addUserFloat("vChi2",jpsiChi2);
				        


				//get Lepton
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
		             L1Px_fit =     fitL1->currentState().kinematicParameters().momentum().x();
			         L1Py_fit =     fitL1->currentState().kinematicParameters().momentum().y();
		             L1Pz_fit =     fitL1->currentState().kinematicParameters().momentum().z();
                }
			    reco::CompositeCandidate recoL1(L1Q_fit, math::XYZTLorentzVector(L1Px_fit, L1Py_fit, L1Pz_fit,
		                                        sqrt(L1M_fit*L1M_fit + L1Px_fit*L1Px_fit + L1Py_fit*L1Py_fit +
		                                        L1Pz_fit*L1Pz_fit)), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
                reco::CompositeCandidate msrdL1(lept1->charge(), math::XYZTLorentzVector(lept1->px(), lept1->py(), lept1->pz(),
                                                sqrt((lept1->mass())*(lept1->mass()) + (lept1->px())*(lept1->px()) + (lept1->py())*(lept1->py()) +
                                                (lept1->pz())*(lept1->pz()))), math::XYZPoint(ZVtxX_fit, ZVtxY_fit, ZVtxZ_fit));
                pat::CompositeCandidate pat_msrdL1(msrdL1);
	
                pat::CompositeCandidate patL1(recoL1);
           		patL1.addUserFloat("Dxy"	,ldxy1);
           		patL1.addUserFloat("Dz"		,ldz1);

                patL1.addUserFloat("Dxy_gsf"        , lept1->gsfTrack()->dxy(PV->position()));
                patL1.addUserFloat("Dz_gsf"         , lept1->gsfTrack()->dz(PV->position()));
                
                //only data, no userFloat for mc
                float corrEt_1;//  lept1->et() * lept1->userFloat("ecalTrkEnergyPostCorr")/lept1->energy();
                float corrfactor_1;// lept1->userFloat("ecalTrkEnergyPostCorr") / lept1->energy();
                try {
                    corrEt_1 = lept1->et() * lept1->userFloat("ecalTrkEnergyPostCorr")/lept1->energy();
                    corrfactor_1 = lept1->userFloat("ecalTrkEnergyPostCorr") / lept1->energy();
                }
                catch (...) {
                    corrEt_1 = 0;
                    corrfactor_1 = 0;
                }
                patL1.addUserFloat("corrEt_", corrEt_1);
                patL1.addUserFloat("corrfactor_", corrfactor_1);
                
                patL1.addUserFloat("dRIso"	,getIsoVar( *lept1 ) );
                //New
                patL1.addUserFloat("dIP3DSig",tkPVdistel1.second.significance());
          
                patL1.addUserFloat("dIP3D"	,tkPVdistel1.second.value());
                patL1.addUserFloat("dIP3DErr"	,tkPVdistel1.second.error());

                patL1.addUserFloat("dPhiInSeed" , lept1->deltaPhiSuperClusterTrackAtVtx());
                patL1.addUserFloat("dEtaInSeed" , getEtaInSeed( *lept1 )) ;
                patL1.addUserFloat("SigmaIEtaIEta" ,lept1->full5x5_sigmaIetaIeta());
                patL1.addUserFloat("SigmaIPhiIPhi" ,lept1->full5x5_sigmaIphiIphi());
                patL1.addUserFloat("HoverE"     , lept1->hcalOverEcal());
                //patL1.addUserInt("ElecMissHits" , lept1->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
                patL1.addUserFloat("ElecMissHits" , ZLe1_ElecMissHits);
                //patL1.addUserInt("ElecMissHits" , lept1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) );
                //Due to changes in function names from 80X to 94X, the function numberOfHits
                // MY Muon ID
				 
                patL1.addUserInt("ZLe1Qid_", ZLe1Qid);
	                        
                //patL1.addUserFloat("ZLe1Qid_TP_", ZLe1Qid_TP);

                patL1.addUserFloat("ZLe1_TrackerLWM_", ZLe1_TrackerLWM);
                patL1.addUserFloat("ZLe1_PixelLWM_", ZLe1_PixelLWM);
                patL1.addUserFloat("ZLe1_ValPixHit_", ZLe1_ValPixHit);
                patL1.addUserFloat("lept1Ele25wpT_", lept1Ele25wpT);
                patL1.addUserFloat("lept1Ele23_12_", lept1Ele23_12);
                patL1.addUserFloat("lept1Mu8DiEle12_", lept1Mu8DiEle12);
                //get Lepton
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
                patL2.addUserFloat("Dxy"	,ldxy2);
           		patL2.addUserFloat("Dz"		,ldz2);

                patL2.addUserFloat("Dxy_gsf"        , lept2->gsfTrack()->dxy(PV->position()));
                patL2.addUserFloat("Dz_gsf"         , lept2->gsfTrack()->dz(PV->position()));
                
                //only data, no userFloat for mc
                float corrEt_2;//lept2->et() * lept2->userFloat("ecalTrkEnergyPostCorr")/lept2->energy();
                float corrfactor_2;//lept2->userFloat("ecalTrkEnergyPostCorr") / lept2->energy();
                try {
                    corrEt_2 = lept2->et() * lept1->userFloat("ecalTrkEnergyPostCorr")/lept2->energy();
                    corrfactor_2 = lept2->userFloat("ecalTrkEnergyPostCorr") / lept2->energy();
                }
                catch (...) {
                    corrEt_2 = 0;
                    corrfactor_2 = 0;
                }
                patL2.addUserFloat("corrEt_", corrEt_2);
                patL2.addUserFloat("corrfactor_", corrfactor_2);
                
                patL2.addUserFloat("dRIso"	,getIsoVar( *lept2 ) );
                patL2.addUserFloat("dIP3DSig",tkPVdistel2.second.significance());
                patL2.addUserFloat("dIP3D"	,tkPVdistel2.second.value());
                patL2.addUserFloat("dIP3DErr"	,tkPVdistel2.second.error());

                patL2.addUserFloat("dPhiInSeed" , lept2->deltaPhiSuperClusterTrackAtVtx());
                patL2.addUserFloat("dEtaInSeed" , getEtaInSeed( *lept2 )) ;
                patL2.addUserFloat("SigmaIEtaIEta" ,lept2->full5x5_sigmaIetaIeta());
                patL2.addUserFloat("SigmaIPhiIPhi" ,lept2->full5x5_sigmaIphiIphi());
                patL2.addUserFloat("HoverE"     , lept2->hcalOverEcal());
                //patL2.addUserInt("ElecMissHits" , lept2->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) );
                patL2.addUserFloat("ElecMissHits" , ZLe2_ElecMissHits);
                //patL2.addUserInt("ElecMissHits" , lept2->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) );
                //Due to changes in function names from 80X to 94X, the function numberOfHits

				// My leoton-electron Q id
                patL2.addUserInt("ZLe2Qid_", ZLe2Qid);

                //patL2.addUserFloat("ZLe2Qid_TP_", ZLe2Qid_TP);
	              
                patL2.addUserFloat("ZLe2_TrackerLWM_", ZLe2_TrackerLWM);
                patL2.addUserFloat("ZLe2_PixelLWM_", ZLe2_PixelLWM);
                patL2.addUserFloat("ZLe2_ValPixHit_", ZLe2_ValPixHit);

                patL2.addUserFloat("lept2Ele25wpT_", lept2Ele25wpT);
                patL2.addUserFloat("lept2Ele23_12_", lept2Ele23_12);
                patL2.addUserFloat("lept2Mu8DiEle12_", lept2Mu8DiEle12);
                
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

				patZ.addDaughter(jpsi,"jpsi");
                patZ.addDaughter(msrd_jpsi, "msrd_jpsi");
                patZ.addDaughter(patL1,"lepton1");
                patZ.addDaughter(patL2,"lepton2");
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
                Event_Cand++;
                
				ZCandColl->push_back(patZ);
                std::cout<< "something has been pushed"	<< std::endl;
			
			//} Z candidate is valid

	    }

   }
  iEvent.put(std::move(ZCandColl),"ZCandidates");
  //std::cout << "is ZCandColl empty ?" << /*ZCandColl->empty() <<*/ std::endl; 
  //std::cout << "jpsiElecKmcFitter is working ok" << std::endl;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
//void
//jpsiElecKmcFitter::beginStream(edm::StreamID)
//{
//}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
//void
//jpsiElecKmcFitter::endStream() {
//}

// ------------ method called when starting to processes a run  ------------
/*
void
jpsiElecKmcFitter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
jpsiElecKmcFitter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
jpsiElecKmcFitter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
jpsiElecKmcFitter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
Float_t jpsiElecKmcFitter::getIso(const pat::Muon& mu){
		Float_t coriso = 99.0; 
		reco::MuonPFIsolation pfR03 = mu.pfIsolationR03();
	        coriso = pfR03.sumChargedHadronPt + std::max(0., pfR03.sumNeutralHadronEt+pfR03.sumPhotonEt-0.5*pfR03.sumPUPt);
return coriso; 
}

Float_t jpsiElecKmcFitter::getIsoVar(const pat::Electron& el){
                reco::GsfElectron::PflowIsolationVariables pfIso1 = el.pfIsolationVariables();
                float coriso1 = pfIso1.sumChargedHadronPt + std::max(0., pfIso1.sumNeutralHadronEt+pfIso1.sumPhotonEt-0.5*pfIso1.sumPUPt);
return coriso1;
}

Float_t jpsiElecKmcFitter::getEtaInSeed(const pat::Electron& el){
        float dEtaInSeed_ =  el.superCluster().isNonnull() && el.superCluster()->seed().isNonnull() ?
                                   el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
return dEtaInSeed_;
}

bool jpsiElecKmcFitter::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool jpsiElecKmcFitter::IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk->eta());
  double DeltaP   = fabs(mu.p()-tk->p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

int jpsiElecKmcFitter::convertBinaryToDecimal(unsigned long long n)
{
    int decimalNumber = 0, i = 0, remainder;
    while (n!=0)
    {
        remainder = n%10;
        n /= 10;
        decimalNumber += remainder*pow(2,i);
        ++i;
    }
    return decimalNumber;
}

void jpsiElecKmcFitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//////////////////////////////////////////////////////////////////////////////////////////
float jpsiElecKmcFitter::ElectronRelIso(const reco::Candidate *cand, float rho) {
  pat::Electron el = *((pat::Electron*)cand);
  float relIsoWithEA = 0;
  const int nEtaBins = 5;
  const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
  const float effectiveAreaValues[nEtaBins] = {0.1013, 0.0988, 0.0572, 0.0842, 0.1530};
  reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
  float etaSC = el.superCluster()->eta();
  // Find eta bin first. If eta>2.5, the last eta bin is used.
  int etaBin = 0;
  while(etaBin < nEtaBins-1 && abs(etaSC) > etaBinLimits[etaBin+1]) ++etaBin;
  float area = effectiveAreaValues[etaBin];
  relIsoWithEA = (float)(pfIso.sumChargedHadronPt+std::max(float(0.0), pfIso.sumNeutralHadronEt+pfIso.sumPhotonEt-rho*area))/el.pt();
  return relIsoWithEA;
}
//////////////////////////////////////////////////////////////////////////////////////////
float jpsiElecKmcFitter::MuonRelIso(const reco::Candidate *cand, float rho) {
  pat::Muon mu = *((pat::Muon*)cand);
  float relIsoWithEA = 0.001;
  const int nEtaBins = 5;
  const float etaBinLimits[nEtaBins+1] = {0.0, 0.8, 1.3, 2.0, 2.2, 2.5};
  const float effectiveAreaValues[nEtaBins] = {0.0913, 0.0765, 0.0546, 0.0728, 0.1177};
  reco::MuonPFIsolation pfIso = mu.pfIsolationR03();
  // Find eta bin first. If eta>2.5, the last eta bin is used.
  int etaBin = 0;
  while(etaBin < nEtaBins-1 && abs(mu.eta()) > etaBinLimits[etaBin+1]) ++etaBin;
  float area = effectiveAreaValues[etaBin];
  relIsoWithEA = (float)(pfIso.sumChargedHadronPt+std::max(float(0.0),pfIso.sumNeutralHadronEt+pfIso.sumPhotonEt-rho*area))/mu.pt();
  return relIsoWithEA;
}
//define this as a plug-in
DEFINE_FWK_MODULE(jpsiElecKmcFitter);
