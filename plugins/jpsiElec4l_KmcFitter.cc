// -*- C++ -*-
//
// Package:    AnalyzeZll/jpsiElec4l_KmcFitter
// Class:      jpsiElec4l_KmcFitter
// 
/**\class jpsiElec4l_KmcFitter jpsiElec4l_KmcFitter.cc AnalyzeZll/jpsiElec4l_KmcFitter/plugins/jpsiElec4l_KmcFitter.cc

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

class jpsiElec4l_KmcFitter : public edm::stream::EDProducer<>{
   public:
      explicit jpsiElec4l_KmcFitter(const edm::ParameterSet&);
      ~jpsiElec4l_KmcFitter() override{};

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
   // virtual void beginStream(edm::StreamID) override;
      Float_t getIso(const pat::Muon& );
      Float_t getIsoVar(const pat::Electron&);       
      Float_t getEtaInSeed(const pat::Electron&);
    
      Float_t ElectronRelIso(const pat::Electron&);

      void printMCtree(const reco::Candidate *, int);
      std::string printName(int);

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
        //edm::EDGetToken electronsMiniAODToken_;
        //edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
    
        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;
    
        //edm::EDGetTokenT<double> fixedGridRhoFastjetAll_;
        edm::Handle<double> rhoH;
        edm::EDGetTokenT<double> rhoHToken;

};


jpsiElec4l_KmcFitter::jpsiElec4l_KmcFitter(const edm::ParameterSet& iConfig){
	dimuon_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter< edm::InputTag>("dimuon"));
	dielec_Label = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("dilepton"));
    //electronsMiniAODToken_ = mayConsume<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronsMiniAOD"));

    //electronsMiniAODToken_ = mayConsume<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronsMiniAOD"));
    //muonToken_ = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
    
	primaryVertices_Label = consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"));
    
    genCands_ = consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"));
    packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
    
    //fixedGridRhoFastjetAll_ = consumes<double> (iConfig.getParameter <edm::InputTag>("fixedGridRhoFastjetAll"));
    rhoHToken = consumes<double>(iConfig.getUntrackedParameter("fixedGridRhoFastjetAll",edm::InputTag("fixedGridRhoFastjetAll")));

   	produces<pat::CompositeCandidateCollection>("ZCandidates"); 
 
}


//jpsiElec4l_KmcFitter::~jpsiElec4l_KmcFitter(){}


//
// member functions
//
//Try to print mc Tree
std::string jpsiElec4l_KmcFitter::printName(int pdgid){
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

void jpsiElec4l_KmcFitter::printMCtree(const reco::Candidate* mother, int indent=0){
    if (mother == NULL){
         std::cout << "end tree" << std::endl;
         return;
    }
    if (mother->numberOfDaughters() > 0){   
        if(indent){
                std::cout << std::setw(indent) << " ";
        }  
        std::cout << printName(mother->pdgId()) <<" has "<< mother->numberOfDaughters() <<" daughters " <<std::endl;
    }
    int extraIndent = 0;    
    for (size_t i=0; i< mother->numberOfDaughters(); i++){
        const reco::Candidate * daughter = mother->daughter(i);
        if (mother->numberOfDaughters() > 0){       
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
//recursively check is a given particle is ancestor
bool jpsiElec4l_KmcFitter::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool jpsiElec4l_KmcFitter::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}


// ------------ method called to produce the data  ------------
void jpsiElec4l_KmcFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::auto_ptr<pat::ElectronCollection > selectedCollection(new pat::ElectronCollection ); //new
  std::unique_ptr<pat::CompositeCandidateCollection> ZCandColl(new pat::CompositeCandidateCollection); 
     
  edm::Handle<pat::CompositeCandidateCollection> dimuons;
  iEvent.getByToken(dimuon_Label,dimuons);
  
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  iEvent.getByToken(dielec_Label,dileptons);

  //edm::Handle<edm::View<pat::Electron> > electrons; //new
  //if( !electrons.isValid() ) iEvent.getByToken(electronsMiniAODToken_,electrons); //new
  
  //edm::Handle<View<pat::Muon>> muons;
  //iEvent.getByToken(muonToken_, muons);
    
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
  
  // Rho values
  //edm::Handle<double> fixedGridRhoFastjetAllH;
  //iEvent.getByToken(fixedGridRhoFastjetAll_, fixedGridRhoFastjetAllH);
  iEvent.getByToken(rhoHToken, rhoH);

  //////////////////////////////////
  //// Select  the best PV      ////
  //////////////////////////////////
  

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
  //NEED TO RE-MAKE TO READ Z --> 4L (IN THIS CASE 2MU 2EL)
  gen_z_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lepton1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_lepton2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_z_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  int tst = 0;
  int n_Z_dau = 0;
  int Event_Cand = 1;
  //std::cout << "test" << std::endl;
  //NEW MC ALV
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
                if (mom->numberOfDaughters() == 1){
                    if (mom->daughter(0)->pdgId() == 23) zfromZ++; 
                }
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
                if (mfromZ >= 2 && efromZ >= 2 && tfromZ == 0){    
                    if (tst > 2) break;
                    std::cout<< "Found Z with mass: "<< mom->mass() <<", print Tree ... " << tst << std::endl;                    
                    printMCtree(mom, 0);
                    for(size_t k=0; k<packed->size(); k++){
                       const reco::Candidate * stable_dau = &(*packed)[k];
                       int stable_id = (*packed)[k].pdgId();
                       if (stable_dau != nullptr && isAncestor(mom,stable_dau)) { 
                           if (stable_id != 13 && stable_id != -13 && stable_id != 11 && stable_id != -11) continue;
                           if(stable_id == 11 && temp_lep_1.M() == 0){ // if muon- && not previiulsy assigned
                                std::cout << " matched 1 "  << " with Pt: ";
                                std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_lep_1.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                            }
                            else if(stable_id == -11 && temp_lep_2.M() == 0){ // if muon+ && not previusly assigned
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
                        gen_lepton1_p4 = temp_lep_1;
                        gen_lepton2_p4 = temp_lep_2;
                        gen_muon1_p4 = temp_mu_1;
                        gen_muon2_p4 = temp_mu_2; 
                       /* 
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
                       }*/ 
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
  if (tst) std::cout << "x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x NEW EVENT -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x" << std::endl; 
  //NEW for muons and psi pairs
  //NEW for muons and psi pairs
  
  
  int nonia = dimuons->size();
  int nmuons = dileptons->size()*2;
  int nPV    = primaryVertices_handle->size();
  if (gen_z_p4.M() != 0 && nonia && nmuons){
  std::cout<< "+X+X+X+X+X+X+X+X+X+X+X+X+X+X START READING DATA X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X" << std::endl;
  std::cout<< "dimuon length: " << nonia << std::endl;
  std::cout<< "dilepton length:  " << nmuons << std::endl;
  for (pat::CompositeCandidateCollection::const_iterator dilepton = dileptons->begin(); dilepton != dileptons->end(); ++dilepton){
  for (pat::CompositeCandidateCollection::const_iterator dimuon = dimuons->begin(); dimuon != dimuons->end() ; ++dimuon){
        //Jpsi Muons
        std::cout << "test 454 " << std::endl;
	    const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon1"));
	    const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(dimuon->daughter("muon2"));

        const pat::Electron* lept1 = dynamic_cast<const pat::Electron*>(dilepton->daughter("lepton1"));
        const pat::Electron* lept2 = dynamic_cast<const pat::Electron*>(dilepton->daughter("lepton2"));
        std::cout << "test 460 "<<std::endl;    
        //reco::TrackRef glbTrack_l1 = lept1->track(); //leads to segmentation fault :(
        //reco::TrackRef glbTrack_l2 = lept2->track();
        reco::TrackRef glbTrack_l1 = lept1->closestCtfTrackRef();
        reco::TrackRef glbTrack_l2 = lept2->closestCtfTrackRef();
        reco::TrackRef glbTrack_m1 = muon1->track();
        reco::TrackRef glbTrack_m2 = muon2->track();
        std::cout << "test 465 " << std::endl;
        if(glbTrack_l1.isNull() || glbTrack_l2.isNull() || glbTrack_m1.isNull() || glbTrack_m2.isNull()) continue;
        std::cout << "test 467 " << std::endl;
        if(!(glbTrack_m1->quality(reco::TrackBase::highPurity))) continue;
        if(!(glbTrack_m2->quality(reco::TrackBase::highPurity))) continue;
        std::cout << "test 469 " << std::endl;
	    float jpsiVprob=0;
	    float jpsiChi2=0;
	    //jpsiVprob = dimuon->userFloat("vProb");
	    //jpsiChi2  = dimuon->userFloat("vNChi2");
      
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
        std::cout << "test 548 " << std::endl;
        int psiM1_TrackerLWM = muon1->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
        int psiM1_PixelLWM   = muon1->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
        int psiM1_ValPixHit  = muon1->muonBestTrack()->hitPattern().numberOfValidPixelHits();
        int psiM2_TrackerLWM = muon2->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
        int psiM2_PixelLWM   = muon2->muonBestTrack()->hitPattern().pixelLayersWithMeasurement();
        int psiM2_ValPixHit  = muon2->muonBestTrack()->hitPattern().numberOfValidPixelHits();
 	    std::cout << "test 555 " << std::endl;
        reco::TransientTrack lept1TT((*theTTBuilder).build(glbTrack_l1));
        reco::TransientTrack lept2TT((*theTTBuilder).build(glbTrack_l2));
        reco::TransientTrack muon1TT((*theTTBuilder).build(glbTrack_m1));
        reco::TransientTrack muon2TT((*theTTBuilder).build(glbTrack_m2));
        std::cout << "test 560 " << std::endl; 
        std::pair<bool, Measurement1D> tkPVdistm1 = IPTools::absoluteImpactParameter3D(muon1TT, *PV);
        std::pair<bool, Measurement1D> tkPVdistm2 = IPTools::absoluteImpactParameter3D(muon2TT, *PV);
        std::cout<< "test 564 " << std::endl;
        std::pair<bool, Measurement1D> tkPVdistl1 = IPTools::absoluteImpactParameter3D(lept1TT, *PV);
        std::pair<bool, Measurement1D> tkPVdistl2 = IPTools::absoluteImpactParameter3D(lept2TT, *PV);
        std:: cout << "test 566 "<< std::endl;
        if(!tkPVdistl1.first || !tkPVdistl2.first || !tkPVdistm1.first || !tkPVdistm2.first) continue;
         
        //cambiar pues ahora son electrones
        int lept1Ele25wpT = 0;
        int lept2Ele25wpT = 0;
        int lept1Ele23_12 = 0;
        int lept2Ele23_12 = 0;
        int lept1Mu8DiEle12 = 0;
        int lept2Mu8DiEle12 = 0;

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
                ZLe1Qid_n += 10000000000;
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
                ZLe2Qid_n += 10000000000;
                //std::cout << "2L is 94X iso V2 wp90 " << std::endl;
            }
        ZLe1Qid = convertBinaryToDecimal(ZLe1Qid_n);
        ZLe2Qid = convertBinaryToDecimal(ZLe2Qid_n);
        std::cout << "test 685  "<< std::endl;  
        int ZLe1_TrackerLWM = lept1->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
        int ZLe1_PixelLWM   = lept1->gsfTrack()->hitPattern().pixelLayersWithMeasurement();
        int ZLe1_ValPixHit  = lept1->gsfTrack()->hitPattern().numberOfValidPixelHits();
 		
        int ZLe2_TrackerLWM = lept2->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
        int ZLe2_PixelLWM   = lept2->gsfTrack()->hitPattern().pixelLayersWithMeasurement();
        int ZLe2_ValPixHit  = lept2->gsfTrack()->hitPattern().numberOfValidPixelHits();

        int ZLe1_ElecMissHits = lept1->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);
        int ZLe2_ElecMissHits = lept2->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);

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
        float mdz1  = muon1->muonBestTrack()->dz(PV->position());
        float mdxy2 = muon2->muonBestTrack()->dxy(PV->position());
        float mdz2  = muon2->muonBestTrack()->dz(PV->position());
        float ldxy1 = lept1->bestTrack()->dxy(PV->position());
        float ldz1  = lept1->bestTrack()->dz(PV->position());
        float ldxy2 = lept2->bestTrack()->dxy(PV->position());
        float ldz2  = lept2->bestTrack()->dz(PV->position());
		
        if ( dRel1mu1 <0.01 || dRel1mu2<0.01 || dRel2mu1 <0.01 || dRel2mu2<0.01 ) continue;
        std::cout << "test 717" << std::endl;
        ///ORIGINAL CODE
        //reco::TransientTrack tt1 = theTTBuilder->build(lept1->gsfTrack());
        //reco::TransientTrack tt2 = theTTBuilder->build(lept2->gsfTrack());
        //std::pair<bool,Measurement1D> tkPVdistel1 = IPTools::absoluteImpactParameter3D(tt1,*PV);
        //std::pair<bool,Measurement1D> tkPVdistel2 = IPTools::absoluteImpactParameter3D(tt2,*PV);
            
        ////NEW CODE
        //reco::TrackRef dilepTk[2]={lept1->gsfTrack(), lept2->gsfTrack()};
        //reco::TrackRef dilepTk[2]={lept1->closestCtfTrackRef(), lept2->closestCtfTrackRef()};
        ///https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideGsfElectronObject#Other_tracks
        //std::vector<reco::TransientTrack> LLTTks;
		        
        /////////////////////////////////////////////////////
        //    S t a r t   t h e   K i n e m a t i c   F i t //
        //////////////////////////////////////////////////////
        const ParticleMass elMass(0.000511);
        const ParticleMass muMass(0.10565837);
        float muSigma = 1.e-6;
        float elSigma = 1.e-6;

        KinematicParticleFactoryFromTransientTrack pFactory;
        std::vector<RefCountedKinematicParticle> ZDaughters;
        
        float chi = 0;
        float ndf = 0;
         
        try{
            ZDaughters.push_back(pFactory.particle(muon1TT, muMass, chi, ndf, muSigma));
            ZDaughters.push_back(pFactory.particle(muon2TT, muMass, chi, ndf, muSigma));
            ZDaughters.push_back(pFactory.particle(lept1TT, elMass, chi, ndf, elSigma));
            ZDaughters.push_back(pFactory.particle(lept2TT, elMass, chi, ndf, elSigma));
        }
        catch (...){
            continue;
        }
        std::cout << "test 753 " << std::endl;
        KinematicParticleVertexFitter ZVertexFitter;
        RefCountedKinematicTree ZTree;
        try{
           ZTree = ZVertexFitter.fit(ZDaughters);
        }
        catch (...){
           continue;
        }
        
        //std::cout << "is ZTree empty? " << ZTree->isEmpty() << std::endl;
        if (ZTree->isEmpty())continue;
         
        ZTree->movePointerToTheTop();
        RefCountedKinematicParticle fitZ = ZTree->currentParticle();
		RefCountedKinematicVertex ZDecayVertex = ZTree->currentDecayVertex();
        std::cout << "test 769 " << std::endl;
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
        std::cout << "test 803 " << std::endl;
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
        //std::cout << "Z Fit diff: " << ZM_fit - gen_z_p4.M() << std::endl;
        //std::cout << "Z msrd diff: "<< theZ.M() - gen_z_p4.M() << std::endl;
        if (theZ.M() < 60.0) continue;
        if (theZ.M() > 150.0) continue;
        reco::CompositeCandidate recoZ(0, math::XYZTLorentzVector(ZPx_fit, ZPy_fit, ZPz_fit,
                       sqrt(ZM_fit*ZM_fit + ZPx_fit*ZPx_fit + ZPy_fit*ZPy_fit +
                       ZPz_fit*ZPz_fit)), math::XYZPoint(ZVtxX_fit,
                       ZVtxY_fit, ZVtxZ_fit), 23);
        reco::CompositeCandidate msrdZ(0, math::XYZTLorentzVector( Z_px, Z_py, Z_pz,
                       sqrt(theZ.M()*theZ.M() + Z_px*Z_px + Z_py*Z_py +
                       Z_pz*Z_pz)), math::XYZPoint(ZVtxX_fit,
	      		       ZVtxY_fit, ZVtxZ_fit), 23);
        pat::CompositeCandidate patMZ(msrdZ);
		//pat::CompositeCandidate patZ(recoZ);
		pat::CompositeCandidate *patZ = new pat::CompositeCandidate(recoZ);
        //New
        std::cout << "test 832 " << std::endl;
        patZ->addUserInt("passFit_", passFit);
        patZ->addUserInt("nonia_", nonia );
        patZ->addUserInt("nmuons_",nmuons);
        patZ->addUserInt("nPV_",   nPV   );
          
        patZ->addUserFloat("vProb",ZVtxP_fit);
        patZ->addUserFloat("vChi2",ZDecayVertex->chiSquared());
        patZ->addUserFloat("ZvtxX",ZVtxX_fit);
        patZ->addUserFloat("ZvtxY",ZVtxY_fit);
        patZ->addUserFloat("ZvtxZ",ZVtxZ_fit);
        patZ->addUserFloat("dRm1m2",dRm1m2);
        patZ->addUserFloat("dRl1l2",dRel1el2);
        patZ->addUserFloat("dRl1m1",dRel1mu1);
        patZ->addUserFloat("dRl1m2",dRel1mu2);
        patZ->addUserFloat("dRl2m1",dRel2mu1);
        patZ->addUserFloat("dRl2m2",dRel2mu2);

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
        std::cout << "test 862 " << std::endl;   
        if (!child){
            //std::cout << "Mu1" << std::endl;
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
        patMu1.addUserFloat("dIP3DSig", tkPVdistm1.second.significance());
          
        patMu1.addUserFloat("dIP3D", tkPVdistm1.second.value());
        patMu1.addUserFloat("dIP3DErr", tkPVdistm1.second.error());
        ////////////////////////////////
        ///////// MY MUON Q ID /////////
        ////////////////////////////////
        //patMu1.addUserFloat("JpM1Qid_", JpM1Qid);
                                
        //patMu1.addUserFloat("muon1Mu8DiEle12_", muon1Mu8DiEle12);
                       
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
        std::cout << "test 923" << std::endl;  
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
        patMu2.addUserFloat("dIP3DSig", tkPVdistm2.second.significance());
          
        patMu2.addUserFloat("dIP3D", tkPVdistm2.second.value());
        patMu2.addUserFloat("dIP3DErr",tkPVdistm2.second.error());
		          
        //patMu2.addUserFloat("JpM2Qid_", JpM2Qid);
        //patMu2.addUserFloat("muon2Mu8DiEle12_", muon2Mu8DiEle12);
		             
        patMu2.addUserInt("JpM2Qid_", ZMu2Qid);
        //patMu2.addUserFloat("JpM2Qid_TP_", JpM2Qid_TP);
        std::cout << "test 963 " << std::endl;
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
        std::cout << "test 1020 " << std::endl;                
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
        std::cout << "test 1034 " << std::endl;        
        patL1.addUserFloat("dRIso"	,getIsoVar( *lept1 ) );
        //New
        patL1.addUserFloat("dIP3DSig", tkPVdistl1.second.significance());
        std::cout << "test 1038 " << std::endl;  
        patL1.addUserFloat("dIP3D"	, tkPVdistl1.second.value());
        patL1.addUserFloat("dIP3DErr"	, tkPVdistl1.second.error());
                
        patL1.addUserFloat("dRIsoEA", ElectronRelIso(*lept1));
        std::cout << "test 1043 " << std::endl;
        //std::cout << " dRisoEA l1 ~" << ElectronRelIso(*lept1) << std::endl;
        patL1.addUserFloat("trackMomentumAtVtx"   , (float)sqrt(lept1->trackMomentumAtVtx().mag2()));
        patL1.addUserFloat("ecalEnergy"           , (float)lept1->ecalEnergy());
        patL1.addUserFloat("full5x5_sigmaIetaIeta", (float)lept1->full5x5_sigmaIetaIeta());
        patL1.addUserFloat("dEtaIn"               , (float)lept1->deltaEtaSuperClusterTrackAtVtx());
        patL1.addUserFloat("dPhiIn"               , (float)lept1->deltaPhiSuperClusterTrackAtVtx());
        patL1.addUserFloat("HoE"                  , (float)lept1->hadronicOverEm());
        patL1.addUserFloat("ooEmooP"              , (float)fabs(1/lept1->ecalEnergy() - 1/sqrt(lept1->trackMomentumAtVtx().mag2())));
        patL1.addUserFloat("passConversionVeto"   , (float)lept1->passConversionVeto());
        std::cout << "test 1053 " << std::endl;
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
        std::cout << "test 1075 " << std::endl;
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
        patL2.addUserFloat("dIP3DSig",tkPVdistl2.second.significance());
        patL2.addUserFloat("dIP3D"	, tkPVdistl2.second.value());
        patL2.addUserFloat("dIP3DErr"	, tkPVdistl2.second.error());
    
        patL2.addUserFloat("dRIsoEA", ElectronRelIso(*lept2));
        patL2.addUserFloat("trackMomentumAtVtx"   , (float)sqrt(lept2->trackMomentumAtVtx().mag2()));
        patL2.addUserFloat("ecalEnergy"           , (float)lept2->ecalEnergy());
        patL2.addUserFloat("full5x5_sigmaIetaIeta", (float)lept2->full5x5_sigmaIetaIeta());
        patL2.addUserFloat("dEtaIn"               , (float)lept2->deltaEtaSuperClusterTrackAtVtx());
        patL2.addUserFloat("dPhiIn"               , (float)lept2->deltaPhiSuperClusterTrackAtVtx());
        patL2.addUserFloat("HoE"                  , (float)lept2->hadronicOverEm());
        patL2.addUserFloat("ooEmooP"              , (float)fabs(1/lept2->ecalEnergy() - 1/sqrt(lept2->trackMomentumAtVtx().mag2())));
        patL2.addUserFloat("passConversionVeto"   , (float)lept2->passConversionVeto());
    
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

        patZ->addDaughter(jpsi,"jpsi");
        patZ->addDaughter(msrd_jpsi, "msrd_jpsi");
        patZ->addDaughter(patL1,"lepton1");
        patZ->addDaughter(patL2,"lepton2");
        patZ->addDaughter(pat_msrdL1, "msrd_lepton1");
        patZ->addDaughter(pat_msrdL2, "msrd_lepton2");
        patZ->addDaughter(patMZ, "patMZ");
                
        patZ->addDaughter(pat_mcZ, "mcZ");
        patZ->addDaughter(pat_mcPsi, "mcPsi");
        patZ->addDaughter(pat_mcM1, "mcM1");
        patZ->addDaughter(pat_mcM2, "mcM2");
        patZ->addDaughter(pat_mcL1, "mcL1");
        patZ->addDaughter(pat_mcL2, "mcL2");
               
        patZ->addUserInt("Event_Cand_", Event_Cand);
        Event_Cand++;
		ZCandColl->push_back(*patZ);
        std::cout << "test 1216, zcandcoll size? , empty ?"<< ZCandColl->size()<< " | " <<ZCandColl->empty() << std::endl;        
        //std::cout<< "something has been pushed"	<< std::endl;
   } //end dimuon
  }// end dilepton
  std::cout << "test 1220 -<-<-<-<-<-<-<-<-< " << std::endl;
  iEvent.put(std::move(ZCandColl),"ZCandidates");
  std::cout << "is ZCandColl empty ?" << ZCandColl->empty() << std::endl;
  std::cout<< "+X+X+X+X+X+X+X+X+X+X+X+X+X+X END READING DATA X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X+X" << std::endl;
  }// end if good MC
  
  //std::cout << "jpsiElec4l_KmcFitter is working ok" << std::endl;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
//void
//jpsiElec4l_KmcFitter::beginStream(edm::StreamID)
//{
//}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
//void
//jpsiElec4l_KmcFitter::endStream() {
//}

// ------------ method called when starting to processes a run  ------------
/*
void
jpsiElec4l_KmcFitter::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
jpsiElec4l_KmcFitter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
jpsiElec4l_KmcFitter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
jpsiElec4l_KmcFitter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
Float_t jpsiElec4l_KmcFitter::getIso(const pat::Muon& mu){
		Float_t coriso = 99.0; 
		reco::MuonPFIsolation pfR03 = mu.pfIsolationR03();
	        coriso = pfR03.sumChargedHadronPt + std::max(0., pfR03.sumNeutralHadronEt+pfR03.sumPhotonEt-0.5*pfR03.sumPUPt);
return coriso; 
}

Float_t jpsiElec4l_KmcFitter::getIsoVar(const pat::Electron& el){
                reco::GsfElectron::PflowIsolationVariables pfIso1 = el.pfIsolationVariables();
                float coriso1 = pfIso1.sumChargedHadronPt + std::max(0., pfIso1.sumNeutralHadronEt+pfIso1.sumPhotonEt-0.5*pfIso1.sumPUPt);
return coriso1;
}

Float_t jpsiElec4l_KmcFitter::getEtaInSeed(const pat::Electron& el){
        float dEtaInSeed_ =  el.superCluster().isNonnull() && el.superCluster()->seed().isNonnull() ?
                                   el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta() : std::numeric_limits<float>::max();
return dEtaInSeed_;
}

bool jpsiElec4l_KmcFitter::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}
bool jpsiElec4l_KmcFitter::IsTheSame2(const reco::TrackRef& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk->eta());
  double DeltaP   = fabs(mu.p()-tk->p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

int jpsiElec4l_KmcFitter::convertBinaryToDecimal(unsigned long long n)
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

void jpsiElec4l_KmcFitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//////////////////////////////////////////////////////////////////////////////////////////
Float_t jpsiElec4l_KmcFitter::ElectronRelIso(const pat::Electron& el)
{
    float relIsoWithEA = 0;
    //pat::Electron el = *((pat::Electron*)cand);
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Spring15_selection_25ns
    // https://indico.cern.ch/event/370507/contribution/1/attachments/1140657/1633761/Rami_eleCB_ID_25ns.pdf
    // Effective areas from https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf
    // UPDATED VERSION : https://indico.cern.ch/event/732971/contributions/3022843/attachments/1658685/2656462/eleIdTuning.pdf
    const int nEtaBins = 7;
    const float etaBinLimits[nEtaBins+1] = {0.0, 1.0, 1.479, 2.0, 2.2, 2.3, 2.4, 2.5};
    //const float effectiveAreaValues[nEtaBins] = {0.1752 , 0.1862, 0.1411, 0.1534 , 0.1903, 0.2243, 0.2687}; //old
    const float effectiveAreaValues[nEtaBins] = {0.1440 , 0.1562, 0.1032, 0.0859 , 0.1116, 0.1321, 0.1654};  //updated

    reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
    float etaSC = el.superCluster()->eta();

    // Find eta bin first. If eta>2.5, the last eta bin is used.
    int etaBin = 0;
    while(etaBin < nEtaBins-1 && abs(etaSC) > etaBinLimits[etaBin+1])  ++etaBin;

    float area = effectiveAreaValues[etaBin];
    relIsoWithEA = (float)( pfIso.sumChargedHadronPt + std::max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - (*rhoH) * area ) )/el.pt();

    return relIsoWithEA;
}

//define this as a plug-in
DEFINE_FWK_MODULE(jpsiElec4l_KmcFitter);
