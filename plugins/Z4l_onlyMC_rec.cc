// -*- C++ -*-
//
// Package:    AnalyzeZphill/Z4l_onlyMC_rec
// Class:      Z4l_onlyMC_rec
// 
/**\class Z4l_onlyMC_rec Z4l_onlyMC_rec.cc AnalyzeZphill/Z4l_onlyMC_rec/plugins/Z4l_onlyMC_rec.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Artemio Ayala Sanchez
//         Created:  Thu, 20 Jun 2019 19:12:04 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

//PAT includes
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

//Trigger includes
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <DataFormats/HepMCCandidate/interface/GenParticle.h>

#include <string>
#include <vector>
#include <unordered_map>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//this type of declaration is terrible for performance, we sould improve to edm::global::EDAnalizer with no shared Resources 
class Z4l_onlyMC_rec : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Z4l_onlyMC_rec(const edm::ParameterSet&);
      ~Z4l_onlyMC_rec();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      void printMCtree(const reco::Candidate *, int);
      void printMCtreeUP(const reco::Candidate *, int);
      std::string printName(int);
      bool    isAncestor(const reco::Candidate*, const reco::Candidate*);
      bool    isAncestor(int, const reco::Candidate*);
      //edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
      edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;

      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TTree *Z_tree;

      int mfromZ;
      int efromZ;
      int tfromZ;
      int zfromZ;

      int goodmuons;
      int goodelecs;

      TLorentzVector gen_z_p4;
      TLorentzVector gen_muon1N;
      TLorentzVector gen_muon1P;
      TLorentzVector gen_muon2N;
      TLorentzVector gen_muon2P;
      TLorentzVector gen_elec1N;
      TLorentzVector gen_elec1P;
      TLorentzVector gen_elec2N;
      TLorentzVector gen_elec2P;

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
Z4l_onlyMC_rec::Z4l_onlyMC_rec(const edm::ParameterSet& iConfig)
{
   genCands_ = consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"));
   packedGenToken_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");

   // If the analyzer does not use TFileService, please remove
   // the template argument to the base class so the class inherits
   // from  edm::one::EDAnalyzer<> and also remove the line from
   // constructor "usesResource("TFileService");"
   // This will improve performance in multithreaded jobs.
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   //Now we assing the member variables to the result TTree
   Z_tree = fs->make < TTree > ("ZTree", "Tree of Z4leptons");

   Z_tree->Branch("mfromZ", &mfromZ, "mfromZ/I");
   Z_tree->Branch("efromZ", &efromZ, "efromZ/I");
   Z_tree->Branch("tfromZ", &tfromZ, "tfromZ/I");
   Z_tree->Branch("zfromZ", &zfromZ, "zfromZ/I");

   Z_tree->Branch("goodmuons", &goodmuons, "goodmuons/I");
   Z_tree->Branch("goodelecs", &goodelecs, "goodelecs/I");

   Z_tree->Branch("gen_z_p4", "TLorentzVector", &gen_z_p4);
   Z_tree->Branch("gen_muon1N",      "TLorentzVector", &gen_muon1N);
   Z_tree->Branch("gen_muon1P",      "TLorentzVector", &gen_muon1P);
   Z_tree->Branch("gen_muon2N",      "TLorentzVector", &gen_muon2N);
   Z_tree->Branch("gen_muon2P",      "TLorentzVector", &gen_muon2P);
   Z_tree->Branch("gen_elec1N",      "TLorentzVector", &gen_elec1N);
   Z_tree->Branch("gen_elec1P",      "TLorentzVector", &gen_elec1P);
   Z_tree->Branch("gen_elec2N",      "TLorentzVector", &gen_elec2N);
   Z_tree->Branch("gen_elec2P",      "TLorentzVector", &gen_elec2P);

//} //end of NotOnlyGen
}//end of constructor 


Z4l_onlyMC_rec::~Z4l_onlyMC_rec()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Z4l_onlyMC_rec::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //MC
   edm::Handle<reco::GenParticleCollection> pruned;
   iEvent.getByToken(genCands_, pruned);
   // Packed particles are all the status 1, so usable to remake jets
   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
   edm::Handle<pat::PackedGenParticleCollection> packed;
   iEvent.getByToken(packedGenToken_,packed);

    gen_z_p4.SetPtEtaPhiM(0.,0.,0.,0.);

    gen_muon1N.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon1P.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon2N.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon2P.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_elec1N.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_elec1P.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_elec2N.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_elec2P.SetPtEtaPhiM(0.,0.,0.,0.);



    int n_Z_dau = 0;
    int Event_Cand = 1;
    //NEW MC ALV
    int tst = 0;
    int z_gen = 0;
    if (pruned.isValid()){ //if mc exist
    //std::cout<< "Valid pruned container" << std::endl;
        for(size_t i=0; i<pruned->size(); i++){ // loop over generated events
            const reco::Candidate *mom = &(*pruned)[i];
            if (std::isnan(mom->mass())) continue;
            if(abs(mom->pdgId()) == 23){ // if generated is Z boson
                TLorentzVector temp_mu1N, temp_mu1P, temp_mu2N, temp_mu2P; //define tempotals
                TLorentzVector temp_el1N, temp_el1P, temp_el2N, temp_el2P; //define tempotals
                temp_mu1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_mu2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el1P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el1N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el2P.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                temp_el2N.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                mfromZ = 0;
                efromZ = 0;
                tfromZ = 0;
                zfromZ = 0;
                goodmuons = 0;
                goodelecs = 0;
                z_gen++;
                if (mom->numberOfDaughters() == 1){
                    if (mom->daughter(0)->pdgId() == 23) zfromZ++;
                }
                if (zfromZ) continue;
                TLorentzVector suma_rancia;
                suma_rancia.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
                for(size_t k=0; k<packed->size(); k++){
                    const reco::Candidate * stable_dau = &(*packed)[k];
                    int stable_id = (*packed)[k].pdgId();
                    if (stable_dau != nullptr /*&& isAncestor(mom,stable_dau)*/) {
                        if (stable_id == 13 || stable_id == -13) { //electros 11 muons 13 as final states
                            mfromZ++;
                            TLorentzVector muon_pedorro;
                            muon_pedorro.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                            suma_rancia = suma_rancia + muon_pedorro;
                            int isAnc = isAncestor(23, stable_dau) ? 1 : 0;
                            std::cout << "comes from Z? " << isAnc <<std::cout;
                        }
                        else if (stable_id == 11 || stable_id == -11){
                            efromZ++;
                        }
                        else if (stable_id == 16 || stable_id == -16){
                            tfromZ++;
                        }
                    }
                }
                if(tfromZ) continue;
                if(mfromZ < 4) continue;
                std::cout << "taus from Z:      " << tfromZ << std::endl;
                std::cout << "muons from Z:     " << mfromZ << std::endl;
                std::cout << "electrons from Z: " << efromZ << std::endl;
                std::cout << "masa invariante en suma rancia de " << mfromZ << " muones es: " << suma_rancia.M() << std::endl;
                //if ((efromZ == 2 && mfromZ == 2) || mfromZ == 4) {
                //if (mfromZ > 2 && tfromZ == 0){
                if(tfromZ == 0){
                    //if (tst > 2) break;
                    std::cout<< "Found Z with mass: "<< mom->mass() <<", print Tree ... " << tst << std::endl;
                    printMCtree(mom, 0);

                    for(size_t k=0; k<packed->size(); k++){
                       const reco::Candidate * stable_dau = &(*packed)[k];
                       int stable_id = (*packed)[k].pdgId();
                       if (stable_dau != nullptr && isAncestor(mom,stable_dau)) {
                           //if (stable_id != 13 && stable_id != -13) continue;
                           if(stable_id == 13 && temp_mu1N.M() == 0){ // if muon- && not previiulsy assigned
                               // std::cout << " matched 1 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_mu1N.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodmuons++;
                            }
                            else if(stable_id == -13 && temp_mu1P.M() == 0){ // if muon+ && not previusly assigned
                                //std::cout << " matched 2 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_mu1P.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodmuons++;
                            }
                            else if(stable_id == 13 && temp_mu2N.M() == 0){ // if muon- && not previusly assigned
                                //std::cout << " matched 3 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_mu2N.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodmuons++;
                            }
                            else if(stable_id == -13 && temp_mu2P.M() == 0){ // if muon+ && not previusly assigned
                                //std::cout << " matched 4 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_mu2P.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodmuons++;
                            }
                            /////ELECTRONS
                            else if(stable_id == 11 && temp_el1N.M() == 0){ // if electron- && not previiulsy assigned
                                //std::cout << " matched 1 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_el1N.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodelecs++;
                            }
                            else if(stable_id == -11 && temp_el1P.M() == 0){ // if electron+ && not previusly assigned
                                //std::cout << " matched 2 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_el1P.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodelecs++;
                            }
                            else if(stable_id == 11 && temp_el2N.M() == 0){ // if electron- && not previusly assigned
                                //std::cout << " matched 3 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_el2N.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodelecs++;
                            }
                            else if(stable_id == -11 && temp_el2P.M() == 0){ // if electron+ && not previusly assigned
                                //std::cout << " matched 4 "  << " with Pt: ";
                                //std::cout << stable_dau->pt() << " | Eta: "<< stable_dau->eta() << " | Phi: "<< stable_dau->phi() << " | mass: "<< stable_dau->mass() << std::endl;
                                temp_el2P.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                                goodelecs++;
                            }
                       }
                    }
                    tst++;
                }

                if (goodelecs > 1 || goodmuons > 1){ //if leptons have been found
                    gen_z_p4.SetPtEtaPhiM(mom->pt(),mom->eta(),mom->phi(),mom->mass());
                    if(!std::isnan(gen_z_p4.M())){
                       if (gen_z_p4.M()!=0){
                          gen_muon1N = temp_mu1N;
                          gen_muon1P = temp_mu1P;
                          gen_muon2N = temp_mu2N;
                          gen_muon2P = temp_mu2P;

                          gen_elec1N = temp_el1N;
                          gen_elec1P = temp_el1P;
                          gen_elec2N = temp_el2N;
                          gen_elec2P = temp_el2P;

                          Z_tree->Fill();
                          //std::cout << " good gen muons: "<< goodmuons << std::endl;
                          //std::cout << " good gen elecs: "<< goodelecs << std::endl;
                       }
                    }
                }// end if generated is good event
            }//end if generated is Z
        }//end loop over pruned
    }//end if pruned
    if (tst) std::cout << "x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x NEW EVENT -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x" << std::endl;


}//end analyze  
//Try to print mc Tree
std::string Z4l_onlyMC_rec::printName(int pdgid){
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

void Z4l_onlyMC_rec::printMCtree(const reco::Candidate* mother, int indent=0){
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

void Z4l_onlyMC_rec::printMCtreeUP(const reco::Candidate* daughter, int indent = 0){
    if (daughter == NULL){
        std::cout << "end tree" << std::endl;
        return;
    }
    int extraIndent = 0;
    for(size_t i = 0; i < daughter->numberOfMothers(); i++){
        const reco::Candidate* mother = daughter->mother(i);
        if(indent){
            std::cout << std::setw(indent) << " ";
        }
        std::cout<< "mother "<< i+1 << ": "<<printName(mother->pdgId()) << std::endl;
        extraIndent+=4;
        if (mother->numberOfMothers()) printMCtreeUP(daughter, indent+extraIndent);
    }
}

//recursively check is a given particle is ancestor
bool Z4l_onlyMC_rec::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool Z4l_onlyMC_rec::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Z4l_onlyMC_rec::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Z4l_onlyMC_rec::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Z4l_onlyMC_rec::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(Z4l_onlyMC_rec);
