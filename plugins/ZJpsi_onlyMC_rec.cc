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
class Zjpsi_onlyMC_rec : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit Zjpsi_onlyMC_rec(const edm::ParameterSet&);
      ~Zjpsi_onlyMC_rec();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      void printMCtree(const reco::Candidate *, int);
      void printMCtreeUP(const reco::Candidate *, int);
      std::string printName(int);
      void analyzeDecay(const reco::Candidate*, TLorentzVector&, TLorentzVector&, TLorentzVector&,
      TLorentzVector&, TLorentzVector&, int);
      bool    isAncestor(const reco::Candidate*, const reco::Candidate*);
      bool    isAncestor(int, const reco::Candidate*);
      //edm::EDGetTokenT<reco::BeamSpot> BSLabel_;
      edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
      edm::EDGetTokenT<pat::PackedGenParticleCollection > packedGenToken_;

      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      TTree *Z_tree;
      UInt_t    run;
      ULong64_t event;

      TLorentzVector gen_z_p4, gen_dimuon_p4,gen_muon1_p4,gen_muon2_p4,gen_lepton1_p4,gen_lepton2_p4, gen_dilepton_p4;
      TLorentzVector gen_z_t, gen_dimuon_t, gen_muon1_t, gen_muon2_t, gen_lepton1_t, gen_lepton2_t, gen_dilepton_t;
      TLorentzVector gen_z_s, gen_dimuon_s, gen_muon1_s, gen_muon2_s, gen_lepton1_s, gen_lepton2_s, gen_dilepton_s;
      TVector3       gen_z_vtx, gen_jpsi_vtx;
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
Zjpsi_onlyMC_rec::Zjpsi_onlyMC_rec(const edm::ParameterSet& iConfig)
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

   Z_tree->Branch("run",      &run,      "run/i");
   Z_tree->Branch("event",    &event,    "event/i");

   Z_tree->Branch("gen_z_p4", "TLorentzVector", &gen_z_p4);
   Z_tree->Branch("gen_muon1_p4",  "TLorentzVector", &gen_muon1_p4);
   Z_tree->Branch("gen_muon2_p4",  "TLorentzVector", &gen_muon2_p4);
   Z_tree->Branch("gen_dimuon_p4", "TLorentzVector", &gen_dimuon_p4);
   Z_tree->Branch("gen_dilepton_p4", "TLorentzVector", &gen_dilepton_p4);
   Z_tree->Branch("gen_lepton1_p4",  "TLorentzVector", &gen_lepton1_p4);
   Z_tree->Branch("gen_lepton2_p4",  "TLorentzVector", &gen_lepton2_p4);

   Z_tree->Branch("gen_z_t", "TLorentzVector", &gen_z_t);
   Z_tree->Branch("gen_muon1_t",  "TLorentzVector", &gen_muon1_t);
   Z_tree->Branch("gen_muon2_t",  "TLorentzVector", &gen_muon2_t);
   Z_tree->Branch("gen_dimuon_t", "TLorentzVector", &gen_dimuon_t);
   Z_tree->Branch("gen_dilepton_t", "TLorentzVector", &gen_dilepton_t);
   Z_tree->Branch("gen_lepton1_t",  "TLorentzVector", &gen_lepton1_t);
   Z_tree->Branch("gen_lepton2_t",  "TLorentzVector", &gen_lepton2_t);
    
   Z_tree->Branch("gen_z_s", "TLorentzVector", &gen_z_s);
   Z_tree->Branch("gen_muon1_s",  "TLorentzVector", &gen_muon1_s);
   Z_tree->Branch("gen_muon2_s",  "TLorentzVector", &gen_muon2_s);
   Z_tree->Branch("gen_dimuon_s", "TLorentzVector", &gen_dimuon_s);
   Z_tree->Branch("gen_dilepton_s", "TLorentzVector", &gen_dilepton_s);
   Z_tree->Branch("gen_lepton1_s",  "TLorentzVector", &gen_lepton1_s);
   Z_tree->Branch("gen_lepton2_s",  "TLorentzVector", &gen_lepton2_s);
//} //end of NotOnlyGen
}//end of constructor 


Zjpsi_onlyMC_rec::~Zjpsi_onlyMC_rec()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Zjpsi_onlyMC_rec::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //MC
   edm::Handle<reco::GenParticleCollection> pruned;
   iEvent.getByToken(genCands_, pruned);
   // Packed particles are all the status 1, so usable to remake jets
   // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
   edm::Handle<pat::PackedGenParticleCollection> packed;
   iEvent.getByToken(packedGenToken_,packed);

    run       = iEvent.id().run();
    event     = iEvent.id().event();

    gen_z_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton2_p4.SetPtEtaPhiM(0.,0.,0.,0.);

    gen_dilepton_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);

    gen_z_t.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_dimuon_t.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon1_t.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon2_t.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton1_t.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton2_t.SetPtEtaPhiM(0.,0.,0.,0.);

    gen_dilepton_t.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_dimuon_t.SetPtEtaPhiM(0.,0.,0.,0.);
    
    gen_z_s.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_dimuon_s.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon1_s.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_muon2_s.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton1_s.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_lepton2_s.SetPtEtaPhiM(0.,0.,0.,0.);

    gen_dilepton_s.SetPtEtaPhiM(0.,0.,0.,0.);
    gen_dimuon_s.SetPtEtaPhiM(0.,0.,0.,0.);
    
    gen_z_vtx.SetXYZ(0.,0.,0.);
    gen_jpsi_vtx.SetXYZ(0.,0.,0.);
    int n_Z_dau = 0;
    int Event_Cand = 1;
    int tst = 0;
    int nm = 0;
    int foundit = 0;
    std::cout << "enters analyzer" << std::endl;
    if (pruned.isValid()){
        std::cout << "pruned valid "<< std::endl;
        for (size_t i=0; i<pruned->size(); i++) {
            const reco::Candidate *cand = &(*pruned)[i];
            if ((abs(cand->pdgId()) == 23)) {
                std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX " << std::endl;
                //printMCtree(cand, 0);
                std::cout << "print Tree " << std::endl;
                gen_z_t.SetPtEtaPhiM(cand->pt(), cand->eta(), cand->phi(), cand->mass());
                analyzeDecay(cand, gen_lepton1_t, gen_lepton2_t, gen_muon1_t, gen_muon2_t, gen_dimuon_t, 0);
                tst++;
                for (size_t k=0; k<packed->size(); k++) {
                   //const reco::Candidate * dauInPrunedColl = (*packed)[k].mother(0);
                   const reco::Candidate * stable_dau = &(*packed)[k];
                   int stable_id = (*packed)[k].pdgId();
                   if (stable_dau != nullptr && isAncestor(443,stable_dau) && isAncestor(23, stable_dau)) {
                      if(stable_id == 13) { //found muon-
                              gen_muon1_p4.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                      }
                      else if(stable_id == -13){ //found muon+
                           gen_muon2_p4.SetPtEtaPhiM(stable_dau->pt(),stable_dau->eta(),stable_dau->phi(),stable_dau->mass());
                      }
                   }
                   else if (stable_dau != nullptr && ! isAncestor(443, stable_dau) && isAncestor(23, stable_dau)){
                       if(stable_id == 13){
                           gen_lepton1_p4.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_id->phi(), stable_dau->mass());
                       }
                       else if(stable_id==-13){
                           gen_lepton2_p4.SetPtEtaPhiM(stable_dau->pt(), stable_dau->eta(), stable_dau->phi(), stable_dau->mass());
                       }
                   }
                }// end loop over stable
                gen_z_p4 = gen_z_t;
                std::cout<<" Z: "<< gen_z_p4.M()<<" | jpsi: "<<gen_dimuon_p4.Pt()<<" | m1: "<<gen_muon1_p4.Pt()<<" | m2: "<<gen_muon2_p4.Pt()<<" | l1: "<<gen_lepton1_p4.Pt()<<
                " | l2: "<<gen_lepton2_p4.Pt()<<std::endl;
                break;
                
            }
        }//end loop over pruned
    }//end if pruned
    
     //std::cout << "test" << std::endl;
    /*
     if ( pruned.isValid() ) {
       //std::cout << "MC ok " << std::endl;

       for (size_t i=0; i<pruned->size(); i++) {

          const reco::Candidate *dau = &(*pruned)[i];
          ///ndau = dau->numberOfDaughters();

          if ( (abs(dau->pdgId()) == 23) ) { //&& (dau->status() == 2) ) { //found Z
             //foundit++;
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
                 //foundit++;
                 gen_jpsi_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
                 gen_dimuon_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
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
             if (nm == 2 && foundit == 2){
                 Z_tree->Fill();
                 tst++;
                 std::cout<<"works ..."<<std::endl;
                 break;
             }
          } //endif found Z
       }//end for

    } //end pruned
    */
    if (tst) std::cout << "x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x NEW EVENT -x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x-x" << std::endl;


}//end analyze  
//Try to print mc Tree
std::string Zjpsi_onlyMC_rec::printName(int pdgid){
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
    umap[333] = "phi";
    umap[443] = "J/psi";
    std::string retstr;
    if (umap.find(pdgid) == umap.end()){
        retstr = std::to_string(pdgid);
    }
    else{
        retstr = umap.at(pdgid);
    }

    return retstr;
}

void Zjpsi_onlyMC_rec::printMCtree(const reco::Candidate* mother, int indent=0){
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

void Zjpsi_onlyMC_rec::printMCtreeUP(const reco::Candidate* daughter, int indent = 0){
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

void Zjpsi_onlyMC_rec::analyzeDecay(const reco::Candidate* mother, TLorentzVector& l1, TLorentzVector& l2, TLorentzVector& m1, TLorentzVector& m2, TLorentzVector& dim, int indent = 0){

    int momID = mother->pdgId();
    if (mother == NULL){
         //std::cout << "end tree" << std::endl;
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
        int dauID = daughter->pdgId();
        if (mother->numberOfDaughters() > 0){
            if(indent){
                std::cout << std::setw(indent) << " ";
            }
            std::cout << " daugter "<< i+1 <<": "<<  printName(daughter->pdgId()) << " with Pt: ";
            std::cout << daughter->pt() << " | Eta: "<< daughter->eta() << " | Phi: "<< daughter->phi() << " | mass: "<< daughter->mass() << std::endl;
            if(indent){
                std::cout << std::setw(indent) << " ";
            }
            std::cout<< " jpsi: "<<dim.Pt()<<" | l1: "<<l1.Pt()<< " | l2: "<<l2.Pt()<<" | m1: "<<m1.Pt()<<" | m2: "<<m2.Pt()<<std::endl;
            extraIndent+=4;
            if(dauID == 443 && dim.M() == 0) {
                dim.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass()); //jpsi
                std::cout << " jpsi " << std::endl;
                }
            else if(dauID == 13 && l1.M() == 0){
                l1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass()); // m1
                std::cout << "lep 1" << std::endl;
            }
            else if(dauID == -13 && l2.M() ==0){
                l2.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
                std::cout << "lep 2" << std::endl;
            }
            else if(dauID == 13 && m1.M() == 0){
                m1.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
                std::cout << "muon 1" << std::endl;
            }
            else if(dauID == -13 && m2.M() ==0){
                m2.SetPtEtaPhiM(daughter->pt(), daughter->eta(), daughter->phi(), daughter->mass());
                std::cout << "muon 2" << std::endl;
            }
        }
            if (daughter->numberOfDaughters()) analyzeDecay(daughter, l1, l2, m1, m2, dim, indent+extraIndent);
    }
}
//recursively check is a given particle is ancestor
bool Zjpsi_onlyMC_rec::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

//recursively check is a given particle has ancestor with given pdg_id
bool Zjpsi_onlyMC_rec::isAncestor(int a_pdgId, const reco::Candidate * particle) {
    if (a_pdgId == particle->pdgId() ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(a_pdgId,particle->mother(i))) return true;
    }
    return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Zjpsi_onlyMC_rec::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Zjpsi_onlyMC_rec::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Zjpsi_onlyMC_rec::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(Zjpsi_onlyMC_rec);
