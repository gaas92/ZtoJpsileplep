// -*- C++ -*-
//
// Package:    AnalyzeZphill/ZjpsiMCTupler
// Class:      ZjpsiMCTupler
// 
/**\class ZjpsiMCTupler ZjpsiMCTupler.cc AnalyzeZphill/ZjpsiMCTupler/plugins/ZjpsiMCTupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Artemio Ayala Sanchez
//         Created:  Thu, 20 Jun 2019 19:12:04 GMT
//         Updated:  v7-12/09/20
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

//Trigger includes
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <string>
#include <vector>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

//this type of declaration is terrible for performance, we sould improve to edm::global::EDAnalizer with no shared Resources 
class ZjpsiMCTupler : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ZjpsiMCTupler(const edm::ParameterSet&);
      ~ZjpsiMCTupler();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void CheckHLTTriggers(const std::vector<std::string>& TrigList);

      // ----------member data ---------------------------
      // to save all triggers in string 
      std::string hlTriggerResults_;
      char triggersL[10000];
      int  nTrgL;
      //ZCand from phiLepLepKFitter
      edm::EDGetTokenT<pat::CompositeCandidateCollection> theZ_;
      edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
      //To perform the primary vertex check, i think is an overkill
      //edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
      //TTree Data
      UInt_t    run;
      ULong64_t event;
      UInt_t    lumiblock;
      UInt_t    nonia;
      UInt_t    nmuons;
      UInt_t    vertex_i;
      //New
      UInt_t    nPV;
      UInt_t    nCands;
    
      TLorentzVector dimuon_p4;
      TLorentzVector muonN_p4;
      TLorentzVector muonP_p4;
      TLorentzVector dilepton_p4;
      TLorentzVector lepton1_p4;
      TLorentzVector lepton2_p4;
      TLorentzVector Z_p4;
    
      //new
      TLorentzVector leading_lepton_p4;
      TLorentzVector trailing_lepton_p4;
      TLorentzVector leading_muon_p4;
      TLorentzVector trailing_muon_p4;
    
      TLorentzVector gen_z_p4;
      TLorentzVector gen_muonN_p4;
      TLorentzVector gen_muonP_p4;
      TLorentzVector gen_dilepton_p4;
      TLorentzVector gen_dimuon_p4;
      TLorentzVector gen_lepton1_p4;
      TLorentzVector gen_lepton2_p4;
      //new
      TLorentzVector gen_leading_lepton_p4;
      TLorentzVector gen_trailing_lepton_p4;
      TLorentzVector gen_leading_muon_p4;
      TLorentzVector gen_trailing_muon_p4;
    
      TVector3 gen_Zvtx, gen_psi_vtx;
    
    
      TLorentzVector msrd_dimuon_p4;
      TLorentzVector msrd_muonN_p4;
      TLorentzVector msrd_muonP_p4;
      TLorentzVector msrd_dilepton_p4;
      TLorentzVector msrd_lepton1_p4;
      TLorentzVector msrd_lepton2_p4;
      TLorentzVector msrd_Z_p4;

      //new
      UInt_t passFit;
      TVector3 Zvtx;
      Float_t ZvtxP;
      Float_t ZvtxC2;
      Float_t pvChi2;
      //int tst;
      //int decaychannel;
      //int pass_match;
      //int z_gen;

      Float_t psiVtxP;
      Float_t psiVtxC2;

      Float_t dxym1;
      Float_t dxym2;
      Float_t dzm1;
      Float_t dzm2;
      Float_t dxyl1;
      Float_t dxyl2;
      Float_t dzl1;
      Float_t dzl2;

      Float_t dRiso_l1;
      Float_t dRiso_l2;
      //new
      Float_t dRiso_m1;
      Float_t dRiso_m2;
      //new
      Float_t rIsoOverPtl1;
      Float_t rIsoOverPtl2;
      Float_t rIsoOverPtm1;
      Float_t rIsoOverPtm2;
    
      Float_t dR_m1_m2;
      Float_t dR_l1_l2;
      Float_t dR_m1_l1;
      Float_t dR_m1_l2;
      Float_t dR_m2_l1;
      Float_t dR_m2_l2;
    
      Float_t dR_m1_m2_;
      Float_t dR_l1_l2_;
      Float_t dR_m1_l1_;
      Float_t dR_m1_l2_;
      Float_t dR_m2_l1_;
      Float_t dR_m2_l2_;
    
      //new
      Float_t ipSm1;
      Float_t ipSm2;
      Float_t ipSl1;
      Float_t ipSl2;
    
      Float_t dipm1;
      Float_t dipm2;
      Float_t dipl1;
      Float_t dipl2;
      Float_t dipm1Err;
      Float_t dipm2Err;
      Float_t dipl1Err;
      Float_t dipl2Err;
      //new
      Float_t ImparSigl1;
      Float_t ImparSigl2;
      Float_t ImparSigm1;
      Float_t ImparSigm2;
      //no phiVprov

      TTree *Z_tree;
      //For MC
      /*
      Int_t mother_pdgId;
      Int_t dimuon_pdgId;
      TLorentzVector gen_dimuon_p4;
      TLorentzVector gen_muonP_p4;
      TLorentzVector gen_muonM_p4;
      */
      int ZLe1Qid;
      int ZLe2Qid;
    
      int JpM1Qid;
      int JpM2Qid;

      int psiM1_TrackerLWM;
      int psiM1_PixelLWM;
      int psiM1_ValPixHit;

      int psiM2_TrackerLWM;
      int psiM2_PixelLWM;
      int psiM2_ValPixHit;

      int ZLe1_TrackerLWM;
      int ZLe1_PixelLWM;
      int ZLe1_ValPixHit;

      int ZLe2_TrackerLWM;
      int ZLe2_PixelLWM;
      int ZLe2_ValPixHit;

      int ZLe1_pass;
      int ZLe2_pass;
    
      int Event_Cand;
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
ZjpsiMCTupler::ZjpsiMCTupler(const edm::ParameterSet& iConfig):
        trigResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults")))

{
   theZ_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("Zcand"));
   hlTriggerResults_ = iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT"));
   //Personaly i think is an overkill perform again a primary vertex check, allready done in LeptonFilter and in phiZKFitter
   //primaryVertices_Label = consumes<reco::VertexCollection>(iConfig.getParameter< edm::InputTag>("primaryVertices"));

   // If the analyzer does not use TFileService, please remove
   // the template argument to the base class so the class inherits
   // from  edm::one::EDAnalyzer<> and also remove the line from
   // constructor "usesResource("TFileService");"
   // This will improve performance in multithreaded jobs.
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   //Now we assing the member variables to the result TTree
   Z_tree = fs->make < TTree > ("ZTree", "Tree of ZphiMuMu");

   Z_tree->Branch("run",      &run,      "run/i");
   Z_tree->Branch("event",    &event,    "event/l");
   Z_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");
   Z_tree->Branch("vertex_i", &vertex_i, "vertex_i/i");

 //  if (!OnlyGen_) {

    //new
    Z_tree->Branch("nonia",    &nonia,    "nonia/i");
    Z_tree->Branch("nmuons",   &nmuons,   "nmuons/i");
    Z_tree->Branch("nPV",      &nPV,      "nPV/i");
    Z_tree->Branch("passFit",  &passFit,  "passFit/i");
    Z_tree->Branch("nPV",      &nPV,      "nPV/i");
    Z_tree->Branch("nCands",   &nCands,   "nCands/i");

    Z_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    Z_tree->Branch("muonN_p4",  "TLorentzVector", &muonN_p4);
    Z_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
    //new
    Z_tree->Branch("dilepton_p4", "TLorentzVector", &dilepton_p4);
    Z_tree->Branch("lepton1_p4",  "TLorentzVector", &lepton1_p4);
    Z_tree->Branch("lepton2_p4",  "TLorentzVector", &lepton2_p4);
    
    //new
    Z_tree->Branch("leading_lepton_p4",  "TLorentzVector", &leading_lepton_p4);
    Z_tree->Branch("trailing_lepton_p4",  "TLorentzVector", &trailing_lepton_p4);
    Z_tree->Branch("leading_muon_p4",  "TLorentzVector", &leading_muon_p4);
    Z_tree->Branch("trailing_muon_p4",  "TLorentzVector", &trailing_muon_p4);
    
    Z_tree->Branch("Z_p4", "TLorentzVector", &Z_p4);

    Z_tree->Branch("gen_z_p4", "TLorentzVector", &gen_z_p4);
    Z_tree->Branch("gen_muonN_p4",  "TLorentzVector", &gen_muonN_p4);
    Z_tree->Branch("gen_muonP_p4",  "TLorentzVector", &gen_muonP_p4);
    Z_tree->Branch("gen_dimuon_p4", "TLorentzVector", &gen_dimuon_p4);
    Z_tree->Branch("gen_dilepton_p4", "TLorentzVector", &gen_dilepton_p4);
    Z_tree->Branch("gen_lepton1_p4",  "TLorentzVector", &gen_lepton1_p4);
    Z_tree->Branch("gen_lepton2_p4",  "TLorentzVector", &gen_lepton2_p4);
    //new
    Z_tree->Branch("gen_leading_lepton_p4",   "TLorentzVector", &gen_leading_lepton_p4);
    Z_tree->Branch("gen_trailing_lepton_p4",  "TLorentzVector", &gen_trailing_lepton_p4);
    Z_tree->Branch("gen_leading_muon_p4",     "TLorentzVector", &gen_leading_muon_p4);
    Z_tree->Branch("gen_trailing_muon_p4",    "TLorentzVector", &gen_trailing_muon_p4);
    Z_tree->Branch("gen_Zvtx", "TVector3", &gen_Zvtx);
    Z_tree->Branch("gen_psi_vtx", "TVector3", &gen_psi_vtx);
    
    Z_tree->Branch("msrd_dimuon_p4", "TLorentzVector", &msrd_dimuon_p4);
    Z_tree->Branch("msrd_muonN_p4",  "TLorentzVector", &msrd_muonN_p4);
    Z_tree->Branch("msrd_muonP_p4",  "TLorentzVector", &msrd_muonP_p4);
    Z_tree->Branch("msrd_lepton1_p4",  "TLorentzVector", &msrd_lepton1_p4);
    Z_tree->Branch("msrd_lepton2_p4",  "TLorentzVector", &msrd_lepton2_p4);
    Z_tree->Branch("msrd_Z_p4", "TLorentzVector", &msrd_Z_p4);

    Z_tree->Branch("pvChi2", &pvChi2, "pvChi2/F");
    Z_tree->Branch("Zvtx", "TVector3", &Zvtx);
    Z_tree->Branch("ZvtxP", &ZvtxP, "ZvtxP/F");
    Z_tree->Branch("ZvtxC2", &ZvtxC2, "ZvtxC2/F");
    //Z_tree->Branch("tst", &tst, "tst/i");
    //Z_tree->Branch("decaychannel", &decaychannel, "decaychannel/i");
    //Z_tree->Branch("pass_match", &pass_match, "pass_match/i");
    //Z_tree->Branch("z_gen", &z_gen, "z_gen/i");

    Z_tree->Branch("psiVtxP", &psiVtxP, "psiVtxP/F");
    Z_tree->Branch("psiVtxC2", &psiVtxC2, "psiVtxC2/F");

    Z_tree->Branch("dxym1", &dxym1, "dxym1/F");
    Z_tree->Branch("dxym2", &dxym2, "dxym2/F");
    Z_tree->Branch("dzm1", &dzm1, "dzm1/F");
    Z_tree->Branch("dzm2", &dzm2, "dzm2/F");

    Z_tree->Branch("dxyl1", &dxyl1, "dxyl1/F");
    Z_tree->Branch("dxyl2", &dxyl2, "dxyl2/F");
    Z_tree->Branch("dzl1", &dzl1, "dzl1/F");
    Z_tree->Branch("dzl2", &dzl2, "dzl2/F");

    Z_tree->Branch("dRiso_l1", &dRiso_l1, "dRiso_l1/F");
    Z_tree->Branch("dRiso_l2", &dRiso_l2, "dRiso_l2/F");
    //new
    Z_tree->Branch("dRiso_m1", &dRiso_m1, "dRiso_m1/F");
    Z_tree->Branch("dRiso_m2", &dRiso_m2, "dRiso_m2/F");
    //new
    Z_tree->Branch("rIsoOverPtl1", &rIsoOverPtl1, "rIsoOverPtl1/F");
    Z_tree->Branch("rIsoOverPtl2", &rIsoOverPtl2, "rIsoOverPtl2/F");
    Z_tree->Branch("rIsoOverPtm1", &rIsoOverPtm1, "rIsoOverPtm1/F");
    Z_tree->Branch("rIsoOverPtm2", &rIsoOverPtm2, "rIsoOverPtm2/F");
    
    Z_tree->Branch("dR_m1_m2_", &dR_m1_m2_, "dR_m1_m2_/F");
    Z_tree->Branch("dR_l1_l2_", &dR_l1_l2_, "dR_l1_l2_/F");
    Z_tree->Branch("dR_m1_l1_", &dR_m1_l1_, "dR_m1_l1_/F");
    Z_tree->Branch("dR_m1_l2_", &dR_m1_l2_, "dR_m1_l2_/F");
    Z_tree->Branch("dR_m2_l1_", &dR_m2_l1_, "dR_m2_l1_/F");
    Z_tree->Branch("dR_m2_l2_", &dR_m2_l2_, "dR_m2_l2_/F");

    Z_tree->Branch("dR_m1_m2", &dR_m1_m2, "dR_m1_m2/F");
    Z_tree->Branch("dR_l1_l2", &dR_l1_l2, "dR_l1_l2/F");
    Z_tree->Branch("dR_m1_l1", &dR_m1_l1, "dR_m1_l1/F");
    Z_tree->Branch("dR_m1_l2", &dR_m1_l2, "dR_m1_l2/F");
    Z_tree->Branch("dR_m2_l1", &dR_m2_l1, "dR_m2_l1/F");
    Z_tree->Branch("dR_m2_l2", &dR_m2_l2, "dR_m2_l2/F");

    Z_tree->Branch("dipm1",     &dipm1,"dipm1/F" );
    Z_tree->Branch("dipm2",     &dipm2,"dipm2/F");
    Z_tree->Branch("dipl1",     &dipl1,"dipl1/F" );
    Z_tree->Branch("dipl2",     &dipl2,"dipl2/F");
    Z_tree->Branch("dipm1Err",  &dipm1Err,"dipm1Err/F") ;
    Z_tree->Branch("dipm2Err",  &dipm2Err,"dipm2Err/F") ;
    Z_tree->Branch("dipl1Err",  &dipl1Err,"dipl1Err/F") ;
    Z_tree->Branch("dipl2Err",  &dipl2Err,"dipl2Err/F") ;
    //new
    Z_tree->Branch("ipSm1",    &ipSm1,"ipSm1/F") ;
    Z_tree->Branch("ipSm2",    &ipSm2,"ipSm2/F") ;
    Z_tree->Branch("ipSl1",    &ipSl1,"ipSl1/F") ;
    Z_tree->Branch("ipSl2",    &ipSl2,"ipSl2/F") ;
    Z_tree->Branch("ImparSigl1",    &ImparSigl1,"ImparSigl1/F") ;
    Z_tree->Branch("ImparSigl2",    &ImparSigl2,"ImparSigl2/F") ;
    Z_tree->Branch("ImparSigm1",    &ImparSigm1,"ImparSigm1/F") ;
    Z_tree->Branch("ImparSigm2",    &ImparSigm2,"ImparSigm2/F") ;
    


    Z_tree->Branch("nTrgL",      &nTrgL,    "nTrgL/I");
    Z_tree->Branch("triggersL",         &triggersL,  "triggersL[nTrgL]/C");

    Z_tree->Branch("ZLe1Qid",         &ZLe1Qid,  "ZLe1Qid/i");
    Z_tree->Branch("ZLe2Qid",         &ZLe2Qid,  "ZLe2Qid/i");
    Z_tree->Branch("JpM1Qid",         &JpM1Qid,  "JpM1Qid/i");
    Z_tree->Branch("JpM2Qid",         &JpM2Qid,  "JpM2Qid/i");

    Z_tree->Branch("psiM1_TrackerLWM",         &psiM1_TrackerLWM,  "psiM1_TrackerLWM/i");
    Z_tree->Branch("psiM1_PixelLWM",           &psiM1_PixelLWM,    "psiM1_PixelLWM/i");
    Z_tree->Branch("psiM1_ValPixHit",          &psiM1_ValPixHit,   "psiM1_ValPixHit/i");
    Z_tree->Branch("psiM2_TrackerLWM",         &psiM2_TrackerLWM,  "psiM2_TrackerLWM/i");
    Z_tree->Branch("psiM2_PixelLWM",           &psiM2_PixelLWM,    "psiM2_PixelLWM/i");
    Z_tree->Branch("psiM2_ValPixHit",          &psiM2_ValPixHit,   "psiM2_ValPixHit/i");

    Z_tree->Branch("ZLe1_TrackerLWM",         &ZLe1_TrackerLWM,  "ZLe1_TrackerLWM/i");
    Z_tree->Branch("ZLe1_PixelLWM",           &ZLe1_PixelLWM,    "ZLe1_PixelLWM/i");
    Z_tree->Branch("ZLe1_ValPixHit",          &ZLe1_ValPixHit,   "ZLe1_ValPixHit/i");
    Z_tree->Branch("ZLe2_TrackerLWM",         &ZLe2_TrackerLWM,  "ZLe2_TrackerLWM/i");
    Z_tree->Branch("ZLe2_PixelLWM",           &ZLe2_PixelLWM,    "ZLe2_PixelLWM/i");
    Z_tree->Branch("ZLe2_ValPixHit",          &ZLe2_ValPixHit,   "ZLe2_ValPixHit/i");

    Z_tree->Branch("ZLe1_pass",          &ZLe1_pass,   "ZLe1_pass/i");
    Z_tree->Branch("ZLe2_pass",          &ZLe2_pass,   "ZLe2_pass/i");
    
    Z_tree->Branch("Event_Cand", &Event_Cand, "Event_Cand_/i");
//} //end of NotOnlyGen
}//end of constructor 


ZjpsiMCTupler::~ZjpsiMCTupler()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZjpsiMCTupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<pat::CompositeCandidateCollection> ZCands;
   iEvent.getByToken(theZ_,ZCands);
   //i think is an overkill
   /* 
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(primaryVertices_Label, vertices);
   numPrimaryVertices = vertices->size();
   reco::VertexCollection::const_iterator PV = vertices->end();
   for (reco::VertexCollection::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx, ++PV) {
     if (!(vtx->isFake()) && vtx->ndof() >= 4. && vtx->position().Rho() < 2.0 && fabs(vtx->position().Z()) < 24.0){
       PV = vtx;
       break;
       }
    }
   if ( PV==vertices->end() ) return;
   */

   pat::CompositeCandidate z_Cand;
   run       = iEvent.id().run();
   event     = iEvent.id().event();
   lumiblock = iEvent.id().luminosityBlock();
   ///////////////////////
   //Jhovanny's stuff
   // ///////////////////

   edm::Handle<edm::TriggerResults> hltresults1;
   try {
     std::string const &trig = std::string("TriggerResults::HLT");//+hlTriggerResults_;
     iEvent.getByLabel(edm::InputTag(trig),hltresults1);
     std::cout << "Handle on HLT ok " << std::endl;
   }
   catch ( ... )
    {
      std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
    }

    std::vector<std::string> TrigTable; TrigTable.clear();
    // Get hold of trigger names - based on TriggerResults object
    const edm::TriggerNames &triggerNames1_ = iEvent.triggerNames(*hltresults1);

    for (unsigned int itrig = 0; itrig < hltresults1->size(); ++itrig){
        if ((*hltresults1)[itrig].accept() == 1){
            std::string trigName1 = triggerNames1_.triggerName(itrig);
            TrigTable.push_back(trigName1);
        }
    }
    CheckHLTTriggers(TrigTable);
    //End of jhovannys (for all triggers in string)
    //if ( ! OnlyGen_ ) { // we will look for dimuons, then for muons
       if (ZCands.isValid() && !ZCands->empty()) {
         std::cout<<"XXXXXXXXXXXXXXXX-------TUPLER-------XXXXXXXXXXXXXXXX"<<std::endl;
         unsigned int csize = ZCands->size();
         //if (bestCandidateOnly_) csize = 1; //not implemented here
         for ( unsigned int i = 0; i < csize; i++ ) {
           z_Cand = (ZCands->at(i));
             
           //new
           nonia  = z_Cand.userInt("nonia_");
           nmuons = z_Cand.userInt("nmuons_");
           nPV    = z_Cand.userInt("nPV_");
           passFit = z_Cand.userInt("passFit_");
           vertex_i = z_Cand.userInt("pvIndex");
           nCands  = csize;
             
           Z_p4.SetPtEtaPhiM(z_Cand.pt(), z_Cand.eta(), z_Cand.phi(), z_Cand.mass());
           dimuon_p4.SetPtEtaPhiM(z_Cand.daughter("psi")->pt(), z_Cand.daughter("psi")->eta(),
                                  z_Cand.daughter("psi")->phi(), z_Cand.daughter("psi")->mass());

           lepton1_p4.SetPtEtaPhiM(z_Cand.daughter("lepton1")->pt(), z_Cand.daughter("lepton1")->eta(),
                      z_Cand.daughter("lepton1")->phi(), z_Cand.daughter("lepton1")->mass());
           lepton2_p4.SetPtEtaPhiM(z_Cand.daughter("lepton2")->pt(), z_Cand.daughter("lepton2")->eta(),
                      z_Cand.daughter("lepton2")->phi(), z_Cand.daughter("lepton2")->mass());
           //new
           reco::Candidate::LorentzVector dilep;
           dilep = z_Cand.daughter("lepton1")->p4() + z_Cand.daughter("lepton2")->p4();
           dilepton_p4.SetPtEtaPhiM(dilep.pt(), dilep.eta(), dilep.phi(), dilep.mass());
             
             reco::Candidate::LorentzVector vP ;//= z_Cand.daughter("psi")->daughter("muon1")->p4();
             reco::Candidate::LorentzVector vM ;//= z_Cand.daughter("psi")->daughter("muon2")->p4();

           if (z_Cand.daughter("psi")->daughter("muon1")->charge() < 0) {
              vP = z_Cand.daughter("psi")->daughter("muon2")->p4();
              vM = z_Cand.daughter("psi")->daughter("muon1")->p4();
           }
           else{
              vP = z_Cand.daughter("psi")->daughter("muon1")->p4();
              vM = z_Cand.daughter("psi")->daughter("muon2")->p4();
           }
           muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
           muonN_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
           //new
           reco::Candidate::LorentzVector lm;
           reco::Candidate::LorentzVector tm;
           reco::Candidate::LorentzVector ll;
           reco::Candidate::LorentzVector tl;
           if (z_Cand.daughter("lepton1")->pt() > z_Cand.daughter("lepton2")->pt()){
               ll = z_Cand.daughter("lepton1")->p4();
               tl = z_Cand.daughter("lepton2")->p4();
               }
           else{
               ll = z_Cand.daughter("lepton2")->p4();
               tl = z_Cand.daughter("lepton1")->p4();
           }
           if (z_Cand.daughter("psi")->daughter("muon1")->pt() > z_Cand.daughter("psi")->daughter("muon2")->pt()){
               lm = z_Cand.daughter("psi")->daughter("muon1")->p4();
               tm = z_Cand.daughter("psi")->daughter("muon2")->p4();
               }
           else{
               lm = z_Cand.daughter("psi")->daughter("muon2")->p4();
               tm = z_Cand.daughter("psi")->daughter("muon1")->p4();
           }
           leading_lepton_p4.SetPtEtaPhiM( ll.pt(), ll.eta(), ll.phi(), ll.mass());
           trailing_lepton_p4.SetPtEtaPhiM(tl.pt(), tl.eta(), tl.phi(), tl.mass());
           leading_muon_p4.SetPtEtaPhiM(   lm.pt(), lm.eta(), lm.phi(), lm.mass());
           trailing_muon_p4.SetPtEtaPhiM(  tm.pt(), tm.eta(), tm.phi(), tm.mass());
             
             
           //Measured values
           msrd_Z_p4.SetPtEtaPhiM(z_Cand.daughter("patMZ")->pt(), z_Cand.daughter("patMZ")->eta(), 
                     z_Cand.daughter("patMZ")->phi(), z_Cand.daughter("patMZ")->mass());
           msrd_dimuon_p4.SetPtEtaPhiM(z_Cand.daughter("msrd_psi")->pt(), z_Cand.daughter("msrd_psi")->eta(),
                          z_Cand.daughter("msrd_psi")->phi(), z_Cand.daughter("msrd_psi")->mass());
           msrd_lepton1_p4.SetPtEtaPhiM(z_Cand.daughter("msrd_lepton1")->pt(), z_Cand.daughter("msrd_lepton1")->eta(), 
                           z_Cand.daughter("msrd_lepton1")->phi(), z_Cand.daughter("msrd_lepton1")->mass());
           msrd_lepton2_p4.SetPtEtaPhiM(z_Cand.daughter("msrd_lepton2")->pt(), z_Cand.daughter("msrd_lepton2")->eta(),
                           z_Cand.daughter("msrd_lepton2")->phi(), z_Cand.daughter("msrd_lepton2")->mass());
           reco::Candidate::LorentzVector msrd_vP ; // = z_Cand.daughter("msrd_psi")->daughter("msrd_kaon1")->p4();
           reco::Candidate::LorentzVector msrd_vM ; // = z_Cand.daughter("msrd_psi")->daughter("msrd_kaon2")->p4();

           if (z_Cand.daughter("msrd_psi")->daughter("msrd_muon1")->charge() < 0) {
              msrd_vP = z_Cand.daughter("msrd_psi")->daughter("msrd_muon2")->p4();
              msrd_vM = z_Cand.daughter("msrd_psi")->daughter("msrd_muon1")->p4();
           }
           else{
              msrd_vP = z_Cand.daughter("msrd_psi")->daughter("msrd_muon1")->p4();
              msrd_vM = z_Cand.daughter("msrd_psi")->daughter("msrd_muon2")->p4();
           }
           msrd_muonP_p4.SetPtEtaPhiM(msrd_vP.pt(), msrd_vP.eta(), msrd_vP.phi(), msrd_vP.mass());
           msrd_muonN_p4.SetPtEtaPhiM(msrd_vM.pt(), msrd_vM.eta(), msrd_vM.phi(), msrd_vM.mass());
           //new
           reco::Candidate::LorentzVector msrd_dilep;
           msrd_dilep = z_Cand.daughter("msrd_lepton1")->p4() + z_Cand.daughter("msrd_lepton2")->p4();
           msrd_dilepton_p4.SetPtEtaPhiM(msrd_dilep.pt(), msrd_dilep.eta(), msrd_dilep.phi(), msrd_dilep.mass());
             


           gen_z_p4.SetPtEtaPhiM(z_Cand.daughter("mcZ")->pt(), z_Cand.daughter("mcZ")->eta(), z_Cand.daughter("mcZ")->phi(), z_Cand.daughter("mcZ")->mass());
           gen_muonP_p4.SetPtEtaPhiM(z_Cand.daughter("mcM2")->pt(), z_Cand.daughter("mcM2")->eta(), z_Cand.daughter("mcM2")->phi(), z_Cand.daughter("mcM2")->mass());
           gen_muonN_p4.SetPtEtaPhiM(z_Cand.daughter("mcM1")->pt(), z_Cand.daughter("mcM1")->eta(), z_Cand.daughter("mcM1")->phi(), z_Cand.daughter("mcM1")->mass());
                //std::cout<< "MK1 = "<<z_Cand.daughter("mcK1")->pt()<<" "<<z_Cand.daughter("mcK1")->eta()<<" "<<z_Cand.daughter("mcK1")->phi()<<" "<< z_Cand.daughter("mcK1")->mass() << std::endl;
           gen_dimuon_p4.SetPtEtaPhiM(z_Cand.daughter("mcPsi")->pt(), z_Cand.daughter("mcPsi")->eta(),
                                            z_Cand.daughter("mcPsi")->phi(), z_Cand.daughter("mcPsi")->mass());
           gen_lepton1_p4.SetPtEtaPhiM(z_Cand.daughter("mcL1")->pt(), z_Cand.daughter("mcL1")->eta(),
                         z_Cand.daughter("mcL1")->phi(), z_Cand.daughter("mcL1")->mass());
           gen_lepton2_p4.SetPtEtaPhiM(z_Cand.daughter("mcL2")->pt(), z_Cand.daughter("mcL2")->eta(),
                         z_Cand.daughter("mcL2")->phi(), z_Cand.daughter("mcL2")->mass());
           gen_Zvtx.SetXYZ(z_Cand.daughter("mcZ")->vx() ,z_Cand.daughter("mcZ")->vy(),z_Cand.daughter("mcZ")->vz()) ;
           gen_psi_vtx.SetXYZ(z_Cand.daughter("mcPsi")->vx() ,z_Cand.daughter("mcPsi")->vy(),z_Cand.daughter("mcPsi")->vz()) ;
             
           //new
           reco::Candidate::LorentzVector gen_dilep;
           gen_dilep = z_Cand.daughter("mcL1")->p4() + z_Cand.daughter("mcL2")->p4();
           gen_dilepton_p4.SetPtEtaPhiM(gen_dilep.pt(), gen_dilep.eta(), gen_dilep.phi(), gen_dilep.mass());
             //new
             reco::Candidate::LorentzVector gen_lm;
             reco::Candidate::LorentzVector gen_tm;
             reco::Candidate::LorentzVector gen_ll;
             reco::Candidate::LorentzVector gen_tl;
             if (z_Cand.daughter("mcL1")->pt() > z_Cand.daughter("mcL2")->pt()){
                 gen_ll = z_Cand.daughter("mcL1")->p4();
                 gen_tl = z_Cand.daughter("mcL2")->p4();
                }
             else{
                 gen_ll = z_Cand.daughter("mcL2")->p4();
                 gen_tl = z_Cand.daughter("mcL1")->p4();
             }
             if (z_Cand.daughter("mcM1")->pt() > z_Cand.daughter("mcM2")->pt()){
                 gen_lm = z_Cand.daughter("mcM1")->p4();
                 gen_tm = z_Cand.daughter("mcM2")->p4();
                }
             else{
                 gen_lm = z_Cand.daughter("mcM2")->p4();
                 gen_tm = z_Cand.daughter("mcM1")->p4();
             }
             gen_leading_lepton_p4.SetPtEtaPhiM( gen_ll.pt(), gen_ll.eta(), gen_ll.phi(), gen_ll.mass());
             gen_trailing_lepton_p4.SetPtEtaPhiM(gen_tl.pt(), gen_tl.eta(), gen_tl.phi(), gen_tl.mass());
             gen_leading_muon_p4.SetPtEtaPhiM(   gen_lm.pt(), gen_lm.eta(), gen_lm.phi(), gen_lm.mass());
             gen_trailing_muon_p4.SetPtEtaPhiM(  gen_tm.pt(), gen_tm.eta(), gen_tm.phi(), gen_tm.mass());
             
           pvChi2 = z_Cand.userFloat("pvChi2_");
           ZvtxP = z_Cand.userFloat("vProb") ;
           ZvtxC2 = z_Cand.userFloat("vChi2") ;
           //tst = z_Cand.userInt("tst_");
           //decaychannel = z_Cand.userInt("decay_");
           //pass_match = z_Cand.userInt("pass_match_");
             //std::cout << "testing shit -z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z-z in tupler" << std::endl;
             //std::cout << "decay channel is : " << decaychannel << std::endl;
             //std::cout << "pass match is : " << pass_match << std::endl;
           //z_gen = z_Cand.userInt("z_gen_");
           Zvtx.SetXYZ(z_Cand.userFloat("ZvtxX") ,z_Cand.userFloat("ZvtxY"),z_Cand.userFloat("ZvtxZ")) ;
           psiVtxP = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")))->userFloat("vProb");
           psiVtxC2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")))->userFloat("vChi2");

           ZLe1Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userInt("ZLe1Qid_");
           ZLe2Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userInt("ZLe2Qid_");
        
           JpM1Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userInt("ZMu1Qid_");
           JpM2Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userInt("ZMu2Qid_");

           ZLe1_pass = 0;//(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("pass");
           ZLe2_pass = 0;//(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("pass");

           psiM1_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("psiM1_TrackerLWM_");
           psiM1_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("psiM1_PixelLWM_");
           psiM1_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("psiM1_ValPixHit_");
           psiM2_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("psiM2_TrackerLWM_");
           psiM2_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("psiM2_PixelLWM_");
           psiM2_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("psiM2_ValPixHit_");

           ZLe1_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ZLe1_TrackerLWM_");
           ZLe1_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ZLe1_PixelLWM_");
           ZLe1_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ZLe1_ValPixHit_");
           ZLe2_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ZLe2_TrackerLWM_");
           ZLe2_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ZLe2_PixelLWM_");
           ZLe2_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ZLe2_ValPixHit_");

           dxym1= (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("Dxy");
           dxym2= (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("Dxy");;
           dzm1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("Dz"); ;
           dzm2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("Dz"); ;

           dxyl1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("Dxy");
           dxyl2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("Dxy");
           dzl1  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("Dz");
           dzl2  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("Dz");


           dipm1   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("dIP3D");
           dipm2   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("dIP3D");
           dipl1   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dIP3D");
           dipl2   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dIP3D");
           dipm1Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("dIP3DErr");
           dipm2Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("dIP3DErr");
           dipl1Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dIP3DErr");
           dipl2Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dIP3DErr");

           //R iso
           dRiso_l1   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dRIso");
           dRiso_l2   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dRIso");
           //new
           dRiso_m1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("dRIso");
           dRiso_m2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("dRIso");
           rIsoOverPtl1 = dRiso_l1/lepton1_p4.Pt();
           rIsoOverPtl2 = dRiso_l2/lepton2_p4.Pt();
           rIsoOverPtm1 = dRiso_m1/muonN_p4.Pt();
           rIsoOverPtm2 = dRiso_m2/muonP_p4.Pt();
           //std::cout<<"\n"<<std::endl;
           //reco::Candidate::LorentzVector _m1_ = z_Cand.daughter("psi")->daughter("muon1")->p4();
           //reco::Candidate::LorentzVector _m2_ = z_Cand.daughter("psi")->daughter("muon2")->p4();
           
           //std::cout<< "\ndRiso M1: "<< dRiso_m1 << std::endl;
           //std::cout<< "pT M1: "<< _m1_.Pt() << std::endl;
           //std::cout<< "dRisoOverpT M1: "<< dRiso_m1/ _m1_.Pt() << std::endl;
           //std::cout<< "dRiso M2: "<< dRiso_m2 << std::endl;
           //std::cout<< "pT M2: "<< _m2_.Pt() << std::endl;
           //std::cout<< "dRisoOverpT M2: "<< dRiso_m2/ _m2_.Pt() << std::endl;
           //new
           ipSm1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon1")))->userFloat("dIP3DSig");
           ipSm2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("psi")->daughter("muon2")))->userFloat("dIP3DSig");
           ipSl1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dIP3DSig");
           ipSl2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dIP3DSig");
           //new
           ImparSigl1 = dipl1/dipl1Err;
           ImparSigl2 = dipl2/dipl2Err;
           ImparSigm1 = dipm1/dipm1Err;
           ImparSigm2 = dipm2/dipm2Err;
             
           dR_m1_m2_ = z_Cand.userFloat("dRm1m2_");
           dR_l1_l2_ = z_Cand.userFloat("dRl1l2_");
           dR_m1_l1_ = z_Cand.userFloat("dRl1m1_");
           dR_m1_l2_ = z_Cand.userFloat("dRl1m2_");
           dR_m2_l1_ = z_Cand.userFloat("dRl2m1_");
           dR_m2_l2_ = z_Cand.userFloat("dRl2m2_");
             
           dR_m1_m2 = z_Cand.userFloat("dRm1m2");
           dR_l1_l2 = z_Cand.userFloat("dRl1l2");
           dR_m1_l1 = z_Cand.userFloat("dRl1m1");
           dR_m1_l2 = z_Cand.userFloat("dRl2m1");
           dR_m2_l1 = z_Cand.userFloat("dRl1m2");
           dR_m2_l2 = z_Cand.userFloat("dRl2m2");

           Event_Cand = z_Cand.userInt("Event_Cand_");
           Event_Cand++;
           Z_tree->Fill();
           //std::cout << "Z-Fill" << std::endl;
         }//end for csize
       }//end if Zcand.isValid and Zcand not empty
    //}// end if not only gen
   //std::cout << "Pass to psi tupler " << std::endl;

}//end analyze  


// ------------ method called once each job just before starting event loop  ------------
void 
ZjpsiMCTupler::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZjpsiMCTupler::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZjpsiMCTupler::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//trigger function from Jhovanny`s
void ZjpsiMCTupler::CheckHLTTriggers(const std::vector<std::string>& TrigList){

    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace helper;


    string AllTrg="";
    string tmptrig;

    int ntrigs=TrigList.size();
    if (ntrigs==0)
        std::cout << "No trigger name given in TriggerResults of the input " << endl;

    for (int itrig=0; itrig< ntrigs; itrig++) {
         //TString trigName = triggerNames_.triggerName(itrig);
         string trigName = TrigList.at(itrig);
         std::cout << "saving ..." << trigName << std::endl;
         tmptrig = (string) trigName; tmptrig +=" ";
         AllTrg += tmptrig;
    }

    //int n = sprintf(triggersL,"%s","");
   int n = sprintf(triggersL,"%s",AllTrg.c_str());
   std::cout<<" INFO: Triggers :  "<<triggersL<<std::endl;
   nTrgL = AllTrg.size();

   return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZjpsiMCTupler);
