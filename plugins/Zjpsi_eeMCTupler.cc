// -*- C++ -*-
//
// Package:    AnalyzeZll/ZRootupler
// Class:      ZRootupler  -> Zjpsi_eeMCTupler
// 
/**\class ZRootupler ZRootupler.cc AnalyzeZll/ZRootupler/plugins/ZRootupler.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rogelio Reyes Almanza
//         Created:  Thu, 09 Aug 2018 10:10:36 GMT
//         Updated:  v7-12/09/20
//

//this is hardcoded from Rogelio's code, it is suposed to be used when que have no trigger selection of match
//we use Jhovannys code for B's lifetime for reference 

// system include files
#include <memory>



// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Gabriel
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "MuonAnalysis/MuonAssociators/interface/PropagateToMuon.h"

#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "HeavyFlavorAnalysis/Onia2MuMu/interface/OniaVtxReProducer.h"


#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
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

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include <string> 
#include <vector> 
//
// class declaration
//

class Zjpsi_eeMCTupler:public edm::EDAnalyzer {
      public:
	explicit Zjpsi_eeMCTupler(const edm::ParameterSet &);
	~Zjpsi_eeMCTupler() override;

	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
        //from Jhovanny's 
        void CheckHLTTriggers(const std::vector<std::string>& TrigList);
        //Gabriel
        void CheckLepHLT(const std::vector<std::string>& TrigList);
        void CheckLepQid(const std::vector<std::string>& TrigList);  //ESA NO ES LA VARIABLE INDICADA WORK IN PROGRES

        UInt_t getTriggerBits(const edm::Event &, std::vector<std::string>);
        bool   isAncestor(const reco::Candidate *, const reco::Candidate *);
        const  reco::Candidate* GetAncestor(const reco::Candidate *);
        UInt_t isTriggerMatchedbyFilter(const pat::CompositeCandidate *);
        UInt_t isTriggerMatchedbyPath(const pat::CompositeCandidate *);
        UInt_t isTriggerMatch(const edm::Event&, const pat::CompositeCandidate*, std::vector<std::string>);

	    //void beginJob() override;
	    void analyze(const edm::Event &, const edm::EventSetup &) override;
	    //void endJob() override;

        //void beginRun(const edm::Run &, const edm::EventSetup &) override;
        //void endRun(edm::Run const &, edm::EventSetup const &) override;
        //void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
        //void endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
        // ----------member data ---------------------------
        std::string file_name;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> theZ_;
        edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  	    edm::EDGetTokenT<edm::TriggerResults> trigResultsToken;
  	    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigObjCollToken;
        //added from Jhovanny's
        std::string hlTriggerResults_;
        //Gabriel
        edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_ ; 
        int  pdgid_;
        std::vector<std::string> FilterNames_;
        //std::vector<std::string> SingleFilterNames_;
	    bool isMC_;
        bool bestCandidateOnly_;
        bool OnlyGen_;
        std::vector<std::string> HLTLastFilters;
        // need to eliminate HLTpaths_ NOTEFORSEARCH
        std::vector<std::string> HLTpaths_;

        //PropagateToMuon prop1_, prop2_; 

	    UInt_t    run;
	    ULong64_t event;
        UInt_t    lumiblock;
        UInt_t    ismatched;
        UInt_t    nonia;
        UInt_t    nmuons;
        UInt_t    nelecs;
        UInt_t    vertex_i;

        UInt_t    trigger;
        UInt_t    triggersingle;
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
        //New
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

        Float_t dxyl1_gsf;
        Float_t dxyl2_gsf;
        Float_t dzl1_gsf;
        Float_t dzl2_gsf;
    
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
    
        //double coriso= 99;
	    TTree *Z_tree;

        Int_t mother_pdgId;
        Int_t dimuon_pdgId;

          
        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;
        edm::EDGetTokenT<pat::PackedGenParticleCollection> packCands_;

        //from Jhovanny
        char triggersL[10000];
        int  nTrgL;
        //Gabriel
        
        //Trigger Match
        //int JpM1Mu8DiEle12;
        //int JpM2Mu8DiEle12;

        int JpM1Qid;

        int JpM1_TrackerLWM;
        int JpM1_PixelLWM;
        int JpM1_ValPixHit;

        int JpM2Qid;

        int JpM2_TrackerLWM;
        int JpM2_PixelLWM;
        int JpM2_ValPixHit;
        //Triger Match
        //int ZLe1Ele25wpT;
        //int ZLe2Ele23_12;
        //int ZLe1Mu8DiEle12;
        //int ZLe2Ele25wpT;
        //int ZLe1Ele23_12;
        //int ZLe2Mu8DiEle12;

        // Electron Q_ID
        int ZLe1Qid;
    
        float ZLe1dRIsoEA;
        float ZLe1trackMomentumAtVtx;
        float ZLe1ecalEnergy;
        float ZLe1full5x5_sigmaIetaIeta;
        float ZLe1dEtaIn;
        float ZLe1dPhiIn;
        float ZLe1HoE;
        float ZLe1ooEmooP;
        float ZLe1passConversionVeto;
    
        int ZLe1_TrackerLWM;
        int ZLe1_PixelLWM;
        int ZLe1_ValPixHit;

        float ZLe1dPhiInSeed;
        float ZLe1dEtaInSeed;
        float ZLe1SigmaIEtaIEta;
        float ZLe1SigmaIPhiIPhi;
        float ZLe1HoverE;
        float ZLe1ElecMissHits;
    
        float l1_ecalEnergyPreCorr;
        float l1_ecalEnergyErrPreCorr;
        float l1_ecalEnergyPostCorr;
        float l1_ecalEnergyErrPostCorr;
        float l1_ecalTrkEnergyPreCorr;
        float l1_ecalTrkEnergyErrPreCorr;
        float l1_ecalTrkEnergyPostCorr;
        float l1_ecalTrkEnergyErrPostCorr;
        float l1_energyScaleValue;
        float l1_energySigmaValue;
    
        float l2_ecalEnergyPreCorr;
        float l2_ecalEnergyErrPreCorr;
        float l2_ecalEnergyPostCorr;
        float l2_ecalEnergyErrPostCorr;
        float l2_ecalTrkEnergyPreCorr;
        float l2_ecalTrkEnergyErrPreCorr;
        float l2_ecalTrkEnergyPostCorr;
        float l2_ecalTrkEnergyErrPostCorr;
        float l2_energyScaleValue;
        float l2_energySigmaValue;
    
        // Electron Q_ID
        int ZLe2Qid;
    
        float ZLe2dRIsoEA;
        float ZLe2trackMomentumAtVtx;
        float ZLe2ecalEnergy;
        float ZLe2full5x5_sigmaIetaIeta;
        float ZLe2dEtaIn;
        float ZLe2dPhiIn;
        float ZLe2HoE;
        float ZLe2ooEmooP;
        float ZLe2passConversionVeto;
    
        float ZLe2dPhiInSeed;
        float ZLe2dEtaInSeed;
        float ZLe2SigmaIEtaIEta;
        float ZLe2SigmaIPhiIPhi;
        float ZLe2HoverE;
        float ZLe2ElecMissHits; //mva Spring16 v1 HZZ wp Loose

        int ZLe2_TrackerLWM;
        int ZLe2_PixelLWM;
        int ZLe2_ValPixHit;
    
        int Event_Cand;

};

//
// constructors and destructor
//

Zjpsi_eeMCTupler::Zjpsi_eeMCTupler(const edm::ParameterSet & iConfig):
	trigResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults")))
{
	theZ_ = consumes<pat::CompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("Zcand")); 
    hlTriggerResults_ = iConfig.getUntrackedParameter<std::string>("HLTriggerResults",std::string("TriggerResults::HLT"));
    edm::Service < TFileService > fs;
    
    Z_tree = fs->make < TTree > ("ZTree", "Tree of Z2jpsielel");

    Z_tree->Branch("run",      &run,      "run/i");
    Z_tree->Branch("event",    &event,    "event/l");
    Z_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");
    Z_tree->Branch("ismatched",&ismatched,"ismatched/i");
    Z_tree->Branch("vertex_i", &vertex_i, "vertex_i/i")

    //if (!OnlyGen_) {
    //new
    Z_tree->Branch("nonia",    &nonia,    "nonia/i");
    Z_tree->Branch("nmuons",   &nmuons,   "nmuons/i");
    Z_tree->Branch("nelecs",   &nelecs,   "nelecs/i");
    Z_tree->Branch("nPV",      &nPV,      "nPV/i");
    Z_tree->Branch("passFit",  &passFit,  "passFit/i");
    Z_tree->Branch("nCands",   &nCands,   "nCands/i");

    Z_tree->Branch("trigger",  &trigger,  "trigger/i");
    Z_tree->Branch("triggersingle",&triggersingle,"triggersingle/i");

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
    Z_tree->Branch("mard_dilepton_p4", "TLorentzVector", &msrd_dilepton_p4);
    Z_tree->Branch("msrd_lepton1_p4",  "TLorentzVector", &msrd_lepton1_p4);
    Z_tree->Branch("msrd_lepton2_p4",  "TLorentzVector", &msrd_lepton2_p4);
    Z_tree->Branch("msrd_Z_p4", "TLorentzVector", &msrd_Z_p4);

    Z_tree->Branch("pvChi2", &pvChi2, "pvChi2/F");
    Z_tree->Branch("Zvtx", "TVector3", &Zvtx);
    Z_tree->Branch("ZvtxP", &ZvtxP, "ZvtxP/F");
    Z_tree->Branch("ZvtxC2", &ZvtxC2, "ZvtxC2/F");

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

    Z_tree->Branch("dxyl1_gsf", &dxyl1_gsf, "dxyl1_gsf/F");
    Z_tree->Branch("dxyl2_gsf", &dxyl2_gsf, "dxyl2_gsf/F");
    Z_tree->Branch("dzl1_gsf", &dzl1_gsf, "dzl1_gsf/F");
    Z_tree->Branch("dzl2_gsf", &dzl2_gsf, "dzl2_gsf/F");                  
                                       
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
    

    Z_tree->Branch("dipm1",	&dipm1,"dipm1/F" );	  
    Z_tree->Branch("dipm2",	&dipm2,"dipm2/F");
    Z_tree->Branch("dipl1",	&dipl1,"dipl1/F" );
    Z_tree->Branch("dipl2",	&dipl2,"dipl2/F");
    Z_tree->Branch("dipm1Err",	&dipm1Err,"dipm1Err/F") ;
    Z_tree->Branch("dipm2Err",	&dipm2Err,"dipm2Err/F")	;
    Z_tree->Branch("dipl1Err",	&dipl1Err,"dipl1Err/F") ;
    Z_tree->Branch("dipl2Err",	&dipl2Err,"dipl2Err/F")	;
    //new
    Z_tree->Branch("ipSm1",    &ipSm1,"ipSm1/F") ;
    Z_tree->Branch("ipSm2",    &ipSm2,"ipSm2/F") ;
    Z_tree->Branch("ipSl1",    &ipSl1,"ipSl1/F") ;
    Z_tree->Branch("ipSl2",    &ipSl2,"ipSl2/F") ;
    Z_tree->Branch("ImparSigl1",    &ImparSigl1,"ImparSigl1/F") ;
    Z_tree->Branch("ImparSigl2",    &ImparSigl2,"ImparSigl2/F") ;
    Z_tree->Branch("ImparSigm1",    &ImparSigm1,"ImparSigm1/F") ;
    Z_tree->Branch("ImparSigm2",    &ImparSigm2,"ImparSigm2/F") ;
    

    //from Jhovanny
    Z_tree->Branch("nTrgL",      &nTrgL,    "nTrgL/I");  
    Z_tree->Branch("triggersL",         &triggersL,  "triggersL[nTrgL]/C");
    //Gabriel
/*    Z_tree->Branch("JpM1Trig",         &JpM1Trig,  "JpM1Trig[nTrgL]/C");
    Z_tree->Branch("JpM2Trig",         &JpM2Trig,  "JpM2Trig[nTrgL]/C");
    Z_tree->Branch("ZLe1Trig",         &ZLe1Trig,  "ZLe1Trig[nTrgL]/C");
    Z_tree->Branch("Zle2Trig",         &ZLe2Trig,  "ZLe2Trig[nTrgL]/C");*/
/*    Z_tree->Branch("JpM1Qid",         &JpM1Qid,  "JpM1Qid/i");
    Z_tree->Branch("JpM2Qid",         &JpM2Qid,  "JpM2Qid/i");
    Z_tree->Branch("ZLe1Qid",         &ZLe1Qid,  "ZLe1Qid/i");
    Z_tree->Branch("ZLe2Qid",         &ZLe2Qid,  "ZLe2Qid/i"); */
    // Jpsi Trigger Match
    // Z_tree->Branch("JpM1Mu8DiEle12",         &JpM1Mu8DiEle12,  "JpM1Mu8DiEle12/i");
    // Z_tree->Branch("JpM2Mu8DiEle12",         &JpM2Mu8DiEle12,  "JpM2Mu8DiEle12/i");
    
    Z_tree->Branch("JpM1Qid",         &JpM1Qid,  "JpM1Qid/i");

    Z_tree->Branch("JpM1_TrackerLWM",         &JpM1_TrackerLWM,  "JpM1_TrackerLWM/i");
    Z_tree->Branch("JpM1_PixelLWM",           &JpM1_PixelLWM,    "JpM1_PixelLWM/i");
    Z_tree->Branch("JpM1_ValPixHit",          &JpM1_ValPixHit,   "JpM1_ValPixHit/i");

    Z_tree->Branch("JpM2Qid",         &JpM2Qid,  "JpM2Qid/i");

    Z_tree->Branch("JpM2_TrackerLWM",         &JpM2_TrackerLWM,  "JpM2_TrackerLWM/i");
    Z_tree->Branch("JpM2_PixelLWM",           &JpM2_PixelLWM,    "JpM2_PixelLWM/i");
    Z_tree->Branch("JpM2_ValPixHit",          &JpM2_ValPixHit,   "JpM2_ValPixHit/i");

    Z_tree->Branch("ZLe1Qid",         &ZLe1Qid,  "ZLe1Qid/i");
    
    Z_tree->Branch("ZLe1dRIsoEA",               &ZLe1dRIsoEA,              "ZLe1dRIsoEA/i");
    Z_tree->Branch("ZLe1trackMomentumAtVtx",    &ZLe1trackMomentumAtVtx,   "ZLe1trackMomentumAtVtx/i");
    Z_tree->Branch("ZLe1ecalEnergy",            &ZLe1ecalEnergy,           "ZLe1ecalEnergy/i");
    Z_tree->Branch("ZLe1full5x5_sigmaIetaIeta", &ZLe1full5x5_sigmaIetaIeta,"ZLe1full5x5_sigmaIetaIeta/i");
    Z_tree->Branch("ZLe1dEtaIn",                &ZLe1dEtaIn,               "ZLe1dEtaIn/i");
    Z_tree->Branch("ZLe1dPhiIn",                &ZLe1dPhiIn,               "ZLe1dPhiIn/i");
    Z_tree->Branch("ZLe1HoE",                   &ZLe1HoE,                  "ZLe1HoE/i");
    Z_tree->Branch("ZLe1ooEmooP",               &ZLe1ooEmooP,              "ZLe1ooEmooP/i");
    Z_tree->Branch("ZLe1passConversionVeto",    &ZLe1passConversionVeto,   "ZLe1passConversionVeto/i");
                   
    Z_tree->Branch("ZLe1_TrackerLWM",         &ZLe1_TrackerLWM,  "ZLe1_TrackerLWM/i");
    Z_tree->Branch("ZLe1_PixelLWM",           &ZLe1_PixelLWM,    "ZLe1_PixelLWM/i");
    Z_tree->Branch("ZLe1_ValPixHit",          &ZLe1_ValPixHit,   "ZLe1_ValPixHit/i");

    //Z_tree->Branch("ZLe1Ele25wpT",         &ZLe1Ele25wpT,  "ZLe1Ele25wpT/i");
    //Z_tree->Branch("ZLe1Ele23_12",         &ZLe1Ele23_12,  "ZLe1Ele23_12/i");
    //Z_tree->Branch("ZLe1Mu8DiEle12",         &ZLe1Mu8DiEle12,  "ZLe1Mu8DiEle12/i");

    Z_tree->Branch("ZLe2Qid",         &ZLe2Qid,  "ZLe2Qid/i");

    Z_tree->Branch("ZLe2dRIsoEA",               &ZLe2dRIsoEA,              "ZLe2dRIsoEA/i");
    Z_tree->Branch("ZLe2trackMomentumAtVtx",    &ZLe2trackMomentumAtVtx,   "ZLe2trackMomentumAtVtx/i");
    Z_tree->Branch("ZLe2ecalEnergy",            &ZLe2ecalEnergy,           "ZLe2ecalEnergy/i");
    Z_tree->Branch("ZLe2full5x5_sigmaIetaIeta", &ZLe2full5x5_sigmaIetaIeta,"ZLe2full5x5_sigmaIetaIeta/i");
    Z_tree->Branch("ZLe2dEtaIn",                &ZLe2dEtaIn,               "ZLe2dEtaIn/i");
    Z_tree->Branch("ZLe2dPhiIn",                &ZLe2dPhiIn,               "ZLe2dPhiIn/i");
    Z_tree->Branch("ZLe2HoE",                   &ZLe2HoE,                  "ZLe2HoE/i");
    Z_tree->Branch("ZLe2ooEmooP",               &ZLe2ooEmooP,              "ZLe2ooEmooP/i");
    Z_tree->Branch("ZLe2passConversionVeto",    &ZLe2passConversionVeto,   "ZLe2passConversionVeto/i");
    
    //Z_tree->Branch("ZLe2Ele25wpT",         &ZLe2Ele25wpT,  "ZLe2Ele25wpT/i");
    //Z_tree->Branch("ZLe2Ele23_12",         &ZLe2Ele23_12,  "ZLe2Ele23_12/i");
    //Z_tree->Branch("ZLe2Mu8DiEle12",         &ZLe2Mu8DiEle12,  "ZLe2Mu8DiEle12/i");

    Z_tree->Branch("ZLe2_TrackerLWM",         &ZLe2_TrackerLWM,  "ZLe2_TrackerLWM/i");
    Z_tree->Branch("ZLe2_PixelLWM",           &ZLe2_PixelLWM,    "ZLe2_PixelLWM/i");
    Z_tree->Branch("ZLe2_ValPixHit",          &ZLe2_ValPixHit,   "ZLe2_ValPixHit/i");

    Z_tree->Branch("ZLe1dPhiInSeed",         &ZLe1dPhiInSeed,     "ZLe1dPhiInSeed/F");
    Z_tree->Branch("ZLe1dEtaInSeed",         &ZLe1dEtaInSeed,     "ZLe1dEtaInSeed/F");
    Z_tree->Branch("ZLe1SigmaIEtaIEta",      &ZLe1SigmaIEtaIEta,  "ZLe1SigmaIEtaIEta/F");
    Z_tree->Branch("ZLe1SigmaIPhiIPhi",      &ZLe1SigmaIPhiIPhi,  "ZLe1SigmaIPhiIPhi/F");
    Z_tree->Branch("ZLe1HoverE",             &ZLe1HoverE,         "ZLe1HoverE/F");
    Z_tree->Branch("ZLe1ElecMissHits",       &ZLe1ElecMissHits,   "ZLe1ElecMissHits/F");

    Z_tree->Branch("ZLe2dPhiInSeed",         &ZLe2dPhiInSeed,     "ZLe2dPhiInSeed/F");
    Z_tree->Branch("ZLe2dEtaInSeed",         &ZLe2dEtaInSeed,     "ZLe2dEtaInSeed/F");
    Z_tree->Branch("ZLe2SigmaIEtaIEta",      &ZLe2SigmaIEtaIEta,  "ZLe2SigmaIEtaIEta/F");
    Z_tree->Branch("ZLe2SigmaIPhiIPhi",      &ZLe2SigmaIPhiIPhi,  "ZLe2SigmaIPhiIPhi/F");
    Z_tree->Branch("ZLe2HoverE",             &ZLe2HoverE,         "ZLe2HoverE/F");
    Z_tree->Branch("ZLe2ElecMissHits",       &ZLe2ElecMissHits,   "ZLe2ElecMissHits/F");
    
    Z_tree->Branch("l1_ecalEnergyPreCorr",        &l1_ecalEnergyPreCorr,        "l1_ecalEnergyPreCorr/F"       );
    Z_tree->Branch("l1_ecalEnergyErrPreCorr",     &l1_ecalEnergyErrPreCorr,     "l1_ecalEnergyErrPreCorr/F"    );
    Z_tree->Branch("l1_ecalEnergyPostCorr",       &l1_ecalEnergyPostCorr,       "l1_ecalEnergyPostCorr/F"      );
    Z_tree->Branch("l1_ecalEnergyErrPostCorr",    &l1_ecalEnergyErrPostCorr,    "l1_ecalEnergyErrPostCorr/F"   );
    Z_tree->Branch("l1_ecalTrkEnergyPreCorr",     &l1_ecalTrkEnergyPreCorr,     "l1_ecalTrkEnergyPreCorr/F"    );
    Z_tree->Branch("l1_ecalTrkEnergyErrPreCorr",  &l1_ecalTrkEnergyErrPreCorr,  "l1_ecalTrkEnergyErrPreCorr/F" );
    Z_tree->Branch("l1_ecalTrkEnergyPostCorr",    &l1_ecalTrkEnergyPostCorr,    "l1_ecalTrkEnergyPostCorr/F"   );
    Z_tree->Branch("l1_ecalTrkEnergyErrPostCorr", &l1_ecalTrkEnergyErrPostCorr, "l1_ecalTrkEnergyErrPostCorr",);
    Z_tree->Branch("l1_energyScaleValue",         &l1_energyScaleValue,         "l1_energyScaleValue",        );
    Z_tree->Branch("l1_energySigmaValue",         &l1_energySigmaValue,         "l1_energySigmaValue",        );
    
    Z_tree->Branch("l2_ecalEnergyPreCorr",        &l2_ecalEnergyPreCorr,        "l2_ecalEnergyPreCorr",       );
    Z_tree->Branch("l2_ecalEnergyErrPreCorr",     &l2_ecalEnergyErrPreCorr,     "l2_ecalEnergyErrPreCorr",    );
    Z_tree->Branch("l2_ecalEnergyPostCorr",       &l2_ecalEnergyPostCorr,       "l2_ecalEnergyPostCorr",      );
    Z_tree->Branch("l2_ecalEnergyErrPostCorr",    &l2_ecalEnergyErrPostCorr,    "l2_ecalEnergyErrPostCorr",   );
    Z_tree->Branch("l2_ecalTrkEnergyPreCorr",     &l2_ecalTrkEnergyPreCorr,     "l2_ecalTrkEnergyPreCorr",    );
    Z_tree->Branch("l2_ecalTrkEnergyErrPreCorr",  &l2_ecalTrkEnergyErrPreCorr,  "l2_ecalTrkEnergyErrPreCorr", );
    Z_tree->Branch("l2_ecalTrkEnergyPostCorr",    &l2_ecalTrkEnergyPostCorr,    "l2_ecalTrkEnergyPostCorr",   );
    Z_tree->Branch("l2_ecalTrkEnergyErrPostCorr", &l2_ecalTrkEnergyErrPostCorr, "l2_ecalTrkEnergyErrPostCorr",);
    Z_tree->Branch("l2_energyScaleValue",         &l2_energyScaleValue,         "l2_energyScaleValue",        );
    Z_tree->Branch("l2_energySigmaValue",         &l2_energySigmaValue,         "l2_energySigmaValue",        );
 
    
    
    
    
    
    Z_tree->Branch("ZLe1CorrEt",         &ZLe1CorrEt,     "ZLe1CorrEt/F");
    Z_tree->Branch("ZLe1CorrFact",       &ZLe1CorrFact,   "ZLe1CorrFact/F");
    
    Z_tree->Branch("ZLe2CorrEt",         &ZLe2CorrEt,     "ZLe2CorrEt/F");
    Z_tree->Branch("ZLe2CorrFact",       &ZLe2CorrFact,   "ZLe2CorrFact/F");
    
    Z_tree->Branch("Event_Cand", &Event_Cand, "Event_Cand_/i");

//  }
/*ç*/
    /*
  if (isMC_ || OnlyGen_) {
     std::cout << "Zjpsi_eeMCTupler::Zjpsi_eeMCTupler: Onia id " << pdgid_ << std::endl;
     Z_tree->Branch("mother_pdgId",  &mother_pdgId,     "mother_pdgId/I");
     Z_tree->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
     Z_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
     Z_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
     Z_tree->Branch("gen_muonN_p4",  "TLorentzVector",  &gen_muonM_p4);
  }
  genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"prunedGenParticles");
  packCands_ = consumes<pat::PackedGenParticleCollection>((edm::InputTag)"packedGenParticles");
*/
}

Zjpsi_eeMCTupler::~Zjpsi_eeMCTupler() {}

//
// member functions
//

//trigger function from Jhovanny`s
void Zjpsi_eeMCTupler::CheckHLTTriggers(const std::vector<std::string>& TrigList){

    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace helper;
    
    
    string AllTrg="";
    string tmptrig;

    int ntrigs=TrigList.size();
    if (ntrigs==0)
        //std::cout << "No trigger name given in TriggerResults of the input " << endl;
    
    for (int itrig=0; itrig< ntrigs; itrig++) {
        //TString trigName = triggerNames_.triggerName(itrig);
        string trigName = TrigList.at(itrig);
         //std::cout << "saving ..." << trigName << std::endl; 
         tmptrig = (string) trigName; tmptrig +=" ";
         AllTrg += tmptrig;
    }

    //int n = sprintf(triggersL,"%s","");
   int n = sprintf(triggersL,"%s",AllTrg.c_str());
   //std::cout<<" INFO: Triggers :  "<<triggersL<<std::endl;
   nTrgL = AllTrg.size();

   return;
}


UInt_t Zjpsi_eeMCTupler::getTriggerBits(const edm::Event& iEvent, std::vector<std::string> TestFilterNames_ ) {
   UInt_t trigger = 0;
   edm::Handle<edm::TriggerResults> triggerResults_handle;
   iEvent.getByToken(trigResultsToken, triggerResults_handle);
   if (triggerResults_handle.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerResults_handle);
      for (unsigned int i = 0; i < TestFilterNames_.size(); i++) {
         for (int version = 1; version < 9; version++) {
            std::stringstream ss;
            ss << TestFilterNames_[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
            if (bit < triggerResults_handle->size() && triggerResults_handle->accept(bit) && !triggerResults_handle->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }
   } else std::cout << "Zjpsi_eeMCTupler::getTriggerBits: *** NO triggerResults found *** " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
   return trigger;
}

// ------------ method called for each event  ------------
void Zjpsi_eeMCTupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle<pat::CompositeCandidateCollection> ZCands;
  iEvent.getByToken(theZ_,ZCands);

  pat::CompositeCandidate z_Cand;  
   

 
  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();
    
  ismatched = 0;

  // try to eliminate NOTEFORSEARCH
  trigger       = getTriggerBits(iEvent,HLTpaths_);
//  triggersingle = getTriggerBits(iEvent,SingleFilterNames_);
 
  dxym1 =0.0;   
  dxym2 =0.0;  
  dzm1  =0.0;  
  dzm2  =0.0;  
         
  dxyl1 =0.0; 
  dxyl2 =0.0;
  dzl1  =0.0;
  dzl2  =0.0;
 
  dxyl1_gsf =0.0;
  dxyl2_gsf =0.0;
  dzl1_gsf  =0.0;
  dzl2_gsf  =0.0;
 
  dRiso_l1   =0.0; 
  dRiso_l2   =0.0;
      
  //Jpsi_vprob =0.0;
  //Jpsi_vchi2 =0.0;

  dR_m1_m2 =0.0;   
  dR_l1_l2 =0.0;
  dR_m1_l1 =0.0;
  dR_m1_l2 =0.0;
  dR_m2_l1 =0.0;  
  dR_m2_l2 =0.0;

   dipm1   =0.0;  
   dipm2   =0.0;
   dipl1   =0.0;
   dipl2   =0.0;
   dipm1Err=0.0;
   dipm2Err=0.0;
   dipl1Err=0.0; 
   dipl2Err=0.0;

  ///////////////////////
  //Jhovanny's stuff
  // ///////////////////
  //chek if not repeated
  
  edm::Handle<edm::TriggerResults> hltresults1;
  try {
    std::string const &trig = std::string("TriggerResults::HLT");//+hlTriggerResults_;
    //std::string const &trig = std::string("TriggerResults::")+trigResultsToken;
    iEvent.getByLabel(edm::InputTag(trig),hltresults1);
    //std::cout << "Handle on HLT ok " << std::endl;
  }
  catch ( ... ) 
    {
      //std::cout << "Couldn't get handle on HLT Trigger!" << std::endl;
    }
    
    //HLTConfigProvider hltrigConfig_;
    //std::cout << "Filling trig table: "<< std::endl; 
    std::vector<std::string> TrigTable; TrigTable.clear();
    // Get hold of trigger names - based on TriggerResults object
    const edm::TriggerNames &triggerNames1_ = iEvent.triggerNames(*hltresults1);
    
    for (unsigned int itrig = 0; itrig < hltresults1->size(); ++itrig){
        if ((*hltresults1)[itrig].accept() == 1){
            std::string trigName1 = triggerNames1_.triggerName(itrig);
            //int trigPrescale = hltConfig_.prescaleValue(itrig, trigName1);
            TrigTable.push_back(trigName1);
           // std::cout << trigName1 << std::endl;
        }
    }
  //std::cout << "Finished trig table" << std::endl;
  CheckHLTTriggers(TrigTable);

  //std::cout << "is Empty ? " << ZCands->empty() << std::endl;
  if (ZCands.isValid() && !ZCands->empty()) {
      std::cout << "ZCands is not empty: tupler working ..." << std::endl;
      unsigned int csize = ZCands->size();
      //if (bestCandidateOnly_) csize = 1;

      for ( unsigned int i = 0; i < csize; i++ ) {
      		
		
         	     
	   z_Cand = (ZCands->at(i));
	   

	  //ismatched = isTriggerMatch(iEvent, z_Cand, HLTpaths_); 	
       //new
       nonia  = z_Cand.userInt("nonia_");
       nmuons = z_Cand.userInt("nmuons_");
       nelecs = z_Cand.userInt("nelecs_");
       nPV    = z_Cand.userInt("nPV_");
       passFit = z_Cand.userInt("passFit_");
       vertex_i = z_Cand.userInt("pvIndex");
          
       nCands  = csize;
          
	   Z_p4.SetPtEtaPhiM(z_Cand.pt(), z_Cand.eta(), z_Cand.phi(), z_Cand.mass());
	   dimuon_p4.SetPtEtaPhiM(z_Cand.daughter("jpsi")->pt(), z_Cand.daughter("jpsi")->eta(),
                                  z_Cand.daughter("jpsi")->phi(), z_Cand.daughter("jpsi")->mass());
       
	   lepton1_p4.SetPtEtaPhiM(z_Cand.daughter("lepton1")->pt(), z_Cand.daughter("lepton1")->eta(), z_Cand.daughter("lepton1")->phi(), z_Cand.daughter("lepton1")->mass());
	   lepton2_p4.SetPtEtaPhiM(z_Cand.daughter("lepton2")->pt(), z_Cand.daughter("lepton2")->eta(), z_Cand.daughter("lepton2")->phi(), z_Cand.daughter("lepton2")->mass());
       //new
       reco::Candidate::LorentzVector dilep;
       dilep = z_Cand.daughter("lepton1")->p4() + z_Cand.daughter("lepton2")->p4();
       dilepton_p4.SetPtEtaPhiM(dilep.pt(), dilep.eta(), dilep.phi(), dilep.mass());
    
       reco::Candidate::LorentzVector vP;//= z_Cand.daughter("jpsi")->daughter("muon1")->p4();
       reco::Candidate::LorentzVector vM;//= z_Cand.daughter("jpsi")->daughter("muon2")->p4();

	   if (z_Cand.daughter("jpsi")->daughter("muon1")->charge() < 0) {
	      vP = z_Cand.daughter("jpsi")->daughter("muon2")->p4();
	      vM = z_Cand.daughter("jpsi")->daughter("muon1")->p4();
	   }
	   else{
	      vP = z_Cand.daughter("jpsi")->daughter("muon1")->p4();
	      vM = z_Cand.daughter("jpsi")->daughter("muon2")->p4();
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
       if (z_Cand.daughter("jpsi")->daughter("muon1")->pt() > z_Cand.daughter("jpsi")->daughter("muon2")->pt()){
           lm = z_Cand.daughter("jpsi")->daughter("muon1")->p4();
           tm = z_Cand.daughter("jpsi")->daughter("muon2")->p4();
          }
       else{
           lm = z_Cand.daughter("jpsi")->daughter("muon2")->p4();
           tm = z_Cand.daughter("jpsi")->daughter("muon1")->p4();
       }
       leading_lepton_p4.SetPtEtaPhiM( ll.pt(), ll.eta(), ll.phi(), ll.mass());
       trailing_lepton_p4.SetPtEtaPhiM(tl.pt(), tl.eta(), tl.phi(), tl.mass());
       leading_muon_p4.SetPtEtaPhiM(   lm.pt(), lm.eta(), lm.phi(), lm.mass());
       trailing_muon_p4.SetPtEtaPhiM(  tm.pt(), tm.eta(), tm.phi(), tm.mass());
          
          
       //Measured values
       msrd_Z_p4.SetPtEtaPhiM(z_Cand.daughter("patMZ")->pt(), z_Cand.daughter("patMZ")->eta(), z_Cand.daughter("patMZ")->phi(), z_Cand.daughter("patMZ")->mass());
       msrd_dimuon_p4.SetPtEtaPhiM(z_Cand.daughter("msrd_jpsi")->pt(), z_Cand.daughter("msrd_jpsi")->eta(), z_Cand.daughter("msrd_jpsi")->phi(), z_Cand.daughter("msrd_jpsi")->mass());
       msrd_lepton1_p4.SetPtEtaPhiM(z_Cand.daughter("msrd_lepton1")->pt(), z_Cand.daughter("msrd_lepton1")->eta(), z_Cand.daughter("msrd_lepton1")->phi(), z_Cand.daughter("msrd_lepton1")->mass());
       msrd_lepton2_p4.SetPtEtaPhiM(z_Cand.daughter("msrd_lepton2")->pt(), z_Cand.daughter("msrd_lepton2")->eta(), z_Cand.daughter("msrd_lepton2")->phi(), z_Cand.daughter("msrd_lepton2")->mass());
	   reco::Candidate::LorentzVector msrd_vP = z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon1")->p4();
	   reco::Candidate::LorentzVector msrd_vM = z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon2")->p4();
       //new
       reco::Candidate::LorentzVector msrd_dilep;
       msrd_dilep = z_Cand.daughter("msrd_lepton1")->p4() + z_Cand.daughter("msrd_lepton2")->p4();
       msrd_dilepton_p4.SetPtEtaPhiM(msrd_dilep.pt(), msrd_dilep.eta(), msrd_dilep.phi(), msrd_dilep.mass());
          
       gen_z_p4.SetPtEtaPhiM(z_Cand.daughter("mcZ")->pt(), z_Cand.daughter("mcZ")->eta(), z_Cand.daughter("mcZ")->phi(), z_Cand.daughter("mcZ")->mass());
       gen_muonP_p4.SetPtEtaPhiM(z_Cand.daughter("mcM2")->pt(), z_Cand.daughter("mcM2")->eta(), z_Cand.daughter("mcM2")->phi(), z_Cand.daughter("mcM2")->mass());
       gen_muonN_p4.SetPtEtaPhiM(z_Cand.daughter("mcM1")->pt(), z_Cand.daughter("mcM1")->eta(), z_Cand.daughter("mcM1")->phi(), z_Cand.daughter("mcM1")->mass());
       gen_dimuon_p4.SetPtEtaPhiM(z_Cand.daughter("mcPsi")->pt(), z_Cand.daughter("mcPsi")->eta(), z_Cand.daughter("mcPsi")->phi(), z_Cand.daughter("mcPsi")->mass());
       gen_lepton1_p4.SetPtEtaPhiM(z_Cand.daughter("mcL1")->pt(), z_Cand.daughter("mcL1")->eta(), z_Cand.daughter("mcL1")->phi(), z_Cand.daughter("mcL1")->mass());
       gen_lepton2_p4.SetPtEtaPhiM(z_Cand.daughter("mcL2")->pt(), z_Cand.daughter("mcL2")->eta(), z_Cand.daughter("mcL2")->phi(), z_Cand.daughter("mcL2")->mass());
       gen_Zvtx.SetXYZ(z_Cand.daughter("mcZ")->vx() ,z_Cand.daughter("mcZ")->vy(),z_Cand.daughter("mcZ")->vz()) ;
       gen_psi_vtx.SetXYZ(z_Cand.daughter("mcPsi")->vx() ,z_Cand.daughter("mcPsi")->vy(),z_Cand.daughter("mcPsi")->vz()) ;
	   if (z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon1")->charge() < 0) {
	      msrd_vP = z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon2")->p4();
	      msrd_vM = z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon1")->p4();
	   }
	   else{
	      msrd_vP = z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon1")->p4();
	      msrd_vM = z_Cand.daughter("msrd_jpsi")->daughter("msrd_muon2")->p4();
           } 
	   msrd_muonP_p4.SetPtEtaPhiM(msrd_vP.pt(), msrd_vP.eta(), msrd_vP.phi(), msrd_vP.mass());
	   msrd_muonN_p4.SetPtEtaPhiM(msrd_vM.pt(), msrd_vM.eta(), msrd_vM.phi(), msrd_vM.mass());
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
          
// 	   ismatched 	= isTriggerMatchedbyPath(Z_);

           
	   //Gabriel
	   //if(z_Cand.daughter("jpsi")->daughter("muon1")->isTightMuon(*PV)) std::cout<< "fbcudbckbdsckjbdBKSDDSBJSDKBJ" <<std::endl;
           //if (z_Cand.daughter("lepton1")->isTightMuon(*PV)) std::cout << "DSNVCOINEFBDVB" << std::endl;
           //std::cout <<  "Provando" << std::endl;

       pvChi2 = z_Cand.userFloat("pvChi2_");

	   ZvtxP = z_Cand.userFloat("vProb") ;
       ZvtxC2 = z_Cand.userFloat("vChi2") ;
	   Zvtx.SetXYZ(z_Cand.userFloat("ZvtxX") ,z_Cand.userFloat("ZvtxY"),z_Cand.userFloat("ZvtxZ")) ;
          
       psiVtxP = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")))->userFloat("vProb");
       psiVtxC2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")))->userFloat("vChi2");
       //Gabriel
       JpM1Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userInt("JpM1Qid_");
       JpM2Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userInt("JpM2Qid_");
       ZLe1Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userInt("ZLe1Qid_");
       ZLe2Qid = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userInt("ZLe2Qid_");
           
       //Trigger Matching
       /*
       JpM1Mu8DiEle12 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("muon1Mu8DiEle12_");
       JpM2Mu8DiEle12 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("muon2Mu8DiEle12_");
       */
       JpM1_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("psiM1_TrackerLWM_");
       JpM1_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("psiM1_PixelLWM_");
       JpM1_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("psiM1_ValPixHit_");

       JpM2_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("psiM2_TrackerLWM_");
       JpM2_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("psiM2_PixelLWM_");
       JpM2_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("psiM2_ValPixHit_");

       ZLe1_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ZLe1_TrackerLWM_");
       ZLe1_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ZLe1_PixelLWM_");
       ZLe1_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ZLe1_ValPixHit_");
       /*
       ZLe1Ele25wpT = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("lept1Ele25wpT_");
       ZLe1Ele23_12 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("lept1Ele23_12_");
       ZLe1Mu8DiEle12 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("lept1Mu8DiEle12_");
       */
       ZLe2_TrackerLWM = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ZLe2_TrackerLWM_");
       ZLe2_PixelLWM   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ZLe2_PixelLWM_");
       ZLe2_ValPixHit  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ZLe2_ValPixHit_");
       /*
       ZLe2Ele25wpT = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("lept2Ele25wpT_");
       ZLe2Ele23_12 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("lept2Ele23_12_");
       ZLe2Mu8DiEle12 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("lept2Mu8DiEle12_");
       */
           //int temp = isTriggerMatchedbyFilter(&z_Cand);
           //const pat::Muon* lept1 = dynamic_cast<const pat::Muon*>(z_Cand.daughter("lepton1"));
           //const pat::Muon* lept2 = dynamic_cast<const pat::Muon*>(z_Cand.daughter("lepton2"));
           /*try {
               const pat::TriggerObjectStandAloneCollection muHLTMatches1_t1 = lept1->triggerObjectMatchesByFilter("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09");
               std::cout << "IT WORKS ... for now" << std::endl;
verE
               //if (muHLTMatches1_t1.size() > 0) std::cout << "IT WORKS!!!!!" << std::endl;
           }
           catch ( ... ){
                std::cout << "Esta madre no jala" << std::endl;
           }*/
       dxym1= (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("Dxy");
       dxym2= (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("Dxy");;
       dzm1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("Dz"); ;
       dzm2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("Dz"); ;
        
       dxyl1_gsf = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("Dxy_gsf");
       dxyl2_gsf = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("Dxy_gsf");
       dzl1_gsf  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("Dz_gsf");
       dzl2_gsf  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("Dz_gsf");
   
       dxyl1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("Dxy");
       dxyl2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("Dxy");
       dzl1  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("Dz");
       dzl2  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("Dz");

	   dipm1   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("dIP3D");
       dipm2   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("dIP3D");
       dipl1   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dIP3D");
       dipl2   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dIP3D");
       dipm1Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("dIP3DErr");
       dipm2Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("dIP3DErr");
       dipl1Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dIP3DErr");
       dipl2Err=(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dIP3DErr");
		
	   //R iso
       dRiso_l1   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dRIso");
       dRiso_l2   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dRIso");
       //new
       dRiso_m1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("dRIso");
       dRiso_m2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("dRIso");
       rIsoOverPtl1 = dRiso_l1/lepton1_p4.Pt();
       rIsoOverPtl2 = dRiso_l2/lepton2_p4.Pt();
       rIsoOverPtm1 = dRiso_m1/muonN_p4.Pt();
       rIsoOverPtm2 = dRiso_m2/muonP_p4.Pt();
       
       //new
       ipSm1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon1")))->userFloat("dIP3DSig");
       ipSm2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")->daughter("muon2")))->userFloat("dIP3DSig");
       ipSl1 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dIP3DSig");
       ipSl2 = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dIP3DSig");
       //new
       ImparSigl1 = dipl1/dipl1Err;
       ImparSigl2 = dipl2/dipl2Err;
       ImparSigm1 = dipm1/dipm1Err;
       ImparSigm2 = dipm2/dipm2Err;
          
       //Jpsi_vprob =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")))->userFloat("vProb");
       //Jpsi_vchi2 =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("jpsi")))->userFloat("vChi2");
           
       ZLe1dPhiInSeed    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dPhiInSeed");
       ZLe1dEtaInSeed    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dEtaInSeed");
       ZLe1SigmaIEtaIEta =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("SigmaIEtaIEta");
       ZLe1SigmaIPhiIPhi =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("SigmaIPhiIPhi");
       ZLe1HoverE        =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("HoverE");
       ZLe1ElecMissHits  =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ElecMissHits");

       ZLe1dRIsoEA               =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dRIsoEA");
       ZLe1trackMomentumAtVtx    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("trackMomentumAtVtx");
       ZLe1ecalEnergy            =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalEnergy");
       ZLe1full5x5_sigmaIetaIeta =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("full5x5_sigmaIetaIeta");
       ZLe1dEtaIn                =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dEtaIn");
       ZLe1dPhiIn                =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("dPhiIn");
       ZLe1HoE                   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("HoE");
       ZLe1ooEmooP               =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ooEmooP");
       ZLe1passConversionVeto    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("passConversionVeto");
          
          
       ZLe2dPhiInSeed    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dPhiInSeed");
       ZLe2dEtaInSeed    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dEtaInSeed");
       ZLe2SigmaIEtaIEta =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("SigmaIEtaIEta");
       ZLe2SigmaIPhiIPhi =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("SigmaIPhiIPhi");
       ZLe2HoverE        =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("HoverE");
       ZLe2ElecMissHits  =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ElecMissHits");

       ZLe2dRIsoEA               =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dRIsoEA");
       ZLe2trackMomentumAtVtx    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("trackMomentumAtVtx");
       ZLe2ecalEnergy            =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalEnergy");
       ZLe2full5x5_sigmaIetaIeta =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("full5x5_sigmaIetaIeta");
       ZLe2dEtaIn                =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dEtaIn");
       ZLe2dPhiIn                =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("dPhiIn");
       ZLe2HoE                   =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("HoE");
       ZLe2ooEmooP               =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ooEmooP");
       ZLe2passConversionVeto    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("passConversionVeto");
          
       ZLe1CorrEt    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("corrEt_");
       ZLe1CorrFact  =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("corrfactor_");
       ZLe2CorrEt    =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("corrEt_");
       ZLe2CorrFact  =(dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("corrfactor_");

       l1_ecalEnergyPreCorr         = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalEnergyPreCorr_");
       l1_ecalEnergyErrPreCorr      = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalEnergyErrPreCorr_");
       l1_ecalEnergyPostCorr        = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalEnergyPostCorr_");
       l1_ecalEnergyErrPostCorr     = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalEnergyErrPostCorr_");
       l1_ecalTrkEnergyPreCorr      = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalTrkEnergyPreCorr_");
       l1_ecalTrkEnergyErrPreCorr   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalTrkEnergyErrPreCorr_");
       l1_ecalTrkEnergyPostCorr     = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalTrkEnergyPostCorr_");
       l1_ecalTrkEnergyErrPostCorr  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("ecalTrkEnergyErrPostCorr_");
       l1_energyScaleValue          = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("energyScaleValue_");
       l1_energySigmaValue          = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton1")))->userFloat("energySigmaValue_");
 
       l2_ecalEnergyPreCorr         = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalEnergyPreCorr_");
       l2_ecalEnergyErrPreCorr      = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalEnergyErrPreCorr_");
       l2_ecalEnergyPostCorr        = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalEnergyPostCorr_");
       l2_ecalEnergyErrPostCorr     = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalEnergyErrPostCorr_");
       l2_ecalTrkEnergyPreCorr      = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalTrkEnergyPreCorr_");
       l2_ecalTrkEnergyErrPreCorr   = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalTrkEnergyErrPreCorr_");
       l2_ecalTrkEnergyPostCorr     = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalTrkEnergyPostCorr_");
       l2_ecalTrkEnergyErrPostCorr  = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("ecalTrkEnergyErrPostCorr_");
       l2_energyScaleValue          = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("energyScaleValue_");
       l2_energySigmaValue          = (dynamic_cast<const pat::CompositeCandidate*>(z_Cand.daughter("lepton2")))->userFloat("energySigmaValue_");
          
       //ONLY WHEN Z->MU MU
       dR_m1_m2 = z_Cand.userFloat("dRm1m2");
       dR_l1_l2 = z_Cand.userFloat("dRl1l2");
       dR_m1_l1 = z_Cand.userFloat("dRl1m1");
       dR_m1_l2 = z_Cand.userFloat("dRl2m1");
       dR_m2_l1 = z_Cand.userFloat("dRl1m2");
       dR_m2_l2 = z_Cand.userFloat("dRl2m2");
          
       Event_Cand = z_Cand.userInt("Event_Cand_");
       Z_tree->Fill();
      } 
    }
 //std::cout << "code works ZElectron " << std::endl;
}

// Look if dimuon candidate is trigger matched
UInt_t Zjpsi_eeMCTupler::isTriggerMatchedbyFilter(const pat::CompositeCandidate *diMuon_cand) {
  //const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  //const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lepton1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lepton2"));
  UInt_t matched = 0;  // if no list is given, is not matched 

// if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<HLTLastFilters.size(); iTr++ ) {
     const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByFilter(HLTLastFilters[iTr]);
     const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByFilter(HLTLastFilters[iTr]);
     if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()){
         matched += (1<<iTr);
         //std::cout << "IT WORKS !!!" << std::endl;
     } 
  }
  return matched;
}
UInt_t Zjpsi_eeMCTupler::isTriggerMatchedbyPath(const pat::CompositeCandidate *diMuon_cand) {
  //const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon1"));
  //const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("muon2"));
  const pat::Muon* muon1 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lepton1"));
  const pat::Muon* muon2 = dynamic_cast<const pat::Muon*>(diMuon_cand->daughter("lepton2"));
  UInt_t matched = 0;  // if no list is given, is not matched 

// if matched a given trigger, set the bit, in the same order as listed
  for (unsigned int iTr = 0; iTr<FilterNames_.size(); iTr++ ) {
        for (int version = 1; version < 9; version++) {
            std::stringstream ss;
            ss << FilterNames_[iTr] << "_v" << version;
     const std::string triggerM = ss.str(); 
     const pat::TriggerObjectStandAloneCollection mu1HLTMatches = muon1->triggerObjectMatchesByPath(triggerM);
     const pat::TriggerObjectStandAloneCollection mu2HLTMatches = muon2->triggerObjectMatchesByPath(triggerM);
     if (!mu1HLTMatches.empty() && !mu2HLTMatches.empty()) matched += (1<<iTr); 
     }
  }
  return matched;
}

UInt_t Zjpsi_eeMCTupler::isTriggerMatch(const edm::Event& iEvent ,const pat::CompositeCandidate* cand, std::vector<std::string> TestFilterNames_){
 
  const pat::Electron* el1 = dynamic_cast<const pat::Electron*>(cand->daughter("lepton1"));
  const pat::Electron* el2 = dynamic_cast<const pat::Electron*>(cand->daughter("lepton2"));
 
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(trigResultsToken,triggerBits);
  const edm::TriggerNames& names = iEvent.triggerNames(*triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(trigObjCollToken,triggerObjects);
  UInt_t matched = 0;
   double dR1 = 0;  
   double DPtRel1 = 0; 
   double dR2 = 0;  
   double DPtRel2 = 0; 
   for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
      obj.unpackPathNames(names);
      if ( not (obj.hasTriggerObjectType(trigger::TriggerElectron))) continue;
      dR1 = deltaR(obj, *(el1->gsfTrack()));
      DPtRel1 = fabs( obj.pt() - el1->pt()) /el1->pt();
      dR2 = deltaR(obj, *(el2->gsfTrack()));
      DPtRel2 = fabs( obj.pt() - el2->pt()) /el2->pt();

      for (unsigned int i = 0; i < TestFilterNames_.size(); i++) {
         for (int version = 1; version < 9; version++) {
            std::stringstream ss;
            ss << TestFilterNames_[i] << "_v" << version;
            bool hasPath = obj.hasPathName(ss.str() ,true,true); 
            if (hasPath){
		std::cout << ss.str() << std::endl;  
  		if( (dR1 < 0.1  || dR2 < 0.1 ) && (DPtRel1 < 0.5 ||  DPtRel2 < 0.5) ) matched += (1<<i); 
	       }
            }
        }
    
   }
   return matched;  

}



// ------------ method called once each job just before starting event loop  ------------
//void ZRootupler::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
//void ZRootupler::endJob() {}

// ------------ method called when starting to processes a run  ------------
//void ZRootupler::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
//    prop1_.init(iSetup);
//    prop2_.init(iSetup);
//}

// ------------ method called when ending the processing of a run  ------------
//void ZRootupler::endRun(edm::Run const &, edm::EventSetup const &) {}

// ------------ method called when starting to processes a luminosity block  ------------
//void ZRootupler::beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method called when ending the processing of a luminosity block  ------------
//void ZRootupler::endLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Zjpsi_eeMCTupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Zjpsi_eeMCTupler);
