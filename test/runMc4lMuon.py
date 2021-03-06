import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootupler")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/mc/RunIIAutumn18MiniAOD/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext1-v2/10000/4E82A14A-00E7-2A48-A9A7-C2DE26FE6F7A.root',
    #'/store/mc/RunIIAutumn18MiniAOD/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15_ext2-v2/20000/EB08DA54-6990-F74F-8D3D-96B606DE386F.root'
     '/store/mc/RunIISummer16MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/769AED45-70C7-E611-B925-001E67504445.root' #stefanos
   )

)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Z4l_mm_MC18.root'),
)



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

#Trigger unpacker para poder hacer el trigger Id por objeto
#process.load("AnalyzeZll.OniaRootupler.mySlimmedMuonsTriggerUnpacker")
#Ignore Trigger Selection

process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('slimmedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("isGlobalMuon")
#process.onia2MuMuPAT.dimuonSelection=cms.string("2.55 < mass && mass < 3.65") ## linea 149
process.onia2MuMuPAT.dimuonSelection=cms.string("0.01 < mass && mass < 100.0") ## linea 149
process.onia2MuMuPAT.addMCTruth = cms.bool(False)


process.dileptonFilter = cms.EDProducer('LeptonMcFilter',
    muons     	    = cms.InputTag("slimmedMuons"),
    primaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices')
)

#Corta con abs(tkPVdist.second.significance())>4. para ambos dileptones
#Corta con dR1, dR2, dR3 y dR4 < 0.2, podemos subir en Onia...cc
#Fit Cinematico
process.Zfitter    = cms.EDProducer("jpsi4LepLepKmcFitter",
        dimuon            = cms.InputTag("onia2MuMuPAT"),
        dilepton          = cms.InputTag("dileptonFilter","dilepton"),
        #muons             = cms.InputTag("slimmedMuons"),
        primaryVertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
        GenParticles      = cms.InputTag("prunedGenParticles"),
        packedGenParticles = cms.InputTag("packedGenParticles")
)

process.oniarootupler = cms.EDAnalyzer('ZjpsiMCTupler',
#process.oniarootupler = cms.EDAnalyzer('ZRootupler',
                          Zcand	 = cms.InputTag("Zfitter","ZCandidates"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
                          isMC 		  = cms.bool(False),
                          OnlyBest 	  = cms.bool(False), 
			  onia_pdgid	  = cms.uint32(443),
                          OnlyGen  	  = cms.bool(False)
           )



process.oniaSequence = cms.Sequence(process.onia2MuMuPAT) 
process.leptonSequence = cms.Sequence(process.dileptonFilter*process.Zfitter)

process.p = cms.Path(process.oniaSequence*process.leptonSequence*process.oniarootupler)

