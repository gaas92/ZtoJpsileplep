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
        '/store/mc/RunIIAutumn18MiniAOD/ZZTo4L_TuneCP5_DoubleScattering_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v5/00000/56984832-8D15-3A4D-8C03-F145A5DD2627.root' ### 2018MC
   )
)

#new sample Dataset for 4l: /ZZTo4L_TuneCP5_DoubleScattering_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v5/MINIAODSIM
process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Z4l_ee_MC18.root'),
)

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
runEnergyCorrections=False, #corrections by default are fine so no need to re-run
era='2018-Prompt')


###jspsi sequence
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('slimmedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.55 < mass && mass < 3.65") ## linea 149
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

###Lepton filter
process.dileptonFilter = cms.EDProducer('LeptonMcFilterElectron',
    electronsMiniAOD        = cms.InputTag("slimmedElectrons"),
    primaryVertices         = cms.InputTag('offlineSlimmedPrimaryVertices')
)

####Kinematic Fit
process.Zfitter    = cms.EDProducer("jpsiElec4l_KmcFitter",
                          dimuon = cms.InputTag("onia2MuMuPAT"),
                          dilepton = cms.InputTag("dileptonFilter","dilepton"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          GenParticles   = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles")
)

#### Z tupler
process.oniarootupler = cms.EDAnalyzer('Zjpsi_eeMCTupler',
                          Zcand     = cms.InputTag("Zfitter","ZCandidates"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT")
           )


process.oniaSequence = cms.Sequence(process.onia2MuMuPAT) ##No trigger matching for 2017yet
process.leptonSequence = cms.Sequence(process.dileptonFilter*process.Zfitter)

process.p = cms.Path(process.egammaPostRecoSeq*process.oniaSequence*process.leptonSequence*process.oniarootupler)
