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
        #'/store/mc/RunIIAutumn18MiniAOD/ZZTo4L_TuneCP5_DoubleScattering_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v5/00000/56984832-8D15-3A4D-8C03-F145A5DD2627.root' ### 2018MC
        '/store/mc/RunIIAutumn18MiniAOD/ZZTo4L_TuneCP5_DoubleScattering_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v5/00000/6072A673-CF57-DB43-8744-6D0D46E0D670.root'
   )
)

#new sample Dataset for 4l: /ZZTo4L_TuneCP5_DoubleScattering_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v5/MINIAODSIM
process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Z4l_ee_MC18.root'),
)

#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
runEnergyCorrections=False, #corrections by default are fine so no need to re-run
era='2018-Prompt')


####Kinematic Fit
process.Zfitter    = cms.EDProducer("jpsiElec4l_KmcFitter",
                          electronsMiniAOD = cms.InputTag("slimmedElectrons"),
                          muons             = cms.InputTag("slimmedMuons"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          GenParticles   = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles")
)

#### Z tupler
process.oniarootupler = cms.EDAnalyzer('Zjpsi_eeMCTupler',
                          Zcand     = cms.InputTag("Zfitter","ZCandidates"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT")
           )


process.leptonSequence = cms.Sequence(process.Zfitter)

process.p = cms.Path(process.egammaPostRecoSeq*process.leptonSequence*process.oniarootupler)
