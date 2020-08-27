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
     '/store/mc/RunIISummer16MiniAODv2/ZToJPsiMuMu_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/363AF3B9-B1DE-E611-BE16-02163E017665.root' #jpsi mumu
   )

)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Zjpsi_only_gen.root'),
)



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))


process.Zfitter = cms.EDAnalyzer('Zjpsi_onlyMC_rec',
                        GenParticles = cms.InputTag("prunedGenParticles"),
                        packedGenParticles = cms.InputTag("packedGenParticles")
           )

process.p = cms.Path(process.Zfitter)

