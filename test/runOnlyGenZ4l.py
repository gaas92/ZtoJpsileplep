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
        fileName = cms.string('Z4l_only_gen.root'),
)



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))


process.oniarootupler = cms.EDAnalyzer('Z4l_onlyMC_rec',
                        GenParticles = cms.InputTag("prunedGenParticles"),
                        packedGenParticles = cms.InputTag("packedGenParticles")
           )

process.p = cms.Path(process.Zfitter)

