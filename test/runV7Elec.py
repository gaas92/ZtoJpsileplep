import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.fixedGridRhoProducerFastjet_cfi import fixedGridRhoFastjetAll

process = cms.Process("Rootupler")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'/store/data/Run2016G/DoubleMuon/MINIAOD/17Jul2018-v1/50000/7A20CF7C-4D8C-E811-99B9-0242AC1C0501.root'
        '/store/data/Run2018C/MuonEG/MINIAOD/17Sep2018-v1/80000/ACED0351-B365-4E4D-8A4A-B23DD48A4081.root' ### 2018MC
        #'/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/1110000/E2C423A9-E3F1-E611-8EA4-047D7BD6DE30.root'
   )

)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Zjpsi_elel_v7.root'),
)
process.SimpleMemoryCheck = cms.Service('SimpleMemoryCheck',
                                        ignoreTotal=cms.untracked.int32(0),
                                        oncePerEventMode = cms.untracked.bool(False)
                                        )

process.Timing = cms.Service("Timing")


#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
runEnergyCorrections=False, #corrections by default are fine so no need to re-run
era='2018-Prompt')

#update for v7
####this apply de BF cuts to dimuon & dilepton 1, 3, 4 & 7
####For electros use PATElectronSelector
process.muonFilter = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuons'),
   cut = cms.string(
                    '(pfIsolationR03().sumChargedHadronPt + max(0., pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt-0.5*pfIsolationR03().sumPUPt))/pt() < 0.6'
                    #' && innerTrack.hitPattern.trackerLayersWithMeasurement > 4'
                    #' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    #' && innerTrack.quality(\"highPurity\") '
                    ' && abs(eta) <= 2.5 && pt >= 1'),
   filter = cms.bool(True)
)
####
process.electonFilter = cms.EDFilter('PATElectronSelector',
   src = cms.InputTag('slimmedElectrons'),
   cut = cms.string(
                    '(pfIsolationR03().sumChargedHadronPt + max(0., pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt-0.5*pfIsolationR03().sumPUPt))/pt() < 0.6'
                    #' && innerTrack.hitPattern.trackerLayersWithMeasurement > 4'
                    #' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    #' && innerTrack.quality(\"highPurity\") '
                    ' && abs(eta) <= 2.5 && pt >= 1'),
   filter = cms.bool(True)
)
###jspsi sequence
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('slimmedMuons')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.0 < mass && mass < 100") ## linea 149
process.onia2MuMuPAT.addMCTruth = cms.bool(False)

###Lepton filter
process.dileptonFilter = cms.EDProducer('LeptonMcFilterElectron',
    electronsMiniAOD        = cms.InputTag("slimmedElectrons"),
    primaryVertices         = cms.InputTag('offlineSlimmedPrimaryVertices')
)

####Kinematic Fit
process.Zfitter    = cms.EDProducer("jpsiElecKmcFitter",
                          dimuon = cms.InputTag("onia2MuMuPAT"),
                          dilepton = cms.InputTag("dileptonFilter","dilepton"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          GenParticles   = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles")
                          
#                          fixedGridRhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll")
)

#### Z tupler
process.oniarootupler = cms.EDAnalyzer('Zjpsi_eeMCTupler',
                          Zcand     = cms.InputTag("Zfitter","ZCandidates"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT")
           )


process.oniaSequence = cms.Sequence(process.onia2MuMuPAT) ##No trigger matching for 2017yet
process.leptonSequence = cms.Sequence(process.dileptonFilter*process.Zfitter)

process.p = cms.Path(process.egammaPostRecoSeq*process.muonFilter*process.electonFilter*process.oniaSequence*process.leptonSequence*process.oniarootupler)
