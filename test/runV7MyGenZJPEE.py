import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootupler")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( open('my_gen-NF_zjpee_18.txt').readlines()
    #fileNames = cms.untracked.vstring(
    #        '/store/mc/RunIISummer16MiniAODv3/ZToJPsiEE_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/40000/542780DC-860A-EA11-9709-0CC47AFF04A4.root',
    #        '/store/mc/RunIISummer16MiniAODv3/ZToJPsiEE_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/40000/88138499-A50A-EA11-8459-0CC47AFF24BA.root',
    #        '/store/mc/RunIISummer16MiniAODv3/ZToJPsiEE_TuneCUEP8M1_13TeV-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/40000/8C2C4791-780A-EA11-93D3-0CC47AF9B1D6.root'
    #
   )

)
"""
process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Zjpsi_only_gen.root'),
)



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))


process.Zfitter = cms.EDAnalyzer('ZJpsi_onlyMC_ee_rec',
                        GenParticles = cms.InputTag("prunedGenParticles"),
                        packedGenParticles = cms.InputTag("packedGenParticles")
           )

process.p = cms.Path(process.Zfitter)

"""
process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Zjpsi_elel_test_v7.root'),
)



#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process, runVID=True, #if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2017-Nov17ReReco')  #era is new to select between 2016 / 2017,  it defaults to 2017

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
                    #' && innerTrack.hitPattern.trackerLayersWithMeasurement > 4'
                    #' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    #' && innerTrack.quality(\"highPurity\") '
                    ' abs(eta) <= 2.6 && pt >= 1'),
   filter = cms.bool(True)
)
###jspsi sequence
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
process.onia2MuMuPAT.muons=cms.InputTag('muonFilter')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.dimuonSelection=cms.string("2.5 < mass && mass < 80.0") ## linea 149
process.onia2MuMuPAT.addMCTruth = cms.bool(False)


#update for v7
####this apply de BF cuts to dimuon & dilepton 5 & 6
process.Zfitter    = cms.EDProducer("jpsiElecKmcFitter",
                          dimuon = cms.InputTag("onia2MuMuPAT"),
                          leptons = cms.InputTag("electonFilter"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          GenParticles   = cms.InputTag("prunedGenParticles"),
                          packedGenParticles = cms.InputTag("packedGenParticles"),
                          isMC4l              = cms.bool(False),
                          ImparSigm           = cms.double(4.5),
                          ImparSigl           = cms.double(4.5),
                          dxym                = cms.double(0.5),
                          dxyl                = cms.double(0.5),
                          dzm                 = cms.double(1.0),
                          dzl                 = cms.double(1.0),
                          trackerLayersWithMeasurement = cms.double(4),
                          pixelLayersWithMeasurement   = cms.double(1)
#                          fixedGridRhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll")
)

#### Z tupler
process.oniarootupler = cms.EDAnalyzer('Zjpsi_eeMCTupler',
                          Zcand     = cms.InputTag("Zfitter","ZCandidates"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT")
           )


process.oniaSequence = cms.Sequence(process.onia2MuMuPAT) ##No trigger matching for 2017yet
process.leptonSequence = cms.Sequence(process.Zfitter)

process.p = cms.Path(process.egammaPostRecoSeq*process.muonFilter*process.electonFilter*process.oniaSequence*process.leptonSequence*process.oniarootupler)
