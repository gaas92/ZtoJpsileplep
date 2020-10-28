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
    fileNames = cms.untracked.vstring(
       '/store/data/Run2016G/DoubleMuon/MINIAOD/17Jul2018-v1/50000/7A20CF7C-4D8C-E811-99B9-0242AC1C0501.root'
       #'/store/mc/RunIIAutumn18MiniAOD/ZToJPsiMuMu_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/250000/C3523B94-D36A-ED4B-AE7C-3B3D14BF9A71.root'   ##mc 18 mumu
        #'/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/1110000/E2C423A9-E3F1-E611-8EA4-047D7BD6DE30.root'
   )

)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Zjpsi_mumu_v7.root'),
)



process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

#Trigger unpacker para poder hacer el trigger Id por objeto
#process.load("AnalyzeZll.OniaRootupler.mySlimmedMuonsTriggerUnpacker")
#Ignore Trigger Selection

#update for v7
####this apply de BF cuts to dimuon & dilepton 1, 3, 4 & 7
####For electros use PATElectronSelector
process.muonFilter = cms.EDFilter('PATMuonSelector',
   src = cms.InputTag('slimmedMuons'),
   cut = cms.string(
                    #for test comment
                    '(pfIsolationR03().sumChargedHadronPt + max(0., pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt-0.5*pfIsolationR03().sumPUPt))/pt() < 0.6'
                    #' && innerTrack.hitPattern.trackerLayersWithMeasurement > 4'
                    #' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    #' && innerTrack.quality(\"highPurity\") '
                    ' && abs(eta) <= 2.5 && pt >= 1'),
   filter = cms.bool(True)
)
process.load("HeavyFlavorAnalysis.Onia2MuMu.onia2MuMuPAT_cfi")
#process.onia2MuMuPAT.muons=cms.InputTag('slimmedMuons')
process.onia2MuMuPAT.muons=cms.InputTag('muonFilter')
process.onia2MuMuPAT.primaryVertexTag=cms.InputTag('offlineSlimmedPrimaryVertices')
process.onia2MuMuPAT.beamSpotTag=cms.InputTag('offlineBeamSpot')
process.onia2MuMuPAT.higherPuritySelection=cms.string("isGlobalMuon")
process.onia2MuMuPAT.lowerPuritySelection=cms.string("isGlobalMuon")
#process.onia2MuMuPAT.dimuonSelection=cms.string("2.55 < mass && mass < 3.65") ## linea 149
process.onia2MuMuPAT.dimuonSelection=cms.string("2.5 < mass && mass < 80.0") ## linea 149
process.onia2MuMuPAT.addMCTruth = cms.bool(False)



#update for v7
####this apply de BF cuts to dimuon & dilepton 5 & 6
process.Zfitter    = cms.EDProducer("jpsiLepLepKmcFitter",
                          dimuon = cms.InputTag("onia2MuMuPAT"),
			                 leptons = cms.InputTag("muonFilter"),
                          primaryVertices     = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          GenParticles        = cms.InputTag("prunedGenParticles"),
                          packedGenParticles  = cms.InputTag("packedGenParticles"),
                          isMC4l              = cms.bool(False),
                          ImparSigm           = cms.double(4.5),
                          ImparSigl           = cms.double(4.5),
                          dxym                = cms.double(0.5),
                          dxyl                = cms.double(0.5),
                          dzm                 = cms.double(1.0),
                          dzl                 = cms.double(1.0),
                          trackerLayersWithMeasurement = cms.double(4),
                          pixelLayersWithMeasurement   = cms.double(1)
)

process.oniarootupler = cms.EDAnalyzer('ZjpsiMCTupler',
#process.oniarootupler = cms.EDAnalyzer('ZRootupler',
                          Zcand	 = cms.InputTag("Zfitter","ZCandidates"),
                          TriggerResults  = cms.InputTag("TriggerResults", "", "HLT")
           )



process.oniaSequence = cms.Sequence(process.onia2MuMuPAT) ##No trigger matching for 2017yet
process.leptonSequence = cms.Sequence(process.Zfitter)

process.p = cms.Path(process.muonFilter*process.oniaSequence*process.leptonSequence*process.oniarootupler)
