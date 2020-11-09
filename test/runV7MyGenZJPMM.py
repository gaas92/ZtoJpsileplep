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
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_1.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_2.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_3.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_4.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_5.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_6.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_7.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_8.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_9.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_10.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_11.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_12.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_13.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_14.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_15.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_16.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_17.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_18.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_19.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_20.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_21.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_22.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_23.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_24.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_25.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_26.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_27.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_28.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_29.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_30.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_31.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_32.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_33.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_34.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_35.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_36.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_37.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_38.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_39.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_40.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_41.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_42.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_43.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_44.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_45.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_46.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_47.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_48.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_49.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_50.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_51.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_52.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_53.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_54.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_55.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_57.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_58.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_59.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_60.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_61.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_62.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_63.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_64.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_65.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_66.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_67.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_68.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_69.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_70.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_71.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_72.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_73.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_74.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_75.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_76.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_77.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_78.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_79.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_80.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_81.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_82.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_83.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_84.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_85.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_86.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_87.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_88.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_89.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_90.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_91.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_92.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_93.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_94.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_95.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_96.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_97.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_98.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_99.root',
        'root://cmseos.fnal.gov//eos/uscms/store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-08-19-21/201108_182128/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_100.root'
       #'/store/mc/RunIIAutumn18MiniAOD/ZToJPsiMuMu_TuneCP5_13TeV-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/250000/C3523B94-D36A-ED4B-AE7C-3B3D14BF9A71.root'   ##mc 18 mumu
       #'/store/data/Run2016D/DoubleMuon/MINIAOD/03Feb2017-v1/1110000/E2C423A9-E3F1-E611-8EA4-047D7BD6DE30.root'
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

"""
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
"""