#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException
from CRABClient.UserUtilities import config

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.requestName = None
        #config.General.workArea = 'ZMuondecay'
        config.General.workArea = 'MyGenZJPMM_2018_v7'
	config.General.transferOutputs = True
	config.General.transferLogs = False

        config.JobType.pluginName = 'Analysis'
	config.JobType.psetName = '/afs/cern.ch/work/g/gayalasa/public/Zll/z_sl7/ZtoJpsill/CMSSW_10_6_3/src/AnalizeZll/ZtoJpsileplep/test/runV7MyGenZJPMM.py' #2018 DT configfile
	config.JobType.allowUndistributedCMSSW = True

        config.Data.inputDataset = None
	config.Data.inputDBS = 'global'
   #     config.Data.splitting = 'Automatic'
        config.Data.splitting = 'LumiBased'
   #     config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 30
   #     config.Data.unitsPerJob = 1
   #     config.Data.totalUnits = 30
	#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON_MuonPhys.txt' #has nosence in Mc
	config.Data.publication = True
        config.Data.outputDatasetTag = None
	config.Data.outLFNDirBase = '/store/user/%s/MyGenZJPMM_2018_v7/' % ("gayalasa")
	config.Site.storageSite = 'T3_US_FNALLPC'
        #config.Site.storageSite = None # Choose your site. 
        #--------------------------------------------------------
        
        # Will submit one task for each of these input datasets.
        inputDatasets = [ 
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_100.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_101.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_102.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_103.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_104.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_105.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_106.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_107.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_108.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_109.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_10.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_110.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_111.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_112.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_113.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_114.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_115.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_116.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_117.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_118.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_119.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_11.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_120.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_121.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_122.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_123.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_124.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_125.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_126.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_127.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_128.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_129.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_12.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_130.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_131.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_132.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_133.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_134.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_135.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_136.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_137.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_138.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_139.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_13.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_140.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_141.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_142.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_143.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_144.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_145.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_146.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_147.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_148.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_149.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_14.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_150.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_151.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_153.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_154.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_155.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_156.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_157.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_158.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_159.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_15.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_160.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_161.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_162.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_163.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_165.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_166.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_167.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_168.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_169.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_16.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_170.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_171.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_172.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_173.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_174.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_175.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_176.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_177.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_178.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_179.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_17.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_180.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_181.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_182.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_183.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_185.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_186.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_187.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_188.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_189.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_18.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_190.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_191.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_192.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_193.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_194.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_195.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_196.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_197.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_198.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_199.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_19.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_1.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_200.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_201.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_202.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_203.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_204.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_206.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_207.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_208.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_209.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_20.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_210.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_211.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_212.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_213.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_214.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_215.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_216.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_217.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_218.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_219.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_220.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_221.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_222.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_223.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_224.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_225.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_226.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_227.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_228.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_229.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_22.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_230.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_231.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_232.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_233.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_234.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_235.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_236.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_237.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_238.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_239.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_23.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_240.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_241.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_242.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_243.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_245.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_246.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_247.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_24.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_250.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_251.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_252.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_253.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_254.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_255.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_256.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_257.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_258.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_259.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_25.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_260.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_261.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_262.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_263.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_264.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_265.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_266.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_267.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_268.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_269.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_26.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_270.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_271.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_272.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_273.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_274.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_275.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_276.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_277.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_278.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_279.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_27.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_280.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_281.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_283.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_284.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_285.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_286.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_287.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_288.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_289.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_28.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_290.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_291.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_292.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_293.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_294.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_295.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_296.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_297.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_298.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_299.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_29.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_2.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_300.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_301.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_302.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_303.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_304.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_305.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_306.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_307.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_308.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_309.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_30.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_310.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_311.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_312.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_313.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_314.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_315.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_316.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_317.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_318.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_319.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_320.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_321.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_323.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_325.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_326.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_327.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_329.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_32.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_330.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_331.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_332.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_333.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_334.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_335.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_336.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_337.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_338.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_339.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_33.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_340.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_341.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_342.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_344.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_345.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_346.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_347.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_348.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_349.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_34.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_350.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_351.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_352.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_353.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_354.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_355.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_356.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_357.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_358.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_359.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_35.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_360.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_362.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_363.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_364.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_365.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_366.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_367.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_368.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_369.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_36.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_370.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_371.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_372.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_374.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_375.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_376.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_377.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_378.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_379.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_37.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_380.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_381.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_382.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_383.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_384.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_385.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_386.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_387.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_388.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_389.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_38.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_390.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_391.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_392.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_393.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_394.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_395.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_396.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_397.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_398.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_399.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_39.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_3.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_400.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_40.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_41.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_42.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_43.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_44.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_45.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_46.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_47.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_48.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_49.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_4.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_50.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_51.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_52.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_53.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_54.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_56.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_57.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_58.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_59.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_5.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_60.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_61.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_62.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_63.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_64.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_65.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_66.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_67.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_68.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_69.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_6.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_70.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_71.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_72.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_73.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_74.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_75.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_76.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_77.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_79.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_7.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_80.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_81.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_82.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_84.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_85.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_86.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_87.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_88.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_8.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_90.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_91.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_92.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_93.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_94.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_95.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_96.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_97.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_98.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_99.root",
        "root://cmsxrootd.fnal.gov//store/user/gabaasz/Zpsi_mm_MINIAODSIM-BPH-FW/PrivateMC-2018-ZtoJpsiMuMu/crab_PrivateMC-2018-ZtoJpsiMuMu-2020-11-10-06-50/201110_055117/0000/step3-MINIAODSIM-ZtoJpsiMuMu_BPH-FW-result_9.root"
                 	]
 
        for inDS in inputDatasets:
             # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            #config.General.requestName = (inDS.split('/')[2]).split('-')[0]+(inDS.split('/')[2]).split('-')[-1]
            config.General.requestName = 'job'+(inDS.split('/')[-1]).split('_')[.1].replace('.root', '')
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
            # Submit.
            try:
                print "Submitting for input dataset %s" % (inDS)
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Submission for input dataset %s failed: %s" % (inDS, hte.headers)
            except ClientException as cle:
                print "Submission for input dataset %s failed: %s" % (inDS, cle)

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print "-"*len(msg)
            print msg
            print "-"*len(msg)
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
            except ClientException as cle:
                print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)


if __name__ == '__main__':
    main()
