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
   #     config.Data.splitting = 'LumiBased'
        config.Data.splitting = 'FileBased'
   #     config.Data.unitsPerJob = 30
        config.Data.unitsPerJob = 1
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
                 	]
 
        for inDS in inputDatasets:
             # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.General.requestName = (inDS.split('/')[2]).split('-')[0]+(inDS.split('/')[2]).split('-')[-1]
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
