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
        config.General.workArea = '2018MC_psi_mm_v7_3'
	config.General.transferOutputs = True
	config.General.transferLogs = False

        config.JobType.pluginName = 'Analysis'
	config.JobType.psetName = '/afs/cern.ch/work/g/gayalasa/public/Zll/z_sl7/ZtoJpsill/CMSSW_10_6_3/src/AnalizeZll/ZtoJpsileplep/test/runV7Muon.py' #2018 DT configfile
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
	config.Data.outLFNDirBase = '/store/user/%s/Zpsi_mm18mc_v7_3/' % ("gayalasa")
	config.Site.storageSite = 'T3_US_FNALLPC'
        #config.Site.storageSite = None # Choose your site. 
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [ 
                          '/ZToJPsiMuMu_TuneCP5_13TeV-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
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
