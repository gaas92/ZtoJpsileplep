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
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 40
	config.Data.publication = True
        config.Data.outputDatasetTag = None
	config.Data.outLFNDirBase = '/store/user/%s/config.General.requestName/' % ("gayalasa")
	config.Site.storageSite = 'T3_US_FNALLPC'
        #config.Site.storageSite = None # Choose your site. 
        #--------------------------------------------------------
        
        # New stuff
        config.Data.userInputFiles = open('my_gen_zjpmm_files.txt').readlines()
        config.General.requestName = "MyGenZJPMM_2018_v7_reqname"
        config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
        config.Data.outputPrimaryDataset = 'crab_'+'config.General.requestName'
        # Submit.
        try:
            print "Submitting for input dataset " 
            crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
        except HTTPException as hte:
            print "Submission for input dataset failed: %s" % (hte.headers)
        except ClientException as cle:
            print "Submission for input dataset failed: %s" % (cle)

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
