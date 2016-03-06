#!/usr/bin/env python
'''
A script to run the InitialStateAnalysis analyzers on ntuples output from FinalStateAnalysis.

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import glob
import pwd
import argparse
import errno
import socket
import signal
import logging
import math

from multiprocessing import Pool

from InitialStateAnalysis.Utilities.utilities import *
from InitialStateAnalysis.Analyzers.AnalyzerZ import AnalyzerZ
from InitialStateAnalysis.Analyzers.AnalyzerWZ import AnalyzerWZ, AnalyzerWZ_ZFakeRate, AnalyzerWZ_TTFakeRate
from InitialStateAnalysis.Analyzers.AnalyzerWZ_W import AnalyzerWZ_WFakeRate
from InitialStateAnalysis.Analyzers.AnalyzerWZ_Dijet import AnalyzerWZ_DijetFakeRate
from InitialStateAnalysis.Analyzers.AnalyzerHpp2l import AnalyzerHpp2l, AnalyzerHpp2l_Z, AnalyzerHpp2l_Charge, AnalyzerHpp2l_TT
from InitialStateAnalysis.Analyzers.AnalyzerHpp3l import AnalyzerHpp3l, AnalyzerHpp3l_WZ, AnalyzerHpp3l_LowMass
from InitialStateAnalysis.Analyzers.AnalyzerHpp4l import AnalyzerHpp4l, AnalyzerHpp4l_ZZ

def run_analyzer(args):
    '''Run the analysis'''
    analysis, channel, sample_name, filelist, outfile, period, metShift, loglevel = args
    analyzerMap = {
        'Z'       : {
                    'Z'       : AnalyzerZ,
                    },
        'Hpp2l'   : {
                    'Hpp2l'   : AnalyzerHpp2l,
                    'Z'       : AnalyzerHpp2l_Z,
                    'Charge'  : AnalyzerHpp2l_Charge,
                    'TT'      : AnalyzerHpp2l_TT,
                    },
        'WZ'      : {
                    'WZ'      : AnalyzerWZ,
                    'FakeRate': AnalyzerWZ_ZFakeRate,
                    'TTFakeRate': AnalyzerWZ_TTFakeRate,
                     },
        'WZ_W'    : {
                    'FakeRate': AnalyzerWZ_WFakeRate,
                    },
        'WZ_Dijet': {
                    'FakeRate': AnalyzerWZ_DijetFakeRate,
                    },
        'Hpp3l'   : {
                    'Hpp3l'   : AnalyzerHpp3l,
                    'WZ'      : AnalyzerHpp3l_WZ,
                    'LowMass' : AnalyzerHpp3l_LowMass,
                    },
        'Hpp4l'   : {
                    'Hpp4l'   : AnalyzerHpp4l,
                    'ZZ'      : AnalyzerHpp4l_ZZ,
                    },
    }
    theAnalyzer = analyzerMap[analysis][channel]
    with theAnalyzer(sample_name,filelist,outfile,period,metShift=metShift,loglevel=loglevel) as analyzer:
        analyzer.analyze()

def get_sample_names(analysis,period,samples,**kwargs):
    '''Get unix sample names'''
    customDir = kwargs.pop('customDir','')
    ntupleDict = {
        8: {
            'Z'          : '2015-06-01-8TeV-2l',
            'Charge'     : '2015-06-01-8TeV-2l',
            'TT'         : '2015-06-01-8TeV-2l',
            'Hpp2l'      : '2015-06-01-8TeV-2l',
            #'WZ'         : '2015-06-01-8TeV', 
            'WZ'         : '2015-10-17-8TeV-3l', # add fakerate no iso variables 
            'WZ_W'       : 'N/A',
            'WZ_Dijet'   : 'N/A',
            #'Hpp3l'      : '2015-06-01-8TeV',
            'Hpp3l'      : '2015-10-17-8TeV-3l', # add fakerate no iso variables
            #'Hpp4l'      : '2015-08-26-8TeV-4l', 
            'Hpp4l'      : '2015-10-17-8TeV-4l', # add fakerate no iso variables
        },
        13: {
            'Z'          : 'N/A',
            'TT'         : 'N/A',
            'Hpp2l'      : 'N/A',
            #'WZ'         : '2015-08-03-13TeV-WZ', # last 50ns with old CBIDs
            #'WZ'         : '2015-08-27-13TeV-WZ', # add new egamma CBIDs
            #'WZ'         : '2015-09-01-13TeV-WZ', # move dr 0.1 to 0.01 from veto definition
            #'WZ'         : '2015-09-13-13TeV-WZ', # corrected isolation to new effective areas
            #'WZ'         : '2015-09-18-13TeV-WZ', # add ht variable
            #'WZ'         : '2015-09-28-13TeV-WZ', # add summed weights
            #'WZ'         : '2015-10-06-13TeV-WZ', # add hzz veto
            #'WZ'         : '2015-10-12-13TeV-WZ', # lower trigger
            #'WZ'         : '2015-10-15-13TeV-WZ', # fixed trigger and add WZ no iso ID
            #'WZ'         : '2015-10-24-13TeV-WZ', # miniaodv2, metfilters, new met uncertainty
            #'WZ'         : '2015-10-25-13TeV-WZ', # all samples and bug fix
            #'WZ'         : '2015-11-06-13TeV-WZ', # latest jec, metfilters, metuncertainty, new samples
            #'WZ'         : '2015-11-07-13TeV-WZ', # update to met uncertainty, also met shifts
            #'WZ'         : '2015-11-12-13TeV-WZ', # add new muon medium ID counts
            #'WZ'         : '2015-11-18-13TeV-WZ', # jet clean based on tight IDs and loose jet id pt>20
            #'WZ'         : '2015-11-27-13TeV-WZ', # fix to event number (no more negatives) and add some met uncertainty stuff
            #'WZ'         : '2016-01-18-13TeV-WZ', # move to medium muon and veryTight electron
            'WZ'         : '2016-01-29-13TeV-WZ', # move to WW ids
            #'WZ_W'       : '2015-08-03-13TeV-2l',
            #'WZ_W'       : '2015-11-19-13TeV-2l', # all the udpates above
            'WZ_W'       : '2016-01-30-13TeV-2l', # move to WW ids
            #'WZ_Dijet'   : '2015-08-17-13TeV-1l',
            #'WZ_Dijet'   : '2015-09-14-13TeV-1l', # updated with WZ changes
            #'WZ_Dijet'   : '2015-09-24-13TeV-1l', # add ht
            #'WZ_Dijet'   : '2015-09-28-13TeV-1l', # add summed weights
            #'WZ_Dijet'   : '2015-10-06-13TeV-1l', # add hzz veto
            #'WZ_Dijet'   : '2015-10-12-13TeV-1l', # lower trigger
            #'WZ_Dijet'   : '2015-10-15-13TeV-1l', # fixed trigger and WZ no iso ID
            #'WZ_Dijet'   : '2015-11-06-13TeV-1l', # latest jec, metfilters, metuncertainty, new samples
            #'WZ_Dijet'   : '2015-11-19-13TeV-1l', # The jet cleaning fixes
            #'WZ_Dijet'   : '2016-01-18-13TeV-1l', # move to medium muon and veryTight electron
            'WZ_Dijet'   : '2016-01-29-13TeV-1l', # move to WW ids
            'Hpp3l'      : '2015-03-30-13TeV-3l',
            'Hpp4l'      : '2015-03-30-13TeV-4l',
        },
    }
    base_dir = '/hdfs/store/user/dntaylor/data'
    root_dir = '%s/%s' % (base_dir, ntupleDict[period][analysis])

    if customDir: # use a custom directory (or path ind data)
        if os.path.isdir(customDir):
            root_dir = customDir
        elif os.path.isdir('{0}/{1}'.format(base_dir,customDir)):
            root_dir = '{0}/{1}'.format(base_dir,customDir)
        else:
            logging.warning('Warning: directory {0} does not exist, using default'.format(customDir))

    sample_names = [os.path.basename(fname)
                    for string in samples
                    for fname in glob.glob("%s/%s" % (root_dir, string))]

    return root_dir, sample_names

def run_ntuples(analysis, channel, period, samples, loglevel, **kwargs):
    '''Run a given analyzer for the analysis'''
    logger = logging.getLogger(__name__)
    test = kwargs.pop('test',False)
    metShift = kwargs.pop('metShift','')
    ntup_dir = './ntuples/%s_%iTeV_%s' % (analysis, period, channel)
    python_mkdir(ntup_dir)
    root_dir, sample_names = get_sample_names(analysis,period,samples,**kwargs)


    filelists = {}
    for sample in sample_names:
        sampledir = '%s/%s' % (root_dir, sample)
        filelists[sample] = ['%s/%s' % (sampledir, x) for x in os.listdir(sampledir)]

    if len(sample_names)==1 or test: # only one, its a test, dont use map
        name = sample_names[0]
        outname =  "%s/%s.root" % (ntup_dir, name)
        if test: outname = 'test.root'
        run_analyzer((analysis, channel, name, filelists[name], outname, period, metShift, loglevel))
        return 0

    p = Pool(8)
    try:
        p.map_async(run_analyzer, [(analysis, channel, name, filelists[name], "%s/%s.root" % (ntup_dir, name), period, metShift, loglevel) for name in sample_names]).get(999999)
    except KeyboardInterrupt:
        p.terminate()
        logger.info('Analyzer cancelled')
        sys.exit(1)
   
    return 0

def submitFwkliteJob(sampledir,args):
    '''
    Submit a job using farmoutAnalysisJobs --fwklite
    '''
    logger = logging.getLogger(__name__)
    jobName = args.jobName
    analysis = args.analysis
    channel = args.channel
    period = args.period
    dryrun = args.dryrun
    metShift = args.metShift
    sample_name = os.path.basename(sampledir)

    if 'uwlogin' in socket.gethostname():
        sample_dir = '/scratch/%s/%s/%s' % (pwd.getpwuid(os.getuid())[0], jobName, sample_name)
    else:
        sample_dir = '/nfs_scratch/%s/%s/%s' % (pwd.getpwuid(os.getuid())[0], jobName, sample_name)

    # create submit dir
    submit_dir = '%s/submit' % (sample_dir)
    if os.path.exists(submit_dir):
        logger.warning('Submission directory exists for %s %s.' % (jobName, sample_name))
        return

    # create dag dir
    dag_dir = '%s/dags/dag' % (sample_dir)
    os.system('mkdir -p %s' % (os.path.dirname(dag_dir)))
    os.system('mkdir -p %s' % (dag_dir+'inputs'))

    # output dir
    output_dir = 'srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/%s/%s/%s/'\
                 % (pwd.getpwuid(os.getuid())[0], jobName, sample_name)

    # create file list
    filelist = ['%s/%s' % (sampledir, x) for x in os.listdir(sampledir)]
    numfiles = len(filelist)
    totalsize = sum([os.path.getsize(f) for f in filelist])
    averagesize = totalsize/numfiles
    filesperjob = int(math.ceil(50000000./averagesize)) # average them to 50MB per job
    input_name = '%s/%s.txt' % (dag_dir+'inputs', sample_name)
    with open(input_name,'w') as file:
        for f in filelist:
            file.write('%s\n' % f.replace('/hdfs','',1))

    # create bash script
    bash_name = '%s/%s_%s_%i_%s.sh' % (dag_dir+'inputs', analysis, channel, period, sample_name)
    bashScript = '#!/bin/bash\npython $CMSSW_BASE/src/InitialStateAnalysis/Analyzers/python/Analyzer%s.py %s %s $INPUT $OUTPUT %i' % (analysis, channel, sample_name, period)
    if metShift: bashScript += ' --metShift {0}'.format(metShift)
    bashScript += '\n'
    with open(bash_name,'w') as file:
        file.write(bashScript)
    os.system('chmod +x %s' % bash_name)

    # create farmout command
    farmoutString = 'farmoutAnalysisJobs --infer-cmssw-path --assume-input-files-exist --fwklite --input-file-list=%s' % (input_name)
    farmoutString += ' --submit-dir=%s --output-dag-file=%s --output-dir=%s' % (submit_dir, dag_dir, output_dir)
    #if period == 8:
    #    farmoutString += ' --input-files-per-job=20 %s %s' % (jobName, bash_name)
    #else:
    #    farmoutString += ' --input-files-per-job=10 %s %s' % (jobName, bash_name)
    farmoutString += ' --input-files-per-job=%i %s %s' % (filesperjob, jobName, bash_name)

    if not args.dryrun:
        logger.info('Submitting %s' % sample_name)
        os.system(farmoutString)
    else:
        print farmoutString

    return

def parse_command_line(argv):
    parser = get_parser("Run the desired analyzer on FSA n-tuples")

    parser.add_argument('sample_names', nargs='+',help='Sample names w/ UNIX wildcards')
    parser.add_argument('-s','--submit',action='store_true',help='Submit jobs to condor')
    parser.add_argument('-dr','--dryrun',action='store_true',help='Create jobs but dont submit')
    parser.add_argument('-t','--test',action='store_true',help='Do a test (output test.root)')
    parser.add_argument('-jn','--jobName',nargs='?',type=str,const='',help='Job Name for condor submission')
    parser.add_argument('-d','--customDir',nargs='?',type=str,const='',help='Custom input directory')
    parser.add_argument('-ms','--metShift',nargs='?',type=str,const='',help='Shift the met')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    if args.period == 7:
        logger.warning("7 TeV not implemented")
    else:
        logger.info("Running %s:%s %i TeV analyzer" %(args.analysis, args.channel, args.period))
        if args.submit:
            root_dir, sample_names = get_sample_names(args.analysis, args.period, args.sample_names, customDir=args.customDir)
            for sample in sample_names:
                sampledir = '%s/%s' % (root_dir, sample)
                submitFwkliteJob(sampledir,args)
        else:
            run_ntuples(args.analysis, args.channel, args.period, args.sample_names, args.log, customDir=args.customDir, test=args.test, metShift=args.metShift)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
