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

from multiprocessing import Pool

from InitialStateAnalysis.Utilities.utilities import *
from InitialStateAnalysis.Analyzers.AnalyzerZ import AnalyzerZ
from InitialStateAnalysis.Analyzers.AnalyzerWZ import AnalyzerWZ, AnalyzerWZ_Z
from InitialStateAnalysis.Analyzers.AnalyzerWZ_W import AnalyzerWZ_W
from InitialStateAnalysis.Analyzers.AnalyzerWZ_FakeRate import AnalyzerWZ_FakeRate
from InitialStateAnalysis.Analyzers.AnalyzerHpp2l import AnalyzerHpp2l, AnalyzerHpp2l_Z, AnalyzerHpp2l_TT
from InitialStateAnalysis.Analyzers.AnalyzerHpp3l import AnalyzerHpp3l, AnalyzerHpp3l_WZ, AnalyzerHpp3l_LowMass
from InitialStateAnalysis.Analyzers.AnalyzerHpp4l import AnalyzerHpp4l

def run_analyzer(args):
    '''Run the analysis'''
    analysis, channel, sample_name, filelist, outfile, period, loglevel = args
    analyzerMap = {
        'Z'       : {
                    'Z'       : AnalyzerZ,
                    },
        'Hpp2l'   : {
                    'Hpp2l'   : AnalyzerHpp2l,
                    'Z'       : AnalyzerHpp2l_Z,
                    'TT'      : AnalyzerHpp2l_TT,
                    },
        'WZ'      : {
                    'WZ'      : AnalyzerWZ,
                    'Z'       : AnalyzerWZ_Z,
                     },
        'WZ_W'    : {
                    'W'       : AnalyzerWZ_W,
                    },
        'WZ_FakeRate' : {
                    'FakeRate': AnalyzerWZ_FakeRate,
                    },
        'Hpp3l'   : {
                    'Hpp3l'   : AnalyzerHpp3l,
                    'WZ'      : AnalyzerHpp3l_WZ,
                    'LowMass' : AnalyzerHpp3l_LowMass,
                    },
        'Hpp4l'   : {
                    'Hpp4l'   : AnalyzerHpp4l,
                    },
    }
    theAnalyzer = analyzerMap[analysis][channel]
    with theAnalyzer(sample_name,filelist,outfile,period,loglevel=loglevel) as analyzer:
        analyzer.analyze()

def get_sample_names(analysis,period,samples):
    '''Get unix sample names'''
    ntupleDict = {
        '8': {
            'Z'          : '2015-06-01-8TeV-2l',
            'TT'         : '2015-06-01-8TeV-2l',
            'Hpp2l'      : '2015-06-01-8TeV-2l',
            'WZ'         : '2015-06-01-8TeV', 
            'WZ_W'       : 'N/A',
            'WZ_FakeRate': 'N/A',
            'Hpp3l'      : '2015-06-01-8TeV',
            'Hpp4l'      : '2015-08-26-8TeV-4l', 
        },
        '13': {
            'Z'          : 'N/A',
            'TT'         : 'N/A',
            'Hpp2l'      : 'N/A',
            #'WZ'         : '2015-08-03-13TeV-WZ', # last 50ns with old CBIDs
            #'WZ'         : '2015-08-27-13TeV-WZ', # add new egamma CBIDs
            'WZ'         : '2015-08-31-13TeV-WZ', # remove dr 0.1 from veto definition
            'WZ_W'       : '2015-08-03-13TeV-2l',
            'WZ_FakeRate': '2015-08-17-13TeV-1l',
            'Hpp3l'      : '2015-03-30-13TeV-3l',
            'Hpp4l'      : '2015-03-30-13TeV-4l',
        },
    }
    root_dir = '/hdfs/store/user/dntaylor/data/%s' % ntupleDict[period][analysis]

    sample_names = [os.path.basename(fname)
                    for string in samples
                    for fname in glob.glob("%s/%s" % (root_dir, string))]

    return root_dir, sample_names

def run_ntuples(analysis, channel, period, samples, loglevel):
    '''Run a given analyzer for the H++ analysis'''
    logger = logging.getLogger(__name__)
    ntup_dir = './ntuples/%s_%sTeV_%s' % (analysis, period, channel)
    python_mkdir(ntup_dir)
    root_dir, sample_names = get_sample_names(analysis,period,samples)


    filelists = {}
    for sample in sample_names:
        sampledir = '%s/%s' % (root_dir, sample)
        filelists[sample] = ['%s/%s' % (sampledir, x) for x in os.listdir(sampledir)]

    if len(sample_names)==1: # only one, its a test, dont use map
        name = sample_names[0]
        run_analyzer((analysis, channel, name, filelists[name], "%s/%s.root" % (ntup_dir, name), period, loglevel))
        return 0

    p = Pool(8)
    try:
        p.map_async(run_analyzer, [(analysis, channel, name, filelists[name], "%s/%s.root" % (ntup_dir, name), period, loglevel) for name in sample_names]).get(999999)
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
    input_name = '%s/%s.txt' % (dag_dir+'inputs', sample_name)
    with open(input_name,'w') as file:
        for f in filelist:
            file.write('%s\n' % f.replace('/hdfs','',1))

    # create bash script
    bash_name = '%s/%s_%s_%s_%s.sh' % (dag_dir+'inputs', analysis, channel, period, sample_name)
    bashScript = '#!/bin/bash\npython $CMSSW_BASE/src/InitialStateAnalysis/Analyzers/python/Analyzer%s.py %s %s $INPUT $OUTPUT %s\n' % (analysis, channel, sample_name, period)
    with open(bash_name,'w') as file:
        file.write(bashScript)
    os.system('chmod +x %s' % bash_name)

    # create farmout command
    farmoutString = 'farmoutAnalysisJobs --infer-cmssw-path --fwklite --input-file-list=%s' % (input_name)
    farmoutString += ' --submit-dir=%s --output-dag-file=%s --output-dir=%s' % (submit_dir, dag_dir, output_dir)
    if period == '8':
        farmoutString += ' --input-files-per-job=20 %s %s' % (jobName, bash_name)
    else:
        farmoutString += ' --input-files-per-job=10 %s %s' % (jobName, bash_name)

    logger.info('Submitting %s' % sample_name)
    os.system(farmoutString)

    return

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Run the desired analyzer on "
                                                 "FSA n-tuples")

    parser.add_argument('analysis', type=str, choices=['Z','WZ','WZ_W','WZ_FakeRate','Hpp2l','Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('channel', type=str, choices=['Z','WZ','W','FakeRate','TT','Hpp2l','Hpp3l','Hpp4l','LowMass'], help='Channel to run for given analysis')
    parser.add_argument('period', type=str, choices=['8','13'], help='Energy (TeV)')
    parser.add_argument('sample_names', nargs='+',help='Sample names w/ UNIX wildcards')
    parser.add_argument('-s','--submit',action='store_true',help='Submit jobs to condor')
    parser.add_argument('-jn','--jobName',nargs='?',type=str,const='',help='Job Name for condor submission')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    if args.period == '7':
        logger.warning("7 TeV not implemented")
    else:
        logger.info("Running %s:%s %s TeV analyzer" %(args.analysis, args.channel, args.period))
        if args.submit:
            root_dir, sample_names = get_sample_names(args.analysis, args.period, args.sample_names)
            for sample in sample_names:
                sampledir = '%s/%s' % (root_dir, sample)
                submitFwkliteJob(sampledir,args)
        else:
            run_ntuples(args.analysis, args.channel, args.period, args.sample_names, args.log)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)