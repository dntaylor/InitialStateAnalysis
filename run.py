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

from multiprocessing import Pool

from utilities.utilities import *
from analyzers.AnalyzerWZ import AnalyzerWZ, AnalyzerTT
from analyzers.AnalyzerHpp3l import AnalyzerHpp3l, AnalyzerFakeRate
from analyzers.AnalyzerHpp4l import AnalyzerHpp4l

def run_analyzer(args):
    '''Run the analysis'''
    analysis, channel, location, outfile, period = args
    analyzerMap = {
        'WZ'      : AnalyzerWZ,
        'TT'      : AnalyzerTT,
        'Hpp3l'   : AnalyzerHpp3l,
        'Hpp4l'   : AnalyzerHpp4l,
        'FakeRate': AnalyzerFakeRate,
    }
    theAnalyzer = analyzerMap[channel]
    with theAnalyzer(location,outfile,period) as analyzer:
        analyzer.analyze()

def get_sample_names(analysis,period,samples):
    '''Get unix sample names'''
    ntupleDict = {
        '7': {
            'WZ'   : 'N/A',
            'Hpp3l': 'N/A',
            'Hpp4l': 'N/A', 
        },
        '8': {
            'WZ'   : '2015-04-19-8TeV', 
            'Hpp3l': '2014-12-04-8TeV',
            'Hpp4l': 'N/A', 
        },
        '13': {
            'WZ'   : '2015-04-28-13TeV',
            'Hpp3l': '2015-03-30-13TeV-3l',
            'Hpp4l': '2015-03-30-13TeV-4l',
        },
    }
    root_dir = '/hdfs/store/user/dntaylor/data/%s' % ntupleDict[period][analysis]
    #root_dir = 'root://cmsxrootd.fnal.gov//store/user/dntaylor/data/%s' % ntupleDict[period][analysis]

    sample_names = [os.path.basename(fname)
                    for string in samples
                    for fname in glob.glob("%s/%s" % (root_dir, string))]

    return root_dir, sample_names

def run_ntuples(analysis, channel, period, samples):
    '''Run a given analyzer for the H++ analysis'''
    ntup_dir = './ntuples%s_%stev_%s' % (analysis, period, channel)
    python_mkdir(ntup_dir)
    root_dir, sample_names = get_sample_names(analysis,period,samples)

    p = Pool(8)

    p.map(run_analyzer, [(analysis, channel, "%s/%s" % (root_dir, name), "%s/%s.root" % (ntup_dir, name), period) for name in sample_names])

    return 0

def submitJob(jobName,runArgs):
    '''
    Submit a job
    jobName: jobname for executable
    runArgs: [analysis, channel, period, sample]
    '''
    userDir = '/nfs_scratch/%s' % pwd.getpwuid(os.getuid())[0]
    submitDir = '%s/%s' % (userDir,jobName)
    if os.path.exists(submitDir):
        print "Submit directory exists, use a different JOBNAME."
        return 0
    python_mkdir(submitDir)

    # copy code
    os.system('tar -zcf %s/userCode.tar.gz analyzers utilities run.py' % submitDir)

    # copy submit script
    submitInfile = open('utilities/condor.submit')
    submitOutfile = open('%s/condor.submit' % submitDir, 'w')

    scriptName = '_'.join(runArgs)
    submitReplacements = {
        'EXECUTABLE' : scriptName,
    }

    for line in submitInfile:
        for src, target in submitReplacements.iteritems():
            line = line.replace(src, target)
        submitOutfile.write(line)

    submitInfile.close()
    submitOutfile.close()

    # copy executable
    executable = '%s/%s.sh' % (submitDir,scriptName)
    infile = open('utilities/doJob.sh')
    outfile = open(executable, 'w')

    replacements = {
        'OPTIONS' : ' '.join(runArgs),
    }

    for line in infile:
        for src, target in replacements.iteritems():
            line = line.replace(src, target)
        outfile.write(line)

    infile.close()
    outfile.close()

    os.system('chmod +x %s' % executable)

    # submit job
    os.system('cd %s; condor_submit condor.submit' % submitDir)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Run the desired analyzer on "
                                                 "FSA n-tuples")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('channel', type=str, choices=['WZ','TT','Hpp3l','Hpp4l','FakeRate'], help='Channel to run for given analysis')
    parser.add_argument('period', type=str, choices=['7','8','13'], help='Energy (TeV)')
    parser.add_argument('sample_names', nargs='+',help='Sample names w/ UNIX wildcards')
    parser.add_argument('-s','--submit',action='store_true',help='Submit jobs to condor')
    parser.add_argument('-jn','--jobName',nargs='?',type=str,const='',help='Job Name for condor submission')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.period == '7':
        print "7 TeV not implemented"
    else:
        print "Running %s:%s %s TeV analyzer" %(args.analysis, args.channel, args.period)
        if args.submit:
            root_dir, sample_names = get_sample_names(args.analysis, args.period, args.sample_names)
            for sample in sample_names:
                print sample
                runArgs = [args.analysis, args.channel, args.period, sample]
                jobName = '%s/%s' % (args.jobName, sample)
                submitJob(jobName,runArgs)
        else:
            run_ntuples(args.analysis, args.channel, args.period, args.sample_names)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
