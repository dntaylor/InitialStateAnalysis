#!/usr/bin/env python
'''
A script to run full CLs on an input datacard

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
import ROOT
import subprocess

_3L_MASSES = [170, 200, 250, 300, 350, 400, 450, 500, 600, 700]
_4L_MASSES = [130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise



def submitFwkliteJob(args):
    '''
    Submit a job using farmoutAnalysisJobs --fwklite
    '''
    logger = logging.getLogger(__name__)
    jobName = args.jobName
    dryrun = args.dryrun
    datacard = args.datacard
    mass = args.mass
    analysis = args.analysis
    bp = args.branchingPoint

    lookup = False
    if not datacard: # no datacard defined, look it up from bp, analysis, mass
        lookup = True
        if analysis in ['3lAP', 'AP'] and mass not in _3L_MASSES: return
        analysisMap = {
           'Combined': 'HppComb_8tev_HppComb',
           'AP': 'HppAP_8tev_HppAP',
           'PP': 'HppPP_8tev_HppPP',
           '3lAP': 'Hpp3l_8tev_Hpp3l',
           '3lPP': 'Hpp3l_8tev_Hpp3l',
           '3lAPandPP': 'Hpp3l_8tev_Hpp3l',
           '4lPP': 'Hpp4l_8tev_Hpp4l',
        }
        endString = {
            '3lPP': '_4l',
            '3lAPandPP': '_APandPP',
        }
        datacard = 'datacards/{0}/{1}/{2}/{1}_comb{3}.txt'.format(analysisMap[analysis],bp,mass,endString[analysis] if analysis in endString else '')

    # get datacard path relative to $CMSSW_BASE
    dfull = os.path.abspath(os.path.join(os.path.dirname(__file__),datacard))
    ddir = os.path.dirname(dfull)
    dname = os.path.basename(dfull)
    cmssw_base = os.environ['CMSSW_BASE']
    dsplit = dfull.split('/')
    srcpos = dsplit.index('src')
    drel = '/'.join(dsplit[srcpos:])
    dreldir = '/'.join(dsplit[srcpos:-1])

    sample_dir = '/nfs_scratch/{0}/{1}/'.format(pwd.getpwuid(os.getuid())[0], jobName)
    if lookup: sample_dir += '{0}/{1}/{2}/'.format(analysis,bp,mass)

    # create submit dir
    submit_dir = '{0}/submit'.format(sample_dir)
    if os.path.exists(submit_dir):
        if lookup:
            logger.warning('Submission directory exists for {0} {1} {2} {3}.'.format(jobName, analysis, bp, mass))
        else:
            logger.warning('Submission directory exists for {0}.'.format(jobName))
        return


    # first, get the approximate bounds from asymptotic
    work = 'working'
    os.system('rm -rf working')
    python_mkdir(work)
    workfull = os.path.join(os.path.dirname(__file__),work)
    command = 'pushd {0}; combine -M Asymptotic {1} -m {2};'.format(workfull,dfull,mass) 
    logger.info('Finding initial bounds with Asymptotic limit: {0}'.format(datacard))
    out = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

    fname = os.path.join(workfull, "higgsCombineTest.Asymptotic.mH{0}.root".format(mass))
    file = ROOT.TFile(fname,"READ")
    tree = file.Get("limit")
    if not tree: 
        logger.warning('Asymptotic presearch failed')
        return
    quartiles = []
    for i, row in enumerate(tree):
        quartiles += [row.limit]

    # setup the job parameters
    rmin = 0.8*min(quartiles)
    rmax = 1.2*max(quartiles)
    num_points = 100
    points_per_job = 10
    toys = 10000

    # create dag dir
    dag_dir = '{0}/dags/dag'.format(sample_dir)
    os.system('mkdir -p {0}'.format(os.path.dirname(dag_dir)))
    os.system('mkdir -p {0}'.format(dag_dir+'inputs'))

    # output dir
    output_dir = 'srm://cmssrm.hep.wisc.edu:8443/srm/v2/server?SFN=/hdfs/store/user/{0}/{1}/'.format(pwd.getpwuid(os.getuid())[0], jobName)
    if lookup: output_dir += '{0}/{1}/{2}/'.format(analysis,bp,mass)

    # create file list
    rlist = [r*(rmax-rmin)/num_points + rmin for r in range(int(num_points/points_per_job))]
    input_name = '{0}/rvalues.txt'.format(dag_dir+'inputs')
    with open(input_name,'w') as file:
        for r in rlist:
            file.write('{0}\n'.format(r))

    # create bash script
    bash_name = '{0}/{1}.sh'.format(dag_dir+'inputs', jobName)
    bashScript = '#!/bin/bash\n'
    bashScript += 'printenv\n'
    bashScript += 'read -r RVAL < $INPUT\n'
    for i in range(points_per_job):
        dr = i*(rmax-rmin)/points_per_job
        bashScript += 'combine $CMSSW_BASE/{0} -M HybridNew --freq -s -1 --singlePoint $(bc -l <<< "$RVAL+{1}") --saveToys --fullBToys --clsAcc 0 --saveHybridResult -m {2} -n Tag -T {3} -i 2\n'.format(drel,dr,mass,toys)
    bashScript += 'hadd $OUTPUT higgsCombineTag.HybridNew.mH{0}.*.root\n'.format(mass)
    bashScript += 'rm higgsCombineTag.HybridNew.mH{0}.*.root\n'.format(mass)
    with open(bash_name,'w') as file:
        file.write(bashScript)
    os.system('chmod +x {0}'.format(bash_name))

    # create farmout command
    farmoutString = 'farmoutAnalysisJobs --infer-cmssw-path --fwklite --input-file-list={0} --assume-input-files-exist'.format(input_name)
    farmoutString += ' --submit-dir={0} --output-dag-file={1} --output-dir={2}'.format(submit_dir, dag_dir, output_dir)
    farmoutString += ' --extra-usercode-files="{0}" {1} {2}'.format(dreldir, jobName, bash_name)

    if not args.dryrun:
        if lookup:
            logger.info('Submitting {0} {1} {2} {3}'.format(jobName, analysis, bp, mass))
        else:
            logger.info('Submitting {0}'.format(jobName))
        os.system(farmoutString)
    else:
        print farmoutString

    return

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit limts to grid')

    parser.add_argument('jobName',type=str,help='Job Name for condor submission')
    parser.add_argument('-d','--datacard', nargs='?', type=str, const='', help='Datacard to process')
    parser.add_argument('-a','--analysis', nargs='?',type=str,const='Combined',choices=['Combined','AP','PP','3lAP','3lPP','3lAPandPP','4lPP'],help='Analysis to process')
    parser.add_argument('-m','--mass', nargs='?',type=str,const='120',help='Mass for Higgs combine')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses')
    parser.add_argument('-aa','--allAnalyses',action='store_true',help='Run over all anlayses')
    parser.add_argument('-dr','--dryrun',action='store_true',help='Create jobs but dont submit')
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

    allowedAnalyses = ['Combined','AP','PP','3lAP','3lPP','3lAPandPP','4lPP'] if args.allAnalyses else [args.analysis]
    allowedBranchingPoints = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'] if args.allBranchingPoints else [args.branchingPoint]
    allowedMasses = _4L_MASSES if args.allMasses else [args.mass]

    for an in allowedAnalyses:
        for bp in allowedBranchingPoints:
            for m in allowedMasses:
                args.analysis = an
                args.branchingPoint = bp
                args.mass = m
                submitFwkliteJob(args)

if __name__ == "__main__":
    status = main()
    sys.exit(status)

