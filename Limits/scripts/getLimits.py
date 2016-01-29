#!/usr/bin/env python
'''
A script to retrieve the limits

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
from multiprocessing import Pool

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

def limitsWrapper(args):
    getLimits(*args)

def getLimits(analysis,bp,mass,outDir):
    '''
    Submit a job using farmoutAnalysisJobs --fwklite
    '''
    logger = logging.getLogger(__name__)

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

    # first, get the approximate bounds from asymptotic
    work = 'working'
    python_mkdir(work)
    workfull = os.path.join(os.path.dirname(__file__),work)
    combineCommand = 'combine -M Asymptotic {0} -m {1}'.format(dfull,mass)
    command = 'pushd {0}; nice {1};'.format(workfull,combineCommand) 
    logger.info('{0}:{1}:{2}: Finding Asymptotic limit: {3}'.format(analysis,bp,mass,datacard))
    logger.debug('{0}:{1}:{2}: {3}'.format(analysis,bp,mass,combineCommand))
    out = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

    fname = os.path.join(workfull, "higgsCombineTest.Asymptotic.mH{0}.root".format(mass))
    file = ROOT.TFile(fname,"READ")
    tree = file.Get("limit")
    if not tree: 
        logger.warning('Asymptotic presearch failed')
        quartiles = [0., 0., 0., 0., 0., 0.]
    else:
        quartiles = []
        for i, row in enumerate(tree):
            quartiles += [row.limit]

    # now get the fullCLs
    args = [
        ['Expected 0.025', '--expectedFromGrid 0.025', 'higgsCombineTest.HybridNew.mH{0}.quant0.025.root'],
        ['Expected 0.160', '--expectedFromGrid 0.160', 'higgsCombineTest.HybridNew.mH{0}.quant0.160.root'],
        ['Expected 0.500', '--expectedFromGrid 0.500', 'higgsCombineTest.HybridNew.mH{0}.quant0.500.root'],
        ['Expected 0.840', '--expectedFromGrid 0.840', 'higgsCombineTest.HybridNew.mH{0}.quant0.840.root'],
        ['Expected 0.975', '--expectedFromGrid 0.975', 'higgsCombineTest.HybridNew.mH{0}.quant0.975.root'],
        ['Observed',       '',                         'higgsCombineTest.HybridNew.mH{0}.root'],
    ]
    
    # merge the output
    gridfile = 'grid_{0}.root'.format(mass)
    sourceDir = '/hdfs/store/user/dntaylor/2016-01-21_allLimits_10KToys_100Points_v1/{0}/{1}/{2}'.format(analysis,bp,mass)
    logger.info('{0}:{1}:{2}: Merging: {3}'.format(analysis,bp,mass,sourceDir))
    haddCommand = 'hadd {0} {1}/*.root'.format(gridfile,sourceDir)
    command = 'pushd {0}; nice {1};'.format(workfull,haddCommand)
    logger.debug('{0}:{1}:{2}: {3}'.format(analysis,bp,mass,haddCommand))
    out = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

    # get CL
    fullQuartiles = []
    for i in range(len(args)):
        logger.info('{0}:{1}:{2}: Calculating: {3}'.format(analysis,bp,mass,args[i][0]))
        combineCommand = 'combine {0} -M HybridNew --freq --grid={1} -m {2} --rAbsAcc 0.001 --rRelAcc 0.001 {3}'.format(dfull, gridfile, mass, args[i][1])
        command = 'pushd {0}; nice {1};'.format(workfull, combineCommand)
        logger.debug('{0}:{1}:{2}: {3}'.format(analysis,bp,mass,combineCommand))
        outfile = '{0}/{1}'.format(workfull,args[i][2].format(mass))
        out = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

        # read the limit
        file = ROOT.TFile(outfile,"READ")
        tree = file.Get("limit")
        if not tree:
            logger.warning('HybridNew failed')
            val = 0.
        else:
            val = 0.
            for i, row in enumerate(tree):
                val = row.limit
        fullQuartiles += [val]

    # save the values
    quartileMap = {
        'asymptotic' : quartiles,
        'fullCLs' : fullQuartiles,
    }
    for name in ['asymptotic','fullCLs']:
        fileDir = '{0}/{1}/{2}/{3}'.format(name,analysis,bp,mass)
        if outDir: fileDir = outDir + '/' + fileDir
        python_mkdir(fileDir)
        fileName = '{0}/limits.txt'.format(fileDir)
        with open(fileName,'w') as f:
            outline = ' '.join([str(x) for x in quartileMap[name]])
            logger.info('{0}:{1}:{2}: Limits: {3} - {4}'.format(analysis,bp,mass,name, outline))
            f.write(outline)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Submit limts to grid')

    parser.add_argument('-d','--directory', nargs='?',type=str,const='',help='Custom out directory')
    parser.add_argument('-a','--analysis', nargs='?',type=str,const='Combined',choices=['Combined','AP','PP','3lAP','3lPP','3lAPandPP','4lPP'],help='Analysis to process')
    parser.add_argument('-m','--mass', nargs='?',type=str,const='120',help='Mass for Higgs combine')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses')
    parser.add_argument('-aa','--allAnalyses',action='store_true',help='Run over all anlayses')
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
            # clean working
            os.system('rm -rf working')
            if len(allowedMasses)==1:
                getLimits(an,bp,allowedMasses[0],args.directory)
            else:
                allArgs = []
                for m in allowedMasses:
                    newArgs = [an,bp,m,args.directory]
                    allArgs += [newArgs]
                p = Pool(6)
                try:
                    p.map_async(limitsWrapper, allArgs).get(999999)
                except KeyboardInterrupt:
                    p.terminate()
                    print 'limits cancelled'
                    sys.exit(1)

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)

