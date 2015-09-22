#!/usr/bin/env python

import os
import sys
import subprocess
import errno
import argparse
import glob
import logging
from InitialStateAnalysis.Plotters.plotUtils import _3L_MASSES, _4L_MASSES, python_mkdir
from InitialStateAnalysis.Utilities.utilities import *

def doDatacards(analysis,region,period,bp,bgMode,do4l):
    '''A function to move into the combined limits folder, run higgs combine tool on all datacards
       produced by the mklimits script, then copy the root files back here.'''
    logging.info('Processing %s with mode %s' % (bp,bgMode))
    pipe = subprocess.PIPE
    # do the combination
    combos = {
        'HppComb': ['Hpp3l', 'Hpp4l'],
        'HppAP': ['Hpp3l'],
        'HppPP': ['Hpp4l'],
    }
    datacardDir = 'datacards/%s_%itev_%s/%s' % (analysis, period, region, bp)
    combineDatacardDir = 'combineWorking/%s_%itev_%s' % (analysis, period, region)
    datacardLimitsDir = 'limitData/%s_%itev_%s/%s' % (analysis, period, region, bp)
    masses = _3L_MASSES if analysis == 'Hpp3l' else _4L_MASSES
    if do4l: masses = _4L_MASSES
    if period==13: masses = [500]
    datacardString = '' if bgMode == "sideband" else "_{0}".format(bgMode)
    if do4l: datacardString += "_4l"
    python_mkdir(combineDatacardDir)
    # merge for combine
    if analysis in combos:
        for mass in masses:
            # merge the inputs
            logging.info('%s: Merging %i' % (bp,mass))
            mergeCommands = []
            mergeCommands += ['pushd datacards/{1}/{0}/{2}; combineCards.py {0}_[em][em][em].txt > {0}_comb.txt'.format(bp,'Hpp3l_8tev_Hpp3l',mass)]
            mergeCommands += ['pushd datacards/{1}/{0}/{2}; combineCards.py {0}_[em][em][em]_4l.txt > {0}_comb_4l.txt'.format(bp,'Hpp3l_8tev_Hpp3l',mass)]
            mergeCommands += ['pushd datacards/{1}/{0}/{2}; combineCards.py {0}_[em][em][em][em].txt > {0}_comb.txt'.format(bp,'Hpp4l_8tev_Hpp4l',mass)]
            for command in mergeCommands:
                out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
            # merge the merged cards
            dirsToCombine = ['datacards/%s_%itev_%s/%s' % (a, period, a, bp) for a in combos[analysis]]
            theCards = ['%s/%i/%s%s.txt' %(x,mass,bp,datacardString) for x in dirsToCombine]
            # manually add Hpp3l 4l
            if analysis in ['HppComb', 'HppPP']: theCards += ['datacards/Hpp3l_%itev_Hpp3l/%s/%i/%s%s_4l.txt' % (period, bp, mass,bp,datacardString)]
            cardsToCombine = [x for x in theCards if os.path.isfile(x)]
            outCard = '%s/%i/%s%s.txt' %(datacardDir,mass,bp,datacardString)
            python_mkdir('%s/%i' %(datacardDir,mass))
            command = 'combineCards.py %s > %s' % (' '.join(cardsToCombine),outCard)
            out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
    else: # merge for individual analyses
        for mass in masses:
            logging.info('%s: Merging %i' % (bp,mass))
            if analysis in ['Hpp3l']:
                if do4l: command = 'pushd {1}/{2}; combineCards.py {0}_[em][em][em]_4l.txt > {0}_comb_4l.txt'.format(bp,datacardDir,mass)
                if not do4l: command = 'pushd {1}/{2}; combineCards.py {0}_[em][em][em].txt > {0}_comb.txt'.format(bp,datacardDir,mass)
            if analysis in ['Hpp4l']:
                command = 'pushd {1}/{2}; combineCards.py {0}_[em][em][em][em].txt > {0}_comb.txt'.format(bp,datacardDir,mass)
            out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
            
    command = 'cp -r %s %s' %(datacardDir, combineDatacardDir)
    out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
    python_mkdir(datacardLimitsDir)
    for mass in masses:
        logging.info('%s: Calculating limit for %i' % (bp,mass))
        command = 'cd %s/%s/%i; combine -m %i -M Asymptotic %s%s.txt' % (combineDatacardDir, bp, mass, mass, bp, datacardString)
        out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
        command = 'cp %s/%s/*/higgsCombineTest.Asymptotic.mH%i.root %s/higgsCombineTest.Asymptotic.mH%i%s.root' % (combineDatacardDir, bp, mass, datacardLimitsDir, mass, datacardString)
        out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]

def parse_command_line(argv):
    parser = get_parser("Produce datacards")

    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='comb',default='comb',choices=['mc','sideband','comb'],help='Choose BG estimation')
    parser.add_argument('-df','--do4l', action='store_true',help='Run the 4l lepton limits')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    branchingPoints = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']

    if args.period == 7:
        print "7 TeV not implemented"
    elif args.allBranchingPoints:
        for bp in branchingPoints:
            doDatacards(args.analysis,args.channel,args.period,bp,args.bgMode,args.do4l)
    else:
        doDatacards(args.analysis,args.channel,args.period,args.branchingPoint,args.bgMode,args.do4l)

    return 0


if __name__ == "__main__":
    main()
              
