#!/usr/bin/env python

import os
import sys
import subprocess
import errno
import argparse
import glob
from plotters.plotUtils import _3L_MASSES, _4L_MASSES, python_mkdir

def doDatacards(analysis,period,combineDir,bp):
    '''A function to move into the combined limits folder, run higgs combine tool on all datacards
       produced by the mklimits script, then copy the root files back here.'''
    datacardDir = 'datacards/%s_%itev/%s' % (analysis, period, bp)
    combineDatacardDir = '%s/%s_%itev' % (combineDir, analysis, period)
    datacardLimitsDir = 'limitData/%s_%itev/%s' % (analysis, period, bp)
    python_mkdir(combineDatacardDir)
    os.system('cp -r %s %s' %(datacardDir, combineDatacardDir))
    masses = _3L_MASSES if analysis == 'Hpp3l' else _4L_MASSES
    if period==13: masses = [500]
    for mass in masses:
        os.system('cd %s/%s/%i; pwd; eval `scramv1 runtime -sh`; combine -M MaxLikelihoodFit -t -1 --expectSignal 0 %s.txt; python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a mlfit.root -g plots.root; combine -M MaxLikelihoodFit -t -1 --expectSignal 0 %s.txt; python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -a mlfit.root -g plots.root' % (combineDatacardDir, bp, mass, bp, bp))

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('region', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[7, 8, 13], help='Run period')
    parser.add_argument('combineDir', type=str, help='Directory of the Higgs Combine CMSSW src directory')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')


    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    branchingPoints = ['ee100','em100','mm100','BP1','BP2','BP3','BP4']

    if args.period == 7:
        print "7 TeV not implemented"
    elif args.allBranchingPoints:
        for bp in branchingPoints:
            doDatacards(args.region,args.period,args.combineDir,bp)
    else:
        doDatacards(args.region,args.period,args.combineDir,args.branchingPoint)

    return 0


if __name__ == "__main__":
    main()
              
