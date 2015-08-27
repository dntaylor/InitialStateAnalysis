#!/usr/bin/env python

import os
import sys
import subprocess
import errno
import argparse
import glob
from InitialStateAnalysis.Plotters.plotUtils import _3L_MASSES, _4L_MASSES, python_mkdir

def doDatacards(analysis,region,period,bp,bgMode,do4l):
    '''A function to move into the combined limits folder, run higgs combine tool on all datacards
       produced by the mklimits script, then copy the root files back here.'''
    # do the combination
    combos = {
        'HppComb': ['Hpp3l', 'Hpp4l'],
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
    # merge Hpp3l
    for mass in masses:
        if not do4l: os.system('pushd {1}/{2}; combineCards.py {0}_[em][em][em].txt > {0}_comb.txt'.format(bp,datacardDir,mass))
        if do4l: os.system('pushd {1}/{2}; combineCards.py {0}_[em][em][em]_4l.txt > {0}_comb_4l.txt'.format(bp,datacardDir,mass))
    # merge for combine
    if analysis in combos:
        for mass in masses:
            dirsToCombine = ['datacards/%s_%itev_%s/%s' % (a, period, a, bp) for a in combos[analysis]]
            theCards = ['%s/%i/%s%s.txt' %(x,mass,bp,datacardString) for x in dirsToCombine]
            cardsToCombine = [x for x in theCards if os.path.isfile(x)]
            outCard = '%s/%i/%s%s.txt' %(datacardDir,mass,bp,datacardString)
            python_mkdir('%s/%i' %(datacardDir,mass))
            print 'Creating combined card mass %i' % mass
            os.system('combineCards.py %s > %s' % (combineDatacardDir,' '.join(cardsToCombine),outCard))
            
    os.system('cp -r %s %s' %(datacardDir, combineDatacardDir))
    python_mkdir(datacardLimitsDir)
    for mass in masses:
        os.system('cd %s/%s/%i; combine -m %i -M Asymptotic %s%s.txt' % (combineDatacardDir, bp, mass, mass, bp, datacardString))
        os.system('cp %s/%s/*/higgsCombineTest.Asymptotic.mH%i.root %s/higgsCombineTest.Asymptotic.mH%i%s.root' % (combineDatacardDir, bp, mass, datacardLimitsDir, mass, datacardString))

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l','HppComb'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','HppComb','WZ'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[7, 8, 13], help='Run period')
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

    branchingPoints = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']

    if args.period == 7:
        print "7 TeV not implemented"
    elif args.allBranchingPoints:
        for bp in branchingPoints:
            doDatacards(args.analysis,args.region,args.period,bp,args.bgMode,args.do4l)
    else:
        doDatacards(args.analysis,args.region,args.period,args.branchingPoint,args.bgMode,args.do4l)

    return 0


if __name__ == "__main__":
    main()
              
