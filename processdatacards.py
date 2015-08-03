#!/usr/bin/env python

import os
import sys
import subprocess
import errno
import argparse
import glob
from plotters.plotUtils import _3L_MASSES, _4L_MASSES, python_mkdir

def doDatacards(analysis,region,period,combineDir,bp,bgMode):
    '''A function to move into the combined limits folder, run higgs combine tool on all datacards
       produced by the mklimits script, then copy the root files back here.'''
    # do the combination
    combos = {
        'HppComb': ['Hpp3l', 'Hpp4l'],
    }
    datacardDir = 'datacards/%s_%itev_%s/%s' % (analysis, period, region, bp)
    combineDatacardDir = '%s/%s_%itev_%s' % (combineDir, analysis, period, region)
    datacardLimitsDir = 'limitData/%s_%itev_%s/%s' % (analysis, period, region, bp)
    masses = _3L_MASSES if analysis == 'Hpp3l' else _4L_MASSES
    datacardString = '' if bgMode == "sideband" else "_{0}".format(bgMode)
    python_mkdir(combineDatacardDir)
    if analysis in combos:
        for mass in masses:
            dirsToCombine = ['datacards/%s_%itev_%s/%s' % (a, period, a, bp) for a in combos[analysis]]
            theCards = ['%s/%i/%s%s.txt' %(x,mass,bp,datacardString) for x in dirsToCombine]
            cardsToCombine = [x for x in theCards if os.path.isfile(x)]
            outCard = '%s/%i/%s%s.txt' %(datacardDir,mass,bp,datacardString)
            python_mkdir('%s/%i' %(datacardDir,mass))
            print 'Creating combined card mass %i' % mass
            os.system('pushd %s;  eval `scramv1 runtime -sh`; popd; combineCards.py %s > %s' % (combineDatacardDir,' '.join(cardsToCombine),outCard))
            
    os.system('cp -r %s %s' %(datacardDir, combineDatacardDir))
    if period==13: masses = [500]
    python_mkdir(datacardLimitsDir)
    for mass in masses:
        os.system('cd %s/%s/%i; pwd; eval `scramv1 runtime -sh`; combine -m %i -M Asymptotic %s%s.txt' % (combineDatacardDir, bp, mass, mass, bp, datacardString))
        os.system('cp %s/%s/*/higgsCombineTest.Asymptotic.mH%i.root %s/higgsCombineTest.Asymptotic.mH%i%s.root' % (combineDatacardDir, bp, mass, datacardLimitsDir, mass, datacardString))

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l','HppComb'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','HppComb','WZ'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[7, 8, 13], help='Run period')
    parser.add_argument('combineDir', type=str, help='Directory of the Higgs Combine CMSSW src directory')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='sideband',default='sideband',choices=['mc','sideband','comb'],help='Choose BG estimation')


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
            doDatacards(args.analysis,args.region,args.period,args.combineDir,bp,args.bgMode)
    else:
        doDatacards(args.analysis,args.region,args.period,args.combineDir,args.branchingPoint,args.bgMode)

    return 0


if __name__ == "__main__":
    main()
              
