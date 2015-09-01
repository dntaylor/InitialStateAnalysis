#!/usr/bin/env python
from InitialStateAnalysis.Plotters.limits import *
import sys
import os
import argparse


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot limits")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l','HppComb'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','HppComb'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[8, 13], help='Run period')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point')
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

    datacardString = '' if args.bgMode == "sideband" else "_{0}".format(args.bgMode)

    if args.period == 7:
        print "7 TeV not implemented"
    elif args.allBranchingPoints:
        for bp in branchingPoints:
            print 'Plotting limit for %s' % bp
            plot_limits(args.analysis,args.region,args.period,'limits_%s_%itev_%s%s'%(args.region,args.period,bp,datacardString),branchingPoint=bp,bgMode=args.bgMode,do4l=args.do4l)
    else:
        print 'Plotting limit for %s' % args.branchingPoint
        plot_limits(args.analysis,args.region,args.period,'limits_%s_%itev_%s%s'%(args.region,args.period,args.branchingPoint,datacardString),branchingPoint=args.branchingPoint,bgMode=args.bgMode,do4l=args.do4l)

    return 0


if __name__ == "__main__":
    main()
