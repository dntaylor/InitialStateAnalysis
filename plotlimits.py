#!/usr/bin/env python
from plotters.limits import *
import sys
import os
import argparse


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot limits")

    parser.add_argument('region', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[7, 8, 13], help='Run period')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='',choices=['ee100','em100','mm100','BP1','BP2','BP3','BP4'],help='Choose branching point')
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
            print 'Plotting limit for %s' % bp
            plot_limits(args.region,args.period,'limits_%s_%itev_%s'%(args.region,args.period,bp),branchingPoint=bp)
    else:
        print 'Plotting limit for %s' % args.branchingPoint
        plot_limits(args.region,args.period,'limits_%s_%itev_%s'%(args.region,args.period,args.branchingPoint),branchingPoint=args.branchingPoint)

    return 0


if __name__ == "__main__":
    main()
