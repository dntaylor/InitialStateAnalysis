#!/usr/bin/env python
from InitialStateAnalysis.Plotters.limits import *
from InitialStateAnalysis.Utilities.utilities import *
import sys
import os
import argparse


def parse_command_line(argv):
    parser = get_parser("Plot limits")

    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='comb',default='comb',choices=['mc','sideband','comb'],help='Choose BG estimation')
    parser.add_argument('-ub','--unblind', action='store_true',help='Unblind the analysis')
    parser.add_argument('-df','--do4l', action='store_true',help='Run the 4l lepton limits')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    branchingPoints = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']

    datacardString = '' if args.bgMode == "sideband" else "_{0}".format(args.bgMode)

    outstring = ''
    if args.period == 7:
        print "7 TeV not implemented"
        return 0
    elif args.allBranchingPoints:
        for bp in branchingPoints:
            print 'Plotting limit for %s' % bp
            limvals = plot_limits(args.analysis,args.channel,args.period,'limits_%s_%itev_%s%s'%(args.channel,args.period,bp,datacardString),branchingPoint=bp,bgMode=args.bgMode,do4l=args.do4l,unblind=args.unblind)
            outstring += '{3}: {0} [+{1},-{2}]\n'.format(limvals[0],limvals[1],limvals[2],bp)
            if args.analysis in ['HppComb']: plot_combined_limits(args.period,'limits_combinedCrossSection_%itev_%s%s'%(args.period,bp,datacardString),branchingPoint=bp,bgMode=args.bgMode,unblind=args.unblind)
    else:
        print 'Plotting limit for %s' % args.branchingPoint
        limvals = plot_limits(args.analysis,args.channel,args.period,'limits_%s_%itev_%s%s'%(args.channel,args.period,args.branchingPoint,datacardString),branchingPoint=args.branchingPoint,bgMode=args.bgMode,do4l=args.do4l,unblind=args.unblind)
        outstring += '{3}: {0} [+{1},-{2}]\n'.format(limvals[0],limvals[1],limvals[2],args.branchingPoint)
        if args.analysis in ['HppComb']: plot_combined_limits(args.period,'limits_combinedCrossSection_%itev_%s%s'%(args.period,args.branchingPoint,datacardString),branchingPoint=args.branchingPoint,bgMode=args.bgMode,unblind=args.unblind)

    savename = 'plots/limits/limits_%s_%itev'%(args.channel,args.period)
    if args.do4l: savename += '_4l'

    with open(savename,'w') as f:
        f.write(outstring)

    return 0


if __name__ == "__main__":
    main()
