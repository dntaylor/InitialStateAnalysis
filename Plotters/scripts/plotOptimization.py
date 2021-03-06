#!/usr/bin/env python
from InitialStateAnalysis.Plotters.OptimizerPlotter import *
from InitialStateAnalysis.Utilities.utilities import *
import sys
import os
import argparse


def parse_command_line(argv):
    parser = get_parser("Plot a given channel and period")

    parser.add_argument('-va','--variable',nargs='?',type=str,const='',choices=['st','met','zmass','zmassUnder','zmassOver','dR','hmass','hmassUnder','hmassOver'],help='Choose variable to plot')
    parser.add_argument('-av','--allVariables',action='store_true',help='Run over all variables')
    parser.add_argument('-nt','--numTaus',nargs='?',type=int,const='',choices=[0,1,2],help='Choose number of taus')
    parser.add_argument('-at','--allTaus',action='store_true',help='Run over all taus')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    variables = [args.variable]
    if args.allVariables: variables = ['st','met','zmass','zmassUnder','zmassOver','dR','hmass','hmassUnder','hmassOver']
    numTaus = [args.numTaus]
    if args.allTaus: numTaus = [0,1,2]
    for t in numTaus:
        for v in variables:
            plotOptimization(args.analysis,args.period,v,t)
    return 0


if __name__ == "__main__":
    main()
