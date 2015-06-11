#!/usr/bin/env python

from plotters.Plotter import Plotter
from plotters.plotUtils import *
import argparse
import itertools
import sys

ZMASS = 91.1876

def unfold(analysis,channel,period):
    '''
    Unfold distributions for a given analysis.
    '''
    # open analysis distributions
    # open gen distributions
    # create 2d histograms
    # create roounfold
    # unfold

    return 0

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Unfold distributions")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    unfold(args.analysis,args.channel,args.period)

    return 0


if __name__ == "__main__":
    main()
