#!/usr/bin/env python
import logging
import sys
import argparse
import numpy as np
import itertools

from InitialStateAnalysis.Plotters.plotUtils import getSigMap, getIntLumiMap, getChannels, getMergeDict, ZMASS, getChannelBackgrounds
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Limits.WZLimits import WZLimits
from InitialStateAnalysis.Limits.limitUtils import *
from InitialStateAnalysis.Utilities.utilities import *
from multiprocessing import Pool

def parse_command_line(argv):
    parser = get_parser("Produce datacards")

    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to limits.')
    parser.add_argument('-ub','--unblind',action='store_true',help='unblind')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')
    parser.add_argument('--outputDirectory',type=str,default='',help='Directory to output datacards')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    wzlimits(args.analysis,args.channel,args.period,scalefactor=args.scaleFactor,cut=args.cut,unblind=args.unblind,outputDirectory=args.outputDirectory)

    return 0


if __name__ == "__main__":
    main()

