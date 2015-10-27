#!/usr/bin/env python

import argparse
import itertools
import sys
import pickle
import json
import logging
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES
from InitialStateAnalysis.Utilities.utilities import *

def getCrossSection(analysis,channel,period,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',False)
    cut = kwargs.pop('cut','1')
    scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    directory = kwargs.pop('directory','')

    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,period,channel)
    saves = '%s_%s_%iTeV' % (analysis,channel,period)
    mergeDict = getMergeDict(period)
    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ'   : 3,
        'WZ_W' : 1,
        'WZ_Dijet': 1,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[analysis]
    finalStates, leptons = getChannels(nl)
    if analysis in ['WZ']: finalStates = ['eee','eem','mme','mmm']
    sigMap = getSigMap(nl)
    channelBackground =  getChannelBackgrounds(period)
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,dataDriven=doDataDriven)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel]])
    plotter.initializeDataSamples([sigMap[period]['data']])
    intLumi = getIntLumiMap()[period]
    plotter.setIntLumi(intLumi)

    # get acceptance and efficiency numbers
    # just WZ
    # zWindow, fiducial
    # 4 channels
    accEff = getAcceptanceEfficiency(analysis)

    
    


def parse_command_line(argv):
    parser = get_parser("Extract cross section")

    parser.add_argument('-dd','--doDataDriven',action='store_true',help='Do datadriven background estimation')
    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to cross section')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for cross section')
    parser.add_argument('-d','--directory',type=str,default='',help='Custom subdirectory (to keep more than one ntuple at a time)')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    getCrossSection(args.analysis,args.channel,args.period,doDataDriven=args.doDataDriven,cut=args.cut,scaleFactor=args.scaleFactor,directory=args.directory)

    return 0


if __name__ == "__main__":
    main()

