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

def getChargeSystematic(analysis,channel,runPeriod,**kwargs):
    '''A function to simplify plotting multiple channels and run periods.'''
    mass = kwargs.pop('mass',500)
    myCut = kwargs.pop('myCut','1')
    scaleFactor = kwargs.pop('scaleFactor','event.pu_weight*event.lep_scale*event.trig_scale')
    loglevel = kwargs.pop('loglevel','INFO')
    for key, value in kwargs.iteritems():
        logger.warning("Unrecognized parameter '" + key + "' = " + str(value))
        return 0

    logging.info("%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass))
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
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = getChannelBackgrounds(runPeriod)

    finalStates, leptons = getChannels(nl)
    mergeDict = getMergeDict(runPeriod)
    cutFlowMap = {}
    cutFlowMap[channel] = defineCutFlowMap(channel,finalStates,mass)

    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel,rootName='charge_sys')
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])

    systMap = {}
    for numElec in range(nl+1):
        logging.info('%s:%s:%iTeV: Num electrons: %i' % (analysis, channel, runPeriod, numElec))
        systMap[numElec] = {}
        states = [x for x in finalStates if x.count('e')==numElec]
        chanCuts = ['channel=="{0}"'.format(x) for x in states]
        chanCut = ' || '.join(chanCuts)
        cut = '{0} && ({1})'.format(myCut,chanCut)
        sig_num = plotter.getSignalEntries(cut,customScale='event.pu_weight*event.lep_scale*event.trig_scale*event.charge_uncertainty')
        bg_num = plotter.getBackgroundEntries(cut,customScale='event.pu_weight*event.lep_scale*event.trig_scale*event.charge_uncertainty')
        data_num = plotter.getDataEntries(cut,customScale='charge_uncertainty')
        sig_den = plotter.getSignalEntries(cut)
        bg_den = plotter.getBackgroundEntries(cut)
        data_den = plotter.getDataEntries(cut)
        sig_syst = sig_num/sig_den if sig_den else 1.
        bg_syst = bg_num/bg_den if bg_den else 1.
        data_syst = data_num/data_den if data_den else 1.
        systMap[numElec]['sig'] = sig_syst
        systMap[numElec]['bg'] = bg_syst
        systMap[numElec]['data'] = data_syst

    return systMap

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    loglevel = getattr(logging,'INFO')
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    systMap = {}
    for analysis in ['Hpp3l','Hpp4l']:
        masses = _3L_MASSES if analysis in ['Hpp3l'] else _4L_MASSES
        systMap[analysis] = {}
        for mass in masses:
            systMap[analysis][mass] = getChargeSystematic(analysis,analysis,8,mass=mass,myCut='select.passTight')

    print json.dumps(systMap,sort_keys=True,indent=4)

    return 0


if __name__ == "__main__":
    main()

