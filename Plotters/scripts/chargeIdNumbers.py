#!/usr/bin/env python

import logging
import json
import math

from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES
from InitialStateAnalysis.Utilities.utilities import *

def getRatio(num,denom):
    val = num[0]/denom[0] if denom[0] else -1
    err = val * math.sqrt((num[1]/num[0])**2 + (denom[1]/denom[0])**2)
    return (val,err)

def getSum(*vals):
    val = sum([x[0] for x in vals])
    err = math.sqrt(sum([x[1]**2 for x in vals]))
    return (val,err)

def main(argv=None):
    loglevel = getattr(logging,'INFO')
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    nl = 2
    mass = 500
    analysis = 'Hpp2l'
    runPeriod = 8
    channel = 'Charge'
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = getChannelBackgrounds(runPeriod)
    finalStates, leptons = getChannels(nl)
    mergeDict = getMergeDict(runPeriod)
    logger.info('Load plotter')
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,rootName='chargeId')
    plotter.initializeBackgroundSamples([sigMap[runPeriod]['Z']])
    #plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotMCData'
    plotMethod = getattr(plotter,plotMode)
    results = {'e': {}, 'm': {}}
    myCut = 'finalstate.met<20.'
    samesign = 'l1.Chg==l2.Chg && {0}'.format(myCut)
    oppsign = 'l1.Chg!=l2.Chg && {0}'.format(myCut)
    logger.info('Calculate values')
    ptBins = [10.,20.,30.,40.,60.,100.,200.]
    etaBins = {
        'e': [0,1.479,2.5],
        'm': [0,1.2,2.4],
    }
    for l in ['e']:
        flv = 'z1Flv=="{0}" && fabs(z1.mass-91.1876)<10.'.format(l+l)
        results[l] = {'ss':{},'os':{},'unc':{},}

        for p in range(len(ptBins)-1):
            logger.info('{0}: pT {1:f}-{2:f}'.format(l,ptBins[p],ptBins[p+1]))
            # tag first lepton, probe second, both in same eta region
            ptcut = 'l2.Pt>={0:f} && l2.Pt<{1:f}'.format(ptBins[p],ptBins[p+1])

            for e in range(len(etaBins[l])-1):
                logger.info('{0}: eta {1:f}-{2:f}'.format(l,etaBins[l][e],etaBins[l][e+1]))
                etacut = 'fabs(l1.Eta)>={0:f} && fabs(l1.Eta)<{1:f} && fabs(l2.Eta)>={0:f} && fabs(l2.Eta)<{1:f}'.format(etaBins[l][e],etaBins[l][e+1])

                mc = 'mc{0}{1}'.format(p,e)
                data = 'data{0}{1}'.format(p,e)
                thisCut = '{0} && {1} && {2} && {3}'.format(samesign,flv,ptcut,etacut)
                results[l]['ss'][mc] = plotter.getBackgroundEntries(thisCut,doError=True)
                results[l]['ss'][data] = plotter.getDataEntries(thisCut,doError=True)
                logger.info('{0}: Same sign: MC: {1:f} +/- {2:f}; Data: {3:f} +/- {4:f}'.format(l,results[l]['ss'][mc][0],results[l]['ss'][mc][1],results[l]['ss'][data][0],results[l]['ss'][data][1]))
                plotMethod('z1.mass', [60,60,120],   'z1Mass_samesign_{0}_{1}_{2}'.format(l,p,e),      yaxis='Events/1.0 GeV', xaxis='M_{\\ell^{+}\\ell^{+}} (GeV)',     legendpos=43,logy=0,cut=thisCut)

                thisCut = '{0} && {1} && {2} && {3}'.format(oppsign,flv,ptcut,etacut)
                results[l]['os'][mc] = plotter.getBackgroundEntries(thisCut,doError=True)
                results[l]['os'][data] = plotter.getDataEntries(thisCut,doError=True)
                logger.info('{0}: Opposite sign: MC: {1:f} +/- {2:f}; Data: {3:f} +/- {4:f}'.format(l,results[l]['os'][mc][0],results[l]['os'][mc][1],results[l]['os'][data][0],results[l]['os'][data][1]))
                plotMethod('z1.mass', [60,60,120],   'z1Mass_oppositesign_{0}_{1}_{2}'.format(l,p,e),  yaxis='Events/1.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=thisCut)

                mcnum = results[l]['ss'][mc]
                mcdenom = getSum(results[l]['ss'][mc],results[l]['os'][mc])
                mcratio = getRatio(mcnum,mcdenom)
                results[l]['unc'][mc] = mcratio
                datanum = results[l]['ss'][data]
                datadenom = getSum(results[l]['ss'][data],results[l]['os'][data])
                dataratio = getRatio(datanum,datadenom)
                results[l]['unc'][data] = dataratio
                ratio = getRatio(results[l]['unc'][data],results[l]['unc'][mc])
                logger.info('{0}: Uncertainty: MC: {1:f} +/- {2:f}; Data: {3:f} +/- {4:f}; Data/MC: {5:f} +/- {6:f}'.format(l,results[l]['unc'][mc][0],results[l]['unc'][mc][1],results[l]['unc'][data][0],results[l]['unc'][data][1],ratio[0],ratio[1]))

    print json.dumps(results,sort_keys=True,indent=4)



if __name__ == "__main__":
    main()
