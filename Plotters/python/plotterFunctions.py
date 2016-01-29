import itertools
import os
import sys
import errno
import argparse

from Plotter import Plotter
from plotUtils import *

def getSignalOverBackground(analysis,region,period,cut,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',True)
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')

    nl = 3
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_sOverB',mergeDict=mergeDict,scaleFactor=scalefactor,datadriven=doDataDriven)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[region+'datadriven' if doDataDriven else region]])
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    sig = 0.
    sigErr2 = 0.
    bg = 0.
    bgErr2 = 0.
    for b in plotter.backgrounds:
        val, err = plotter.getNumEntries(cut,b,doError=True)
        if b=='WZJets':
            sig += val
            sigErr2 += err**2
        else:
            bg += val
            bgErr2 += err**2
    sigErr = sigErr2**0.5
    bgErr = bgErr2**0.5

    return sig/bg if bg else -1.

def getPerChannelYields(analysis,region,period,cut,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',True)
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    unblind = kwargs.pop('unblind',False)
    tightW = kwargs.pop('tightW',False)

    nl = 3
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    bgs = channelBackground[region+'datadriven' if doDataDriven else region]

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_statUnc',mergeDict=mergeDict,scaleFactor=scalefactor,datadriven=doDataDriven,tightW=tightW)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in bgs])
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    yields = {}
    errs = {}
    for chan in channels:
        yields[chan] = {}
        errs[chan] = {}
        chanCut = '{0} && channel=="{1}"'.format(cut,chan)
        for b in bgs + ['datadriven']:
            val, err = plotter.getNumEntries(chanCut,sigMap[period][b],doError=True)
            yields[chan][b] = val
            errs[chan][b] = err

        data, dataErr = plotter.getDataEntries(chanCut,doError=True)
        yields[chan]['data'] = data
        errs[chan]['data'] = dataErr

    yields['wz'] = {}
    errs['wz'] = {}
    for i in bgs + ['datadriven'] + ['data']:
        tot = sum([yields[chan][i] for chan in ['eee','eem','mme','mmm']])
        totErr = sum([errs[chan][i]**2 for chan in ['eee','eem','mme','mmm']])**0.5
        yields['wz'][i] = tot
        errs['wz'][i] = totErr

    if not unblind:
        for c in channels + ['wz']:
            yields[c]['data'] = 0.
            errs[c]['data'] = 0.

    return yields, errs


def getYieldsErrors(analysis,region,period,cut,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',True)
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')

    nl = 3
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_statUnc',mergeDict=mergeDict,scaleFactor=scalefactor,datadriven=doDataDriven)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[region+'datadriven' if doDataDriven else region]])
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    yields = {}
    for chan in ['eee','eem','mme','mmm']:
        yields[chan] = {}
        chanCut = '{0} && channel=="{1}"'.format(cut,chan)
        sig = 0.
        sigErr2 = 0.
        bg = 0.
        bgErr2 = 0.
        for b in plotter.backgrounds:
            val, err = plotter.getNumEntries(chanCut,b,doError=True)
            if b=='WZJets':
                sig += val
                sigErr2 += err**2
            else:
                bg += val
                bgErr2 += err**2
        sigErr = sigErr2**0.5
        bgErr = bgErr2**0.5

        data, dataErr = plotter.getDataEntries(chanCut,doError=True)
        dataErr2 = dataErr**2

        yields[chan]['sig'] = [sig, sigErr]
        yields[chan]['bg'] = [bg, bgErr]
        yields[chan]['data'] = [data, dataErr]

    yields['wz'] = {}
    for i in ['sig','bg','data']:
        tot = sum([yields[chan][i][0] for chan in ['eee','eem','mme','mmm']])
        totErr = sum([yields[chan][i][1]**2 for chan in ['eee','eem','mme','mmm']])**0.5
        yields['wz'][i] = [tot, totErr]

    return yields

def getBackgroundEstimation(analysis,region,period,cut,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',True)
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')

    nl = 3
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    bgChannels = [sigMap[period][x] for x in channelBackground[region+'datadriven' if doDataDriven else region]]

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_bgEstimation',mergeDict=mergeDict,scaleFactor=scalefactor,datadriven=doDataDriven)
    plotter.initializeBackgroundSamples(bgChannels)
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    estimates = {}
    for chan in channels:
        thisCut = '{0} && channel=="{1}"'.format(cut,chan)
        estimates[chan] = {}
        obs, obsErr = plotter.getDataEntries(thisCut,doError=True)
        obsErr2 = obsErr**2
        dd = 0.
        ddErr2 = 0.
        mc = 0.
        mcErr2 = 0.
        for bg in bgChannels + ['datadriven']:
            val, err = plotter.getNumEntries(thisCut,bg,doError=True)
            if bg in ['datadriven']:
                dd += val
                ddErr2 += err**2
            else:
                mc += val
                mcErr2 += err**2
        fromMC = [obs - mc, (obsErr2 + mcErr2)**0.5]
        fromDD = [dd, ddErr2**0.5]
        estimates[chan]['mc'] = fromMC
        estimates[chan]['datadriven'] = fromDD

    estimates['total'] = {}
    estimates['total']['mc'] = [sum([estimates[x]['mc'][0] for x in channels]), sum([estimates[x]['mc'][1]**2 for x in channels])**0.5]
    estimates['total']['datadriven'] = [sum([estimates[x]['datadriven'][0] for x in channels]), sum([estimates[x]['datadriven'][1]**2 for x in channels])**0.5]

    return estimates

def getCrossSections(plotter,cut,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',True)
    sigMap = getSigMap(3)
    analysis = plotter.getAnalysis()
    channel = analysis
    period = plotter.getPeriod()
    channelBackground = getChannelBackgrounds(period)
    if doDataDriven:
        backgroundMC = [x for x in channelBackground[channel] if x not in ['Z','ZG','TT','T','WZ']] + ['datadriven']
    else:
        backgroundMC = [x for x in channelBackground[channel] if x not in ['WZ']]

    yields = {}
    for chan in ['eee','eem','mme','mmm']:
        thisCut = '{0} && channel=="{1}"'.format(cut,chan)
        yields[chan] = {}
        # get signal yield
        wz, wzErr = plotter.getNumEntries(thisCut,sigMap[period]['WZ'],doError=True)
        wzErr2 = wzErr**2
        # get data yield
        data, dataErr = plotter.getDataEntries(thisCut,doError=True)
        dataErr2 = dataErr**2
        # get background
        bg = 0
        bgErr2 = 0
        for b in backgroundMC:
            bgyield = plotter.getNumEntries(thisCut,sigMap[period][b],doError=True)
            bg += bgyield[0]
            bgErr2 += bgyield[1]**2
        # subtract background from data
        extractedYield = data - bg
        extractedErr2 = dataErr2 + bgErr2
        # fill the dict
        yields[chan]['mc'] = [wz, wzErr]
        yields[chan]['data'] = [extractedYield, extractedErr2**0.5]

    return yields


