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

def getDataDrivenCompontents(analysis,region,period,cut,**kwargs):
    baseScaleFactor = kwargs.pop('baseScaleFactor','event.gen_weight*event.pu_weight*event.trig_scale')
    nameMap = {
        0: 'z1.PassMedium1',
        1: 'z1.PassMedium2',
        2: 'w1.PassTight1',
    }
    scaleMap = {
        0: 'z1.LepScaleMedium1',
        1: 'z1.LepScaleMedium2',
        2: 'w1.LepScaleTight1',
    }
    scaleLooseMap = {
        0: 'z1.LepScaleLoose1',
        1: 'z1.LepScaleLoose2',
        2: 'w1.LepScaleLoose1',
    }

    lepScaleMap = {
        'P' : scaleMap,
        'F' : scaleLooseMap,
    }

    nl = 3
    fullScale = '{0}*{1}'.format(baseScaleFactor,'*'.join([scaleMap[i] for i in range(nl)]))
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    if analysis=='WZ': channels = ['eee','eem','mme','mmm']
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    bgs = channelBackground[region+'datadriven']

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_datadrivenComponents',mergeDict=mergeDict,scaleFactor=fullScale,datadriven=True,tightW=True)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in bgs])
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    yields = {}
    errs = {}
    modes = ['PPF','PFP','FPP','PFF','FPF','FFP','FFF']
    for mode in modes:
        yields[mode] = {}
        errs[mode] = {}
        for chan in channels:
            yields[mode][chan] = {}
            errs[mode][chan] = {}
            selection = '{0} && channel=="{1}"'.format(cut,chan)
            mcScale = '{0}*{1}*{2}*{3}'.format(baseScaleFactor,lepScaleMap[mode[0]][0],lepScaleMap[mode[1]][1],lepScaleMap[mode[2]][2])
            for bg in bgs + ['data']:
                if bg in ['datadriven']: continue
                scale = '1' if bg=='data' else mcScale
                singBack = 'data' if bg=='data' else sigMap[period][bg]
                val, err = plotter.getNumEntries(selection,'datadriven',singleComponent=mode,singleBackground=singBack,customScale=scale,doError=True)
                yields[mode][chan][bg] = val
                errs[mode][chan][bg] = err

    for m in modes:
        yields[m]['wz'] = {}
        errs[m]['wz'] = {}
        for i in bgs + ['data']:
            if i=='datadriven': continue
            tot = sum([yields[m][chan][i] for chan in ['eee','eem','mme','mmm']])
            totErr = sum([errs[m][chan][i]**2 for chan in ['eee','eem','mme','mmm']])**0.5
            yields[m]['wz'][i] = tot
            errs[m]['wz'][i] = totErr

    return yields, errs


def getPerChannelYields(analysis,region,period,cut,**kwargs):
    doDataDriven = kwargs.pop('doDataDriven',True)
    tightW = kwargs.pop('tightW',True)
    fakeMode = kwargs.pop('fakeMode','fakerate')
    baseScaleFactor = kwargs.pop('baseScaleFactor','event.gen_weight*event.pu_weight*event.trig_scale')
    skipLepton = kwargs.pop('skipLepton',False)
    skipDataDriven = kwargs.pop('skipDataDriven',False)

    nameMap = {
        0: 'z1.PassTight1',
        1: 'z1.PassTight2',
        2: 'w1.PassTight1',
    }
    scaleMap = {
        0: 'z1.LepScaleTight1',
        1: 'z1.LepScaleTight2',
        2: 'w1.LepScaleTight1',
    }
    scaleLooseMap = {
        0: 'z1.LepScaleLoose1',
        1: 'z1.LepScaleLoose2',
        2: 'w1.LepScaleLoose1',
    }
    promptMap = {
        0: 'z1.GenIsPrompt1',
        1: 'z1.GenIsPrompt2',
        2: 'w1.GenIsPrompt1',
    }
    if tightW:
        nameMap = {
            0: 'z1.PassMedium1',
            1: 'z1.PassMedium2',
            2: 'w1.PassTight1',
        }
        scaleMap = {
            0: 'z1.LepScaleMedium1',
            1: 'z1.LepScaleMedium2',
            2: 'w1.LepScaleTight1',
        }
        scaleLooseMap = {
            0: 'z1.LepScaleLoose1',
            1: 'z1.LepScaleLoose2',
            2: 'w1.LepScaleLoose1',
        }

    lepScaleMap = {
        'P' : scaleMap,
        'F' : scaleLooseMap,
    }

    nl = 3
    fullScale = '{0}*{1}'.format(baseScaleFactor,'*'.join([scaleMap[i] for i in range(nl)]))
    if skipLepton: fullScale = baseScaleFactor
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    if analysis=='WZ': channels = ['eee','eem','mme','mmm']
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    bgs = channelBackground[region+'datadriven' if doDataDriven else region]

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_statUnc',mergeDict=mergeDict,scaleFactor=fullScale,datadriven=doDataDriven,tightW=tightW,fakeMode=fakeMode)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in bgs])
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    yields = {}
    errs = {}
    systs = {}
    for chan in channels:
        yields[chan] = {}
        errs[chan] = {}
        systs[chan] = {}
        chanCut = '{0} && channel=="{1}"'.format(cut,chan)
        for b in bgs + ['datadriven']:
            if skipDataDriven and b=='datadriven':
                val,err,syst = 0,0,0
            else:
                val, err, syst = plotter.getNumEntries(chanCut,sigMap[period][b],doError=True,doSyst=True,baseScaleFactor=baseScaleFactor)
            yields[chan][b] = val
            errs[chan][b] = err
            systs[chan][b] = syst

        data, dataErr = plotter.getDataEntries(chanCut,doError=True)
        yields[chan]['data'] = data
        errs[chan]['data'] = dataErr
        systs[chan]['data'] = 0.

    yields['wz'] = {}
    errs['wz'] = {}
    systs['wz'] = {}
    for i in bgs + ['datadriven'] + ['data']:
        tot = sum([yields[chan][i] for chan in ['eee','eem','mme','mmm']])
        totErr = sum([errs[chan][i]**2 for chan in ['eee','eem','mme','mmm']])**0.5
        totSyst = sum([systs[chan][i]**2 for chan in ['eee','eem','mme','mmm']])**0.5
        yields['wz'][i] = tot
        errs['wz'][i] = totErr
        systs['wz'][i] = totSyst

    return yields, errs, systs

def getPerChannelControlYields(analysis,region,period,cut,**kwargs):
    tightW = kwargs.pop('tightW',True)

    baseScaleFactor = kwargs.pop('baseScaleFactor','event.gen_weight*event.pu_weight*event.trig_scale')

    nameMap = {
        0: 'z1.PassTight1',
        1: 'z1.PassTight2',
        2: 'w1.PassTight1',
    }
    scaleMap = {
        0: 'z1.LepScaleTight1',
        1: 'z1.LepScaleTight2',
        2: 'w1.LepScaleTight1',
    }
    scaleLooseMap = {
        0: 'z1.LepScaleLoose1',
        1: 'z1.LepScaleLoose2',
        2: 'w1.LepScaleLoose1',
    }
    promptMap = {
        0: 'z1.GenIsPrompt1',
        1: 'z1.GenIsPrompt2',
        2: 'w1.GenIsPrompt1',
    }
    if tightW:
        nameMap = {
            0: 'z1.PassMedium1',
            1: 'z1.PassMedium2',
            2: 'w1.PassTight1',
        }
        scaleMap = {
            0: 'z1.LepScaleMedium1',
            1: 'z1.LepScaleMedium2',
            2: 'w1.LepScaleTight1',
        }
        scaleLooseMap = {
            0: 'z1.LepScaleLoose1',
            1: 'z1.LepScaleLoose2',
            2: 'w1.LepScaleLoose1',
        }

    lepScaleMap = {
        'P' : scaleMap,
        'F' : scaleLooseMap,
    }


    nl = 3
    fullScale = '{0}*{1}'.format(baseScaleFactor,'*'.join([scaleMap[i] for i in range(nl)]))
    sigMap = getSigMap(nl)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(period)
    channelBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl)
    if analysis=='WZ': channels = ['eee','eem','mme','mmm']
    saves = '%s_%s_%iTeV' % (analysis,region,period)
    ntuples = getNtupleDirectory(analysis,region,period)

    bgs = channelBackground[region+'datadriven']

    plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=period,rootName='plots_perChannelYields',mergeDict=mergeDict,scaleFactor=fullScale,tightW=tightW)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in bgs])
    plotter.initializeDataSamples([sigMap[period]['data']])
    plotter.setIntLumi(intLumiMap[period])

    yields = {}
    errs = {}
    controls = ['PPP','PPF','PFP','FPP','PFF','FPF','FFP','FFF']
    mcCut = ' && '.join(['fabs(l{0}.GenPatPdgId)<100'.format(l+1) for l in range(nl)])
    for control in controls:
        yields[control] = {}
        errs[control] = {}
        fakechan = 'fakeChannel_tightW' if tightW else 'fakeChannel'
        scalefactor = '{0}*{1}*{2}*{3}'.format(baseScaleFactor,lepScaleMap[control[0]][0],lepScaleMap[control[1]][1],lepScaleMap[control[2]][2])
        for chan in channels:
            yields[control][chan] = {}
            errs[control][chan] = {}
            chanCut = '{0} && channel=="{1}"'.format(cut,chan)
            controlCut = '{0} && {1}=="{2}"'.format(chanCut,fakechan,control)
            mcFullCut = '{0} && {1}'.format(controlCut,mcCut)
            for b in bgs:
                val, err = plotter.getNumEntries(mcFullCut,sigMap[period][b],customScale=scalefactor,doError=True)
                yields[control][chan][b] = val
                errs[control][chan][b] = err
            yields[control][chan]['datadriven'] = 0.
            errs[control][chan]['datadriven'] = 0.

            data, dataErr = plotter.getDataEntries(controlCut,doError=True)
            yields[control][chan]['data'] = data
            errs[control][chan]['data'] = dataErr

    for control in controls:
        yields[control]['wz'] = {}
        errs[control]['wz'] = {}
        for i in bgs + ['datadriven'] + ['data']:
            tot = sum([yields[control][chan][i] for chan in ['eee','eem','mme','mmm']])
            totErr = sum([errs[control][chan][i]**2 for chan in ['eee','eem','mme','mmm']])**0.5
            yields[control]['wz'][i] = tot
            errs[control]['wz'][i] = totErr

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


