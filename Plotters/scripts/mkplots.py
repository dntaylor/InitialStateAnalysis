#!/usr/bin/env python

import argparse
import itertools
import sys
import pickle
import json
import logging
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.ShapePlotter import ShapePlotter
from InitialStateAnalysis.Plotters.CutFlowPlotter import CutFlowPlotter
from InitialStateAnalysis.Plotters.EfficiencyPlotter import EfficiencyPlotter
from InitialStateAnalysis.Plotters.FakeRatePlotter import FakeRatePlotter
from InitialStateAnalysis.Plotters.CorrelationPlotter import CorrelationPlotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES
from InitialStateAnalysis.Utilities.utilities import *


def plotRegion(analysis,channel,runPeriod,**kwargs):
    '''A function to simplify plotting multiple channels and run periods.'''
    logger = logging.getLogger(__name__)
    blind = kwargs.pop('blind',True)
    mass = kwargs.pop('mass',500)
    runTau = kwargs.pop('runTau',False)
    myCut = kwargs.pop('myCut','1')
    plotFinalStates = kwargs.pop('plotFinalStates',False)
    plotJetBins = kwargs.pop('plotJetBins',False)
    plotOverlay = kwargs.pop('plotOverlay',False)
    plotShapes = kwargs.pop('plotShapes',False)
    plotCutFlow = kwargs.pop('plotCutFlow',False)
    plotFakeRegions = kwargs.pop('plotFakeRegions',False)
    finalStatesToPlot = kwargs.pop('finalStates','all')
    nostack = kwargs.pop('nostack',False)
    normalize = kwargs.pop('normalize',False)
    doDetailed = kwargs.pop('doDetailed',False)
    scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    useSignal = analysis in ['Hpp3l','Hpp4l']
    loglevel = kwargs.pop('loglevel','INFO')
    for key, value in kwargs.iteritems():
        logger.warning("Unrecognized parameter '" + key + "' = " + str(value))
        return 0

    if useSignal: logger.info("%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass))
    isControl = analysis != channel
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

    finalStates, leptons = getChannels(nl,runTau=runTau)
    if analysis in ['WZ']: finalStates = ['eee','eem','mme','mmm']
    if finalStatesToPlot=='all':
        fsToPlot = finalStates
    else:
        fsToPlot = finalStatesToPlot.split(',')
    logger.info('%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut))
    dataplot = (isControl or not blind)
    mergeDict = getMergeDict(runPeriod)
    cutFlowMap = {}
    cutFlowMap[channel] = defineCutFlowMap(channel,finalStates,mass)

    genChannels = {
        'Hpp3l' : {
            'ee': ['eee','eem','eet'],
            'em': ['eme','emm','emt'],
            'et': ['ete','etm','ett'],
            'mm': ['mme','mmm','mmt'],
            'mt': ['mte','mtm','mtt'],
            'tt': ['tte','ttm','ttt'],
        },
        'Hpp4l' : {
            'ee': ['eeee'],
            'em': ['emem'],
            'et': ['etet'],
            'mm': ['mmmm'],
            'mt': ['mtmt'],
            'tt': ['tttt'],
        },
    }

    hppChannels = ['ee','em','et','mm','mt','tt']
    hpChannels = ['e','m','t']

    recoChannels = {
        'Hpp3l' : {
            'ee': ['eee','eem'],
            'em': ['eme','emm','mee','mem'],
            'mm': ['mme','mmm'],
            'et': ['eee','eme','eem','emm','mee','mem'],
            'mt': ['mee','mem','mme','mmm','eme','emm'],
            'tt': ['eee','eem','eme','emm','mee','mem','mme','mmm'],
        },
        'Hpp4l' : {
            'ee': ['eeee'],
            'em': ['emem','emme','meem','meme'],
            'mm': ['mmmm'],
            'et': ['eeee','eeem','eeme','emee','emem','emme','meee','meem','meme'],
            'mt': ['emem','emme','emmm','meem','meme','memm','mmem','mmme','mmmm'],
            'tt': ['eeee','eeem','eeme','eemm','emee','emem','emme','emmm','meee','meem','meme','memm','mmee','mmem','mmme','mmmm'],
        },
    }

    customFinalStates = {
        'Hpp3l' : {},
        'Hpp4l' : {},
    }

    numTauCuts = {
        'Hpp3l' : {},
        'Hpp4l' : {},
    }

    for nt in range(3):
        theCut = '(' + ' || '.join(['genChannel=="%s"' %gChan for hChan in hppChannels for gChan in genChannels['Hpp3l'][hChan] if hChan.count('t')==nt]) + ')'
        numTauCuts['Hpp3l'][nt] = theCut
    numTauCuts['Hpp4l'][0] = '(' + ' || '.join(['genChannel=="%s"' % gChan for gChan in ['eeee','eeem','eemm','emee','emem','emmm','mmee','mmem','mmmm']]) + ')'
    numTauCuts['Hpp4l'][1] = '(' + ' || '.join(['genChannel=="%s"' % gChan for gChan in ['etet','etmt','mtet','mtmt']]) + ')'
    numTauCuts['Hpp4l'][2] = '(' + ' || '.join(['genChannel=="%s"' % gChan for gChan in ['tttt']]) + ')'

    for a in ['Hpp3l','Hpp4l']:
        for c in hppChannels:
            genCut = '(' + ' || '.join(['genChannel=="%s"'%x for x in genChannels[a][c] + ['aaa']]) + ')'
            recoCut = '(' + ' || '.join(['channel=="%s"'%x for x in recoChannels[a][c]]) + ')'
            customFinalStates[a][c] = genCut + ' && ' + recoCut

    # plot efficiencies
    if analysis in ['Hpp3l','Hpp4l']:
        logger.info("%s:%s:%iTeV: Plotting efficiency" % (analysis, channel, runPeriod))
        plotter = EfficiencyPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_efficiency',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
        plotMethod = getattr(plotter,plotMode)
        plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'efficiency',labels=cutFlowMap[channel]['labels'],lumitext=33)
        if plotFinalStates:
            for c in fsToPlot:
                logger.info("%s:%s:%iTeV: Plotting efficiency  %s" % (analysis, channel, runPeriod, c))
                plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/efficiency'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)

        ## plot cut flows overlays
        #plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflow_overlay',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        #plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel] if x not in ['WZ']])
        #plotter.initializeSignalSamples([sigMap[runPeriod]['WZ']])
        #if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        #plotter.setIntLumi(intLumiMap[runPeriod])
        #plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
        #plotMethod = getattr(plotter,plotMode)
        #plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'cutFlow_overlay',labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)
        #if plotFinalStates:
        #    for c in fsToPlot:
        #        logger.info("%s:%s:%iTeV: Plotting cut flow overlay  %s" % (analysis, channel, runPeriod, c))
        #        plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow_overlay'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)

    # plotting correlation
    if analysis in ['Hpp3l', 'Hpp4l']:
        plotter = CorrelationPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName='plots_correlation',loglevel=loglevel)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])

        logger.info('%s:%s:%iTeV: Plotting correlation' % (analysis, channel, runPeriod))
        plotMethod = getattr(plotter,'plotCorrelation')
        plotMethod(cutFlowMap[channel]['cuts'][1:], 'correlation/mc', cut=myCut, labels=cutFlowMap[channel]['labels'][1:], plottype='mc')
        if useSignal: plotMethod(cutFlowMap[channel]['cuts'][1:], 'correlation/sig', cut=myCut, labels=cutFlowMap[channel]['labels'][1:], plottype='sig')

    # do variables on same plot
    if useSignal:
        logger.info("%s:%s:%iTeV: Plotting signal" % (analysis, channel, runPeriod))
        plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_signal',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        masses = _3L_MASSES if nl==3 else _4L_MASSES
        plotter.initializeSignalSamples([sigMap[runPeriod][x] for x in masses])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotSignal'
        plotMethod = getattr(plotter,plotMode)
        plotMethod('h1.mass',       [1000,0,1000],'signal/hppMass',          yaxis='A.U.',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     legendpos=43,cut=myCut,logy=0, normalize=1)
        plotMethod('h2.mass',       [1000,0,1000],'signal/hmmMass',          yaxis='A.U.',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     legendpos=43,cut=myCut,logy=0, normalize=1)
        plotMethod('h1.dPhi',       [100,0,5],    'signal/hppDphi',          yaxis='A.U.',xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',lumitext=33,logy=0,cut=myCut,normalize=1)
        plotMethod('h1.dR',         [100,0,6.28], 'signal/hppDR',            yaxis='A.U.',xaxis='\\Delta R_{\\ell^{\\pm}\\ell^{\\pm}}',         lumitext=33,logy=0,cut=myCut,normalize=1)
        plotMethod('finalstate.sT', [250,0,2000], 'signal/sT',               yaxis='A.U.',xaxis='S_{T} (GeV/c^{2})',                            legendpos=43,logy=0,cut=myCut,overflow=True,normalize=1)
        plotMethod('z1.mass',       [250,0,1000], 'signal/z1Mass_fullWindow',yaxis='A.U.',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             legendpos=43,logy=0,cut=myCut,overflow=True,normalize=1)
        plotMethod('finalstate.met',[200,0,1000], 'signal/met',              yaxis='A.U.',xaxis='E_{T}^{miss} (GeV/c^{2})',                     legendpos=43,logy=0,cut=myCut,overflow=True,normalize=1)
        plotMethod('abs(h1.mass-h2.mass)', [200,0,400],'signal/massdiff',    yaxis='A.U.',xaxis='|M_{\\ell^{+}\\ell^{+}}-M_{\\ell^{-}\\ell^{-}}| (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,normalize=1)
        for nt in range(3):
            if analysis not in ['Hpp3l','Hpp4l']: continue
            theCut = numTauCuts[analysis][nt] + ' && ' + myCut
            plotMethod('h1.mass',       [1000,0,1000],'signal/hppMass_%iTau'%nt,          yaxis='A.U.',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     legendpos=43,cut=theCut,logy=0, normalize=1)
            plotMethod('h2.mass',       [1000,0,1000],'signal/hmmMass_%iTau'%nt,          yaxis='A.U.',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     legendpos=43,cut=theCut,logy=0, normalize=1)
            plotMethod('h1.dPhi',       [100,0,5],    'signal/hppDphi_%iTau'%nt,          yaxis='A.U.',xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',lumitext=33,logy=0,cut=theCut,normalize=1)
            plotMethod('h1.dR',         [100,0,6.28], 'signal/hppDR_%iTau'%nt,            yaxis='A.U.',xaxis='\\Delta R_{\\ell^{\\pm}\\ell^{\\pm}}',         lumitext=33,logy=0,cut=theCut,normalize=1)
            plotMethod('finalstate.sT', [500,0,2000], 'signal/sT_%iTau'%nt,               yaxis='A.U.',xaxis='S_{T} (GeV/c^{2})',                            legendpos=43,logy=0,cut=theCut,overflow=True,normalize=1)
            plotMethod('z1.mass',       [250,0,1000], 'signal/z1Mass_fullWindow_%iTau'%nt,yaxis='A.U.',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             legendpos=43,logy=0,cut=theCut,overflow=True,normalize=1)
            plotMethod('finalstate.met',[200,0,1000], 'signal/met_%iTau'%nt,              yaxis='A.U.',xaxis='E_{T}^{miss} (GeV/c^{2})',                     legendpos=43,logy=0,cut=theCut,overflow=True,normalize=1)
            plotMethod('abs(h1.mass-h2.mass)', [200,0,400],'signal/massdiff_%iTau'%nt,   yaxis='A.U.',xaxis='|M_{\\ell^{+}\\ell^{+}}-M_{\\ell^{-}\\ell^{-}}| (GeV/c^{2})',lumitext=33,logy=0,cut=theCut,overflow=True,normalize=1)
        #plotter.initializeSignalSamples([allSigMap[runPeriod][mass]])
        #plotter.setIntLumi(intLumiMap[runPeriod])
        #plotMode = 'plotSignal'
        #plotMethod = getattr(plotter,plotMode)
        #plotMethod('h1.mass',[400,300,700],'signal/hppMass_cut',yaxis='Events/0.1 GeV/c^{2}',xaxis='M_{\\ell^{+}\\ell^{+}} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,boxes=[[450,550,2]])
        names = {
            'e': 'Elec',
            'm': 'Mu',
            't': 'Tau',
        }
        tex = {
            'e': 'e',
            'm': '\\mu',
            't': '\\tau',
        }
        for l in ['e','m']:
            name = names[l]
            t = tex[l]
            cuts = ['%s & %s' %(myCut,'l%iFlv=="%s"' %((x+1),l)) for x in range(nl)]
            plotMethod(['l%i.Pt'  %(x+1) for x in range(nl)], [100,0,1000], 'signal/%sPt'%name, yaxis='A.U.', xaxis='p_{T}^{%s} (GeV)' %t, legendpos=43, logy=0, cut=cuts, overflow=True, normalize=1)
            for nt in range(3):
                if analysis not in ['Hpp3l','Hpp4l']: continue
                theCuts = [numTauCuts[analysis][nt] + ' && ' + c for c in cuts]
                plotMethod(['l%i.Pt'  %(x+1) for x in range(nl)], [100,0,1000], 'signal/%sPt_%iTau'%(name,nt), yaxis='A.U.', xaxis='p_{T}^{%s} (GeV)' %t, legendpos=43, logy=0, cut=theCuts, overflow=True, normalize=1)




    # Plotting discriminating variables
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName='plots_2d',loglevel=loglevel)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])


    if analysis in ['Hpp3l'] and not blind:
        plotter.plotMCDataSignalRatio2D('z1.mass',       'finalstate.mass',[70,0,140], [100,0,200],'m3l_v_z1_mc',     xaxis='M_{\\ell^{+}\\ell^{-}}',                 yaxis='M_{3\\ell}',        cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('z1.mass',       'finalstate.mass',[70,0,140], [100,0,200],'m3l_v_z1_data',   xaxis='M_{\\ell^{+}\\ell^{-}}',                 yaxis='M_{3\\ell}',        cut=myCut, plotdata=1, plotmc=0, plotsig=0)
    if analysis in ['Hpp3l']:
        #plotter.setScaleFactor('event.lep_scale*event.trig_scale')
        plotter.plotMCDataSignalRatio2D('h1.mass',       'event.nvtx',[24,0,600], [50,0,50],'hppMass_v_pu_sig',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('z1.mass',       'event.nvtx',[80,0,240], [50,0,50],'z1Mass_v_pu_sig', xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('h1.dPhi',       'event.nvtx',[32,0,3.2], [50,0,50],'hppDphi_v_pu_sig',xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('finalstate.sT', 'event.nvtx',[40,0,1000],[50,0,50],'sT_v_pu_sig',     xaxis='S_{T} (GeV/c^{2})',                            yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('finalstate.met','event.nvtx',[40,0,200], [50,0,50],'met_v_pu_sig',    xaxis='E_{T}^{miss} (GeV)',                           yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('h1.mass',       'event.nvtx',[24,0,600], [50,0,50],'hppMass_v_pu_mc', xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('z1.mass',       'event.nvtx',[80,0,240], [50,0,50],'z1Mass_v_pu_mc',  xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('h1.dPhi',       'event.nvtx',[32,0,3.2], [50,0,50],'hppDphi_v_pu_mc', xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('finalstate.sT', 'event.nvtx',[40,0,1000],[50,0,50],'sT_v_pu_mc',      xaxis='S_{T} (GeV/c^{2})',                            yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('finalstate.met','event.nvtx',[40,0,200], [50,0,50],'met_v_pu_mc',     xaxis='E_{T}^{miss} (GeV)',                           yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)

    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])

    plotMode = 'plotMCDataRatio' if dataplot else 'plotMC'
    plotMethod = getattr(plotter,plotMode)
    logger.info("%s:%s:%iTeV: Plotting discriminating variables" % (analysis,channel, runPeriod))
    plotDistributions(plotMethod,myCut,nl,isControl,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)

    # each channel
    if plotFinalStates:
        logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
        for c in fsToPlot:
            logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
            plotDistributions(plotMethod,myCut+'&&channel=="%s"'%c,nl,isControl,savedir=c,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
               sel = customFinalStates[analysis][c]
               logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
               plotDistributions(plotMethod,'%s & %s'%(myCut,sel),nl,isControl,savedir=c,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)

    if plotFakeRegions:
        logger.info("%s:%s:%iTeV: Plotting control regions" % (analysis, channel, runPeriod))
        labels = {0:'F',1:'T'}
        cuts = {0:'!{0}.PassTight{1}',1:'{0}.PassTight{1}'}
        for z1 in [0,1]:
            for z2 in [0,1]:
                for w1 in [0,1]:
                    thisCut = '{0} && '.format(myCut) + ' && '.join([cuts[z1].format('z1','1'),cuts[z2].format('z1','2'),cuts[w1].format('w1','1')])
                    thisName = ''.join([labels[z1],labels[z2],labels[w1]])
                    logger.info("%s:%s:%iTeV: Plotting control regions - %s" % (analysis, channel, runPeriod, thisName))
                    plotDistributions(plotMethod,thisCut,nl,isControl,savedir=thisName,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed,doMinimal=True)
                    if plotFinalStates:
                        for c in fsToPlot:
                            logger.info("%s:%s:%iTeV: Plotting control regions - %s - channel %s" % (analysis, channel, runPeriod, thisName, c))
                            newCut = '{0} && channel=="{1}"'.format(thisCut,c)
                            plotDistributions(plotMethod,newCut,nl,isControl,savedir='{0}/{1}'.format(thisName,c),analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed,doMinimal=True)

    # some partially blind plots for h++
    if runPeriod==8 and not dataplot and analysis in ['Hpp3l']:
        logger.info("%s:%s:%iTeV: Plotting partially blinded variables" % (analysis, channel, runPeriod))
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotModeUnblind = 'plotMCDataRatio'
        plotMethodUnblind = getattr(plotter,plotModeUnblind)
        plotMethodUnblind('h1.mass',      [24,0,600], 'hppMass_unblind',          yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',lumitext=33,logy=1,cut=myCut,overflow=True,blinder=[150,99999])
        plotMethodUnblind('finalstate.sT',[50,0,500], 'sT_unblind',               yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',                       lumitext=33,logy=0,cut=myCut,overflow=True,blinder=[200,99999])
        plotMethodUnblind('z1.mass',      [42,70,112],'z1Mass_unblind',           yaxis='Events/1.0 GeV',       xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',        legendpos=43,logy=0,cut=myCut,blinder=[112,99999])
        plotMethodUnblind('z1.mass',      [80,0,240], 'z1Mass_fullWindow_unblind',yaxis='Events/3.0 GeV',       xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',        legendpos=43,logy=0,cut=myCut,blinder=[112,99999])
        if plotFinalStates:
            logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
            for c in fsToPlot:
                logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                plotMethodUnblind('h1.mass',      [24,0,600], c+'/hppMass_unblind',          yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',lumitext=33,logy=1,cut=myCut+'&&channel=="%s"'%c,overflow=True,blinder=[150,99999])
                plotMethodUnblind('finalstate.sT',[50,0,500], c+'/sT_unblind',               yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',                       lumitext=33,logy=0,cut=myCut+'&&channel=="%s"'%c,overflow=True,blinder=[200,99999])
                plotMethodUnblind('z1.mass',      [42,70,112],c+'/z1Mass_unblind',           yaxis='Events/1.0 GeV',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',               legendpos=43,logy=0,cut=myCut+'&&channel=="%s"'%c,blinder=[112,99999])
                plotMethodUnblind('z1.mass',      [80,0,240], c+'/z1Mass_fullWindow_unblind',yaxis='Events/3.0 GeV',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',               legendpos=43,logy=0,cut=myCut+'&&channel=="%s"'%c,blinder=[112,99999])

    # setup signal overlay plots
    if useSignal:
        plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_overlay',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotMCDataSignalRatio' if dataplot else 'plotMCSignalRatio'
        plotMethod = getattr(plotter,plotMode)
        if analysis in ['Hpp3l']:
            plotter.plotMCDataSignalRatio2D('h1.mass','h1.dPhi', [80,0,800], [32,0,3.2], 'h1mass_v_h1dphi_mc',  plotdata=0, plotmc=1, plotsig=0, cut=myCut, xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}}',yaxis='#Delta#phi_{\\ell^{\\pm}\\ell^{\\pm}}')
            plotter.plotMCDataSignalRatio2D('h1.mass','h1.dPhi', [80,0,800], [32,0,3.2], 'h1mass_v_h1dphi_sig', plotdata=0, plotmc=0, plotsig=1, cut=myCut, xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}}',yaxis='#Delta#phi_{\\ell^{\\pm}\\ell^{\\pm}}')
    if plotOverlay and useSignal:
        # plot the signal overlay
        logger.info("%s:%s:%iTeV: Plotting signal overlay discriminating variables" % (analysis, channel, runPeriod))
        plotDistributions(plotMethod,myCut,nl,isControl,savedir='overlay',signalscale=100,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)


    # plot shapes
    if plotShapes:
        logger.info("%s:%s:%iTeV: Plotting shapes" % (analysis, channel, runPeriod))
        plotter = ShapePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_shapes',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotter.plotMC('z1.mass',['channel=="mmm"','channel=="emm"'],[42,70,112],'zMass_mc_mm',yaxis='Normalized',xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,cutNames=['mmm','emm'])
        if dataplot: plotter.plotData('z1.mass',['channel=="mmm"','channel=="emm"'],[42,70,112],'zMass_data_mm',yaxis='Normalized',xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,cutNames=['mmm','emm'])
        plotter.plotMC('z1.mass',['channel=="eee"','channel=="eme"'],[42,70,112],'zMass_mc_ee',yaxis='Normalized',xaxis='M(l^{+}l^{-}) (GeV)',logy=0,cut=myCut,cutNames=['eee','eme'])
        if dataplot: plotter.plotData('z1.mass',['channel=="eee"','channel=="eme"'],[42,70,112],'zMass_data_ee',yaxis='Normalized',xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,cutNames=['eee','eme'])

    # plot cut flows (each cut)
    logger.info("%s:%s:%iTeV: Plotting cut flow" % (analysis, channel, runPeriod))
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutFlowSelections',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotMCDataRatio' if dataplot else 'plotMC'
    plotMethod = getattr(plotter,plotMode)
    if plotCutFlow:
        for i in range(len(cutFlowMap[channel]['cuts'])):
            logger.info('%s:%s:%iTeV: Plotting cut flow selections %s' % (analysis, channel, runPeriod, cutFlowMap[channel]['labels_simple'][i]))
            thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
            plotDistributions(plotMethod,'%s&%s'%(myCut,thisCut),nl,isControl,savedir='cutflow/%s'%cutFlowMap[channel]['labels_simple'][i],analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)
            if plotFinalStates:
                logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
                for c in fsToPlot:
                    logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                    plotDistributions(plotMethod,'%s&channel=="%s"&%s'%(myCut,c,thisCut),nl,isControl,savedir='cutflow/%s/%s' %(cutFlowMap[channel]['labels_simple'][i],c),analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)
            thisCut = cutFlowMap[channel]['cuts'][i]
            plotDistributions(plotMethod,'%s&%s'%(myCut,thisCut),nl,isControl,savedir='cutflow/%s_only'%cutFlowMap[channel]['labels_simple'][i],analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)
            if plotFinalStates:
                logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
                for c in fsToPlot:
                    logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                    plotDistributions(plotMethod,'%s&channel=="%s"&%s'%(myCut,c,thisCut),nl,isControl,savedir='cutflow/%s_only/%s' %(cutFlowMap[channel]['labels_simple'][i],c),analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass,doDetailed=doDetailed)

    # plot cut flows on same plot
    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflow',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    if useSignal: plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
    plotMethod = getattr(plotter,plotMode)
    plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'cutFlow',labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)
    if plotFinalStates:
        for c in fsToPlot:
            logger.info("%s:%s:%iTeV: Plotting cut flow  %s" % (analysis, channel, runPeriod, c))
            plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
               sel = customFinalStates[analysis][c]
               logger.info("%s:%s:%iTeV: Plotting cut flow  %s" % (analysis, channel, runPeriod, c))
               plotMethod(['%s&& %s &&%s' %(x,sel,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)


    # setup individual channel cuts on same plot
    plotChannelStrings, plotChannelCuts = getChannelStringsCuts(channel,finalStates)
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    if useSignal: plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
    plotMethod = getattr(plotter,plotMode)
    plotMethod([myCut]+['%s&&%s' %(x,myCut) for x in plotChannelCuts],'individualChannels',labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0,signalscale=1000,numcol=2)
    if plotFinalStates:
        logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
        if analysis in customFinalStates:
            plotGenChannelStrings, plotGenChannelCuts = getGenChannelStringsCuts(channel,customFinalStates[analysis])
            for c in fsToPlot:
                logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                plotMethod(['%s&&channel=="%s"'%(myCut,c)]+['%s&&%s&&channel=="%s"' %(x,myCut,c) for x in plotGenChannelCuts],'%s/individualGenChannels'%c,labels=['Total']+plotGenChannelStrings,nosum=True,lumitext=33,logy=0,signalscale=1000,numcol=2)
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
               sel = customFinalStates[analysis][c]
               logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
               plotMethod(['%s&&%s' %(myCut,sel)]+['%s&&%s&&%s' %(x,myCut,sel) for x in plotChannelCuts],'%s/individualChannels'%c,labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0,signalscale=1000,numcol=2)


    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflowSelectionsChannels',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    plotMethod = getattr(plotter,plotMode)
    for i in range(len(cutFlowMap[channel]['cuts'])):
        thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
        plotMethod(['%s&&%s'%(myCut,thisCut)]+['%s&&%s&&%s' %(x,myCut,thisCut) for x in plotChannelCuts],'cutflow/%s/individualChannels'%cutFlowMap[channel]['labels_simple'][i],labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0,numcol=2)


def plotFakeRate(analysis,channel,runPeriod,**kwargs):
    '''Plot fake rate for an analysis.'''
    loglevel = kwargs.pop('loglevel','INFO')
    logger = logging.getLogger(__name__)
    blind = kwargs.pop('blind',True)
    mass = kwargs.pop('mass',500)
    runTau = kwargs.pop('runTau',False)
    myCut = kwargs.pop('myCut','1')
    plotFinalStates = kwargs.pop('plotFinalStates',False)
    plotJetBins = kwargs.pop('plotJetBins',False)
    plotOverlay = kwargs.pop('plotOverlay',False)
    plotShapes = kwargs.pop('plotShapes',False)
    plotCutFlow = kwargs.pop('plotCutFlow',False)
    finalStatesToPlot = kwargs.pop('finalStates','all')
    nostack = kwargs.pop('nostack',False)
    normalize = kwargs.pop('normalize',False)
    scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    useSignal = analysis in ['Hpp3l','Hpp4l']
    doDetailed = kwargs.pop('doDetailed',False)
    for key, value in kwargs.iteritems():
        logger.warning("Unrecognized parameter '" + key + "' = " + str(value))
        return 0

    dataDriven = not blind
    if useSignal: logger.info("%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass))
    isControl = analysis != channel
    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ_W' : 2,
        'WZ_Dijet' : 1,
        'WZ'   : 3,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[analysis]
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = getChannelBackgrounds(runPeriod)

    finalStates, leptons = getChannels(nl,runTau=runTau)
    if finalStatesToPlot=='all':
        fsToPlot = finalStates
    else:
        fsToPlot = finalStatesToPlot.split(',')
    logger.info('%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut))
    dataplot = (isControl or not blind)
    mergeDict = getMergeDict(runPeriod)

    # TODO: load from pickle file
    # define fake regions
    lepName = {'e': 'Elec', 'm': 'Muon', 't': 'Tau'}
    fakeRegions = {}
    fakeRegions['WZ'] = {}
    for f in ['e', 'm']:
        #for p in ['Loose', 'Tight']:
        for p in ['Tight']:
            # select leading Z pt, Z window [60,120], tight (or loose) Z, low met, m3l>100, w1 mass < 30
            if analysis in ['WZ']:
                #for z in ['Loose', 'Tight']:
                for z in ['Tight']:
                    fakeRegion = 'Z{0}Probe{1}{2}'.format(z,lepName[f],p)
                    #denom = 'z1.Pass{0}1 && z1.Pass{0}2 && finalstate.mass>100. && w1.mass<20. && w1.dR1_z1_1>0.1 && w1.dR1_z1_2>0.1 && w1Flv=="{1}"'.format(z,f)
                    denom = 'z1.Pass{0}1 && z1.Pass{0}2 && w1.dR1_z1_1>0.02 && w1.dR1_z1_2>0.02 && w1Flv=="{1}"'.format(z,f)
                    #denom = 'z1.Pass{0}1 && z1.Pass{0}2 && w1Flv=="{1}"'.format(z,f)
                    numer = '{0} & w1.Pass{1}1'.format(denom,p)
                    fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
                    if p=='Tight' and channel not in ['HZZFakeRate']:
                       fakeRegion += '_LooseProbe'
                       denom += ' && w1.PassLoose1'
                       numer += ' && w1.PassLoose1'
                       fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
            # select w lepton pt, z veto, met
            #'W' : 'w1.Pt1>20. & (z1.mass<60. | z1.mass>120.) & finalstate.met>30. & w1.mass>30.',
            if analysis in ['WZ_W']:
                #for w in ['Loose','Tight']:
                for w in ['Tight']:
                    fakeRegion = 'W{0}Probe{1}{2}'.format(w,lepName[f],p)
                    denom = 'w1.Pt1>20. && w1.mass>30. && finalstate.met>30. && (z1.mass<60. || z1.mass>120.) && l1.Chg==l2.Chg && z1.dR>0.1 && w1.Pass{0}1 && w2Flv=="{1}"'.format(w,f)
                    numer = '{0} && w2.Pass{1}1'.format(denom,p)
                    fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
                    if p=='Tight' and channel not in ['HZZFakeRate']:
                       fakeRegion += '_LooseProbe'
                       denom += ' && w2.PassLoose1'
                       numer += ' && w2.PassLoose1'
                       fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
            # ntuple cuts: zVeto 60-120, met vet 20, w veto 20, jet pt > 20, jet dr > 1.0
            if analysis in ['WZ_Dijet']:
               fakeRegion = 'FakeRateProbe{0}{1}'.format(lepName[f],p)
               denom = 'l1Flv=="{0}"'.format(f)
               numer = '{0} && w1.Pass{1}1'.format(denom,p)
               fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
               if p=='Tight' and channel not in ['HZZFakeRate']:
                  fakeRegion += '_LooseProbe'
                  denom += ' && w1.PassLoose1'
                  numer += ' && w1.PassLoose1'
                  fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}



    # setup selections
    ptBins = [0,5,7,10,20,30,40,50,80]
    etaBins = {
        'e': [0,1.479,2.5],
        #'m': [0,1.2,2.4],
        'm': [0,1.479,2.4],
    }
    for fakeRegion in fakeRegions['WZ']:
        logger.info("%s:%s:%iTeV: Fake Region: %s" % (analysis,channel, runPeriod, fakeRegion))
        denom = fakeRegions['WZ'][fakeRegion]['denom']
        numer = fakeRegions['WZ'][fakeRegion]['numer']
        probe = fakeRegions['WZ'][fakeRegion]['probe']
        ptvar = fakeRegions['WZ'][fakeRegion]['ptVar']
        etavar = fakeRegions['WZ'][fakeRegion]['etaVar']
        ptcut = '{0} >= {1} && {0} < {2}'
        etacut = 'abs({0}) >= {1} && abs({0}) < {2}'

        #plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,rootName='{0}_fakeplots'.format(fakeRegion),scaleFactor='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale*event.trig_prescale',dataScaleFactor='event.trig_prescale')
        #plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        #plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        #plotter.setIntLumi(intLumiMap[runPeriod])
        #plotMode = 'plotMCData'
        #plotMethod = getattr(plotter,plotMode)

        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: All Probes" % (analysis,channel, runPeriod))
        #plotDistributions(plotMethod,denom,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_all'.format(fakeRegion),doDetailed=doDetailed)
        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: All Probes - barrel" % (analysis,channel, runPeriod))
        #thiscut = '{0} && {1}'.format(denom,etacut.format(etavar,etaBins[probe][0],etaBins[probe][1]))
        #plotDistributions(plotMethod,thiscut,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_all/barrel'.format(fakeRegion),doDetailed=doDetailed)
        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: All Probes - endcap" % (analysis,channel, runPeriod))
        #thiscut = '{0} && {1}'.format(denom,etacut.format(etavar,etaBins[probe][1],etaBins[probe][2]))
        #plotDistributions(plotMethod,thiscut,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_all/endcap'.format(fakeRegion),doDetailed=doDetailed)

        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: Passing" % (analysis,channel, runPeriod))
        #plotDistributions(plotMethod,numer,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_prompts'.format(fakeRegion),doDetailed=doDetailed)
        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: Passing - barrel" % (analysis,channel, runPeriod))
        #thiscut = '{0} && {1}'.format(numer,etacut.format(etavar,etaBins[probe][0],etaBins[probe][1]))
        #plotDistributions(plotMethod,thiscut,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_prompts/barrel'.format(fakeRegion),doDetailed=doDetailed)
        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: Passing - endcap" % (analysis,channel, runPeriod))
        #thiscut = '{0} && {1}'.format(numer,etacut.format(etavar,etaBins[probe][1],etaBins[probe][2]))
        #plotDistributions(plotMethod,thiscut,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_prompts/endcap'.format(fakeRegion),doDetailed=doDetailed)

        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: Failing" % (analysis,channel, runPeriod))
        #plotDistributions(plotMethod,'{0} & !({1})'.format(denom,numer),nl,isControl,analysis=analysis,savedir='fakeRate/{0}_fakes'.format(fakeRegion),doDetailed=doDetailed)
        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: Failing - barrel" % (analysis,channel, runPeriod))
        #thiscut = '{0} && {1}'.format('{0} & !({1})'.format(denom,numer),etacut.format(etavar,etaBins[probe][0],etaBins[probe][1]))
        #plotDistributions(plotMethod,thiscut,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_fakes/barrel'.format(fakeRegion),doDetailed=doDetailed)
        #logger.info("%s:%s:%iTeV: Plotting discriminating variables: Failing - endcap" % (analysis,channel, runPeriod))
        #thiscut = '{0} && {1}'.format('{0} & !({1})'.format(denom,numer),etacut.format(etavar,etaBins[probe][1],etaBins[probe][2]))
        #plotDistributions(plotMethod,thiscut,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_fakes/endcap'.format(fakeRegion),doDetailed=doDetailed)

        # now plot the fake rates
        logger.info("%s:%s:%iTeV: Computing fake rates" % (analysis,channel, runPeriod))
        plotter = FakeRatePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='{0}_fakerates'.format(fakeRegion),mergeDict=mergeDict,scaleFactor='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale*event.trig_prescale',dataScaleFactor='event.trig_prescale')
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotter.plotFakeRate(numer, denom, 'fakeRate/{0}_fakerate'.format(fakeRegion), ptBins=ptBins, etaBins=etaBins[probe], logx=1, ptVar=ptvar, etaVar=etavar, dataDriven=dataDriven)
        plotter.plotFakeRateProjection(numer, denom, 'fakeRate/{0}_pt_fakerate'.format(fakeRegion), 'pt', ptBins=ptBins, etaBins=etaBins[probe], logx=0, ptVar=ptvar, etaVar=etavar, dataDriven=dataDriven, xaxis='p_{T} (GeV)')
        plotter.plotFakeRateProjection(numer, denom, 'fakeRate/{0}_eta_fakerate'.format(fakeRegion), 'eta', ptBins=ptBins, etaBins=etaBins[probe], logx=0, ptVar=ptvar, etaVar=etavar, dataDriven=dataDriven, xaxis='\eta')

def parse_command_line(argv):
    parser = get_parser("Plot a given channel and period")

    parser.add_argument('-pf','--plotFinalStates',action='store_true',help='Plot individual final states')
    parser.add_argument('-pj','--plotJetBins',action='store_true',help='Plot jet bins')
    parser.add_argument('-po','--plotOverlay',action='store_true',help='Plot overlay')
    parser.add_argument('-ps','--plotShapes',action='store_true',help='Plot shapes')
    parser.add_argument('-pcf','--plotCutFlow',action='store_true',help='Plot cutflow distributions')
    parser.add_argument('-pfr','--plotFakeRegions',action='store_true',help='Plot fake regions')
    parser.add_argument('-rt','--runTau',action='store_true',help='Run Tau finalStates (not implemented)')
    parser.add_argument('-ub','--unblind',action='store_false',help='Unblind signal channel')
    parser.add_argument('-ns','--nostack',action='store_true',help='Plot histograms unstacked')
    parser.add_argument('-no','--normalize',action='store_true',help='Plot histograms normalized to 1')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500)
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses')
    parser.add_argument('-ac','--allControls',action='store_true',help='Run over all controls for a given analysis (3l, 4l)')
    parser.add_argument('-fr','--doFakeRate',action='store_true',help='Make fake rate plots and output fake rate histograms')
    parser.add_argument('-dd','--doDetailed',action='store_true',help='Do detailed lepton plots')
    parser.add_argument('-fs','--finalStates',type=str,default='all',help='Only run given channels (ie: "eee,emm")')
    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to plots.')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for plots.')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    massLists = {
        13 : {
            'Hpp3l' : [500],
            'Hpp4l' : [500]
        }, 
        8 : {
            '3l' : [170, 200, 250, 300, 350, 400, 450, 500, 600, 700],
            '4l' : [110, 130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]
        }
    }

    controlList = {
        'Hpp3l' : ['Hpp3l', 'WZ'],
        'Hpp4l' : ['Hpp4l']
    }

    if args.period == 7:
        logger.warning("7 TeV not implemented")
    elif args.doFakeRate:
        plotFakeRate(args.analysis,args.channel,args.period,mass=args.mass,loglevel=args.log,blind=args.unblind,doDetailed=args.doDetailed)
    else:
        plotRegion(args.analysis,args.channel,args.period,plotFinalStates=args.plotFinalStates,runTau=args.runTau,blind=args.unblind,mass=args.mass,plotJetBins=args.plotJetBins,plotOverlay=args.plotOverlay,plotShapes=args.plotShapes,plotCutFlow=args.plotCutFlow,myCut=args.cut,finalStates=args.finalStates,nostack=args.nostack,normalize=args.normalize,scaleFactor=args.scaleFactor,loglevel=args.log,doDetailed=args.doDetailed,plotFakeRegions=args.plotFakeRegions)

    return 0


if __name__ == "__main__":
    main()
