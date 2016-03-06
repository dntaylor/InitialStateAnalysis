#!/usr/bin/env python

import argparse
import itertools
import sys
import pickle
import json
import logging
from copy import deepcopy
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
    directory = kwargs.pop('directory','')
    outputDir = kwargs.pop('outputDir','')
    skipCutflow = kwargs.pop('skipCutflow',False)
    blind = kwargs.pop('blind',True)
    mass = kwargs.pop('mass',500)
    runTau = kwargs.pop('runTau',False)
    myCut = kwargs.pop('myCut','1')
    skipDistributions = kwargs.pop('skipDistributions',False)
    plotFinalStates = kwargs.pop('plotFinalStates',False)
    plotJetBins = kwargs.pop('plotJetBins',False)
    plotOverlay = kwargs.pop('plotOverlay',False)
    plotCutFlow = kwargs.pop('plotCutFlow',False)
    plotNMinusOne = kwargs.pop('plotNMinusOne',False)
    plotSignal = kwargs.pop('plotSignal',False)
    plotFakeRegions = kwargs.pop('plotFakeRegions',False)
    finalStatesToPlot = kwargs.pop('finalStates','all')
    doDetailed = kwargs.pop('doDetailed',False)
    doDataDriven = kwargs.pop('doDataDriven',False)
    scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    useSignal = analysis in ['Hpp3l','Hpp4l']
    allowedPlots = kwargs.pop('allowedPlots',['all'])
    tightW = kwargs.pop('tightW',False)
    loglevel = kwargs.pop('loglevel','INFO')
    for key, value in kwargs.iteritems():
        logger.warning("Unrecognized parameter '" + key + "' = " + str(value))
        return 0

    if useSignal: logger.info("%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass))
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
    ntuples = getNtupleDirectory(analysis,channel,runPeriod)
    if directory: ntuples = directory
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    if directory: saves = '{0}/{1}'.format(saves,directory)
    if outputDir: saves = '{0}/{1}'.format(saves,outputDir)
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
    dataplot = not blind
    mergeDict = getMergeDict(runPeriod)
    cutFlowMap = {}
    cutFlowMap[channel] = defineCutFlowMap(channel,finalStates,mass)



    ####################################
    ### Setup the plotter parameters ###
    ####################################

    #########################
    ### Cut flow plotters ###
    #########################
    cutflowPlotParams = {
        # common plotter parameters
        'type' : 'CutFlowPlotter',         # plotter type: Plotter, CutFlowPlotter
        'channel' : channel,               # channel to pass to plotter (name of tree)
        'plotterArgs' : {                  # optional arguments to pass to plotter (kwargs)
            'ntupleDir'     : ntuples,     # ntuple directory
            'saveDir'       : saves,       # directory to save plots
            'period'        : runPeriod,   # run period: 8,13
            'mergeDict'     : mergeDict,   # sample merge dictionary
            'scaleFactor'   : scaleFactor, # scale factor ro apply to mc
            'loglevel'      : loglevel,    # logging module level
            'rootName'      : 'cutflow',   # root file for working
            'tightW'        : tightW,      # tight W for WZ   
            'datadriven'    : doDataDriven,# run data driven
            'baseSelection' : '',          # base selection to skim ntuple (save compute time) NB: you CANNOT skim for data driven plotting
        },
        # sample to add to ploter
        'backgroundSamples' : [sigMap[runPeriod][x] for x in channelBackground[channel] + ['Sig']] if useSignal else [sigMap[runPeriod][x] for x in channelBackground[channel]],
        'signalSamples' : [],
        'dataSamples' : [sigMap[runPeriod]['data']],
        # integrated luminosity to scale the MC
        'intLumi' : intLumiMap[runPeriod],
        # type of plot (data, ratio, etc)
        'plotMode' : 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC',
        # parameters for the 
        'finalStates' : finalStates,
        'savedir'     : '',
        'allPlots' : {},
    }
    if doDataDriven: cutflowPlotParams['backgroundSamples'] = [sigMap[runPeriod][x] for x in channelBackground[channel+'datadriven'] + ['Sig']] if useSignal else [sigMap[runPeriod][x] for x in channelBackground[channel+'datadriven']]

    cutflowParams = {}


    if doDataDriven:
        cutflowParams['datadriven'] = deepcopy(cutflowPlotParams)
        cutflowParams['datadriven'].update({
            'cut' : myCut,
            'savedir' : 'datadriven',
        })
        cutflowParams['datadriven']['plotterArgs'].update({
            'rootName' : 'cutflow_datadriven',
        })
    else:
        cutflowParams['cutflow'] = deepcopy(cutflowPlotParams)
        cutflowParams['cutflow'].update({
            'cut' : myCut,
        })

    # now add the plots
    for param in cutflowParams:
        channel      = cutflowParams[param]['channel']
        finalStates  = cutflowParams[param]['finalStates']
        savedir      = cutflowParams[param]['savedir']
        cut          = cutflowParams[param]['cut']
        allCutflowPlots = getCutflowParams(channel,finalStates,cut,savedir)
        cutflowParams[param]['allPlots'] = allCutflowPlots

    # generate the plotters (serial for now)
    plotTypes = {
        'Plotter'        : Plotter,
        'CutFlowPlotter' : CutFlowPlotter,
    }
    for param in cutflowParams:
        if skipCutflow: break
        plotType          = cutflowParams[param]['type']
        plotMode          = cutflowParams[param]['plotMode']
        channel           = cutflowParams[param]['channel']
        plotterArgs       = cutflowParams[param]['plotterArgs']
        backgroundSamples = cutflowParams[param]['backgroundSamples']
        signalSamples     = cutflowParams[param]['signalSamples']
        dataSamples       = cutflowParams[param]['dataSamples']
        intLumi           = cutflowParams[param]['intLumi']

        logger.info('%s:%s:%iTeV: Generating Cutflow Plots for %s' % (analysis, channel, runPeriod, param))
        plotter = plotTypes[plotType](channel,**plotterArgs)
        plotter.initializeBackgroundSamples(backgroundSamples)
        if signalSamples: plotter.initializeSignalSamples(signalSamples)
        if dataSamples: plotter.initializeDataSamples(dataSamples)
        plotter.setIntLumi(intLumi)
        plotMethod = getattr(plotter,plotMode)

        allPlots = cutflowParams[param]['allPlots']
        for plot in allPlots:
            logger.info('%s:%s:%iTeV: %s - %s' % (analysis, channel, runPeriod, param, plot))
            args   = cutflowParams[param]['allPlots'][plot]['args']
            kwargs = cutflowParams[param]['allPlots'][plot]['kwargs']
            plotMethod(*args,**kwargs)



    ###############################
    ### multiplot distributions ###
    ###############################
    # each of these will create a single plotter and produce a lot of plots with plotDistributions
    basicPlotParams = {
        # common plotter parameters
        'type' : 'Plotter',                # plotter type: Plotter, CutFlowPlotter
        'channel' : channel,               # channel to pass to plotter (name of tree)
        'plotterArgs' : {                  # optional arguments to pass to plotter (kwargs)
            'ntupleDir'     : ntuples,     # ntuple directory
            'saveDir'       : saves,       # directory to save plots
            'period'        : runPeriod,   # run period: 8,13
            'mergeDict'     : mergeDict,   # sample merge dictionary
            'scaleFactor'   : scaleFactor, # scale factor ro apply to mc
            'loglevel'      : loglevel,    # logging module level
            'rootName'      : 'plots',     # root file for working
            'tightW'        : tightW,      # tight W for WZ   
            'datadriven'    : doDataDriven,# run data driven
            'baseSelection' : '',          # base selection to skim ntuple (save compute time) NB: you CANNOT skim for data driven plotting
        },
        # sample to add to ploter
        'backgroundSamples' : [sigMap[runPeriod][x] for x in channelBackground[channel] + ['Sig']] if useSignal else [sigMap[runPeriod][x] for x in channelBackground[channel]],
        'signalSamples' : [],
        'dataSamples' : [sigMap[runPeriod]['data']],
        # integrated luminosity to scale the MC
        'intLumi' : intLumiMap[runPeriod],
        # type of plot (data, ratio, etc)
        'plotMode' : 'plotMCDataRatio' if dataplot else 'plotMC',
        # common distribution parameters
        'cut' : myCut,                    # cut to apply to plot
        'numLeps' : nl,                   # number of leptons in final state
        'distributionArgs' : {            # plotDistributions optional arguments (kwargs)
            'analysis': analysis,         # analysis (Hpp3l, WZ, ...)
            'region' : channel,           # region in analysis
            'savedir' : '',               # subdirectory to save this set of plots
            'mass' : mass,                # mass of h++, only effects that
            'doDetailed' : doDetailed,    # plot a whole bunch of additional plots
        },
        'allPlots' : {},                  # all plots to be plotted
    }
    if doDataDriven: basicPlotParams['backgroundSamples'] = [sigMap[runPeriod][x] for x in channelBackground[channel+'datadriven'] + ['Sig']] if useSignal else [sigMap[runPeriod][x] for x in channelBackground[channel+'datadriven']]

    def addFinalStates(setname,params):
        newParams = {}
        baseSel = deepcopy(params['plotterArgs']['baseSelection'])
        savedir = deepcopy(params['distributionArgs']['savedir'])
        for c in fsToPlot:
            name = '{0}_{1}'.format(setname,c)
            newParams[name] = deepcopy(params)
            newParams[name].update({
                'cut' : '{0} && channel=="{1}"'.format(params['cut'],c),
            })
            newParams[name]['plotterArgs'].update({
                'baseSelection' : '{0} && channel=="{1}"'.format(baseSel,c) if baseSel else 'channel=="{0}"'.format(c),
                'rootName' : 'plots_{0}_{1}'.format(setname,c),
            })
            newParams[name]['distributionArgs'].update({
                'savedir' : '{0}/{1}'.format(savedir,c) if savedir else c,
            })
        #if analysis in combinedFinalStates:
        #    for c in combinedFinalStatesOrdered[analysis]:
        #        name = '{0}_combo'.format(c)
        #        sel = combinedFinalStates[analysis][c]['cut']
        #        newParams[name] = deepcopy(params)
        #        newParams[name].update({
        #            'cut' : '{0} && {1}'.format(params['cut'],sel),
        #        })
        #        newParams[name]['plotterArgs'].update({
        #            'baseSelection' : '{0} && {1}'.format(baseSel,sel) if baseSel else sel,
        #            'rootName' : 'plots_{0}_{1}'.format(setname,name),
        #        })
        #        newParams[name]['distributionArgs'].update({
        #            'savedir' : '{0}/{1}'.format(savedir,name) if savedir else name,
        #        })
        #if analysis in customFinalStates:
        #    for c in customFinalStates[analysis]:
        #        sel = customFinalStates[analysis][c]
        #        name = '{0}_{1}'.format(setname,c)
        #        newParams[name] = deepcopy(params)
        #        newParams[name].update({
        #            'cut' : '{0} && {1}'.format(params['cut'],sel),
        #        })
        #        newParams[name]['plotterArgs'].update({
        #            'baseSelection' : '{0} && {1}'.format(baseSel,sel) if baseSel else sel,
        #            'rootName' : 'plots_{0}_{1}'.format(setname,c),
        #        })
        #        newParams[name]['distributionArgs'].update({
        #            'savedir' : '{0}/{1}'.format(savedir,c) if savedir else c,
        #        })
        return newParams


    # parameters for plotDistributions
    plotParams = {}
    if doDataDriven:
        if not skipDistributions:
            plotParams['datadrivenDistributions'] = deepcopy(basicPlotParams)
            plotParams['datadrivenDistributions'].update({
                'cut' : myCut,
            })
            plotParams['datadrivenDistributions']['plotterArgs'].update({
                'rootName' : 'plots_datadriven',
            })
            plotParams['datadrivenDistributions']['distributionArgs'].update({
                'savedir' : 'datadriven',
            })
            if plotFinalStates:
                plotParams.update(addFinalStates('datadrivenDistributions',plotParams['datadrivenDistributions']))
    else:
        if not skipDistributions:
            plotParams['mcDistributions'] = deepcopy(basicPlotParams)
            plotParams['mcDistributions'].update({
               'cut' : myCut,
            })
            plotParams['mcDistributions']['plotterArgs'].update({
                'baseSelection' : myCut,
                'rootName' : 'plots_distributions',
            })
            if plotFinalStates:
                plotParams.update(addFinalStates('mcDistributions',plotParams['mcDistributions']))
        if plotFakeRegions:
            tightCuts = {
                0: 'z1.PassTight1',
                1: 'z1.PassTight2',
                2: 'w1.PassTight1',
            }
            lepScale = {
                0: 'z1.LepScaleTight1',
                1: 'z1.LepScaleTight2',
                2: 'w1.LepScaleTight1',
            }
            lepScaleLoose = {
                0: 'z1.LepScaleLoose1',
                1: 'z1.LepScaleLoose2',
                2: 'w1.LepScaleLoose1',
            }
            if tightW:
                tightCuts = {
                    0: 'z1.PassMedium1',
                    1: 'z1.PassMedium2',
                    2: 'w1.PassTight1',
                }
                lepScale = {
                    0: 'z1.LepScaleMedium1',
                    1: 'z1.LepScaleMedium2',
                    2: 'w1.LepScaleTight1',
                }
            scaleMap = {
                'P' : lepScale,
                'F' : lepScaleLoose,
            }
            labels = {0:'F',1:'P'}
            allCuts = getPassTightDefinition(analysis,channel,runPeriod)
            if analysis in ['Hpp3l']: allCuts = '1'
            fakeCut = myCut
            fakeCut = fakeCut.replace('select.passTight',allCuts)
            for i in range(nl):
                fakeCut = fakeCut.replace(tightCuts[i],'1')
            for f in itertools.product('PF',repeat=nl):
                thisName = ''.join(f)
                name = 'fakeregion_{0}'.format(thisName)
                fakechan = 'fakeChannel'
                if tightW: fakechan = 'fakeChannel_tightW'
                thisCut = '{0} && {1}=="{2}"'.format(fakeCut,fakechan,thisName)
                baseScale = 'event.gen_weight*event.pu_weight*event.trig_scale'
                scale = '{0}*{1}'.format(baseScale,'*'.join([scaleMap[f[i]][i] for i in range(nl)]))
                plotParams[name] = deepcopy(basicPlotParams)
                plotParams[name].update({
                    'cut' : thisCut,
                })
                plotParams[name]['plotterArgs'].update({
                    'baseSelection' : thisCut,
                    'scaleFactor' : scale,
                    'rootName' : 'plots_fakeregions_{0}'.format(thisName),
                })
                plotParams[name]['distributionArgs'].update({
                        'savedir' : thisName,
                })
                if plotFinalStates:
                    plotParams.update(addFinalStates(name,plotParams[name]))
    if plotOverlay:
        plotParams['overlay'] = deepcopy(basicPlotParams)
        plotParams['overlay'].update({
            'backgroundSamples' : [sigMap[runPeriod][x] for x in channelBackground[channel]],
            'signalSamples' : [sigMap[runPeriod][x] for x in ['Sig']],
            'plotMode' : 'plotMCDataSignalRatio' if dataplot else 'plotMCSignalRatio',
            'cut' : myCut,
        })
        plotParams['overlay']['plotterArgs'].update({
            'baseSelection' : '' if doDataDriven else myCut,
            'rootName' : 'plots_overlay',
        })
        plotParams['overlay']['distributionArgs'].update({
            'savedir' : 'datadriven/overlay' if doDataDriven else 'overlay',
            'signalscale' : 100,
        })
        if plotFinalStates:
            plotParams.update(addFinalStates('overlay',plotParams['overlay']))
    if plotCutFlow or plotNMinusOne:
        baseCut = 'z1.PassTight1==1 && z1.PassTight2==1 && w1.PassTight1==1'
        if tightW: baseCut = 'z1.PassMedium1==1 && z1.PassMedium2==1 && w1.PassTight1==1'
        if analysis in ['Hpp3l']: baseCut = 'h1.PassTight1==1 && h1.PassTight2==1 && h2.PassTight1==1'
        for i in range(len(cutFlowMap[channel]['cuts'])):
            # regular cutflow
            thisName = cutFlowMap[channel]['labels_simple'][i]
            if plotCutFlow:
                thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
                fullCut = '{0} && {1}'.format(baseCut,thisCut)
                name = 'cutflow_{0}'.format(thisName)
                if doDataDriven: name += '_datadriven'
                plotParams[name] = deepcopy(basicPlotParams)
                plotParams[name].update({
                    'cut' : fullCut,
                })
                plotParams[name]['plotterArgs'].update({
                    'baseSelection' : '' if doDataDriven else fullCut,
                    'rootName' : 'plots_{0}'.format(name),
                })
                plotParams[name]['distributionArgs'].update({
                    'savedir' : 'datadriven/cutflow/{0}'.format(thisName) if doDataDriven else 'cutflow/{0}'.format(thisName),
                })
                if plotFinalStates:
                    plotParams.update(addFinalStates(name,plotParams[name]))
            # n-1 plots
            if plotNMinusOne:
                nMinusOneCut = ' && '.join([cutFlowMap[channel]['cuts'][x] for x in range(len(cutFlowMap[channel]['cuts'])) if x!=i])
                fullNMinusOneCut = '{0} && {1}'.format(baseCut,nMinusOneCut)
                nMinusOneName = 'nMinusOne_{0}'.format(thisName)
                if doDataDriven: nMinusOneName += '_datadriven'
                plotParams[nMinusOneName] = deepcopy(basicPlotParams)
                plotParams[nMinusOneName].update({
                    'cut' : fullNMinusOneCut,
                })
                plotParams[nMinusOneName]['plotterArgs'].update({
                    'baseSelection' : '' if doDataDriven else fullNMinusOneCut,
                    'rootName' : 'plots_{0}'.format(nMinusOneName),
                })
                plotParams[nMinusOneName]['distributionArgs'].update({
                    'savedir' : 'datadriven/nMinusOne/{0}'.format(thisName) if doDataDriven else 'nMinusOne/{0}'.format(thisName),
                })
                if plotFinalStates:
                    plotParams.update(addFinalStates(nMinusOneName,plotParams[nMinusOneName]))


    # now add the plots
    for param in plotParams:
        plotMode = plotParams[param]['plotMode']
        thisCut  = plotParams[param]['cut']
        numLeps  = plotParams[param]['numLeps']
        distributionArgs = plotParams[param]['distributionArgs']
        indivPlotParams = getPlotParams(plotMode,thisCut,numLeps,allowedPlots,**distributionArgs)
        plotParams[param]['allPlots'] = indivPlotParams

    # generate the plotters (serial for now)
    plotTypes = {
        'Plotter'        : Plotter,
        'CutFlowPlotter' : CutFlowPlotter,
    }
    for param in plotParams:
        plotType          = plotParams[param]['type']
        channel           = plotParams[param]['channel']
        plotterArgs       = plotParams[param]['plotterArgs']
        backgroundSamples = plotParams[param]['backgroundSamples']
        signalSamples     = plotParams[param]['signalSamples']
        dataSamples       = plotParams[param]['dataSamples']
        intLumi           = plotParams[param]['intLumi']

        logger.info('%s:%s:%iTeV: Generating Plots for %s' % (analysis, channel, runPeriod, param))
        plotter = plotTypes[plotType](channel,**plotterArgs)
        plotter.initializeBackgroundSamples(backgroundSamples)
        if signalSamples: plotter.initializeSignalSamples(signalSamples)
        if dataSamples: plotter.initializeDataSamples(dataSamples)
        plotter.setIntLumi(intLumi)

        allPlots = plotParams[param]['allPlots']
        for plot in allPlots:
            logger.info('%s:%s:%iTeV: %s - %s' % (analysis, channel, runPeriod, param, plot))
            plotMode = plotParams[param]['allPlots'][plot]['plotMode']
            variable = plotParams[param]['allPlots'][plot]['variable']
            binning  = plotParams[param]['allPlots'][plot]['binning']
            savename = plotParams[param]['allPlots'][plot]['savename']
            plotArgs = plotParams[param]['allPlots'][plot]['plotArgs']
            plotMethod = getattr(plotter,plotMode)
            plotMethod(variable,binning,savename,**plotArgs)

    logger.info('%s:%s:%iTeV: Finished' % (analysis, channel, runPeriod))

def parse_command_line(argv):
    parser = get_parser("Plot a given channel and period")

    parser.add_argument('-sd','--skipDistributions',action='store_true',help='Skip the basic distributions')
    parser.add_argument('-pf','--plotFinalStates',action='store_true',help='Plot individual final states')
    parser.add_argument('-pj','--plotJetBins',action='store_true',help='Plot jet bins')
    parser.add_argument('-pfr','--plotFakeRegions',action='store_true',help='Plot fake regions')
    parser.add_argument('-po','--plotOverlay',action='store_true',help='Plot overlay')
    parser.add_argument('-pcf','--plotCutFlow',action='store_true',help='Plot cutflow distributions')
    parser.add_argument('-scf','--skipCutFlowPlotter',action='store_true',help='Skip the cutflow plotter plots')
    parser.add_argument('-pnmo','--plotNMinusOne',action='store_true',help='Plot n minus one distributions')
    parser.add_argument('-dd','--datadriven',action='store_true',help='Plot data driven')
    parser.add_argument('--tightW',action='store_true',help='TightW for WZ')
    parser.add_argument('-rt','--runTau',action='store_true',help='Run Tau finalStates (not implemented)')
    parser.add_argument('-ub','--unblind',action='store_false',help='Unblind signal channel')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500)
    parser.add_argument('--detailed',action='store_true',help='Do detailed lepton plots')
    parser.add_argument('-fs','--finalStates',type=str,default='all',help='Only run given channels (ie: "eee,emm")')
    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to plots.')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for plots.')
    parser.add_argument('-p','--plots',nargs='+',type=str,default=['all'])
    parser.add_argument('-d','--directory',type=str,default='',help='Custom subdirectory (to keep more than one ntuple at a time)')
    parser.add_argument('-o','--outputDirectory',type=str,default='',help='Custom ouput directory (to keep more than one set of plots at a time)')
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
    else:
        plotRegion(args.analysis,
                   args.channel,
                   args.period,
                   skipDistributions=args.skipDistributions,
                   plotFinalStates=args.plotFinalStates,
                   runTau=args.runTau,
                   blind=args.unblind,
                   mass=args.mass,
                   plotJetBins=args.plotJetBins,
                   plotFakeRegions=args.plotFakeRegions,
                   plotOverlay=args.plotOverlay,
                   plotCutFlow=args.plotCutFlow,
                   plotNMinusOne=args.plotNMinusOne,
                   myCut=args.cut,
                   finalStates=args.finalStates,
                   scaleFactor=args.scaleFactor,
                   loglevel=args.log,
                   doDetailed=args.detailed,
                   doDataDriven=args.datadriven,
                   tightW=args.tightW,
                   allowedPlots=args.plots,
                   directory=args.directory,
                   skipCutflow=args.skipCutFlowPlotter,
                   outputDir=args.outputDirectory)

    return 0


if __name__ == "__main__":
    main()
