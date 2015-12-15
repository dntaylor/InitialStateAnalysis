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
    directory = kwargs.pop('directory','')
    blind = kwargs.pop('blind',True)
    mass = kwargs.pop('mass',500)
    runTau = kwargs.pop('runTau',False)
    myCut = kwargs.pop('myCut','1')
    plotFinalStates = kwargs.pop('plotFinalStates',False)
    plotJetBins = kwargs.pop('plotJetBins',False)
    plotEfficiency = kwargs.pop('plotEfficiency',False)
    plotOverlay = kwargs.pop('plotOverlay',False)
    plotShapes = kwargs.pop('plotShapes',False)
    plotCutFlow = kwargs.pop('plotCutFlow',False)
    plotCorrelation = kwargs.pop('plotCorrelation',False)
    plotSignal = kwargs.pop('plotSignal',False)
    plotFakeRegions = kwargs.pop('plotFakeRegions',False)
    finalStatesToPlot = kwargs.pop('finalStates','all')
    nostack = kwargs.pop('nostack',False)
    normalize = kwargs.pop('normalize',False)
    doDetailed = kwargs.pop('doDetailed',False)
    doDataDriven = kwargs.pop('doDataDriven',False)
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
    ntuples = getNtupleDirectory(analysis,channel,runPeriod)
    if directory: ntuples = directory
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    if directory: saves = '{0}/{1}'.format(saves,directory)
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



    ####################################
    ### Setup the plotter parameters ###
    ####################################

    ####################
    ### single plots ###
    ####################
    # each of these will create a signle plotter and then produce the listed plots
    # TODO

    ###############################
    ### multiplot distributions ###
    ###############################
    # each of these will create a single plotter and produce a lot of plots with plotDistributions
    basicPlotParams = {
        # common plotter parameters
        'type' : 'Plotter',               # plotter type: Plotter, CutFlowPlotter
        'channel' : channel,              # channel to pass to plotter (name of tree)
        'plotterArgs' : {                 # optional arguments to pass to plotter (kwargs)
            'ntupleDir': ntuples,         # ntuple directory
            'saveDir' : saves,            # directory to save plots
            'period' : runPeriod,         # run period: 8,13
            'mergeDict' : mergeDict,      # sample merge dictionary
            'scaleFactor' : scaleFactor,  # scale factor ro apply to mc
            'loglevel' : loglevel,        # logging module level
            'rootName' : 'plots',         # root file for working
            'baseSelection' : '',         # base selection to skim ntuple (save compute time) NB: you CANNOT skim for data driven plotting
        },
        # sample to add to ploter
        'backgroundSamples' : [sigMap[runPeriod][x] for x in channelBackground[channel] + ['Sig']] if useSignal else [sigMap[runPeriod][x] for x in channelBackground[channel]],
        'dataSamples' : [sigMap[runPeriod]['data']] if dataplot else [],
        # integrated luminosity to scale the MC
        'intLumi' : intLumiMap[runPeriod],
        # type of plot (data, ratio, etc)
        'plotMode' : 'plotMCDataRatio' if dataplot else 'plotMC',
        # common distribution parameters
        'cut' : myCut,                    # cut to apply to plot
        'numLeps' : nl,                   # number of leptons in final state
        'isControl' : isControl,          # archaic... DELETE!
        'distributionArgs' : {            # plotDistributions optional arguments (kwargs)
            'analysis': analysis,         # analysis (Hpp3l, WZ, ...)
            'region' : channel,           # region in analysis
            'savedir' : '',               # subdirectory to save this set of plots
            'nostack' : nostack,          # dont stack the plots... just dont use
            'normalize' : normalize,      # normalize all to 1... just dont use
            'mass' : mass,                # mass of h++, only effects that
            'doDetailed' : doDetailed,    # plot a whole bunch of additional plots
        },
    }

    def addFinalStates(setname,params):
        newParams = {}
        baseSel = params['plotterArgs']['baseSelection']
        savedir = params['distributionArgs']['savedir']
        for c in fsToPlot:
            name = '{0}_{1}'.format(setname,c)
            newParams[name] = params
            newParams[name].update({
                'plotterArgs' : {
                    'baseSelection' : '{0} && channel=="{1}"'.format(baseSel,c) if baseSel else 'channel=="{0}"'.format(c),
                    'rootName' : 'plots_{0}_{1}'.format(setname,c),
                },
                'cut' : '{0} && channel=="{1}"'.format(params['cut'],c),
                'distributionArgs' : {
                    'savedir' : '{0}/{1}'.format(savedir,c) if savedir else c,
                },
            })
        if analysis in combinedFinalStates:
            for c in combinedFinalStatesOrdered[analysis]:
                name = '{0}_combo'.format(c)
                sel = combinedFinalStates[analysis][c]['cut']
                newParams[name] = params
                newParams[name].update({
                    'plotterArgs' : {
                        'baseSelection' : '{0} && {1}'.format(baseSel,sel) if baseSel else sel,
                        'rootName' : 'plots_{0}_{1}'.format(setname,name),
                    },
                    'cut' : '{0} && {1}'.format(params['cut'],sel),
                    'distributionArgs' : {
                        'savedir' : '{0}/{1}'.format(savedir,name) if savedir else name,
                    },
                })
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
                sel = customFinalStates[analysis][c]
                name = '{0}_{1}'.format(setname,c)
                newParams[name] = params
                newParams[name].update({
                    'plotterArgs' : {
                        'baseSelection' : '{0} && {1}'.format(baseSel,sel) if baseSel else sel,
                        'rootName' : 'plots_{0}_{1}'.format(setname,c),
                    },
                    'cut' : '{0} && {1}'.format(params['cut'],sel),
                    'distributionArgs' : {
                        'savedir' : '{0}/{1}'.format(savedir,c) if savedir else c,
                    },
                })
        return newParams


    # parameters for plotDistributions
    plotParams = {}
    plotParams['distributions'] = basicPlotParams
    plotParams['distributions'].update({
        'plotterArgs' : {
            'baseSelection' : myCut,
            'rootName' : 'plots_distributions',
        },
       'cut' : myCut,
    })
    if plotFinalStates:
        plotParams.update(addFinalStates('distributions',plotParams['distributions']))
    if plotFakeRegions:
        labels = {0:'F',1:'P'}
        allCuts = getPassTightDefinition(analysis,channel,runPeriod)
        if analysis in ['Hpp3l']: allCuts = '1'
        fakeCut = myCut.replace('select.passTight',allCuts)
        for f in itertools.product('PF',repeat=nl):
            thisName = ''.join(f)
            name = 'fakeregion_{0}'.format(thisName)
            thisCut = '{0} && fakeChannel=="{1}"'.format(fakeCut,thisName)
            plotParams[name] = basicPlotParams
            plotParams[name].update({
                'plotterArgs' : {
                    'baseSelection' : thisCut,
                    'rootName' : 'plots_fakeregions_{0}'.format(thisName),
                },
                'cut' : thisCut,
                'distributionArgs' : {
                    'savedir' : thisName,
                },
            })
            if plotFinalStates:
                plotParams.update(addFinalStates(name,plotParams[name]))
    if doDataDriven:
        plotParams['datadriven'] = basicPlotParams
        plotParams['datadriven'].update({
            'plotterArgs' : {
                'rootName' : 'plots_datadriven',
            },
            'cut' : myCut,
            'distributionArgs' : {
                'savedir' : 'datadriven',
            },
        })
        if plotFinalStates:
            plotParams.update(addFinalStates('datadriven',plotParams['datadriven']))
    if plotOverlay:
        plotParams['overlay'] = basicPlotParams
        plotParams['overlay'].update({
            'plotterArgs' : {
                'baseSelection' : myCut,
                'rootName' : 'plots_overlay',
            },
            'backgroundSamples' : [sigMap[runPeriod][x] for x in channelBackground[channel]],
            'signalSamples' : [sigMap[runPeriod][x] for x in ['Sig']],
            'plotMode' : 'plotMCDataSignalRatio' if dataplot else 'plotMCSignalRatio',
            'cut' : myCut,
            'distributionArgs' : {
                'savedir' : 'overlay',
                'signalscale' : 100,
            },

        })
        if plotFinalStates:
            plotParams.update(addFinalStates('overlay',plotParams['overlay']))
    if plotCutFlow:
        baseCut = 'z1.PassTight1==1 && z1.PassTight2==1 && wz.PassTight1==1'
        if analysis in ['Hpp3l']: baseCut = 'h1.PassTight1==1 && h1.PassTight2==1 && h2.PassTight1==1'
        for i in range(len(cutFlowMap[channel]['cuts'])):
            # regular cutflow
            thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
            fullCut = '{0} && {1}'.format(baseCut,thisCut)
            thisName = cutFlowMap[channel]['labels_simple'][i]
            name = 'cutflow_{0}'.format(thisName)
            plotParams[name] = basicPlotParams
            plotParams[name].update({
                'plotterArgs' : {
                    'baseSelection' : fullCut,
                    'rootName' : 'plots_{0}'.format(name),
                },
                'cut' : fullCut,
                'distributionArgs' : {
                    'savedir' : 'cutflow/{0}'.format(thisName),
                },
            })
            # n-1 plots
            nMinusOneCut = ' && '.join([cutFlowMap[channel]['cuts'][x] for x in range(len(cutFlowMap[channel]['cuts'])) if x!=i])
            fullNMinusOneCut = '{0} && {1}'.format(baseCut,nMinusOneCut)
            nMinusOneName = 'nMinusOne_{0}'.format(thisName)
            plotParams[nMinusOneName] = basicPlotParams
            plotParams[nMinusOneName.update({
                'plotterArgs' : {
                    'baseSelection' : fullNMinusOneCut,
                    'rootName' : 'plots_{0}'.format(nMinusOneName),
                },
                'cut' : fullNMinusOneCut,
                'distributionArgs' : {
                    'savedir' : 'nMinusOne/{0}'.format(thisName),
                },
            })
            if plotFinalStates:
                plotParams.update(addFinalStates(name,plotParams[name]))
                plotParams.update(addFinalStates(nMinusOneName,plotParams[nMinusOneName]))
            if doDataDriven:
                # cutflow
                ddname = 'datadriven_{0}'.format(name)
                plotParams[ddname] = basicPlotParams
                plotParams[ddname].update({
                    'plotterArgs' : {
                        'rootName' : 'plots_{0}'.format(ddname),
                    },
                    'cut' : fullCut,
                    'distributionArgs' : {
                        'savedir' : 'datadriven/cutflow/{0}'.format(thisName),
                    },
                })
                # n-1
                nMinusOneDdname = 'datadriven_{0}'.format(nMinusOneName)
                plotParams[nMinusOneDdname] = basicPlotParams
                plotParams[nMinusOneDdname].update({
                    'plotterArgs' : {
                         'rootName' : 'plots_{0}'.format(nMinusOneDdname),
                    },
                    'cut' : fullNMinusOneCut,
                    'distributionArgs' : {
                        'savedir' : 'datadriven/nMinusOne/{0}'.format(thisName),
                    },
                })
                if plotFinalStates:
                    plotParams.update(addFinalStates(ddname,plotParams[ddname]))
                    plotParams.update(addFinalStates(nMinusOneDdname,plotParams[nMinusOneDdname]))






def parse_command_line(argv):
    parser = get_parser("Plot a given channel and period")

    parser.add_argument('-pf','--plotFinalStates',action='store_true',help='Plot individual final states')
    parser.add_argument('-pj','--plotJetBins',action='store_true',help='Plot jet bins')
    parser.add_argument('-pe','--plotEfficiency',action='store_true',help='Plot efficiency')
    parser.add_argument('-po','--plotOverlay',action='store_true',help='Plot overlay')
    parser.add_argument('-ps','--plotShapes',action='store_true',help='Plot shapes')
    parser.add_argument('-pcf','--plotCutFlow',action='store_true',help='Plot cutflow distributions')
    parser.add_argument('-pfr','--plotFakeRegions',action='store_true',help='Plot fake regions')
    parser.add_argument('-pd','--plotDataDriven',action='store_true',help='Plot data driven')
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
        plotRegion(args.analysis,args.channel,args.period,plotFinalStates=args.plotFinalStates,runTau=args.runTau,blind=args.unblind,mass=args.mass,plotJetBins=args.plotJetBins,plotOverlay=args.plotOverlay,plotShapes=args.plotShapes,plotCutFlow=args.plotCutFlow,myCut=args.cut,finalStates=args.finalStates,nostack=args.nostack,normalize=args.normalize,scaleFactor=args.scaleFactor,loglevel=args.log,doDetailed=args.doDetailed,plotFakeRegions=args.plotFakeRegions,doDataDriven=args.plotDataDriven,plotEfficiency=args.plotEfficiency,directory=args.directory)

    return 0


if __name__ == "__main__":
    main()
