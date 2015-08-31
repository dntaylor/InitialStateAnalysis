#!/usr/bin/env python

import os
import logging
import sys
import argparse
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.Optimizer import Optimizer
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES

def initializePlotter(analysis, period, plotName, nl, runTau):
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,period,analysis)
    saves = '%s_%s_%iTeV' % (analysis,analysis,period)
    sigMap = getSigMap(nl,500)
    intLumiMap = getIntLumiMap()
    regionBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl,runTau=runTau)
    mergeDict = getMergeDict(period)
    scaleFactor = 'event.pu_weight*event.lep_scale*event.trig_scale'
    masses = _3L_MASSES if nl==3 else _4L_MASSES
    plotter = Plotter(analysis,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName=plotName)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in regionBackground[analysis]])
    plotter.initializeSignalSamples([sigMap[period][x] for x in masses])
    plotter.setIntLumi(intLumiMap[period])
    return plotter

def initializeOptimizer(analysis, period, plotName, preselection, sigSelection, numTaus):
    nl = 3
    plotter = initializePlotter(analysis, period, plotName, nl, True)
    optimizer = Optimizer(plotter)
    optimizer.setSelection(preselection,signalSelection=sigSelection,numTaus=numTaus)
    return optimizer

def optimize(analysis, period):
    tauFlavor = {
        0: ['ee','em','mm'],
        1: ['et','mt'],
        2: ['tt']
    }

    genChannels = {
        'ee': ['eee','eem','eet'],
        'em': ['eme','emm','emt'],
        'et': ['ete','etm','ett'],
        'mm': ['mme','mmm','mmt'],
        'mt': ['mte','mtm','mtt'],
        'tt': ['tte','ttm','ttt'],
    }

    recoChannels = {
        'ee': ['eee','eem'],
        'em': ['eme','emm','mee','mem'],
        'mm': ['mme','mmm'],
        'et': ['eee','eme','eem','emm','mee','mem'],
        'mt': ['mee','mem','mme','mmm','eme','emm'],
        'tt': ['eee','eem','eme','emm','mee','mem','mme','mmm'],
    }

    for nTaus in range(3):
        genSelect = '(' + ' || '.join(['genChannel=="%s"'%x for y in tauFlavor[nTaus] for x in genChannels[y]]) + ')'
        recoSelect = '(' + ' || '.join(['channel=="%s"'%x for y in tauFlavor[nTaus] for x in recoChannels[y]]) + ')'

        optimizer = initializeOptimizer(analysis, period, 'plots_optimize_%iTau' %(nTaus), 'select.PassTight & %s' % recoSelect, genSelect, nTaus)
        
        optimizer.addCut('st',         'finalstate.sT >',            100., 1500., 10.)
        optimizer.addCut('met',        'finalstate.met >',             0.,  500.,  5.)
        optimizer.addCut('dR',         'h1.dR <',                      1.,    5.,  0.05)
        optimizer.addCut('zmass',      'fabs(z1.mass-%f) >' % ZMASS,   0.,  200.,  5.)
        optimizer.addCut('zmassUnder', '%f-z1.mass >' % ZMASS,         0.,  200.,  5.)
        optimizer.addCut('zmassOver',  'z1.mass-%f >' % ZMASS,         0.,  200.,  5.)
        optimizer.addCut('hmass',      'fabs(h1.mass-MASS) <',         0.,  500.,  5.)
        optimizer.addCut('hmassUnder', 'MASS-h1.mass >',               0.,  500.,  5.)
        optimizer.addCut('hmassOver',  'h1.mass-MASS >',               0.,  500.,  5.)
        #optimizer.addCut('st',         'finalstate.sT >',            100., 1500., 200.)
        #optimizer.addCut('dR',         'h1.dR <',                      1.,    5.,  1.)
        #optimizer.addCut('zmass',      'fabs(z1.mass-%f) >' % ZMASS,   0.,  100.,  20.)
        #optimizer.addCut('hmass',      'fabs(h1.mass-MASS) <',         0.,  400.,  100.)
        #optimizer.addCut('hmassUnder', 'h1.mass <',                    0., 1000., 200.)
        #optimizer.addCut('hmassOver',  'h1.mass >',                    0., 1000., 200.)

        print 'Optimizing'
        optimizer.optimize(masses=_3L_MASSES)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','WZ'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[8, 13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500,help='Mass for signal')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses for signal')
    parser.add_argument('-ub','--unblind',action='store_true',help='unblind')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    optimize(args.analysis, args.period)

if __name__ == "__main__":
    main()
