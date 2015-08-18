#!/usr/bin/env python

import os
import logging
import sys
import argparse
from InitialStateAnalyis.Plotters.Plotter import Plotter
from InitialStateAnalyis.Plotters.Optimizer import Optimizer
from InitialStateAnalyis.Plotters.plotUtils import *
from InitialStateAnalyis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES

def initializePlotter(analysis, period, mass, plotName, nl, runTau):
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,period,analysis)
    saves = '%s_%s_%iTeV' % (analysis,analysis,period)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    regionBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl,runTau=runTau)
    mergeDict = getMergeDict(period)
    scaleFactor = 'event.pu_weight*event.lep_scale*event.trig_scale'
    plotter = Plotter(analysis,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName=plotName)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in regionBackground[analysis]])
    plotter.initializeSignalSamples([sigMap[period]['Sig']])
    plotter.setIntLumi(intLumiMap[period])
    return plotter

def initializeOptimizer(analysis, period, mass, plotName, preselection, sigSelection):
    nl = 3
    plotter = initializePlotter(analysis, period, mass, plotName, nl, True)
    optimizer = Optimizer(plotter)
    optimizer.setSelection(preselection,signalSelection=sigSelection)
    return optimizer

def optimize(analysis, period, mass):
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

    nTaus = 0

    genSelect = '(' + ' | '.join(['genChannel=="%s"'%x for y in tauFlavor[nTaus] for x in genChannels[y]]) + ')'
    recoSelect = '(' + ' | '.join(['channel=="%s"'%x for y in tauFlavor[nTaus] for x in recoChannels[y]]) + ')'

    optimizer = initializeOptimizer(analysis, period, mass, 'optimize.root', 'select.PassTight & %s' % recoSelect, genSelect)
    
    #optimizer.addCut('st',    'finalstate.sT >',            100., 1500., 10.)
    #optimizer.addCut('dR',    'h1.dR <',                    1.,   5.,    0.1)
    #optimizer.addCut('zmass', 'fabs(z1.mass-%f) >' % ZMASS, 0.,   100.,  5.)
    #optimizer.addCut('hmass', 'fabs(h1.mass-%f) <' % mass,  0.,   200.,  5.)
    optimizer.addCut('st',    'finalstate.sT >',            100., 1500., 200.)
    optimizer.addCut('dR',    'h1.dR <',                    1.,   5.,    1.)
    optimizer.addCut('zmass', 'fabs(z1.mass-%f) >' % ZMASS, 0.,   100.,  20.)
    optimizer.addCut('hmass', 'fabs(h1.mass-%f) <' % mass,  0.,   200.,  40.)

    print 'Optimizing %i' % mass
    optimizer.optimize()

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

    optimize(args.analysis, args.period, args.mass)

if __name__ == "__main__":
    main()
