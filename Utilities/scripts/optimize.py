#!/usr/bin/env python

import os
import logging
import sys
import argparse
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.Optimizer import Optimizer
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES
from InitialStateAnalysis.Utilities.utilities import *

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
    nl = 3 if analysis in ['Hpp3l'] else 4
    plotter = initializePlotter(analysis, period, plotName, nl, True)
    optimizer = Optimizer(plotter,analysis,period)
    optimizer.setSelection(preselection,signalSelection=sigSelection,numTaus=numTaus)
    return optimizer

def optimize(analysis, period):
    #tauFlavor = {
    #    0: ['ee','em','mm'],
    #    1: ['et','mt'],
    #    2: ['tt']
    #}

    #genChannels = {
    #    'ee': ['eee','eem','eet'],
    #    'em': ['eme','emm','emt'],
    #    'et': ['ete','etm','ett'],
    #    'mm': ['mme','mmm','mmt'],
    #    'mt': ['mte','mtm','mtt'],
    #    'tt': ['tte','ttm','ttt'],
    #}

    #recoChannels = {
    #    'ee': ['eee','eem'],
    #    'em': ['eme','emm','mee','mem'],
    #    'mm': ['mme','mmm'],
    #    'et': ['eee','eme','eem','emm','mee','mem'],
    #    'mt': ['mee','mem','mme','mmm','eme','emm'],
    #    'tt': ['eee','eem','eme','emm','mee','mem','mme','mmm'],
    #}

    genChannels = {
        'Hpp3l': {
            0: ['eee','eem','eet','eme','emm','emt','mme','mmm','mmt'],
            1: ['ete','etm','ett','mte','mtm','mtt'],
            2: ['tte','ttm','ttt'],
        },
        'Hpp4l': {
            0: ['eeee','eeem','eemm','emee','emem','emmm','mmee','mmem','mmmm'],
            1: ['etet','etmt','mtet','mtmt'],
            2: ['tttt'],
        },
    }

    recoChannels = {
        'Hpp3l': {
            0: ['eee','eem','eme','emm','mee','mem','mme','mmm'],
            1: ['eee','eem','eme','emm','mee','mem','mme','mmm'],
            2: ['eee','eem','eme','emm','mee','mem','mme','mmm'],
        },
        'Hpp4l': {
            0: ['eeee','eeem','eeme','eemm','emee','emem','emme','emmm','meee','meem','meme','memm','mmee','mmem','mmme','mmmm'],
            1: ['eeee','eeem','eeme','eemm','emee','emem','emme','emmm','meee','meem','meme','memm','mmee','mmem','mmme','mmmm'],
            2: ['eeee','eeem','eeme','eemm','emee','emem','emme','emmm','meee','meem','meme','memm','mmee','mmem','mmme','mmmm'],
        },
    }

    for nTaus in range(3):
        genSelect = '(' + ' || '.join(['genChannel=="%s"'%x for  x in genChannels[analysis][nTaus]]) + ')'
        recoSelect = '(' + ' || '.join(['channel=="%s"'%x for x in recoChannels[analysis][nTaus]]) + ')'

        optimizer = initializeOptimizer(analysis, period, 'plots_optimize_%iTau' %(nTaus), 'select.PassTight & %s' % recoSelect, genSelect, nTaus)
        
        optimizer.addCut('st',         ['finalstate.sT >'],            100., 1500., 10.)
        optimizer.addCut('met',        ['finalstate.met >'],             0.,  500.,  5.)
        if analysis in ['Hpp3l']:
            optimizer.addCut('dR',         ['h1.dR <'],                      1.,    5.,  0.05)
        else:
            optimizer.addCut('dR',         ['h1.dR <','h2.dR <'],                                         1.,    5.,  0.05)
        if analysis in ['Hpp3l']:
            optimizer.addCut('zmass',      ['fabs(z1.mass-%f) >' % ZMASS],   0.,  200.,  5.)
            optimizer.addCut('zmassUnder', ['%f-z1.mass >' % ZMASS],         0.,  200.,  5.)
            optimizer.addCut('zmassOver',  ['z1.mass-%f >' % ZMASS],         0.,  200.,  5.)
        else:
            optimizer.addCut('zmass',      ['fabs(z1.mass-%f) >' % ZMASS,'fabs(z2.mass-%f) >' % ZMASS],   0.,  200.,  5.)
            optimizer.addCut('zmassUnder', ['%f-z1.mass >' % ZMASS,'%f-z2.mass >' % ZMASS],               0.,  200.,  5.)
            optimizer.addCut('zmassOver',  ['z1.mass-%f >' % ZMASS,'z2.mass-%f >' % ZMASS],               0.,  200.,  5.)
        if analysis in ['Hpp3l']:
            optimizer.addCut('hmass',      ['fabs(h1.mass-MASS) <'],         0.,  500.,  5.)
            optimizer.addCut('hmassUnder', ['MASS-h1.mass >'],               0.,  500.,  5.)
            optimizer.addCut('hmassOver',  ['h1.mass-MASS >'],               0.,  500.,  5.)
        else:
            optimizer.addCut('hmass',      ['fabs(h1.mass-MASS) <','fabs(h2.mass-MASS) <'],               0.,  500.,  5.)
            optimizer.addCut('hmassUnder', ['MASS-h1.mass >','MASS-h2.mass >'],                           0.,  500.,  5.)
            optimizer.addCut('hmassOver',  ['h1.mass-MASS >','h2.mass-MASS >'],                           0.,  500.,  5.)
        #optimizer.addCut('st',         ['finalstate.sT >'],            100., 1500., 200.)
        #optimizer.addCut('dR',         ['h1.dR <'],                      1.,    5.,  1.)
        #optimizer.addCut('zmass',      ['fabs(z1.mass-%f) >' % ZMASS],   0.,  100.,  20.)
        #optimizer.addCut('hmass',      ['fabs(h1.mass-MASS) <'],         0.,  400.,  100.)
        #optimizer.addCut('hmassUnder', ['h1.mass <'],                    0., 1000., 200.)
        #optimizer.addCut('hmassOver',  ['h1.mass >'],                    0., 1000., 200.)

        print 'Optimizing'
        masses = _3L_MASSES if analysis in ['Hpp3l'] else _4L_MASSES
        optimizer.optimize(masses=masses)

def parse_command_line(argv):
    parser = get_parser("Merge the output ISA ntuples")

    args = parser.parse_args(argv)
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500,help='Mass for signal')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses for signal')
    parser.add_argument('-ub','--unblind',action='store_true',help='unblind')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    optimize(args.analysis, args.period)

if __name__ == "__main__":
    main()
