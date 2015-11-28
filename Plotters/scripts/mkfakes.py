#!/usr/bin/env python

from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.ShapePlotter import ShapePlotter
from InitialStateAnalysis.Plotters.CutFlowPlotter import CutFlowPlotter
from InitialStateAnalysis.Plotters.EfficiencyPlotter import EfficiencyPlotter
from InitialStateAnalysis.Plotters.FakeRatePlotter import FakeRatePlotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS
from InitialStateAnalysis.Utilities.utilities import *
import argparse
import itertools
import sys
import pickle
import json

def makeFakes(analysis,channel,runPeriod,**kwargs):
    '''Plot fake rate for an analysis.'''
    myCut = kwargs.pop('myCut','1')
    scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    mass = kwargs.pop('mass',500)
    runTau = kwargs.pop('runTau',0)
    for key, value in kwargs.iteritems():
        print "Unrecognized parameter '" + key + "' = " + str(value)
        return 0

    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ'   : 3,
        'WZ_W' : 2,
        'WZ_Dijet' : 1,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[analysis]
    ntuples = getNtupleDirectory(analysis,channel,runPeriod)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = getChannelBackgrounds(runPeriod)

    finalStates, leptons = getChannels(nl,runTau=runTau)
    print 'MKFAKES:%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut)
    dataplot = runPeriod in [7,8]
    mergeDict = getMergeDict(runPeriod)

    # define fake regions
    fakeRegions, ptBins, etaBins = getFakeParams(analysis)
    fakes = {}
    for fakeRegion in fakeRegions[analysis]:
        logger.info("%s:%s:%iTeV: Fake Region: %s" % (analysis,channel, runPeriod, fakeRegion))
        denom = fakeRegions[analysis][fakeRegion]['denom']
        numer = fakeRegions[analysis][fakeRegion]['numer']
        probe = fakeRegions[analysis][fakeRegion]['probe']
        ptvar = fakeRegions[analysis][fakeRegion]['ptVar']
        etavar= fakeRegions[analysis][fakeRegion]['etaVar']
        ptcut = '{0} >= {1} && {0} < {2}'
        etacut = 'abs({0}) >= {1} && abs({0}) < {2}'

        logger.info("%s:%s:%iTeV: Computing fake rates" % (analysis,channel, runPeriod))
        plotter = FakeRatePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='{0}_fakerates'.format(fakeRegion),mergeDict=mergeDict,scaleFactor='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale*event.trig_prescale',dataScaleFactor='event.trig_prescale')
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        hist = plotter.getFakeRate(numer, denom, ptBins, etaBins[probe], ptvar, etavar, fakeRegion, dataDriven=True, subtractSamples=['WZJets','ZZJets','TTVJets','VVVJets'])
        fakes[fakeRegion] = []
        for p,pt in enumerate(ptBins[:-1]):
            for e,eta in enumerate(etaBins[probe][:-1]):
                ptBin = p+1
                etaBin = e+1
                ptLow = pt
                etaLow = eta
                ptHigh = ptBins[p+1]
                etaHigh = etaBins[probe][e+1]
                fakekey = 'pt_{0:.3f}-{1:.3f}_eta_{2:.3f}-{3:.3f}'.format(ptLow,ptHigh,etaLow,etaHigh)
                fakerate = hist.GetBinContent(ptBin,etaBin)
                fakeerror = hist.GetBinError(ptBin,etaBin)
                fakes[fakeRegion] += [{'fakerate':fakerate,'error':fakeerror,'pt_low':ptLow,'pt_high':ptHigh,'eta_low':etaLow,'eta_high':etaHigh}]

    print json.dumps(fakes,sort_keys=True,indent=4)
    # save pickle file
    with open('fakes.pkl','wb') as f:
        pickle.dump(fakes,f)
    # save json file
    with open('fakes.json','w') as f:
        json.dump(fakes,f,sort_keys=True,indent=4)
                
               

def parse_command_line(argv):
    parser = get_parser("Plot a given channel and period")

    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to plots.')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for plots.')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    makeFakes(args.analysis,args.channel,args.period)

    return 0


if __name__ == "__main__":
    main()
