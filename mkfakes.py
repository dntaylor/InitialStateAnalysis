#!/usr/bin/env python

from plotters.Plotter import Plotter
from plotters.ShapePlotter import ShapePlotter
from plotters.CutFlowPlotter import CutFlowPlotter
from plotters.EfficiencyPlotter import EfficiencyPlotter
from plotters.FakeRatePlotter import FakeRatePlotter
from plotters.plotUtils import *
from plotters.plotUtils import ZMASS
import argparse
import itertools
import sys
import pickle
import json

def makeFakes(analysis,channel,runPeriod,**kwargs):
    '''Plot fake rate for an analysis.'''
    myCut = kwargs.pop('myCut','1')
    scaleFactor = kwargs.pop('scaleFactor','event.pu_weight*event.lep_scale*event.trig_scale')
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
        'WZ_W' : 3,
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
    print 'MKFAKES:%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut)
    dataplot = runPeriod in [7,8]
    mergeDict = getMergeDict(runPeriod)

    # define fake regions
    lepName = {'e': 'Elec', 'm': 'Muon', 't': 'Tau'}
    fakeRegions = {}
    fakeRegions['WZ'] = {}
    for f in ['e', 'm']:
        for p in ['Loose', 'Tight']:
            # select leading Z pt, Z window [60,120], tight (or loose) Z, low met, m3l>100, w1 mass < 30
            if analysis in ['WZ']:
                for z in ['Loose', 'Tight']:
                    fakeRegion = 'Z{0}Probe{1}{2}'.format(z,lepName[f],p)
                    denom = 'z1.Pt1>20. & z1.mass>60. & z1.mass<120. & z1.Pass{0}1 & z1.Pass{0}2 & finalstate.met<20. & finalstate.mass>100. & w1.mass<20. & w1.dR1_z1_1>0.1 & w1.dR1_z1_2>0.1 & w1Flv=="{1}"'.format(z,f)
                    numer = '{0} & w1.Pass{1}1'.format(denom,p)
                    fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
                    if p=='Tight':
                       fakeRegion += '_LooseProbe'
                       denom += ' & w1.PassLoose1'
                       numer += ' & w1.PassLoose1'
                       fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
            # select w lepton pt, z veto, met
            #'W' : 'w1.Pt1>20. & (z1.mass<60. | z1.mass>120.) & finalstate.met>30. & w1.mass>30.',
            if analysis in ['WZ_W']:
                for w in ['Loose','Tight']:
                    fakeRegion = 'W{0}Probe{1}{2}'.format(w,lepName[f],p)
                    denom = 'w1.Pt1>20. & w1.mass>30. & finalstate.met>30. & (z1.mass<60. | z1.mass>120.) & l1.Chg==l2.Chg & z1.dR>0.1 & w1.Pass{0}1 & w2Flv=="{1}"'.format(w,f)
                    numer = '{0} & w2.Pass{1}1'.format(denom,p)
                    fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
                    if p=='Tight':
                       fakeRegion += '_LooseProbe'
                       denom += ' & w2.PassLoose1'
                       numer += ' & w2.PassLoose1'
                       fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}

            # veto Z, low met, w mass veto
            #'QCD' : '(z1.mass<60. | z1.mass>120.) & finalstate.met<20. & w1.mass<20.',

    # setup selections
    ptBins = [10,20,40,100,1000]
    etaBins = {
        'e': [0,1.479,2.5],
        'm': [0,1.2,2.4],
    }
    fakes = {}
    for fakeRegion in fakeRegions['WZ']:
        print "MKPLOTS:%s:%s:%iTeV: Fake Region: %s" % (analysis,channel, runPeriod, fakeRegion)
        denom = fakeRegions['WZ'][fakeRegion]['denom']
        numer = fakeRegions['WZ'][fakeRegion]['numer']
        probe = fakeRegions['WZ'][fakeRegion]['probe']
        ptvar = fakeRegions['WZ'][fakeRegion]['ptVar']
        etavar = fakeRegions['WZ'][fakeRegion]['etaVar']

        # now plot the fake rates
        print "MKFAKES:%s:%s:%iTeV: Computing fake rates" % (analysis,channel, runPeriod)
        plotter = FakeRatePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='{0}_fakerates'.format(fakeRegion),mergeDict=mergeDict)
        # this should be done in data... using mc until we get some!
        # subtract WZ, ZZ, ttbar contributions from data
        # only initialize Z... in MC?
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        #plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        hist = plotter.getFakeRate(numer, denom, ptBins, etaBins[probe], ptvar, etavar, fakeRegion)
        fakes[fakeRegion] = {}
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
                fakes[fakeRegion][fakekey] = {'fakerate':fakerate,'error':fakeerror}

    print json.dumps(fakes,sort_keys=True,indent=4)
    # save pickle file
    pickle.dump(fakes,open('fakes.pkl','wb'))
                
               

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['Z','WZ','WZ_W','Hpp2l','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['Z','WZ','W','TT','Hpp2l','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to plots.')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for plots.')
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
