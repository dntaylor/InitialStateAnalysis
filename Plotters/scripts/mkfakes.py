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
    fakeRegions, ptBins, etaBins = getFakeParams(analysis)
    #lepName = {'e': 'Elec', 'm': 'Muon', 't': 'Tau'}
    #fakeRegions = {}
    #fakeRegions['WZ'] = {}
    #for f in ['e', 'm']:
    #    #for p in ['Loose', 'Tight']:
    #    for p in ['Tight']:
    #        # select leading Z pt, Z window [60,120], tight (or loose) Z, low met, m3l>100, w1 mass < 30
    #        if analysis in ['WZ']:
    #            #for z in ['Loose', 'Tight']:
    #            for z in ['Tight']:
    #                fakeRegion = 'Z{0}Probe{1}{2}'.format(z,lepName[f],p)
    #                #denom = 'z1.Pass{0}1 && z1.Pass{0}2 && finalstate.mass>100. && w1.mass<20. && w1.dR1_z1_1>0.1 && w1.dR1_z1_2>0.1 && w1Flv=="{1}"'.format(z,f)
    #                #denom = 'z1.Pass{0}1 && z1.Pass{0}2 && w1.dR1_z1_1>0.02 && w1.dR1_z1_2>0.02 && w1.mass<25 && w1Flv=="{1}"'.format(z,f)
    #                denom = 'z1.Pass{0}1 && z1.Pass{0}2 && w1.dR1_z1_1>0.02 && w1.dR1_z1_2>0.02 && finalstate.met<25. && w1Flv=="{1}"'.format(z,f)
    #                #denom = 'z1.Pass{0}1 && z1.Pass{0}2 && w1Flv=="{1}"'.format(z,f)
    #                numer = '{0} && w1.Pass{1}1'.format(denom,p)
    #                fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
    #                #if p=='Tight' and channel not in ['HZZFakeRate']:
    #                #   fakeRegion += '_LooseProbe'
    #                #   denom += ' && w1.PassLoose1'
    #                #   numer += ' && w1.PassLoose1'
    #                #   fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
    #        # select w lepton pt, z veto, met
    #        #'W' : 'w1.Pt1>20. & (z1.mass<60. | z1.mass>120.) & finalstate.met>30. & w1.mass>30.',
    #        if analysis in ['WZ_W']:
    #            #for w in ['Loose','Tight']:
    #            for w in ['Tight']:
    #                fakeRegion = 'W{0}Probe{1}{2}'.format(w,lepName[f],p)
    #                denom = 'w1.Pt1>20. && w1.mass>30. && finalstate.met>30. && (z1.mass<60. || z1.mass>120.) && l1.Chg==l2.Chg && z1.dR>0.1 && w1.Pass{0}1 && w2Flv=="{1}"'.format(w,f)
    #                numer = '{0} && w2.Pass{1}1'.format(denom,p)
    #                fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
    #                #if p=='Tight' and channel not in ['HZZFakeRate']:
    #                #   fakeRegion += '_LooseProbe'
    #                #   denom += ' && w2.PassLoose1'
    #                #   numer += ' && w2.PassLoose1'
    #                #   fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
    #        # ntuple cuts: zVeto 60-120, met vet 20, w veto 20, jet pt > 20, jet dr > 1.0
    #        if analysis in ['WZ_Dijet']:
    #           fakeRegion = 'FakeRateProbe{0}{1}'.format(lepName[f],p)
    #           denom = 'l1Flv=="{0}"'.format(f)
    #           numer = '{0} && w1.Pass{1}1'.format(denom,p)
    #           fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
    #           #if p=='Tight' and channel not in ['HZZFakeRate']:
    #           #   fakeRegion += '_LooseProbe'
    #           #   denom += ' && w1.PassLoose1'
    #           #   numer += ' && w1.PassLoose1'
    #           #   fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}




    ## setup selections
    #ptBins = [0,10,20,30,200]
    #if analysis in ['WZ_Dijet']: ptBins = [0,10,20,30,40,60,80,200]
    #etaBins = {
    #    'e': [0.,1.47,2.5],
    #    'm': [0.,1.2,2.4],
    #    #'e': [0.,2.5],
    #    #'m': [0.,2.4],
    #}
    fakes = {}
    for fakeRegion in fakeRegions['WZ']:
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
        hist = plotter.getFakeRate(numer, denom, ptBins, etaBins[probe], ptvar, etavar, fakeRegion, dataDriven=True, subtractSamples=['WZJets','ZZJets'])
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
