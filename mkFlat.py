#!/usr/bin/env python

# A script to make flat histograms given a cut

from plotters.Plotter import Plotter
from plotters.CutFlowPlotter import CutFlowPlotter
from plotters.plotUtils import *
from plotters.plotUtils import ZMASS
import glob
import argparse
import itertools
import sys
import os
import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

def generate(analysis, channel, period, **kwargs):
    '''
    Generate root files
    Files are of the form:
      AnalysisChannel
        Channel0
          Variable0
            sample0
            sample1
            ...
            data
          Variable1
          ...
    '''
    cut = kwargs.pop('cut','1')
    filename = kwargs.pop('filename','variables')
    scaleFactor = kwargs.pop('scaleFactor','event.pu_weight*event.lep_scale*event.trig_scale')
    print 'MKFLAT:%s:%s:%iTeV: Running with selection %s' % (analysis, channel, period, cut)
    rootPath = 'rootfiles/%s_%s_%iTeV' % (analysis,channel, period)
    python_mkdir(rootPath)
    savefilename = '%s/%s.root' %(rootPath,filename)
    print 'MKFLAT:%s:%s:%iTeV: Will save in %s' % (analysis, channel, period, savefilename)
    # potentially hash the cut, and the input ntuples so it doesnt get recalculated???
    savefile = ROOT.TFile(savefilename,'recreate')
    selectionDict = {
        # name             : (variable,                        binning,     selection),
        'sT'               : (['finalstate.sT'],               [40,0,1000], cut),
        'numJets30'        : (['finalstate.jetVeto30'],        [8,0,8],     cut),
        'elecVetoLoose'    : (['finalstate.elecVetoLoose'],    [8,0,8],     cut),
        'muonVetoLoose'    : (['finalstate.muonVetoLoose'],    [8,0,8],     cut),
        'elecVetoTight'    : (['finalstate.elecVetoTight'],    [8,0,8],     cut),
        'muonVetoTight'    : (['finalstate.muonVetoTight'],    [8,0,8],     cut),
        'bjetVeto30Loose'  : (['finalstate.bjetVeto30Loose'],  [8,0,8],     cut),
        'bjetVeto30Medium' : (['finalstate.bjetVeto30Medium'], [8,0,8],     cut),
        'bjetVeto30Tight'  : (['finalstate.bjetVeto30Tight'],  [8,0,8],     cut),
        'met'              : (['finalstate.met'],              [40,0,200],  cut),
        'mass'             : (['finalstate.mass'],             [40,0,400],  cut),
        'puVertices'       : (['event.nvtx'],                  [50,0,50],   cut),
        
    }
    if analysis in ['Z', 'Hpp3l', 'Hpp4l', 'WZ'] or region in ['Z', 'TT']:
        selectionDict['zMass']     = (['z1.mass'], [42,70,112], cut)
        selectionDict['zMassFull'] = (['z1.mass'], [80,0,240],  cut)
        selectionDict['zPt']       = (['z1.Pt'],   [40,0,400],  cut)
        selectionDict['zPt1']      = (['z1.Pt1'],  [40,0,200],  cut)
        selectionDict['zPt2']      = (['z1.Pt2'],  [40,0,200],  cut)
        selectionDict['zIso1']     = (['z1.Iso1'], [50,0,0.5],  cut)
        selectionDict['zIso2']     = (['z1.Iso2'], [50,0,0.5],  cut)
        selectionDict['zDR']       = (['z1.dR'],   [60,0,6],    cut)
    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ'   : 3,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[channel]
    channels, leptons = getChannels(nl)
    ntuples = 'ntuples%s_%stev_%s' % (analysis,period,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,period)
    intLumiMap = getIntLumiMap()
    finalStates, leptons = getChannels(nl)
    mergeDict = getMergeDict(runPeriod)
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,scaleFactor=scaleFactor)
    allSamples = [os.path.basename(fname).rstrip('.root') for fname in glob.glob('%s/*'%ntuples)]
    bgSamples = [x for x in allSamples if 'data' not in x]
    dataSamples = [x for x in allSamples if 'data' in x]
    plotter.initializeBackgroundSamples(bgSamples)
    if dataSamples: plotter.initializeDataSamples(dataSamples)
    plotter.setIntLumi(intLumiMap[period])
    histNames = bgSamples + ['data']
    adir = savefile.mkdir(channel)
    adir.cd()
    print 'MKFLAT:%s:%s:%iTeV: Creating analysis directory' % (analysis, channel, period)
    for c in channels:
        cdir = adir.mkdir(c)
        cdir.cd()
        print 'MKFLAT:%s:%s:%iTeV: Channel %s' % (analysis, channel, period, c)
        for name, (variable, binning, selection) in selectionDict.iteritems():
            vdir = cdir.mkdir(name)
            vdir.cd()
            print 'MKFLAT:%s:%s:%iTeV: Variable %s' % (analysis, channel, period, name)
            chanSelection = selection + ' & channel=="%s"' %c
            for histName in histNames:
                hist = plotter.getData(variable,binning,chanSelection,True) if histName == 'data' else\
                       plotter.getHist(histName,variable,binning,chanSelection,True)
                if not hist: continue
                hist.SetName(histName)
                hist.Write()
            cdir.cd()
        adir.cd()
    print 'MKFLAT:%s:%s:%iTeV: Finished' % (analysis, channel, period)
    savefile.Close()
    return 0


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-c','--cut',type=str,default='select.passTight',help='Cut to be applied to plots (default = "select.passTight").')
    parser.add_argument('-fn','--filename',type=str,default='variables',help='Filename of output (default = "variables").')

    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    generate(args.analysis, args.channel, args.period, cut=args.cut, filename=args.filename)

    return 0

if __name__ == "__main__":
    main()
