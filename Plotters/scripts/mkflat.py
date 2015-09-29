#!/usr/bin/env python

# A script to make flat histograms given a cut

from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.CutFlowPlotter import CutFlowPlotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS
from InitialStateAnalysis.Utilities.utilities import *
import glob
import argparse
import itertools
import sys
import os
import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")

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
    logger = logging.getLogger(__name__)
    cut = kwargs.pop('cut','1')
    scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    force = kwargs.pop('force','False')

    if force: logger.info('%s:%s:%iTeV: Forcing reprocessing' % (analysis, channel, period))

    # get hashes
    cut, cuthash = hashcut(cut)
    scaleFactor, scalefactorhash = hashscalefactor(scaleFactor)

    logger.info('%s:%s:%iTeV: Running with selection %s' % (analysis, channel, period, cut))
    rootPath = 'rootfiles/%s_%s_%iTeV/%s/%s' % (analysis,channel, period,cuthash,scalefactorhash)
    python_mkdir(rootPath)

    selectionDict = {
        # name             : (variable,                        binning,     selection),
        'sT'               : (['finalstate.sT'],               [40,0,1000], cut),
        'numJets30'        : (['finalstate.jetVeto30'],        [8,0,8],     cut),
        'elecVetoLoose'    : (['finalstate.elecVetoLoose'],    [8,0,8],     cut),
        'muonVetoLoose'    : (['finalstate.muonVetoLoose'],    [8,0,8],     cut),
        'elecVetoTight'    : (['finalstate.elecVetoTight'],    [8,0,8],     cut),
        'muonVetoTight'    : (['finalstate.muonVetoTight'],    [8,0,8],     cut),
        'bjetVeto30Medium' : (['finalstate.bjetVeto30Medium'], [8,0,8],     cut),
        'met'              : (['finalstate.met'],              [40,0,200],  cut),
        'mass'             : (['finalstate.mass'],             [40,0,400],  cut),
        'puVertices'       : (['event.nvtx'],                  [50,0,50],   cut),
        
    }
    if analysis in ['Hpp3l','Hpp4l'] or region in ['Hpp2l']:
        selectionDict['hppMass']   = (['h1.mass'], [24,0,600],  cut)
        selectionDict['hppDPhi']   = (['h1.dPhi'], [32,0,3.2],  cut)
        selectionDict['hppPt']     = (['h1.Pt'],   [40,0,400],  cut)
        selectionDict['hppPt1']    = (['h1.Pt1'],  [40,0,200],  cut)
        selectionDict['hppPt2']    = (['h1.Pt2'],  [40,0,200],  cut)
        selectionDict['hppIso1']   = (['h1.Iso1'], [50,0,0.5],  cut)
        selectionDict['hppIso2']   = (['h1.Iso2'], [50,0,0.5],  cut)
        selectionDict['hppDR']     = (['h1.dR'],   [60,0,6],    cut)
    if analysis in ['Z', 'Hpp3l', 'Hpp4l', 'WZ'] or region in ['Z', 'TT']:
        selectionDict['zMass']     = (['z1.mass'], [60,60,120], cut)
        selectionDict['zMassFull'] = (['z1.mass'], [80,0,240],  cut)
        selectionDict['zPt']       = (['z1.Pt'],   [40,0,400],  cut)
        selectionDict['zPt1']      = (['z1.Pt1'],  [40,0,200],  cut)
        selectionDict['zPt2']      = (['z1.Pt2'],  [40,0,200],  cut)
        selectionDict['zIso1']     = (['z1.Iso1'], [50,0,0.5],  cut)
        selectionDict['zIso2']     = (['z1.Iso2'], [50,0,0.5],  cut)
        selectionDict['zDR']       = (['z1.dR'],   [60,0,6],    cut)
    if analysis in ['Hpp3l','WZ']:
        selectionDict['wPt']       = (['w1.Pt'],   [40,0,400],  cut)
        selectionDict['wPt1']      = (['w1.Pt1'],  [40,0,200],  cut)
        selectionDict['wIso1']     = (['w1.Iso1'], [50,0,0.5],  cut)
        selectionDict['wMass']     = (['w1.mass'], [40,0,200],  cut)
        selectionDict['wDPhi']     = (['w1.dPhi'], [32,0,3.2],  cut)

    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ'   : 3,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[channel]

    for l in range(nl):
        name = 'l%i' % (l+1)
        for f in ['','e','m']:
            thisCut = cut
            if f:
                thisCut = cut + ' & %sFlv=="%s"' % (name,f)
            selectionDict['%sPt%s'%(name,f.upper())]  = (['%s.Pt' %name],  [40,0,200],            thisCut)
            selectionDict['%sEta%s'%(name,f.upper())] = (['%s.Eta' %name], [30,-3.0,3.0],         thisCut)
            selectionDict['%sPhi%s'%(name,f.upper())] = (['%s.Phi' %name], [30,-3.14159,3.14159], thisCut)
            selectionDict['%sIso%s'%(name,f.upper())] = (['%s.Iso' %name], [50,0,.5],             thisCut)

    # setup plotter
    channels, leptons = getChannels(nl)
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,period,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,period)
    intLumiMap = getIntLumiMap()
    finalStates, leptons = getChannels(nl)
    mergeDict = getMergeDict(period)
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,scaleFactor=scaleFactor)
    allSamples = [os.path.basename(fname).rstrip('.root') for fname in glob.glob('%s/*'%ntuples)]
    bgSamples = [x for x in allSamples if 'data' not in x]
    dataSamples = [x for x in allSamples if 'data' in x]
    plotter.initializeBackgroundSamples(bgSamples)
    if dataSamples: plotter.initializeDataSamples(dataSamples)
    plotter.setIntLumi(intLumiMap[period])
    histNames = bgSamples
    if dataSamples: histNames += ['data']

    for sample in allSamples:
        logger.info('%s:%s:%iTeV: Sample %s' % (analysis, channel, period, sample))
        infilename = '%s/%s.root' % (ntuples,sample)
        filehash = hashfile(infilename)
        savedir = '%s/%s' % (rootPath,sample)
        python_mkdir(savedir)
        savefilename = '%s/%s.root' %(savedir,filehash)
        if os.path.isfile(savefilename) and not force: continue # already exists, don't need to do anything
        savefile = ROOT.TFile(savefilename,'recreate')
        adir = savefile.mkdir(channel)
        adir.cd()
        for name, (variable, binning, selection) in selectionDict.iteritems():
            logger.debug('%s:%s:%iTeV: Variable %s' % (analysis, channel, period, name))
            hist = plotter.getHist(sample,variable,binning,selection,True)
            if not hist: continue
            hist.SetName(name)
            hist.Write()
            adir.cd()
    logger.info('%s:%s:%iTeV: Finished' % (analysis, channel, period))
    savefile.Close()
    return 0




def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-c','--cut',type=str,default='select.passTight',help='Cut to be applied to plots (default = "select.passTight").')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')
    parser.add_argument('-f','--force',action='store_true',help='Force reprocessing')

    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    generate(args.analysis, args.channel, args.period, cut=args.cut, force=args.force)

    return 0

if __name__ == "__main__":
    main()
