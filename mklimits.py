#!/usr/bin/env python
from limits.Limits import Limits
import logging
import sys
import argparse
import numpy as np
from plotters.plotUtils import _3L_MASSES, _4L_MASSES, getSigMap, getIntLumiMap, getChannels, getMergeDict, ZMASS
from plotters.Plotter import Plotter

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class Scales(object):
    def __init__(self, br_ee, br_em, br_et, br_mm, br_mt, br_tt):
        self.a_3l = np.array([br_ee, br_em, br_et, br_mm, br_mt, br_tt], dtype=float)
        self.m_4l = np.outer(self.a_3l, self.a_3l)
        self.index = {"ee": 0, "em": 1, "et": 2, "mm": 3, "mt": 4, "tt": 5}
    def scale_Hpp4l(self, hpp, hmm):
        i = self.index[hpp]
        j = self.index[hmm]
        return self.m_4l[i,j] * 36.0
    def scale_Hpp3l(self, hpp, hm='a'):
        i = self.index[hpp]
        return self.a_3l[i] * 6.0
        

def limit(analysis,period,mass,**kwargs):
    doChannels = kwargs.pop('doChannels',False)
    name = kwargs.pop('name','card')
    scale = kwargs.pop('scale',[1.0])
    directory = kwargs.pop('directory','')
    chans = kwargs.pop('channels',['1'])
    mode = kwargs.pop('mode','mc')
    scalefactor = kwargs.pop('scalefactor','event.pu_weight*event.lep_scale*event.trig_scale')
    logger.info("Processing mass-point %i" % mass)

    channels, leptons = getChannels(3 if analysis=='Hpp3l' or analysis=='WZ' else 4)
    #cutMap = defineCutFlowMap(analysis,channels,mass)
    cutMap = {
        'Hpp3l' : {
             'cuts' : ['1',\
                       'finalstate.sT>1.1*%f+60.' %mass,\
                       'fabs(z1.mass-%f)>80.' %ZMASS,\
                       'h1.dPhi<%f/600.+1.95' %mass,\
                       'h1.mass>0.9*%f & h1.mass<1.1*%f' %(mass,mass)],
             'labels' : ['Preselection','s_{T}','Z Veto','#Delta#phi','Mass window']
        },
        'Hpp4l' : {
             'cuts' : ['1',\
                       'finalstate.sT>0.6*%f+130.' %mass,\
                       'h1.mass>0.9*%f & h1.mass<1.1*%f' %(mass,mass)],
             'labels' : ['Preselection','s_{T}','Mass window']
        }
    }

    nl = 3 if analysis=='Hpp3l' or analysis=='WZ' else 4
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()

    cuts = '&&'.join(cutMap[analysis]['cuts'])
    #if chans:
    #    chanCut = '('+'||'.join(['channel=="%s"'%c for c in chans])+')'

    limits = Limits(analysis, period, cuts, './ntuples%s_%itev_%s' % (analysis, period, analysis),
                    './datacards/%s_%itev/%s/%s' % (analysis, period, directory, mass),
                    channels=['dblh%s' % analysis], lumi=intLumiMap[period],
                    blinded=True, bgMode=mode, scalefactor=scalefactor)

    signal =  sigMap[period]['Sig']
    if mode=='mc':
        add_systematics_mc(limits,mass,signal,name,chans,scale,period)
    elif mode=='sideband':
        add_systematics_sideband(limits,mass,signal,name,chans,scale,period)
    elif mode=='fakerate':
        add_systematics_fakerate(limits,mass,signal,name,chans,scale,period)
    else:
        return 0

def BP(analysis,period,mass,bp,**kwargs):
    if bp == 'ee100':
        s = Scales(1., 0., 0., 0., 0., 0.)
    elif bp == 'em100':
        s = Scales(0., 1., 0., 0., 0., 0.)
    elif bp == 'mm100':
        s = Scales(0., 0., 0., 1., 0., 0.)
    elif bp == 'BP1':
        s = Scales(0, 0.1, 0.1, 0.3, 0.38, 0.3)
    elif bp == 'BP2':
        s = Scales(0.5, 0, 0, 0.125, 0.25, 0.125)
    elif bp == 'BP3':
        s = Scales(0.34, 0, 0, 0.33, 0, 0.33)
    elif bp == 'BP4':
        s = Scales(1./6., 1./6., 1./6., 1./6., 1./6., 1./6.)
    else:
        print 'Unknown branching point: %s' %bp
    logger.info("Processing branching point %s" % bp)
    sf = getattr(s,'scale_%s'%analysis)
    chanMap = {
        'Hpp3l': {
             'names': ['ee','em','mm'],
             'ee'   : ['eee','eem'],
             'em'   : ['eme','emm','mee','mem'],
             'mm'   : ['mme','mmm'],
        },
        'Hpp4l': {
             'names': ['eeee','eeem','eemm','emem','emmm','mmmm'],
             'eeee' : ['eeee'],
             'eeem' : ['eeem','eeme','emee','meee'],
             'eemm' : ['eemm','mmee'],
             'emem' : ['emem','emme','meem','meme'],
             'emmm' : ['emmm','memm','mmem','mmme'],
             'mmmm' : ['mmmm'],
        },
    }
    #for c in chanMap[analysis]['names']:
    #    thisScale = sf(c[:2],c[3:])
    #    if thisScale==0: continue
    #    limit(analysis,period,mass,name=c,directory=bp,channels=chanMap[analysis][c],scale=thisScale,**kwargs)

    chanCuts = []
    chanScales = []
    for c in chanMap[analysis]['names']:
        thisScale = sf(c[:2],c[2:])
        if thisScale==0: continue
        chanCut = '('+'||'.join(['channel=="%s"'%x for x in chanMap[analysis][c]])+')'
        chanCuts += [chanCut]
        chanScales += [thisScale]
    limit(analysis,period,mass,name=bp,directory=bp,channels=chanCuts,scale=chanScales,**kwargs)

def add_systematics_mc(limits,mass,signal,name,chans,sigscale,period):
    limits.add_group("hpp%i" % mass, signal, scale=sigscale, isSignal=True)
    if period=='8': limits.add_group("dyjets", "Z*j*")
    if period=='13': limits.add_group("dyjets", "DY*")
    limits.add_group("zz", "ZZJ*")
    limits.add_group("wz", "WZJ*")
    if period=='8': limits.add_group("ww", "WWJ*")
    if period=='8': limits.add_group("zzz", "ZZZ*")
    if period=='8': limits.add_group("wwz", "WWZ*")
    if period=='8': limits.add_group("www", "WWW*")
    limits.add_group("top", "T[(B|b)ar]_*")
    limits.add_group("tt", "TTJ*")
    limits.add_group("ttz", "TTZJ*")
    limits.add_group("ttw", "TTWJ*")
    if period=='8': limits.add_group("data", "data_R*", isData=True)

    lumi = {
        'hpp%i' % mass: 1.026,
        'dyjets':       1.026,
        'zz':           1.026,
        'wz':           1.026,
        'ww':           1.026,
        'zzz':          1.026,
        'wwz':          1.026,
        'www':          1.026,
        'tt':           1.026,
        'ttz':          1.026,
        'ttw':          1.026,
        'top':          1.026
    }
    limits.add_systematics("lumi", "lnN", **lumi)

    eid = {
        'hpp%i' % mass: 1.01,
        'dyjets':       1.01,
        'zz':           1.01,
        'wz':           1.01,
        'ww':           1.01,
        'zzz':          1.01,
        'wwz':          1.01,
        'www':          1.01,
        'tt':           1.01,
        'ttz':          1.01,
        'ttw':          1.01,
        'top':          1.01
    }
    limits.add_systematics("eid", "lnN", **eid)

    muid = {
        'hpp%i' % mass: 1.005,
        'dyjets':       1.005,
        'zz':           1.005,
        'wz':           1.005,
        'ww':           1.005,
        'zzz':          1.005,
        'wwz':          1.005,
        'www':          1.005,
        'tt':           1.005,
        'ttz':          1.005,
        'ttw':          1.005,
        'top':          1.005
    }
    limits.add_systematics("muid", "lnN", **muid)

    muiso = {
        'hpp%i' % mass: 1.002,
        'dyjets':       1.002,
        'zz':           1.002,
        'wz':           1.002,
        'ww':           1.002,
        'zzz':          1.002,
        'wwz':          1.002,
        'www':          1.002,
        'tt':           1.002,
        'ttz':          1.002,
        'ttw':          1.002,
        'top':          1.002
    }
    limits.add_systematics("muiso", "lnN", **muiso)

    sigmc = {'hpp%i' % mass: 1.15}
    limits.add_systematics("sigmc", "lnN", **sigmc)

    mcdata = {
        'zz':           1.105,
        'wz':           1.056,
        'ww':           1.041,
        'tt':           1.024
    }
    limits.add_systematics("mcdata", "lnN", **mcdata)

    pdf = {
        'zzz':          1.026,
        'wwz':          1.051,
        'www':          1.043,
        'ttz':          1.105,
        'ttw':          1.289
    }
    limits.add_systematics("pdf", "lnN", **pdf)

    limits.gen_card("%s.txt" % name,mass=mass,cuts=chans)

def calculateLeptonSystematic(mass,bp):
    analysis = 'Hpp3l'
    region = 'Hpp3l'
    runPeriod = 8
    nl = 3
    ntuples = 'ntuples%s_%itev_%s' % (analysis,runPeriod,region)
    saves = '%s_%s_%sTeV' % (analysis,region,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    mergeDict = getMergeDict(runPeriod)
    regionBackground = {
        'Hpp3l' : ['T','TT', 'TTV','W','Z','VVV','WW','ZZ','WZ'],
        'Hpp4l' : ['TT','Z','DB']
    }
    channels, leptons = getChannels(nl)

    plotter = Plotter(analysis,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='systematics',mergeDict=mergeDict)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in regionBackground[analysis]])
    plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])

    chanMap = {
        'Hpp3l': {
             'names': ['ee','em','mm'],
             'ee'   : ['eee','eem'],
             'em'   : ['eme','emm','mee','mem'],
             'mm'   : ['mme','mmm'],
        },
        'Hpp4l': {
             'names': ['eeee','eeem','eemm','emem','emmm','mmmm'],
             'eeee' : ['eeee'],
             'eeem' : ['eeem','eeme','emee','meee'],
             'eemm' : ['eemm','mmee'],
             'emem' : ['emem','emme','meem','meme'],
             'emmm' : ['emmm','memm','mmem','mmme'],
             'mmmm' : ['mmmm'],
        },
    }
    if bp == 'ee100':
        s = Scales(1., 0., 0., 0., 0., 0.)
    elif bp == 'em100':
        s = Scales(0., 1., 0., 0., 0., 0.)
    elif bp == 'mm100':
        s = Scales(0., 0., 0., 1., 0., 0.)
    elif bp == 'BP1':
        s = Scales(0, 0.1, 0.1, 0.3, 0.38, 0.3)
    elif bp == 'BP2':
        s = Scales(0.5, 0, 0, 0.125, 0.25, 0.125)
    elif bp == 'BP3':
        s = Scales(0.34, 0, 0, 0.33, 0, 0.33)
    elif bp == 'BP4':
        s = Scales(1./6., 1./6., 1./6., 1./6., 1./6., 1./6.)
    else:
        print 'Unknown branching point: %s' %bp
    sf = getattr(s,'scale_%s'%analysis)
    chanCuts = []
    chanScales = []
    for c in chanMap[analysis]['names']:
        thisScale = sf(c[:2],c[3:])
        if thisScale==0: continue
        chanCut = '('+'||'.join(['channel=="%s"'%x for x in chanMap[analysis][c]])+')'
        chanCuts += [chanCut]
        chanScales += [thisScale]

    fullCut = 'finalstate.mass>100&finalstate.sT>1.1*%f+60. & fabs(z1.mass-%f)>80. & h1.dPhi<%f/600.+1.95' %(mass,ZMASS,mass)
    finalSRCut = 'h1.mass>0.9*%f & h1.mass<1.1*%f' %(mass,mass)

    totBG = 0
    totBG_scaled = 0
    for c,s in zip(chanCuts,chanScales):
        chanBG = s*plotter.getNumEntries('%s&%s&%s' %(c,fullCut,finalSRCut),plotter.signal[0],scaleup=False)
        chanBG_scaled = s*plotter.getNumEntries('%s&%s&%s' %(c,fullCut,finalSRCut),plotter.signal[0],scaleup=True)
        totBG += chanBG
        totBG_scaled += chanBG_scaled
    sigSelSys = (totBG_scaled-totBG)/totBG

    return sigSelSys+1


def add_systematics_sideband(limits,mass,signal,name,chans,sigscale,period):
    limits.add_group("hpp%i" % mass, signal, scale=sigscale, isSignal=True)
    limits.add_group("bg", "bg")
    limits.add_group("data", "data_R*", isData=True)

    lumi = {'hpp%i' % mass: 1.026}
    limits.add_systematics("lumi", "lnN", **lumi)

    idSys = calculateLeptonSystematic(mass,name)

    allid = {'hpp%i' % mass: "%0.3f" %idSys}
    limits.add_systematics("id", "lnN", **allid)

    sigmc = {'hpp%i' % mass: 1.15}
    limits.add_systematics("sigmc", "lnN", **sigmc)

    alpha_pdf = {'bg': 1.05}
    limits.add_systematics("alpha_pdf", "lnN", **alpha_pdf)

    limits.gen_card("%s.txt" % name,mass=mass,cuts=chans)

def add_systematics_fakerate(limits,mass,signal,name,chans,sigscale,period):
    limits.add_group("hpp%i" % mass, signal, scale=sigscale, isSignal=True)
    limits.add_group("data", "data_R*", isData=True)

    lumi = {'hpp%i' % mass: 1.026,}
    limits.add_systematics("lumi", "lnN", **lumi)

    muid = {'hpp%i' % mass: 1.005,}
    limits.add_systematics("muid", "lnN", **muid)

    muiso = {'hpp%i' % mass: 1.002,}
    limits.add_systematics("muiso", "lnN", **muiso)

    sigmc = {'hpp%i' % mass: 1.15}
    limits.add_systematics("sigmc", "lnN", **sigmc)

    limits.gen_card("%s.txt" % name,mass=mass,cuts=chans)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[8, 13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500,help='Mass for signal')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses for signal')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','BP1','BP2','BP3','BP4'],help='Choose branching point for H++')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points for H++')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='mc',default='mc',choices=['mc','sideband','fakerate'],help='Choose BG estimation')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    branchingPoints = ['ee100','em100','mm100','BP1','BP2','BP3','BP4']
    masses = _3L_MASSES if args.region=='Hpp3l' else _4L_MASSES

    if args.period == 7:
        print "7 TeV not implemented"
    elif args.allMasses and args.allBranchingPoints:
        for mass in masses:
            for bp in branchingPoints:
                BP(args.region,args.period,mass,bp,mode=args.bgMode,scalefactor=args.scaleFactor)
    elif args.allMasses:
        for mass in masses:
            BP(args.region,args.period,mass,args.branchingPoint,mode=args.bgMode,scalefactor=args.scaleFactor)
    elif args.allBranchingPoints:
        for bp in branchingPoints:
            BP(args.region,args.period,args.mass,bp,mode=args.bgMode,scalefactor=args.scaleFactor)
    else:
        BP(args.region,args.period,args.mass,args.branchingPoint,mode=args.bgMode,scalefactor=args.scaleFactor)

    return 0


if __name__ == "__main__":
    main()
