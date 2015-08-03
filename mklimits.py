#!/usr/bin/env python
from limits.Limits import Limits
import logging
import sys
import argparse
import numpy as np
from plotters.plotUtils import _3L_MASSES, _4L_MASSES, getSigMap, getIntLumiMap, getChannels, getMergeDict, ZMASS
from plotters.Plotter import Plotter
from multiprocessing import Pool

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
        scale = 9./2
        if hpp in ['ee','mm','tt']: scale = 9.
        return self.a_3l[i] * scale
        

def limit(analysis,region,period,mass,**kwargs):
    doChannels = kwargs.pop('doChannels',False)
    doIndividualChannel = kwargs.pop('doIndividualChannel',False)
    doAlphaTest = kwargs.pop('doAlphaTest',False)
    unblind = kwargs.pop('unblind',False)
    name = kwargs.pop('name','card')
    bp = kwargs.pop('bp','')
    scale = kwargs.pop('scale',[1.0])
    directory = kwargs.pop('directory','')
    chans = kwargs.pop('channels',['1'])
    mode = kwargs.pop('mode','sideband')
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
        },
        'AlphaTest' : {
             'cuts' : ['1',\
                       'finalstate.sT<400. & finalstate.sT>150.',\
                       'z1.mass<110.'],
             'labels' : ['Preselection','s_{T}','Z']
        },
    }

    nl = 3 if analysis in ['Hpp3l', 'WZ'] else 4
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()

    cuts = '&&'.join(cutMap[analysis]['cuts'])
    if doAlphaTest: cuts = ' & '.join(cutMap['AlphaTest']['cuts'])

    if doAlphaTest: unblind = True
    limits = Limits(analysis,region, period, cuts, './ntuples/%s_%iTeV_%s' % (analysis, period, region),
                    './datacards/%s_%itev_%s/%s/%s' % (analysis, period, region, directory, mass),
                    channels=['dblh%s' % analysis], lumi=intLumiMap[period],
                    blinded=not unblind, bgMode=mode, scalefactor=scalefactor)

    signal =  sigMap[period]['Sig']
    if mode=='mc':
        add_systematics_mc(limits,mass,signal,name,chans,scale,period,bp,doAlphaTest,doIndividualChannel)
    elif mode=='sideband':
        add_systematics_sideband(limits,mass,signal,name,chans,scale,period,bp,doAlphaTest,doIndividualChannel)
    else:
        return 0

def BPWrapper(args):
    analysis = args[0]
    region = args[1]
    period = args[2]
    mass = args[3]
    bp = args[4]
    bgMode = args[5]
    scaleFactor = args[6]
    doAlphaTest = args[7]
    unblind = args[8]
    BP(analysis,region,period,mass,bp,mode=bgMode,scalefactor=scaleFactor,doAlphaTest=doAlphaTest,unblind=unblind)

def BP(analysis,region,period,mass,bp,**kwargs):
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
             'eegen': ['eee','eem','eet'],
             'emgen': ['eme','emm','emt'],
             'mmgen': ['mme','mmm','mmt'],
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
        chanCut = '('+' | '.join(['channel=="%s"'%x for x in chanMap[analysis][c]])+')'
        chanCut += ' & (' + ' | '.join(['genChannel=="%s"'%x for x in chanMap[analysis][c+'gen']+['aaa']]) + ')'
        limit(analysis,region,period,mass,bp=bp,name='%s_%s'%(bp,c),directory=bp,channels=[chanCut],scale=[thisScale],**kwargs)
        # now do individual channels
        for theChan in chanMap[analysis][c]:
            thisChanCut = 'channel=="%s"' % theChan
            thisChanCut += ' & (' + ' | '.join(['genChannel=="%s"'%x for x in chanMap[analysis][c+'gen']+['aaa']]) + ')'
            limit(analysis,region,period,mass,bp=bp,name='%s_%s'%(bp,theChan),directory=bp,channels=[thisChanCut],scale=[thisScale],doIndividualChannel=True,**kwargs)
        chanCuts += [chanCut]
        chanScales += [thisScale]
    limit(analysis,region,period,mass,bp=bp,name=bp,directory=bp,channels=chanCuts,scale=chanScales,**kwargs)

def add_systematics_mc(limits,mass,signal,name,chans,sigscale,period,bp,doAlphaTest,doIndividualChannel):
    limits.add_group("hpp%i" % mass, signal, scale=sigscale, isSignal=True)
    if period==8: limits.add_group("dyjets", "Z*j*")
    if period==13: limits.add_group("dyjets", "DY*")
    limits.add_group("zz", "ZZJ*")
    limits.add_group("wz", "WZJ*")
    if period==8: limits.add_group("ww", "WWJ*")
    if period==8: limits.add_group("zzz", "ZZZ*")
    if period==8: limits.add_group("wzz", "WZZ*")
    if period==8: limits.add_group("wwz", "WWZ*")
    if period==8: limits.add_group("www", "WWW*")
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
        'wzz':          1.026,
        'wwz':          1.026,
        'www':          1.026,
        'tt':           1.026,
        'ttz':          1.026,
        'ttw':          1.026,
        'top':          1.026
    }

    limits.add_systematics("lumi", "lnN", **lumi)

    if doIndividualChannel:
        chanNames = [name[-3:]]
        scaleMap = calculateChannelLeptonSystematic(mass,chanNames)
        esys = "%0.3f" %scaleMap[chanNames[0]]['e']
        eid = {
            'hpp%i' % mass: esys,
            'dyjets':       esys,
            'zz':           esys,
            'wz':           esys,
            'ww':           esys,
            'zzz':          esys,
            'wzz':          esys,
            'wwz':          esys,
            'www':          esys,
            'tt':           esys,
            'ttz':          esys,
            'ttw':          esys,
            'top':          esys
        }
        msys = "%0.3f" %scaleMap[chanNames[0]]['m']
        mid = {
            'hpp%i' % mass: msys,
            'dyjets':       msys,
            'zz':           msys,
            'wz':           msys,
            'ww':           msys,
            'zzz':          msys,
            'wzz':          msys,
            'wwz':          msys,
            'www':          msys,
            'tt':           msys,
            'ttz':          msys,
            'ttw':          msys,
            'top':          msys
        }
        limits.add_systematics('eid', 'lnN', **eid)
        limits.add_systematics('mid', 'lnN', **mid)
    else:

        idSys = "%0.3f" %calculateLeptonSystematic(mass,chans,sigscale)

        lepid = {
            'hpp%i' % mass: idSys,
            'dyjets':       idSys,
            'zz':           idSys,
            'wz':           idSys,
            'ww':           idSys,
            'zzz':          idSys,
            'wzz':          idSys,
            'wwz':          idSys,
            'www':          idSys,
            'tt':           idSys,
            'ttz':          idSys,
            'ttw':          idSys,
            'top':          idSys
        }
        limits.add_systematics("lepid", "lnN", **lepid)

    sigmc = {'hpp%i' % mass: 1.15}
    limits.add_systematics("sigmc", "lnN", **sigmc)

    # uncertainties taken from relevant smp papers
    mcdata = {
        'zz':           1.105,
        'wz':           1.056,
        'ww':           1.041,
        'tt':           1.024
    }
    limits.add_systematics("mcdata", "lnN", **mcdata)

    # all uncertainties taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
    pdf = {
        'zzz':          1.027,
        'wzz':          1.060,
        'wwz':          1.056,
        'www':          1.047,
        'ttz':          1.117,
        'ttw':          1.289
    }
    limits.add_systematics("pdf", "lnN", **pdf)

    limits.gen_card("%s_mc.txt" % name,mass=mass,cuts=chans,doAlphaTest=doAlphaTest)

def calculateLeptonSystematic(mass,chanCuts,chanScales):
    analysis = 'Hpp3l'
    region = 'Hpp3l'
    runPeriod = 8
    nl = 3
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,runPeriod,region)
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

def calculateChannelLeptonSystematic(mass,chans):
    analysis = 'Hpp3l'
    region = 'Hpp3l'
    runPeriod = 8
    nl = 3
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,runPeriod,region)
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

    fullCut = 'finalstate.mass>100&finalstate.sT>1.1*%f+60. & fabs(z1.mass-%f)>80. & h1.dPhi<%f/600.+1.95' %(mass,ZMASS,mass)
    finalSRCut = 'h1.mass>0.9*%f & h1.mass<1.1*%f' %(mass,mass)

    scaleStrings = {
        0: 'h1.LepScaleTight1',
        1: 'h1.LepScaleTight2',
        2: 'h2.LepScaleTight1',
    }

    scaleMap = {}
    #print 'Generating scale map'
    for c in chans:
        #print 'Channel: %s' %c
        scaleMap[c] = {}
        for l in ['e','m']:
            #print 'Scale for %s' %l
            # get bg
            scaleFactorBase = "event.trig_scale*event.pu_weight"
            individualScales = "*".join([scaleStrings[x[0]] for x in enumerate(c)])
            theScale = "*".join([scaleFactorBase]+[individualScales])
            #print 'Scale: %s' % theScale
            plotter.setScaleFactor(theScale)
            chanBG = plotter.getNumEntries('channel=="%s"&%s&%s' %(c,fullCut,finalSRCut),plotter.signal[0])
            individualScales = "*".join([scaleStrings[x[0]] for x in enumerate(c) if x[1]!=l])
            individualScales_up = "*".join([scaleStrings[x[0]]+'_up' for x in enumerate(c) if x[1]==l])
            theScale = scaleFactorBase
            if individualScales: theScale += '*'+individualScales
            if individualScales_up: theScale += '*'+individualScales_up
            #print 'Scaleup: %s' % theScale
            plotter.setScaleFactor(theScale)
            chanBG_scaled = plotter.getNumEntries('channel=="%s"&%s&%s' %(c,fullCut,finalSRCut),plotter.signal[0])
            #print chanBG, chanBG_scaled
            scaleMap[c][l] = (chanBG_scaled-chanBG)/chanBG + 1

    #print scaleMap
    return scaleMap



def add_systematics_sideband(limits,mass,signal,name,chans,sigscale,period,bp,doAlphaTest,doIndividualChannel):
    limits.add_group("hpp%i" % mass, signal, scale=sigscale, isSignal=True)
    limits.add_group("bg", "bg")
    limits.add_group("data", "data_R*", isData=True)

    lumi = {'hpp%i' % mass: 1.026}
    limits.add_systematics("lumi", "lnN", **lumi)

    if doIndividualChannel:
        chanNames = [name[-3:]]
        #print 'Do chans'
        #print chanNames
        scaleMap = calculateChannelLeptonSystematic(mass,chanNames)
        eid = {'hpp%i' % mass: "%0.3f" %scaleMap[chanNames[0]]['e']}
        mid = {'hpp%i' % mass: "%0.3f" %scaleMap[chanNames[0]]['m']}
        limits.add_systematics('eid', 'lnN', **eid)
        limits.add_systematics('mid', 'lnN', **mid)
    else:
        idSys = calculateLeptonSystematic(mass,chans,sigscale)

        allid = {'hpp%i' % mass: "%0.3f" %idSys}
        limits.add_systematics("id", "lnN", **allid)

    sigmc = {'hpp%i' % mass: 1.15}
    limits.add_systematics("sig_mc_err", "lnN", **sigmc)

    alpha_pdf = {'bg': 1.1}
    limits.add_systematics("alpha_pdf", "lnN", **alpha_pdf)

    limits.gen_card("%s.txt" % name, mass=mass, cuts=chans, doAlphaTest=doAlphaTest)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','WZ'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[8, 13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500,help='Mass for signal')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses for signal')
    parser.add_argument('-da','--doAlphaTest',action='store_true',help='Run the alpha test')
    parser.add_argument('-ub','--unblind',action='store_true',help='unblind')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','BP1','BP2','BP3','BP4'],help='Choose branching point for H++')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points for H++')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='mc',default='sideband',choices=['mc','sideband'],help='Choose BG estimation')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    branchingPoints = ['ee100','em100','mm100','BP1','BP2','BP3','BP4']
    masses = _3L_MASSES if args.analysis=='Hpp3l' else _4L_MASSES

    if not args.allMasses: masses = [args.mass]
    if not args.allBranchingPoints: branchingPoints = [args.branchingPoint]

    poolArgs = [[m,b] for m in masses for b in branchingPoints]

    p = Pool(8)
    try:
        p.map_async(BPWrapper, [(args.analysis,args.region,args.period,job[0],job[1],args.bgMode,args.scaleFactor,args.doAlphaTest,args.unblind) for job in poolArgs]).get(999999)
    except KeyboardInterrupt:
        p.terminate()
        print 'limits cancelled'
        sys.exit(1)
    

    #for mass in masses:
    #    for bp in branchingPoints:
    #        BP(args.analysis,args.region,args.period,mass,bp,mode=args.bgMode,scalefactor=args.scaleFactor,doAlphaTest=args.doAlphaTest,unblind=args.unblind)

    return 0


if __name__ == "__main__":
    main()
