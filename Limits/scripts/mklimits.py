#!/usr/bin/env python
import logging
import sys
import argparse
import numpy as np
import itertools

from InitialStateAnalysis.Plotters.plotUtils import _3L_MASSES, _4L_MASSES, getSigMap, getIntLumiMap, getChannels, getMergeDict, ZMASS, getChannelBackgrounds
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Limits.Limits import Limits
from InitialStateAnalysis.Limits.limitUtils import *
from multiprocessing import Pool


def limit(analysis,region,period,mass,**kwargs):
    doChannels = kwargs.pop('doChannels',False)
    doAlphaTest = kwargs.pop('doAlphaTest',False)
    unblind = kwargs.pop('unblind',False)
    name = kwargs.pop('name','card')
    bp = kwargs.pop('bp','')
    directory = kwargs.pop('directory','')
    channelScales = kwargs.pop('channelScales',[1.0])
    channelCuts = kwargs.pop('channelCuts',['1'])
    recoChannels = kwargs.pop('recoChannels',['1'])
    genChannels = kwargs.pop('genChannels',['1'])
    mode = kwargs.pop('mode','sideband')
    scalefactor = kwargs.pop('scalefactor','event.pu_weight*event.lep_scale*event.trig_scale')
    datacardDir = kwargs.pop('datacardDir','./datacards')
    do4l = kwargs.pop('do4l',False)
    logging.info("Processing BP %s; mass-point %i; card name %s" % (bp,mass,name))

    # get the cut maps for each gen channel
    genCutMap = {}
    alphaGenCuts = {}
    for genChan in genChannels:
        genCutMap[genChan] = getChannelCutFlowMap(analysis,genChan,mass=mass)
        alphaGenCuts[genChan] = getChannelSidebandSignalRegion(analysis,genChan,mass=mass)

    # and the reco channel
    recoCutMap = {}
    alphaRecoCuts = {}
    for recoChan in recoChannels:
        recoCutMap[recoChan] = getChannelCutFlowMap(analysis,recoChan,mass=mass)
        alphaRecoCuts[recoChan] = getChannelSidebandSignalRegion(analysis,recoChan,mass=mass)

    # if we are testing a tau hypothesis, explicitly use a wider mass window and the tau cuts
    # otherwise we need to just use reco
    if bp in ['et100','mt100','tt100']:
        theCutMap = getChannelCutFlowMap(analysis,bp[:2]+bp[:2],mass=mass)
        theAlphaCuts = getChannelSidebandSignalRegion(analysis,bp[:2]+bp[:2],mass=mass)
    else:
        theCutMap = recoCutMap[recoChannels[0]]       # only support 1 reco channel at a time
        theAlphaCuts = alphaRecoCuts[recoChannels[0]]

    fullCut = ' && '.join(theCutMap['cuts'])

    # setup sideband stuff
    chanCuts = '(' + ' || '.join(channelCuts) + ')' # gen cuts for h++ ORed with the default for BG 'aaa'

    sbCut = theAlphaCuts['sbcut']
    srCut = theAlphaCuts['srcut']
    finalSRCut = theAlphaCuts['srcut']

    channels, leptons = getChannels(3 if analysis=='Hpp3l' or analysis=='WZ' else 4)

    nl = 3 if analysis in ['Hpp3l', 'WZ'] else 4
    sigMap = getSigMap(4,mass) if do4l else getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()


    # setup for final selection
    myCut = '1'

    sbcut = '%s && %s && %s' %(myCut,sbCut,chanCuts)
    srcut = '%s && %s && %s && %s' %(myCut,srCut,fullCut, chanCuts)
    base_selections = '%s && %s && %s' %(myCut, srCut, fullCut)
    if region=='WZ': srcut = '%s && %s && %s' %(myCut,srCut, chanCuts)
    if doAlphaTest:
        sbcut = '%s && finalstate.sT<150. && z1.mass<110. && h1.mass<130.' %(chanCuts)
        srcut = '%s && finalstate.sT<400. && finalstate.sT>150. && z1.mass<110. && h1.mass<130.' %(chanCuts)
    if doAlphaTest: unblind = True

    logging.debug('Sideband cut: %s' % sbcut)
    logging.debug('Signal region cut: %s' % srcut)
    logging.debug('Base cut: %s' % base_selections)

    # create the limits object
    # base_selections is the cut applied on top of individual channel cuts
    # sbcut is the sideband selection for alpha calculation
    # srcut is the signal ragion selection for alpha calculation
    limits = Limits(analysis,region, period, base_selections, './ntuples/%s_%iTeV_%s' % (analysis, period, region),
                    '%s/%s_%itev_%s/%s/%s' % (datacardDir, analysis, period, region, directory, mass),
                    channels=['dblh%s' % analysis], lumi=intLumiMap[period],
                    blinded=not unblind, bgMode=mode, scalefactor=scalefactor,
                    sbcut=sbcut, srcut=srcut)

    signal =  sigMap[period]['Sig']
    if mode=='sideband':
        # add groups, signal scales must be list of floats with same length as cuts list in gen card
        limits.add_group("hpp%i" % mass, signal, scale=channelScales, isSignal=True)
        limits.add_group("bg", "bg")
        limits.add_group("data", "data_R*", isData=True)

        # luminosity systematics, 2.6 % for mc
        lumi = {'hpp%i' % mass: 1.026}
        limits.add_systematics("lumi", "lnN", **lumi)
    
        # lepton systematics, electron and muon separately for mc
        chanNames = recoChannels # only one supported for now
        scaleMap = calculateChannelLeptonSystematic(mass,chanNames,do4l=do4l)
        eid = {'hpp%i' % mass: "%0.3f" %scaleMap[chanNames[0]]['e']}
        mid = {'hpp%i' % mass: "%0.3f" %scaleMap[chanNames[0]]['m']}
        limits.add_systematics('eid', 'lnN', **eid)
        limits.add_systematics('mid', 'lnN', **mid)
    
        # signal mc uncertainty
        sigmc = {'hpp%i' % mass: 1.15}
        limits.add_systematics("sig_mc_err", "lnN", **sigmc)
    
        # uncertainty on bg estimation
        alpha_str = recoChannels[0] # only one supported for now
        alpha_pdf = {'bg': 1.1}
        limits.add_systematics("alpha_%s" %alpha_str, "lnN", **alpha_pdf)
    
        # generate the card, passing the cuts to be applied for each gen channel
        limits.gen_card("%s.txt" % name, mass=mass, cuts=channelCuts, doAlphaTest=doAlphaTest)

    elif mode=='mc':
        logging.error('MC needs to be reimplemented')
        #add_systematics_mc(limits,mass,signal,name,channelCuts,scale,period,bp,doAlphaTest,do4l)
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
    do4l = args[9]
    BP(analysis,region,period,mass,bp,mode=bgMode,scalefactor=scaleFactor,doAlphaTest=doAlphaTest,unblind=unblind,do4l=do4l)

def BP(analysis,region,period,mass,bp,**kwargs):
    do4l = kwargs.pop('do4l',False)
    s = getScales(bp)
    genLeps = 4 if analysis=='Hpp4l' or do4l else 3
    recoLeps = 4 if analysis=='Hpp4l' else 3
    channelMap = getChannelMap(bp, genLeps, recoLeps)
    sf = getattr(s,'scale_%s'%analysis)
    if do4l: sf = getattr(s,'scale_Hpp4l')

    higgsChannels = channelMap['names']
    genChannelsMap = channelMap['genmap']
    recoChannelsMap = channelMap['recomap']
    allRecoChannels = channelMap['allreco']

    for r in allRecoChannels: # one card per reco channel
        recoCut = 'channel=="%s"' % r
        cardName = '%s_%s' % (bp, r)
        if do4l: cardName += '_4l'
        genChannels = []
        genCuts = []
        genScales = []
        logging.debug('Adding reco channel %s' %r)
        for h in genChannelsMap:
            if r not in recoChannelsMap[h]: continue
            for g in genChannelsMap[h]:
                scale = sf(g[:2],g[2:])
                if not scale: continue
                logging.debug('Adding gen channel %s with scale %f' %(g,scale))
                genChannels += [g]
                genCuts += ['(%s && genChannel=="%s")' % (recoCut,g)]
                genScales += [scale]
        genChannels += ['aaa']
        genCuts += ['(%s && genChannel=="aaa")' % (recoCut)]
        genScales += [1.]
        limit(analysis,region,period,mass,bp=bp,name=cardName,directory=bp,channelCuts=genCuts,channelScales=genScales,genChannels=genChannels,recoChannels=[r],do4l=do4l,**kwargs)

def add_systematics_mc(limits,mass,signal,name,chans,sigscale,period,bp,doAlphaTest,doIndividualChannel,do4l):
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
        scaleMap = calculateChannelLeptonSystematic(mass,chanNames,do4l=do4l)
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

        idSys = "%0.3f" %calculateLeptonSystematic(mass,chans,sigscale,do4l=do4l)

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

def calculateLeptonSystematic(mass,chanCuts,chanScales,**kwargs):
    do4l = kwargs.pop('do4l',False)
    analysis = 'Hpp3l'
    region = 'Hpp3l'
    runPeriod = 8
    nl = 3
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,runPeriod,region)
    saves = '%s_%s_%sTeV' % (analysis,region,runPeriod)
    sigMap = getSigMap(4,mass) if do4l else getSigMap(nl,mass)
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

    fullCut = 'finalstate.mass>100 && finalstate.sT>1.1*%f+60. && fabs(z1.mass-%f)>80. && h1.dPhi<%f/600.+1.95' %(mass,ZMASS,mass)
    finalSRCut = 'h1.mass>0.9*%f && h1.mass<1.1*%f' %(mass,mass)

    totBG = 0
    totBG_scaled = 0
    for c,s in zip(chanCuts,chanScales):
        chanBG = s*plotter.getNumEntries('%s && %s && %s' %(c,fullCut,finalSRCut),plotter.signal[0],scaleup=False)
        chanBG_scaled = s*plotter.getNumEntries('%s && %s && %s' %(c,fullCut,finalSRCut),plotter.signal[0],scaleup=True)
        totBG += chanBG
        totBG_scaled += chanBG_scaled
    sigSelSys = (totBG_scaled-totBG)/totBG

    return sigSelSys+1

def calculateChannelLeptonSystematic(mass,chans,**kwargs):
    do4l = kwargs.pop('do4l',False)
    analysis = 'Hpp3l'
    region = 'Hpp3l'
    runPeriod = 8
    nl = 3
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,runPeriod,region)
    saves = '%s_%s_%sTeV' % (analysis,region,runPeriod)
    sigMap = getSigMap(4,mass) if do4l else getSigMap(nl,mass)
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

    fullCut = 'finalstate.mass>100 && finalstate.sT>1.1*%f+60. && fabs(z1.mass-%f)>80. && h1.dPhi<%f/600.+1.95' %(mass,ZMASS,mass)
    finalSRCut = 'h1.mass>0.9*%f && h1.mass<1.1*%f' %(mass,mass)

    scaleStrings = {
        0: 'h1.LepScaleTight1',
        1: 'h1.LepScaleTight2',
        2: 'h2.LepScaleTight1',
        3: 'h2.LepScaleTight2',
    }

    scaleMap = {}
    for c in chans:
        scaleMap[c] = {}
        for l in ['e','m']:
            # get bg
            scaleFactorBase = "event.trig_scale*event.pu_weight"
            individualScales = "*".join([scaleStrings[x[0]] for x in enumerate(c)])
            theScale = "*".join([scaleFactorBase]+[individualScales])
            plotter.setScaleFactor(theScale)
            chanBG = plotter.getNumEntries('channel=="%s" && %s && %s' %(c,fullCut,finalSRCut),plotter.signal[0])
            individualScales = "*".join([scaleStrings[x[0]] for x in enumerate(c) if x[1]!=l])
            individualScales_up = "*".join([scaleStrings[x[0]]+'_up' for x in enumerate(c) if x[1]==l])
            theScale = scaleFactorBase
            if individualScales: theScale += '*'+individualScales
            if individualScales_up: theScale += '*'+individualScales_up
            plotter.setScaleFactor(theScale)
            chanBG_scaled = plotter.getNumEntries('channel=="%s" && %s && %s' %(c,fullCut,finalSRCut),plotter.signal[0])
            scaleMap[c][l] = (chanBG_scaled-chanBG)/chanBG + 1

    return scaleMap



def add_systematics_sideband(limits,mass,signal,name,chans,sigscale,period,bp,doAlphaTest,do4l):
    limits.add_group("hpp%i" % mass, signal, scale=sigscale, isSignal=True)
    limits.add_group("bg", "bg")
    limits.add_group("data", "data_R*", isData=True)

    lumi = {'hpp%i' % mass: 1.026}
    limits.add_systematics("lumi", "lnN", **lumi)

    if doIndividualChannel:
        chanNames = [name[-3:]]
        if do4l: chanNames = [name[-6:-3]]
        scaleMap = calculateChannelLeptonSystematic(mass,chanNames,do4l=do4l)
        eid = {'hpp%i' % mass: "%0.3f" %scaleMap[chanNames[0]]['e']}
        mid = {'hpp%i' % mass: "%0.3f" %scaleMap[chanNames[0]]['m']}
        limits.add_systematics('eid', 'lnN', **eid)
        limits.add_systematics('mid', 'lnN', **mid)
    else:
        idSys = calculateLeptonSystematic(mass,chans,sigscale,do4l=do4l)

        allid = {'hpp%i' % mass: "%0.3f" %idSys}
        limits.add_systematics("id", "lnN", **allid)

    sigmc = {'hpp%i' % mass: 1.15}
    limits.add_systematics("sig_mc_err", "lnN", **sigmc)

    if doIndividualChannel:
        alpha_str = name[-3:]
        if do4l: alpha_str = name[-6:-3]
    else:
        alpha_str = 'pdf'

    alpha_pdf = {'bg': 1.1}
    limits.add_systematics("alpha_%s" %alpha_str, "lnN", **alpha_pdf)

    limits.gen_card("%s.txt" % name, mass=mass, cuts=chans, doAlphaTest=doAlphaTest)



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','WZ'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[8, 13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500,help='Mass for signal')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses for signal')
    parser.add_argument('-da','--doAlphaTest',action='store_true',help='Run the alpha test')
    parser.add_argument('-df','--do4l', action='store_true',help='Run the 4l lepton limits')
    parser.add_argument('-ub','--unblind',action='store_true',help='unblind')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='BP4',default='BP4',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4'],help='Choose branching point for H++')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points for H++')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='sideband',default='sideband',choices=['mc','sideband'],help='Choose BG estimation')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    branchingPoints = ['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4']
    masses = _3L_MASSES if args.analysis=='Hpp3l' else _4L_MASSES
    if args.do4l: masses = _4L_MASSES

    if not args.allMasses: masses = [args.mass]
    if not args.allBranchingPoints: branchingPoints = [args.branchingPoint]

    poolArgs = [[m,b] for m in masses for b in branchingPoints]

    if len(poolArgs)==1:
        job = poolArgs[0]
        BPWrapper((args.analysis,args.region,args.period,job[0],job[1],args.bgMode,args.scaleFactor,args.doAlphaTest,args.unblind,args.do4l))
    else:
        p = Pool(8)
        try:
            p.map_async(BPWrapper, [(args.analysis,args.region,args.period,job[0],job[1],args.bgMode,args.scaleFactor,args.doAlphaTest,args.unblind,args.do4l) for job in poolArgs]).get(999999)
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
