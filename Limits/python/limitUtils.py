import itertools
import numpy as np
import logging
import subprocess
import glob
from multiprocessing import Pool

from InitialStateAnalysis.Plotters.plotUtils import getSigMap, getIntLumiMap, getChannels, getMergeDict, ZMASS, getChannelBackgrounds
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Limits.WZLimits import WZLimits
from InitialStateAnalysis.Utilities.utilities import *


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

def getScales(bp):
    if bp == 'ee100':
        s = Scales(1., 0., 0., 0., 0., 0.)
    elif bp == 'em100':
        s = Scales(0., 1., 0., 0., 0., 0.)
    elif bp == 'et100':
        s = Scales(0., 0., 1., 0., 0., 0.)
    elif bp == 'mm100':
        s = Scales(0., 0., 0., 1., 0., 0.)
    elif bp == 'mt100':
        s = Scales(0., 0., 0., 0., 1., 0.)
    elif bp == 'tt100':
        s = Scales(0., 0., 0., 0., 0., 1.)
    elif bp == 'BP1':
        s = Scales(0, 0.01, 0.01, 0.3, 0.38, 0.3)
    elif bp == 'BP2':
        s = Scales(1./2., 0, 0, 1./8., 1./4., 1./8.)
    elif bp == 'BP3':
        s = Scales(1./3., 0, 0, 1./3., 0, 1./3.)
    elif bp == 'BP4':
        s = Scales(1./6., 1./6., 1./6., 1./6., 1./6., 1./6.)
    else:
        logging.error('Unknown branching point: %s' %bp)
        s = Scales(0., 0., 0., 0., 0., 0.)
    return s

def getAllowedHiggsChannels(bp,runTau):
    if bp == 'ee100':
        higgsChannels = ['ee']
    elif bp == 'em100':
        higgsChannels = ['em']
    elif bp == 'et100':
        higgsChannels = ['et']
    elif bp == 'mm100':
        higgsChannels = ['mm']
    elif bp == 'mt100':
        higgsChannels = ['mt']
    elif bp == 'tt100':
        higgsChannels = ['tt']
    elif bp == 'BP1':
        higgsChannels = ['em','et','mm','mt','tt']
    elif bp == 'BP2':
        higgsChannels = ['ee','mm','mt','tt']
    elif bp == 'BP3':
        higgsChannels = ['ee','mm','tt']
    elif bp == 'BP4':
        higgsChannels = ['ee','em','et','mm','mt','tt']
    else:
        logging.error('Unknown branching point: %s' %bp)
        higgsChannels = []
    if not runTau: higgsChannels = [x for x in higgsChannels if 't' not in x]
    return higgsChannels

def unorder(chars):
    theSet = set()
    if len(chars) == 2:
        theSet.add(chars[0]+chars[1])
        theSet.add(chars[1]+chars[0])
    else:
        theSet.add(chars)
    return theSet

def recoFromGen(genSet):
    recoSet = set()
    for gen in genSet:
        if len(gen)>2:
            logging.error('Only 2 flavors allowed in H++ decay, received %s' % gen)
        allowedRecoStates = []
        for c in gen:
            if c=='t':
                allowedLeps = ['e','m']
            else:
                allowedLeps = [c]
            allowedRecoStates += [allowedLeps]
        if len(gen)==1: [recoSet.add(x) for x in allowedRecoStates[0]]
        if len(gen)==2: [recoSet.add(x+y) for x in allowedRecoStates[0] for y in allowedRecoStates[1]]
    return recoSet

def getChannelMap(bp,genLeps,recoLeps,**kwargs):
    runTau = kwargs.pop('runTau',True)
    logging.debug('Getting channel map for %s %i leptons' % (bp,genLeps))
    higgsChannels = getAllowedHiggsChannels(bp,runTau)
    logging.debug('Higgs channels: %s' % str(higgsChannels))

    # gets the gen channels allowed for a bp
    allowedGenChannels = set()
    genChannels = {}
    if genLeps==4:
        for hpp in higgsChannels:
            for hmm in higgsChannels:
                thisName = ''.join(sorted([hpp,hmm]))
                allowedGenChannels.add(hpp+hmm)
                if thisName not in genChannels: genChannels[thisName] = []
                genChannels[thisName] += [hpp+hmm]
    else:
        for hpp in higgsChannels:
            for hm in ['e','m','t']: # always okay to have tau here
                allowedGenChannels.add(hpp+hm)
                if hpp not in genChannels: genChannels[hpp] = []
                genChannels[hpp] += [hpp+hm]
    logging.debug('Gen channels map: %s' % str(genChannels))

    # for each allowed gen channel, gets the associated reco channels
    allowedRecoChannels = set()
    recoChannels = {}
    for namedChannel in genChannels:
        recoChannels[namedChannel] = []
        for genChannel in genChannels[namedChannel]:
            h1 = genChannel[:2]
            h2 = genChannel[2:]
            h1GenAllowedOrdered = unorder(h1)
            h2GenAllowedOrdered = unorder(h2)
            h1RecoAllowed = recoFromGen(h1GenAllowedOrdered)
            h2RecoAllowed = recoFromGen(h2GenAllowedOrdered)
            for reco in itertools.product(h1RecoAllowed,h2RecoAllowed):
                if genLeps==4 and recoLeps==3:
                    # special case for processing h++h-- in h++h- framework
                    thisRecoChan = ''.join([reco[0],reco[1][0]])
                else:
                    thisRecoChan = ''.join(reco)
                allowedRecoChannels.add(thisRecoChan)
                recoChannels[namedChannel] += [thisRecoChan]
    logging.debug('Reco channels map: %s' % str(recoChannels))
    logging.debug('All reco channels to process: %s' % str(allowedRecoChannels))

    channelMap = {
        'names'   : higgsChannels,
        'genmap'  : genChannels,
        'recomap' : recoChannels,
        'allreco' : list(allowedRecoChannels),
    }

    return channelMap

# WZ limit utils
def wzlimit(analysis,region,period,chan,**kwargs):
    cut = kwargs.pop('cut','1')
    name = kwargs.pop('name','card')
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    datacardDir = kwargs.pop('datacardDir','./datacards')
    logging.info("Processing card name {0}".format(name))

    chanCut = '{0} && channel=="{1}"'.format(cut,chan)

    limits = WZLimits(analysis,region, period, chanCut, './ntuples/%s_%iTeV_%s' % (analysis, period, region),
                    '%s/%s_%itev_%s' % (datacardDir, analysis, period, region), scalefactor=scalefactor)

    # add systematic
    bgnames = ['datadriven','ZZ','WW','TTV','VVV','WG']
    signames = ['WZ']
    mcnames = ['WZ','ZZ','WW','TTV','VVV','WG']

    # lumi
    # current recommendation: 12%
    lumi = {}
    for b in mcnames: lumi[b] = 1.12
    limits.add_systematics("lumi", "lnN", **lumi)

    # datadriven
    # assume 40%
    fake = {'datadriven' : 1.4}
    limits.add_systematics('fake_rate_unc','lnN',**fake)

    # eff uncertainty
    # take scale factors for leptons, propagate up and down based on statistical uncertainty
    lepvals = {
        'eee' : 1.018,
        'eem' : 1.013,
        'mme' : 1.006,
        'mmm' : 1.002,
    }
    lep = {}
    for m in mcnames: lep[m] = lepvals[chan]
    limits.add_systematics('lep_eff_unc','lnN',**lep)

    # pu
    # assume 10% uncertainty on min bias cross section, scale up and down, take largest difference in wz yield
    puvals = {
        'eee' : 1.0042,
        'eem' : 1.0015,
        'mme' : 1.0039,
        'mmm' : 1.0024,
    }
    pu = {}
    for m in mcnames: pu[m] = puvals[chan]
    limits.add_systematics('PU_unc','lnN',**pu)

    # met
    # scale all components up and down independently, add in quadrature the largest
    metvals = {
        'eee' : 1.0146, # placeholder from 8 tev
        'eem' : 1.0150,
        'mme' : 1.0159,
        'mmm' : 1.0117,
    }
    met = {}
    for m in mcnames: met[m] = metvals[chan]
    limits.add_systematics('met_unc','lnN',**met)

    # pdf
    # propagate pdf ucnertainties through the selection, scale up and down, take largest
    #pdfvals = {
    #    'eee' : 1.01407, # for now just taking gen, figure it out later after new fsa
    #    'eem' : 1.01394,
    #    'mme' : 1.01399,
    #    'mmm' : 1.01395,
    #}
    #pdf = {}
    #for s in signames: pdf[s] = pdfvals[chan]
    #limits.add_systematics('pdf_unc','lnN',**pdf)

    # scale
    # propagate scale uncertainties through the selection, scale up and down, take largest
    #scalevals = {
    #    'eee' : 1.04296, # again, now just taking gen, fix fsa later
    #    'eem' : 1.04298,
    #    'mme' : 1.04285,
    #    'mmm' : 1.04298,
    #}
    #scale = {}
    #for s in signames: scale[s] = scalevals[chan]
    #limits.add_systematics('scale_unc','lnN',**scale)

    # gen card
    limits.gen_card("{0}.txt".format(name))

def wzLimitWrapper(args):
    analysis = args[0]
    region = args[1]
    period = args[2]
    chan = args[3]
    name = args[4]
    cut = args[5]
    scalefactor = args[6]
    datacardDir = args[7]
    wzlimit(analysis,region,period,chan,name=name,cut=cut,scalefactor=scalefactor,datacardDir=datacardDir)
    

def wzlimits(analysis,region,period,**kwargs):
    cut = kwargs.pop('cut','1')
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    datacardDir = kwargs.pop('datacardDir','./datacards')

    poolArgs = []
    for chan in ['eee','eem','mme','mmm']:
        poolArgs += [(analysis,region,period,chan,chan,cut,scalefactor,datacardDir)]

    if len(poolArgs)==1:
        job = poolArgs[0]
        wzLimitWrapper(job)
    else:
        p = Pool(8)
        try:
            p.map_async(wzLimitWrapper, poolArgs).get(999999)
        except KeyboardInterrupt:
            p.terminate()
            print 'limits cancelled'
            sys.exit(1)

    return 0


def getSignalStrength(cut,**kwargs):
    '''Get WZ signal strength'''
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    datacardDir = kwargs.pop('datacardDir','datacards_temp')
    # produce limits
    wzlimits('WZ','WZ',13,cut=cut,scalefactor=scalefactor,datacardDir=datacardDir)
    # run combine
    combineDir = '/cms/dntaylor/HIGGSCOMBINE_71X/CMSSW_7_1_5/src'
    sigStrengths = {}
    for chan in ['eee','eem','mme','mmm']:
        # cp card
        command = 'cp {0}/WZ_13tev_WZ/{1}.txt {2};'.format(datacardDir,chan,combineDir)
        # run tool, grepping for signal strength
        command += 'pushd {0}; eval `scramv1 runtime -sh`;'.format(combineDir)
        command += 'combine -M MaxLikelihoodFit {0}.txt 2>&1'.format(chan)
        outString = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]
        sigString = [x for x in outString.split('\n') if 'Best' in x][0]
        sigStrength = sigString.split()[3]
        sigErrors = sigString.split()[4].split('/')
        sigStrengths[chan] = [float(sigStrength),float(sigErrors[1]),float(sigErrors[0])]
    return sigStrengths
