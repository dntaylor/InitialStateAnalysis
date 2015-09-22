import itertools
import numpy as np
import logging

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
        s = Scales(0, 0.1, 0.1, 0.3, 0.38, 0.3)
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
        higgsChannels = ['ee','em','mt','tt']
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
