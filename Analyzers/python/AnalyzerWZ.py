#!/usr/bin/env python
'''
The WZ analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerWZ(AnalyzerBase):
    '''
    An implementation of the AnalyzerBase class for use in WZ analysis.

    Objects:
        z: 2 leptons same flavor opposite sign
        w: 1 lepton and met
    Minimization function:
        Z mass difference
    Selection:
        Trigger: double lepton triggers
        Fiducial cut
        Lepton id: tight muon, cbid medium electron + triggering MVA
        Isolation: 0.12 (0.15) for muon (electron) with deltabeta (rho) corrections
        Invariant mass > 100. GeV
        Z selection: within 20. GeV of Z mass, leading lepton pt > 20.
        W selection: met > 30. lepton pt > 20.
    '''

    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        if not hasattr(self,'channel'): self.channel = 'WZ'
        self.period = period
        self.final_states = ['eee','eem','emm','mmm']
        self.initial_states = ['z1','w1'] # in order of leptons returned in choose_objects
        self.object_definitions = {
            'w1': ['em','n'],
            'z1': ['em','em'],
        }
        self.tightW = False
        #self.alternateIds, self.alternateIdMap = self.defineAlternateIds(period)
        #self.doVBF = (period==13)
        #self.doMetUnc = (period==13)

        super(AnalyzerWZ, self).__init__(sample_name, file_list, out_file, period, **kwargs)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select leptons that best fit the WZ selection.
        The first two leptons are the Z and the third is the W.
        Z are then ordered in pt.
        We select combinatorics by closest to zmass.
        '''

        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            OS1 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) < 0.5 # select opposite sign
            mass = self.getObject(rtrow, 'mass', l[0], l[1])
            massdiff = abs(ZMASS-mass)
            pt2 = self.getObject(rtrow,'pt', l[2])

            ordList = [l[1], l[0], l[2]] if self.getObject(rtrow,'pt',l[0]) < self.getObject(rtrow,'pt', l[1]) else [l[0], l[1], l[2]]

            ## make sure they are all separated (no split tracks)
            #o01 = ordered(l[0],l[1])
            #o02 = ordered(l[0],l[2])
            #o12 = ordered(l[1],l[2])
            #dr01 = self.getObject(rtrow,'dr',o01[0],o01[1]) > 0.02
            #dr02 = self.getObject(rtrow,'dr',o02[0],o02[1]) > 0.02
            #dr12 = self.getObject(rtrow,'dr',o12[0],o12[1]) > 0.02
            #dr = dr01 and dr02 and dr12 
            dr = True
      
            veto = self.veto(rtrow)

            if OS1 and l[0][0]==l[1][0] and dr and veto:
                cands.append((massdiff, -pt2, mass, ordList))

        if not len(cands): return ([],[])

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        massdiff, negpt2, mass, leps = cands[0]

        return ([massdiff,negpt2], leps)

    # overide good_to_store
    # will store via veto
    #@staticmethod
    def good_to_store(self,rtrow, cand1, cand2):
        '''
        Veto on 4th lepton
        '''
        # first select best, (closest to Z and highest pt W)
        good = True
        if len(cand2) < 1: # the case that our previous bestcand didnt exist
            good = True
        elif len(cand1) < 1: # case that our previous best cand existed but this one doesnt
            good = False
        else: # case when we do an actual comparison
            for min1, min2 in zip(cand1, cand2):
                if min1 < min2:
                    good = True
                    break
                if min1 > min2:
                    good = False
                    break

        return good


    #def defineAlternateIds(self,period):
    #    if period==8:
    #        return [], {}
    #    elecIds = ['WZLooseTrigIso','Medium','Tight']
    #    muonIds = ['WZMediumTrigIso']
    #    elecIsos = [0.]
    #    muonIsos = [0.,0.4, 0.12]
    #    idList = []
    #    idMap = {}
    #    for id in elecIds:
    #        for iso in elecIsos:
    #            idName = 'elec%s%0.2f' % (id, iso)
    #            idName = idName.replace('.','p')
    #            idList += [idName]
    #            idMap[idName] = {
    #                'idDef' : {
    #                    'e': id
    #                }
    #            }
    #            if iso:
    #                idMap[idName]['isoCut'] = {
    #                    'e': iso
    #                }
    #    for id in muonIds:
    #        for iso in muonIsos:
    #            idName = 'muon%s%0.2f' % (id, iso)
    #            idName = idName.replace('.','p')
    #            idList += [idName]
    #            idMap[idName] = {
    #                'idDef' : {
    #                    'm': id
    #                }
    #            }
    #            if iso:
    #                idMap[idName]['isoCut'] = {
    #                    'm': iso
    #                }
    #    return idList, idMap

    def getGenChannel(self,rtrow):
        '''Dummy return gen channel string'''
        if 'WZJets' not in self.file_name and 'WZTo3LNu' not in self.file_name: return 'a'
        flav = {'e':'E','m':'Mu','t':'Tau'}
        for z in ['e','m','t']:
            for w in ['e','m','t']:
                zll = getattr(rtrow,'GenDecayZ{0}{0}'.format(flav[z]))
                wln = getattr(rtrow,'GenDecayW{0}Nu'.format(flav[w]))
                if zll and wln: return '{0}{0}{1}'.format(z,w)
        return 'a'

    def getFakeChannel(self,rtrow,**kwargs):
        chan = ''
        for i,obj in enumerate(self.objCand):
            if self.tightW and i==2 and obj[0]=='e':
                chan += 'P' if self.ID(rtrow,obj,**self.getIdArgs('VeryTight')) and self.ID(rtrow,obj,**self.getIdArgs('Tight')) else 'F'
            else:
                chan += 'P' if self.ID(rtrow,obj,**self.getIdArgs('Tight')) else 'F'
        return chan
        

    ###########################
    ### Define preselection ###
    ###########################
    def cutTreeSelections(self):
        cutTree = CutTree()
        cutTree.add(self.returnTrue,'topology')
        cutTree.add(self.trigger,'trigger')
        cutTree.add(self.fiducial,'fiducial')
        cutTree.add(self.ID_loose,'looseID')
        cutTree.add(self.ID_tight,'tightID')
        if self.tightW: cutTree.add(self.ID_tightW,'tightWID')
        cutTree.add(self.mass3l,'mass3l')
        cutTree.add(self.zWindow,'zWindow')
        cutTree.add(self.zLeadPt,'zLeadPt')
        cutTree.add(self.wPt,'wPt')
        cutTree.add(self.wMll,'wMll')
        cutTree.add(self.met,'met')
        cutTree.add(self.bjetVeto,'bjetVeto')
        cutTree.add(self.veto,'veto4thLepton')
        return cutTree

    def veto(self,rtrow):
        # first veto on 4th tight lepton
        veto = (rtrow.eVetoTight + rtrow.muVetoMedium == 0) if self.period==13 else\
               (rtrow.elecVetoWZTight + rtrow.muonVetoWZTight == 0)
        if self.tightW:
            veto = (rtrow.eVetoMedium + rtrow.muVetoMedium == 0) if self.period==13 else\
                   (rtrow.elecVetoWZTight + rtrow.muonVetoWZTight == 0)
        if self.metShift and not self.isData: # only shift MC
            vetoMap = {
                'ees+' : (rtrow.eVetoTight_eesUp + rtrow.muVetoMedium == 0),
                'ees-' : (rtrow.eVetoTight_eesDown + rtrow.muVetoMedium == 0),
                'mes+' : (rtrow.eVetoTight + rtrow.muVetoMedium_mesUp == 0),
                'mes-' : (rtrow.eVetoTight + rtrow.muVetoMedium_mesDown == 0),
            }
            if self.tightW:
                vetoMap = {
                    'ees+' : (rtrow.eVetoMedium_eesUp + rtrow.muVetoMedium == 0),
                    'ees-' : (rtrow.eVetoMedium_eesDown + rtrow.muVetoMedium == 0),
                    'mes+' : (rtrow.eVetoMedium + rtrow.muVetoMedium_mesUp == 0),
                    'mes-' : (rtrow.eVetoMedium + rtrow.muVetoMedium_mesDown == 0),
                }
            if self.metShift in vetoMap:
                veto = vetoMap[self.metShift]

        return veto


    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        #if self.period=='13': cuts.add(self.passAnyId)
        #cuts.add(self.ID_tight)
        cuts.add(self.ID_loose)
        cuts.add(self.veto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        if self.tightW:
            cuts.add(self.ID_tightW)
        else:
            cuts.add(self.ID_tight)
        cuts.add(self.mass3l)
        cuts.add(self.zSelection)
        cuts.add(self.wSelection)
        cuts.add(self.bjetVeto)
        cuts.add(self.veto)
        return cuts

    def getIdArgs(self,type):
        kwargs = {}
        if type=='VeryTight':
            kwargs['idDef'] = {
                'e':'WWTight',
                'm':'WWMedium',
            }
            kwargs['isoCut'] = {
                'e':0.15,
                'm':0.12,
            }
            if self.period==8:
                kwargs['idDef']['e'] = 'WZTight'
                kwargs['idDef']['m'] = 'WZTight'
            if self.period==13:
                kwargs['isoCut']['e'] = 0. # baked into ID
                kwargs['isoCut']['m'] = 0. # baked into ID
        if type=='Tight':
            kwargs['idDef'] = {
                'e':'WWTight',
                'm':'WWMedium',
            }
            kwargs['isoCut'] = {
                'e':0.15,
                'm':0.12,
            }
            if self.period==8:
                kwargs['idDef']['e'] = 'WZTight'
                kwargs['idDef']['m'] = 'WZTight'
            if self.period==13:
                if self.tightW:
                    kwargs['idDef']['e'] = 'WWMedium'
                kwargs['isoCut']['e'] = 0. # baked into ID
                kwargs['isoCut']['m'] = 0. # baked into ID
        if type=='Loose':
            kwargs['idDef'] = {
                'e':'WWLoose',
                'm':'WWLoose',
            }
            kwargs['isoCut'] = {
                'e':0.4,
                'm':1.0,
            }
            if self.period==8:
                kwargs['idDef']['e'] = 'WZLoose'
                kwargs['idDef']['m'] = 'WZLoose'
            if self.period==13:
                kwargs['isoCut']['e'] = 0. # baked into ID
                kwargs['isoCut']['m'] = 0. # baked into ID
        if hasattr(self,'alternateIds'):
            if type in self.alternateIds:
                kwargs = self.alternateIdMap[type]
        kwargs['metShift'] = self.metShift
        return kwargs

    def trigger(self, rtrow):
        if self.period == 8:
            triggers = ["mu17ele8isoPass", "mu8ele17isoPass",
                        "doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == 13:
            triggers = ['singleMuSingleEPass', 'doubleMuPass', 'doubleEPass', 'singleESingleMuPass',
                        'singleIsoMu20', 'singleIsoTkMu20', 'singleIsoMu27', 'singleE23WPLoose']

        good = False
        for t in triggers:
            if getattr(rtrow,t)>0:
                good = True
                break
        return good

    def fiducial(self, rtrow):
        for l in self.objects:
            if l[0]=='e':
                ptcut = 10.0
                etacut = 2.5
            if l[0]=='m':
                ptcut = 10.0
                etacut = 2.4
            pt = self.getObject(rtrow, 'pt', l)
            eta = abs(self.getObject(rtrow, 'eta', l))
            if self.getObject(rtrow, 'pt', l) < ptcut:
                return False
            if abs(self.getObject(rtrow, 'eta', l)) > etacut:
                return False
        return True

    def passAnyId(self,rtrow):
        '''Check to make sure the leptons pass at least 1 ID'''
        passCheck = {'e': False, 'm': False}
        for altId in self.alternateIds:
           if passCheck[altId[0]]: continue
           if self.ID(rtrow,*self.objects,**self.getIdArgs(altId)): passCheck[altId[0]] = True
        return passCheck['e'] and passCheck['m']

    def ID_loose(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Loose'))

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))

    def ID_tightW(self, rtrow):
        return self.ID(rtrow,*self.objects[:-1],**self.getIdArgs('Tight')) and self.ID(rtrow,self.objects[-1],**self.getIdArgs('veryTight'))

    def ID_veryTight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('veryTight'))

    def mass3l(self,rtrow):
        return self.getObject(rtrow,'mass',*self.objCand) > 100

    def zSelection(self,rtrow):
        return self.zWindow(rtrow) and self.zLeadPt(rtrow)

    def zWindow(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = self.getObject(rtrow,'mass',o[0],o[1])
        return abs(m1-ZMASS) < 15.

    def zLeadPt(self,rtrow):
        leps = self.objCand
        l0Pt = self.getObject(rtrow,'pt',leps[0])
        return l0Pt>20.

    def wSelection(self,rtrow):
        return self.wPt(rtrow) and self.wMll(rtrow) and self.met(rtrow)

    def wPt(self,rtrow):
        leps = self.objCand
        return self.getObject(rtrow, 'pt', leps[2])>20.

    def wMll(self,rtrow):
        leps = self.objCand
        for l in leps[:2]:
            o = ordered(l,leps[2])
            mll = self.getObject(rtrow, 'mass', o[0],o[1])
            if mll < 4.: return False
        return True

    def met(self,rtrow):
        return self.getObject(rtrow,'pt','met') > 30.

    def bjetVeto(self, rtrow):
        return rtrow.bjetCISVVeto30Tight==0

class AnalyzerWZ_NoVeto(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_NoVeto, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'NoVeto'

    # reoveride good_to_store
    #@staticmethod
    def good_to_store(self,rtrow, cand1, cand2):
        '''
        Veto on 4th lepton
        '''
        # first select best, (closest to Z and highest pt W)
        good = False
        for min1, min2 in zip(cand1, cand2):
            if min1 < min2:
                good = True
                break
            if min1 > min2:
                good = False
                break

        return good


class AnalyzerWZ_ZFakeRate(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_ZFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'FakeRate'

    def cutTreeSelections(self):
        cutTree = CutTree()
        cutTree.add(self.returnTrue,'topology')
        cutTree.add(self.trigger,'trigger')
        cutTree.add(self.fiducial,'fiducial')
        cutTree.add(self.ID_tight_Z,'looseProbe')
        cutTree.add(self.ID_tight,'tightProbe')
        cutTree.add(self.ID_veryTight,'veryTightProbe')
        cutTree.add(self.zSelection,'zSelection')
        cutTree.add(self.veto,'veto4thLepton')
        return cutTree

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.ID_tight_Z)
        cuts.add(self.zSelection)
        #cuts.add(self.metveto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        #cuts.add(self.metveto)
        return cuts

    def ID_tight_Z(self, rtrow):
        return self.ID(rtrow,*self.objCand[:2],**self.getIdArgs('Tight'))

    def metveto(self,rtrow):
        if self.getObject(rtrow,'pt','met') > 25.: return False
        return True

    def trigger(self, rtrow):
        if self.period == 8:
            triggers = ["doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == 13:
            triggers = ['doubleMuPass', 'doubleEPass']

        nottriggers = []
        # for data, to merge datasets
        if False:
            if 'DoubleMuon' in self.file_name:
                nottriggers = [] # allow all triggers
            elif 'DoubleEG' in self.file_name:
                nottriggers = ['doubleMuPass'] # dont allow doublemuon
            elif 'MuonEG' in self.file_name:
                nottriggers = ['doubleMuPass', 'doubleEPass'] # don't allow
            else:
                nottriggers = []

        good = False
        for t in triggers:
            if getattr(rtrow,t)>0:
                good = True
                break
        for t in nottriggers:
            if getattr(rtrow,t)>0:
                good = False
                break
        return good

    def zSelection(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = self.getObject(rtrow,'mass',o[0],o[1])
        l0Pt = self.getObject(rtrow,'pt',leps[0])
        wl = leps[2]
        o0 = ordered(leps[0],wl)
        o1 = ordered(leps[1],wl)
        dr0 = self.getObject(rtrow,'dr',o0[0],o0[1])
        dr1 = self.getObject(rtrow,'dr',o1[0],o1[1])
        return abs(m1-ZMASS)<10. and l0Pt>20. and dr0>0.02 and dr1>0.02

    def good_to_store(self, rtrow, cand1, cand2):
        '''
        Iterate through minimizing variables.
        '''
        good = False
        for min1, min2 in zip(cand1, cand2):
            if min1 < min2:
                good = True
                break
            if min1 > min2:
                good = False
                break
        if self.objCand[0][0]=='e' and self.objCand[1][0]=='e':
            ename = 'MatchesDoubleE' if self.period == 13 else 'MatchesDoubleEPath'
            match_0 = getattr(rtrow,'%s%s' %(self.objCand[0],ename))
            match_1 = getattr(rtrow,'%s%s' %(self.objCand[1],ename))
        elif self.objCand[0][0]=='m' and self.objCand[1][0]=='m':
            if self.period==13:
                match_0 = getattr(rtrow,'%sMatchesDoubleMu' %self.objCand[0])
                match_1 = getattr(rtrow,'%sMatchesDoubleMu' %self.objCand[1])
            else:
                match_0 = getattr(rtrow,'%sMatchesMu17Mu8Path' %self.objCand[0])>0 + getattr(rtrow,'%sMatchesMu17TrkMu8Path' %self.objCand[0])
                match_1 = getattr(rtrow,'%sMatchesMu17Mu8Path' %self.objCand[1])>0 + getattr(rtrow,'%sMatchesMu17TrkMu8Path' %self.objCand[1])
        else:
            match_0 = 0
            match_1 = 0
        passTrig = match_0 > 0.5 and match_1 > 0.5

        #veto = (rtrow.eVetoTrigIso + rtrow.muVetoTightTrigIso == 0) if self.period==13 else (rtrow.elecVetoWZLoose + rtrow.muonVetoWZLoose == 0)
        veto = (rtrow.eVetoTrigIso + rtrow.muVetoMediumTrigIso == 0) if self.period==13 else (rtrow.elecVetoWZLoose + rtrow.muonVetoWZLoose == 0)

        return good and passTrig and veto


class AnalyzerWZ_TTFakeRate(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_TTFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'TTFakeRate'

    def cutTreeSelections(self):
        cutTree = CutTree()
        cutTree.add(self.returnTrue,'topology')
        cutTree.add(self.trigger,'trigger')
        cutTree.add(self.fiducial,'fiducial')
        cutTree.add(self.ID_tight_TT,'looseProbe')
        cutTree.add(self.ID_tight,'tightProbe')
        cutTree.add(self.ID_veryTight,'veryTightProbe')
        cutTree.add(self.zVeto,'zVeto')
        cutTree.add(self.jetCut,'jetCut')
        cutTree.add(self.bjetCut,'bjetCut')
        cutTree.add(self.metCut,'metCut')
        cutTree.add(self.veto,'veto4thLepton')
        return cutTree

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.ID_tight_TT)
        cuts.add(self.zVeto)
        cuts.add(self.jetCut)
        #cuts.add(self.bjetCut)
        #cuts.add(self.metCut)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.zVeto)
        cuts.add(self.jetCut)
        #cuts.add(self.bjetCut)
        #cuts.add(self.metCut)
        return cuts

    def choose_objects(self, rtrow):
        '''
        Select leptons that best fit the WZ selection.
        The first two leptons are the Z and the third is the W.
        Z are then ordered in pt.
        We select combinatorics by closest to zmass.
        '''
        # first veto on 4th tight lepton
        #veto = (rtrow.eVetoTight + rtrow.muVetoTight == 0) if self.period==13 else\
        veto = (rtrow.eVetoTight + rtrow.muVetoMedium == 0) if self.period==13 else\
               (rtrow.elecVetoWZTight + rtrow.muonVetoWZTight == 0)
        if self.metShift and not self.isData: # only shift MC
            vetoMap = {
                #'ees+' : (rtrow.eVetoTight_eesUp + rtrow.muVetoTight == 0),
                #'ees-' : (rtrow.eVetoTight_eesDown + rtrow.muVetoTight == 0),
                #'mes+' : (rtrow.eVetoTight + rtrow.muVetoTight_mesUp == 0),
                #'mes-' : (rtrow.eVetoTight + rtrow.muVetoTight_mesDown == 0),
                'ees+' : (rtrow.eVetoTight_eesUp + rtrow.muVetoMedium == 0),
                'ees-' : (rtrow.eVetoTight_eesDown + rtrow.muVetoMedium == 0),
                'mes+' : (rtrow.eVetoTight + rtrow.muVetoMedium_mesUp == 0),
                'mes-' : (rtrow.eVetoTight + rtrow.muVetoMedium_mesDown == 0),
            }
            if self.metShift in vetoMap:
                veto = vetoMap[self.metShift]

        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            OS01 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) < 0.5 # select opposite sign
            o02 = ordered(l[0],l[2])
            OS02 = getattr(rtrow, "%s_%s_SS" % (o02[0], o02[1])) < 0.5 # select opposite sign
            o12 = ordered(l[1],l[2])
            OS12 = getattr(rtrow, "%s_%s_SS" % (o12[0], o12[1])) < 0.5 # select opposite sign
            st = self.getObject(rtrow,'pt', l[0]) + self.getObject(rtrow,'pt', l[1]) + self.getObject(rtrow,'pt', l[2])

            ordList = [l[1], l[0], l[2]] if self.getObject(rtrow,'pt',l[0]) < self.getObject(rtrow,'pt', l[1]) else [l[0], l[1], l[2]]

            SF01 = l[0][0] == l[1][0]
            SF02 = l[0][0] == l[2][0]
            SF12 = l[1][0] == l[2][0]

            if OS01 and SF01: continue
            if OS02 and SF02: continue
            if OS12 and SF12: continue

            if OS01 and veto:
                cands.append((-st, ordList))

        if not len(cands): return ([],[])

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        negst, leps = cands[0]

        return ([negst], leps)

    def ID_tight_TT(self, rtrow):
        return self.ID(rtrow,*self.objCand[:2],**self.getIdArgs('Tight'))

    def zVeto(self,rtrow):
        leps = self.objCand
        zWindow = 10
        minPt = 20
        minDr = 0.02
        # check tags against probes for Z
        for l in leps[:2]:
            o = ordered(l,leps[2])
            m = self.getObject(rtrow,'mass',o[0],o[1])
            pt = self.getObject(rtrow,'pt',l)
            dr = self.getObject(rtrow,'dr',o[0],o[1])
            ss = getattr(rtrow, "%s_%s_SS" % (o[0], o[1]))
            if o[0][0]==o[1][0] and ss < 0.5:
                if abs(ZMASS-m)<zWindow: return False
            if pt < minPt: return False
            if dr < minDr: return False
        # check tags for Z
        o = ordered(leps[0],leps[1])
        m = self.getObject(rtrow,'mass',o[0],o[1])
        ss = getattr(rtrow, "%s_%s_SS" % (o[0], o[1]))
        if o[0][0]==o[1][0] and ss < 0.5:
            if abs(ZMASS-m)<zWindow: return False

        return True


    def bjetCut(self,rtrow):
        return rtrow.bjetCISVVeto30Medium > 0

    def jetCut(self,rtrow):
        return rtrow.jetVeto30 > 0

    def metCut(self,rtrow):
        if self.getObject(rtrow,'pt','met') < 20.: return False
        return True

    def trigger(self, rtrow):
        if self.period == 8:
            triggers = ["mu17ele8isoPass", "mu8ele17isoPass",
                        "doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == 13:
            triggers = ['singleMuSingleEPass', 'doubleMuPass', 'doubleEPass', 'singleESingleMuPass']

        nottriggers = []

        good = False
        for t in triggers:
            if getattr(rtrow,t)>0:
                good = True
                break
        for t in nottriggers:
            if getattr(rtrow,t)>0:
                good = False
                break
        return good

    def good_to_store(self, rtrow, cand1, cand2):
        '''
        Iterate through minimizing variables.
        '''
        good = False
        for min1, min2 in zip(cand1, cand2):
            if min1 < min2:
                good = True
                break
            if min1 > min2:
                good = False
                break
        if self.objCand[0][0]=='e' and self.objCand[1][0]=='e':
            lname = 'MatchesDoubleE'
            match_0 = getattr(rtrow,'%s%s' %(self.objCand[0],lname))
            match_1 = getattr(rtrow,'%s%s' %(self.objCand[1],lname))
        elif self.objCand[0][0]=='m' and self.objCand[1][0]=='m':
            lname = 'MatchesDoubleMu'
            match_0 = getattr(rtrow,'%s%s' %(self.objCand[0],lname))
            match_1 = getattr(rtrow,'%s%s' %(self.objCand[1],lname))
        else:
            emname = 'MatchesSingleESingleMu'
            mename = 'MatchesSingleMuSingleE'
            match_0 = getattr(rtrow,'%s%s' %(self.objCand[0],emname)) + getattr(rtrow,'%s%s' %(self.objCand[0],mename))
            match_1 = getattr(rtrow,'%s%s' %(self.objCand[1],emname)) + getattr(rtrow,'%s%s' %(self.objCand[1],mename))
        passTrig = (match_0 > 0.5 and match_1 > 0.5)

        #veto = (rtrow.eVetoTrigIso + rtrow.muVetoTightTrigIso == 0) if self.period==13 else (rtrow.elecVetoWZLoose + rtrow.muonVetoWZLoose == 0)
        veto = (rtrow.eVetoTrigIso + rtrow.muVetoMediumTrigIso == 0) if self.period==13 else (rtrow.elecVetoWZLoose + rtrow.muonVetoWZLoose == 0)

        return good and passTrig and veto


class AnalyzerWZ_HZZFakeRate(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_HZZFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'HZZFakeRate'

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.ID_tight_Z)
        cuts.add(self.crossCleaning)
        cuts.add(self.zSelection)
        cuts.add(self.metveto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.crossCleaning)
        cuts.add(self.zSelection)
        cuts.add(self.metveto)
        return cuts

    def ID_tight_Z(self, rtrow):
        return self.ID(rtrow,*self.objCand[:2],**self.getIdArgs('Tight'))

    def metveto(self,rtrow):
        if rtrow.type1_pfMetEt > 25.: return False
        return True

    def getIdArgs(self,type):
        kwargs = {}
        if type=='Tight':
            kwargs['idDef'] = {
                'e':'ZZTight',
                'm':'ZZTight',
            }
            kwargs['isoCut'] = {
                'e':0.5,
                'm':0.4
            }
        if type=='Loose':
            kwargs['idDef'] = {
                'e':'ZZLoose',
                'm':'ZZLoose',
            }
            kwargs['isoCut'] = {
                'e': 999., # no iso
                'm': 999., # no iso
            }
        if hasattr(self,'alternateIds'):
            if type in self.alternateIds:
                kwargs = self.alternateIdMap[type]
        return kwargs

    def trigger(self, rtrow):
        if self.period == 8:
            triggers = ["doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == 13:
            triggers = ['doubleMuPass', 'doubleEPass']

        for t in triggers:
            if getattr(rtrow,t)>0:
                return True
        return False

    def crossCleaning(self,rtrow):
        leps = self.objCand
        for l0 in leps:
           for l1 in leps:
               if l0[0]=='e' and l1[0]=='m':
                   if self.ID(rtrow,l1,**self.getIdArgs('Tight')):
                       o = ordered(l0,l1)
                       dr = getattr(rtrow,'%s_%s_DR' % (o[0],o[1]))
                       if dr > 0.05: return False
        return True

    def zSelection(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m = getattr(rtrow,'%s_%s_Mass' % (o[0],o[1]))
        l0Pt = getattr(rtrow,'%sPt' %leps[0])
        wl = leps[2]
        o0 = ordered(leps[0],wl)
        o1 = ordered(leps[1],wl)
        dr0 = getattr(rtrow,'%s_%s_DR' % (o0[0],o0[1]))
        dr1 = getattr(rtrow,'%s_%s_DR' % (o1[0],o1[1]))
        m0 = getattr(rtrow,'%s_%s_Mass' % (o0[0],o0[1]))
        m1 = getattr(rtrow,'%s_%s_Mass' % (o1[0],o1[1]))
        return abs(m-ZMASS)<10. and l0Pt>20. and dr0>0.02 and dr1>0.02 and m0>4. and m1>4.

    def good_to_store(self, rtrow, cand1, cand2):
        '''
        Iterate through minimizing variables.
        '''
        good = False
        for min1, min2 in zip([cand1[0]], [cand2[0]]):
            if min1 < min2:
                good = True
                break
            if min1 > min2:
                good = False
                break
        match_0 = getattr(rtrow,'%sMatchesDoubleE' %self.objCand[0])>0 if self.objCand[0][0]=='e' else getattr(rtrow,'%sMatchesDoubleMu' %self.objCand[0])>0
        match_1 = getattr(rtrow,'%sMatchesDoubleE' %self.objCand[1])>0 if self.objCand[1][0]=='e' else getattr(rtrow,'%sMatchesDoubleMu' %self.objCand[1])>0
        passTrig = match_0 > 0.5 and match_1 > 0.5

        veto = (rtrow.eVetoHZZ + rtrow.muVetoHZZ == 0)

        return good and passTrig and veto


##########################
###### Command line ######
##########################
def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('analyzer', type=str)
    parser.add_argument('sample_name', type=str)
    parser.add_argument('file_list', type=str)
    parser.add_argument('out_file', type=str)
    parser.add_argument('period', type=int)
    parser.add_argument('-ms','--metShift',nargs='?',type=str,const='',help='Shift the met')

    args = parser.parse_args(argv)
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.analyzer == 'WZ': analyzer = AnalyzerWZ(args.sample_name,args.file_list,args.out_file,args.period,metShift=args.metShift)
    if args.analyzer == 'NoVeto': analyzer = AnalyzerWZ_NoVeto(args.sample_name,args.file_list,args.out_file,args.period,metShift=args.metShift)
    if args.analyzer == 'FakeRate': analyzer = AnalyzerWZ_ZFakeRate(args.sample_name,args.file_list,args.out_file,args.period,metShift=args.metShift)
    if args.analyzer == 'TTFakeRate': analyzer = AnalyzerWZ_TTFakeRate(args.sample_name,args.file_list,args.out_file,args.period,metShift=args.metShift)
    if args.analyzer == 'HZZFakeRate': analyzer = AnalyzerWZ_HZZFakeRate(args.sample_name,args.file_list,args.out_file,args.period,metShift=args.metShift)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
