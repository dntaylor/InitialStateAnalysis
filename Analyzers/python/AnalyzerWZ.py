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
        self.lepargs = {'tight':True}
        self.cutflow_labels = ['Trigger','Fiducial','ID','Z Selection','W Selection']
        #self.alternateIds, self.alternateIdMap = self.defineAlternateIds(period)
        self.doVBF = (period=='13')
        super(AnalyzerWZ, self).__init__(sample_name, file_list, out_file, period)

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
            mass = getattr(rtrow, "%s_%s_Mass" % (l[0], l[1]))
            massdiff = abs(ZMASS-mass)

            ordList = [l[1], l[0], l[2]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else [l[0], l[1], l[2]]

            if OS1 and l[0][0]==l[1][0]:
                cands.append((massdiff, ordList))

        if not len(cands): return 0

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        massdiff, leps = cands[0]

        return ([massdiff], leps)

    # overide good_to_store
    # will store via veto
    #@staticmethod
    def good_to_store(self,rtrow, cand1, cand2):
        '''
        Veto on 4th lepton
        '''
        return (rtrow.eVetoTight + rtrow.muVetoTight == 0) if self.period=='13' else\
               (rtrow.elecVetoWZTight + rtrow.muonVetoWZTight == 0)

    def defineAlternateIds(self,period):
        if period=='8':
            return [], {}
        elecIds = ['Loose', 'Medium', 'Tight']
        muonIds = ['Loose', 'Tight']
        elecIsos = [0.5, 0.2, 0.15]
        muonIsos = [0.4, 0.2, 0.12]
        idList = []
        idMap = {}
        for id in elecIds:
            for iso in elecIsos:
                idName = 'elec%s%0.2f' % (id, iso) if iso else 'elec%sNoIso' % id
                idName = idName.replace('.','p')
                idList += [idName]
                idMap[idName] = {
                    'idDef' : {
                        'e': id
                    }
                }
                if iso:
                    idMap[idName]['isoCut'] = {
                        'e': iso
                    }
        for id in muonIds:
            for iso in muonIsos:
                idName = 'muon%s%0.2f' % (id, iso) if iso else 'muon%sNoIso' % id
                idName = idName.replace('.','p')
                idList += [idName]
                idMap[idName] = {
                    'idDef' : {
                        'm': id
                    }
                }
                if iso:
                    idMap[idName]['isoCut'] = {
                        'm': iso
                    }
        return idList, idMap
                

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        #if self.period=='13': cuts.add(self.passAnyId)
        #cuts.add(self.ID_tight)
        #cuts.add(self.ID_loose)
        cuts.add(self.ID_veto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.mass3l)
        cuts.add(self.zSelection)
        cuts.add(self.wSelection)
        return cuts

    def getIdArgs(self,type):
        kwargs = {}
        if type=='Tight':
            kwargs['idDef'] = {
                'e':'Medium',
                'm':'Tight',
                't':'Medium'
            }
            kwargs['isoCut'] = {
                'e':0.15,
                'm':0.12
            }
            if self.period=='8':
                kwargs['idDef']['e'] = 'WZTight'
                kwargs['idDef']['m'] = 'WZTight'
            if self.period=='13':
                kwargs['isoCut']['e'] = 9999.
        if type=='Loose':
            kwargs['idDef'] = {
                'e':'Loose',
                'm':'Loose',
                't':'Loose'
            }
            kwargs['isoCut'] = {
                'e':0.2,
                'm':0.2
            }
            if self.period=='8':
                kwargs['idDef']['e'] = 'WZLoose'
                kwargs['idDef']['m'] = 'WZLoose'
            if self.period=='13':
                kwargs['isoCut']['e'] = 9999.
        if type=='Veto':
            kwargs['idDef'] = {
                'e':'Veto',
                'm':'ZZLoose',
                't':'Loose'
            }
            kwargs['isoCut'] = {
                'e':0.4,
                'm':0.4
            }
            if self.period=='13':
                kwargs['isoCut']['e'] = 9999.
        if hasattr(self,'alternateIds'):
            if type in self.alternateIds:
                kwargs = self.alternateIdMap[type]
        return kwargs

    def trigger(self, rtrow):
        if self.period == '8':
            triggers = ["mu17ele8isoPass", "mu8ele17isoPass",
                        "doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == '13':
            triggers = ['singleMuSingleEPass', 'doubleMuPass', 'doubleEPass', 'singleESingleMuPass']

        for t in triggers:
            if getattr(rtrow,t)>0:
                return True
        return False

    def fiducial(self, rtrow):
        for l in self.objects:
            if l[0]=='e':
                ptcut = 10.0
                etacut = 2.5
            if l[0]=='m':
                ptcut = 10.0
                etacut = 2.4
            if l[0]=='t':
                ptcut = 20.0
                etacut = 2.3
            if getattr(rtrow, '%sPt' % l) < ptcut:
                return False
            if getattr(rtrow, '%sAbsEta' % l) > etacut:
                return False
        return True

    def passAnyId(self,rtrow):
        '''Check to make sure the leptons pass at least 1 ID'''
        passCheck = {'e': False, 'm': False}
        for altId in self.alternateIds:
           if passCheck[altId[0]]: continue
           if self.ID(rtrow,*self.objects,**self.getIdArgs(altId)): passCheck[altId[0]] = True
        return passCheck['e'] and passCheck['m']

    def ID_veto(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Veto'))

    def ID_loose(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Loose'))

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))

    def mass3l(self,rtrow):
        return rtrow.Mass > 100.

    def zSelection(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o[0],o[1]))
        l0Pt = getattr(rtrow,'%sPt' %leps[0])
        #return abs(m1-ZMASS)<20. and l0Pt>20.
        return (m1>=60. and m1<=120. and l0Pt>20.)

    def wSelection(self,rtrow):
        leps = self.objCand
        if getattr(rtrow, '%sPt' %leps[2])<20.: return False
        if self.period=='8':
            if rtrow.type1_pfMetEt < 30.: return False
        else:
            if rtrow.pfMetEt < 30.: return False
        for l in leps[:2]:
            o = ordered(l,leps[2])
            dr = getattr(rtrow, '%s_%s_DR' % (o[0],o[1]))
            if dr < 0.1: return False
        return True

class AnalyzerWZ_DataDriven(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_DataDriven, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'DataDriven'

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        return cuts

class AnalyzerWZ_Z(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_Z, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'Z'

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.ID_tight_Z)
        cuts.add(self.zSelection)
        cuts.add(self.metveto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        cuts.add(self.metveto)
        return cuts

    def ID_tight_Z(self, rtrow):
        return self.ID(rtrow,*self.objCand[:2],**self.getIdArgs('Tight'))

    def metveto(self,rtrow):
        if self.period=='8':
            if rtrow.type1_pfMetEt > 20.: return False
        else:
            if rtrow.pfMetEt > 20.: return False
        return True

class AnalyzerWZ_QCD(AnalyzerWZ):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_QCD, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'QCD'

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.metveto)
        cuts.add(self.wveto)
        cuts.add(self.zveto)
        return cuts

    def trigger(self, rtrow):
        triggers = ["singleEPass", "singleMuPass"]

        for t in triggers:
            if getattr(rtrow,t)>0:
                return True
        return False

    def metveto(self,rtrow):
        if self.period=='8':
            if rtrow.type1_pfMetEt > 20.: return False
        else:
            if rtrow.pfMetEt > 20.: return False
        return True

    def wveto(self,rtrow):
        mtVar = 'PFMET' if self.period=='13' else 'PfMet_Ty1'
        obj = self.objCand[2]
        if obj[0] == 'm':
            return getattr(rtrow, "%sMtTo%s" % (obj, mtVar)) < 20.
        return True

    def zveto(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o[0],o[1]))
        if o[0][0] == 'e':
            return abs(m1-ZMASS)>30. and m1>20.
        if o[0][0] == 'm':
            return abs(m1-ZMASS)>15. and m1>20.


##########################
###### Command line ######
##########################
def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('analyzer', type=str)
    parser.add_argument('sample_name', type=str)
    parser.add_argument('file_list', type=str)
    parser.add_argument('out_file', type=str)
    parser.add_argument('period', type=str)

    args = parser.parse_args(argv)
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.analyzer == 'WZ': analyzer = AnalyzerWZ(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer == 'DataDriven': analyzer = AnalyzerWZ_DataDriven(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer == 'Z': analyzer = AnalyzerWZ_Z(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer == 'QCD': analyzer = AnalyzerWZ_QCD(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
