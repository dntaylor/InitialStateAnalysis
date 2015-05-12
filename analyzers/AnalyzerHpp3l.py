#!/usr/bin/env python
'''
The doubly charged higgs associated production analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerHpp3l(AnalyzerBase):
    '''
    The primary analyzer for the doubly charged higgs associated production channel.

    Objects:
        doubly charged higgs and signly charged higgs
        h1: H++ = 2 same sign leptons (not necessarily same flavor)
        h2: H- = 1 lepton and 1 neutrino
    Minimization function:
        None, 4th lepton veto. Override good_to_store
    Selection:
        Trigger: double lepton
        Fiducial
        Trig threshold
        ID: muon tight, cbid tight mva trig
        Isolation: 0.15 (0.12) elec (muon)
        QCD suppression: M(ll) > 12
    '''

    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        runTau = kwargs.pop('runTau',False)
        #runTau=True
        self.channel = 'Hpp3l'
        self.final_states = ['eee','eem','emm','mmm'] # no tau
        if runTau: self.final_states = ['eee','eem','eet','emm','emt','ett','mmm','mmt','mtt','ttt']
        self.initial_states = ['h1','h2']
        self.other_states = [['z1', 'w1']]
        self.object_definitions = {
            'h1': ['em','em'],
            'h2': ['em','n'],
            'z1': ['em','em'],
            'w1': ['em','n'],
        }
        if runTau:
            self.object_definitions['h1'] = ['emt', 'emt']
            self.object_definitions['h2'] = ['emt', 'n']
            self.object_definitions['z1'] = ['emt', 'emt']
            self.object_definitions['w1'] = ['emt', 'n']
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','QCD Suppression']
        super(AnalyzerHpp3l, self).__init__(sample_name, file_list, out_file, period)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select candidate objects
        '''
        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            SS1 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) > 0 # select same sign
            OS = getattr(rtrow, "%sCharge" % l[0]) != getattr(rtrow, "%sCharge" % l[2]) # select opposite sign

            if SS1 and OS:
                cands.append([[0],list(l)]) # minimization is by veto, not variable

        if not len(cands): return 0

        return (cands[0])

    # override choose_alternative_objects
    def choose_alternative_objects(self, rtrow, state):
        '''
        Select alternative candidate.
        '''
        # WZ
        if state == ['z1', 'w1']:
            bestZDiff = float('inf')
            bestLeptons = []

            for l in permutations(self.objects):
                if lep_order(l[0],l[1]):
                    continue

                os1 = getattr(rtrow,'%s_%s_SS' % (l[0], l[1])) < 0.5
                m1 = getattr(rtrow,'%s_%s_Mass' % (l[0], l[1]))

                if l[0][0] == l[1][0] and os1 and abs(m1-ZMASS) < bestZDiff:
                    bestZDiff = abs(m1-ZMASS)
                    bestLeptons = l

            return bestLeptons

    # overide good_to_store
    @staticmethod
    def good_to_store(rtrow, cand1, cand2):
        '''
        Veto on 4th lepton (considered in 4l analysis)
        '''
        return (rtrow.eVetoMVAIsoVtx + rtrow.muVetoPt5IsoIdVtx == 0)

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_loose)
        cuts.add(self.qcd_rejection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.qcd_rejection)
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
        return kwargs

    def trigger(self, rtrow):
        triggers = ["mu17ele8isoPass", "mu8ele17isoPass",
                    "doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == '13':
            triggers = ['muEPass', 'eMuPass', 'doubleMuPass',
                        'doubleEPass', 'tripleEPass']

        for t in triggers:
            if getattr(rtrow,t)>0:
                return True
        return False

    def fiducial(self, rtrow):
        for l in self.objects:
            if l[0]=='e':
                ptcut = 10.0 # CBID: 20, MVA trig: 10, MVA nontrig: 5
                etacut = 2.5
            if l[0]=='m':
                ptcut = 10.0 # ???
                etacut = 2.4
            if l[0]=='t':
                ptcut = 20.0 # 20
                etacut = 2.3
            if getattr(rtrow, '%sPt' % l) < ptcut:
                return False
            if getattr(rtrow, '%sAbsEta' % l) > etacut:
                return False
        return True

    def ID_loose(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Loose'))

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))

    def trigger_threshold(self, rtrow):
        pts = [getattr(rtrow, "%sPt" % l) for l in self.objects]
        pts.sort(reverse=True)
        return pts[0] > 20.0 and pts[1] > 10.0

    def qcd_rejection(self, rtrow):
        qcd_pass = [getattr(rtrow, "%s_%s_Mass" % (l[0], l[1])) > 12.0
                    for l in combinations(self.objects, 2)]
        return all(qcd_pass)

#######################
### Control regions ###
#######################
class AnalyzerHpp3l_WZ(AnalyzerHpp3l):
    '''
    WZ control region for Hpp3l
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerHpp3l_WZ, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'WZ'
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','Mass 3l','Z selection','W selection']

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_loose)
        cuts.add(self.mass3l)
        cuts.add(self.zSelection)
        cuts.add(self.wSelection)

        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.mass3l)
        cuts.add(self.zSelection)
        cuts.add(self.wSelection)
        return cuts

    def mass3l(self,rtrow):
        return rtrow.Mass > 100.

    def zSelection(self,rtrow):
        leps = self.choose_alternative_objects(rtrow, ['z1','w1'])
        if not leps: return False
        o = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o[0],o[1]))
        l0Pt = getattr(rtrow,'%sPt' %leps[0])
        return abs(m1-ZMASS)<20. and l0Pt>20.

    def wSelection(self,rtrow):
        leps = self.choose_alternative_objects(rtrow, ['z1','w1'])
        if not leps: return False
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



#######################
###### Fake rate ######
#######################
class AnalyzerHpp3l_FakeRate(AnalyzerHpp3l):
    '''
    A class to produce ntuples to calculate the fakerate for the leptons.
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerHpp3l_FakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'FakeRate'
        self.final_states = ['emm','mmm','mmt'] 
        self.initial_states = ['z1','f1']
        self.other_states = []
        self.object_definitions = {
            'z1': ['m','m'],
            'f1': ['emt'],
        }
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','QCD Suppression','Z Selection']

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select candidate objects
        '''
        cands = []
        bestZDiff = float('inf')
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            # first two must be the Z candidate
            OS1 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) < 0.5 # select opposite sign
            SF1 = l[0][0]==l[1][0] # select same flavor
            massdiff = abs(getattr(rtrow,'%s_%s_Mass' % (l[0], l[1]))-ZMASS)

            if OS1 and SF1:
                cands.append([massdiff,list(l)])

        if not len(cands): return 0

        cands.sort(key=lambda x: x[0])
        massdiff, leps = cands[0]

        return ([massdiff], leps)

    # reoveride
    def choose_alternative_objects(self, rtorw, state):
        return []

    @staticmethod
    def good_to_store(rtrow, cand1, cand2):
        '''
        Iterate through minimizing variables.
        '''
        for min1, min2 in zip(cand1, cand2):
            if min1 < min2: return True
            if min1 > min2: return False
        return False

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_loose)
        cuts.add(self.qcd_rejection)
        cuts.add(self.z_selection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.qcd_rejection)
        cuts.add(self.z_selection)
        return cuts

    def z_selection(self,rtrow):
        '''Select Z candidate'''
        m1 = getattr(rtrow,'%s_%s_Mass' % (self.objects[0], self.objects[1]))
        return abs(m1-ZMASS)<20.




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

    if args.analyzer=='Hpp3l': analyzer = AnalyzerHpp3l(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer=='WZ': analyzer = AnalyzerHpp3l_WZ(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer=='FakeRate': analyzer = AnalyzerHpp3l_FakeRate(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
