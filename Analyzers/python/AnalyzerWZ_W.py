#!/usr/bin/env python
'''
The WZ analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerWZ_WFakeRate(AnalyzerBase):
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
        if not hasattr(self,'channel'): self.channel = 'FakeRate'
        self.period = period
        self.final_states = ['ee','em','mm']
        self.initial_states = ['w1','w2'] # in order of leptons returned in choose_objects
        self.other_states = [['z1']]
        self.object_definitions = {
            'w1': ['em','n'],
            'w2': ['em','n'],
            'z1': ['em','em'],
        }
        self.lepargs = {'tight':True}
        self.cutflow_labels = []
        self.doVBF = (period==13)
        self.tightW = False
        super(AnalyzerWZ_WFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select two highest pt leptons
        '''
        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            pt0 = getattr(rtrow,'%sPt' % l[0])
            pt1 = getattr(rtrow,'%sPt' % l[1])
            ordList = [l[1], l[0]] if  pt0 < pt1 else [l[0], l[1]]

            cands.append((-1.*(pt0+pt1), ordList))

        if not len(cands): return 0

        # Sort by highest st difference
        cands.sort(key=lambda x: x[0])
        st, leps = cands[0]

        return ([st], leps)

    # override choose_alternative_objects
    def choose_alternative_objects(self, rtrow, state):
        '''
        Select alternative candidate.
        '''
        # Z
        if state == ['z1']:
            bestZDiff = float('inf')
            bestLeptons = []

            for l in permutations(self.objects):
                if lep_order(l[0],l[1]):
                    continue

                m1 = getattr(rtrow,'%s_%s_Mass' % (l[0], l[1]))

                if abs(m1-ZMASS) < bestZDiff:
                    bestZDiff = abs(m1-ZMASS)
                    ordList = [l[1], l[0]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l
                    bestLeptons = ordList

            return bestLeptons

    def getFakeChannel(self,rtrow,**kwargs):
        tightW = kwargs.pop('tightW',False)
        allMedium = kwargs.pop('allMedium',False)
        chan = ''
        for i,obj in enumerate(self.objCand):
            if (tightW and i!=2) or allMedium:
                chan += 'P' if self.ID(rtrow,obj,**self.getIdArgs('Medium')) else 'F'
            else:
                chan += 'P' if self.ID(rtrow,obj,**self.getIdArgs('Tight')) else 'F'
        return chan

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        #cuts.add(self.ID_tight)
        cuts.add(self.ID_loose)
        #cuts.add(self.zVeto)
        #cuts.add(self.wSelection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        #cuts.add(self.zVeto)
        #cuts.add(self.wSelection)
        return cuts

    def getIdArgs(self,type):
        kwargs = {}
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
                kwargs['isoCut']['e'] = 0. # baked into ID
                kwargs['isoCut']['m'] = 0. # baked into ID
        if type=='Medium':
            kwargs['idDef'] = {
                'e':'WWMedium',
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
                        'singleIsoMu20Pass', 'singleIsoTkMu20Pass', 'singleIsoMu27Pass', 'singleE23WPLoosePass']

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

    def ID_veto(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Veto'))

    def ID_loose(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Loose'))

    def ID_medium(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Medium'))

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))

    def mass3l(self,rtrow):
        return rtrow.Mass > 100.

    def zVeto(self,rtrow):
        leps = self.choose_alternative_objects(rtrow, ['z1'])
        o = ordered(leps[0], leps[1])
        m1 = self.getObject(rtrow,'mass',o[0],o[1])
        os = getattr(rtrow,'%s_%s_SS' % (o[0], o[1])) < 0.5
        sf = o[0][0]==o[1][0]
        if os and sf:
            if abs(m1-ZMASS) < 30.: return False
        return True

    def wSelection(self,rtrow):
        leps = self.objCand
        if self.getObject(rtrow, 'pt', leps[0])<20.: return False
        if self.getObject(rtrow,'pt','met') < 40.: return False
        if self.getObject(rtrow, "mt", leps[0], 'met') < 40.: return False
        o = ordered(*leps)
        if self.getObject(rtrow, 'mass', o[0],o[1]) < 4: return False
        return True

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

    args = parser.parse_args(argv)
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.analyzer == 'FakeRate': analyzer = AnalyzerWZ_WFakeRate(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
