#!/usr/bin/env python
'''
The WZ analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerWZ_QCD(AnalyzerBase):
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
        if not hasattr(self,'channel'): self.channel = 'QCD'
        self.period = period
        self.final_states = ['ej','mj']
        self.initial_states = ['w1'] # in order of leptons returned in choose_objects
        self.object_definitions = {
            'w1': ['em','n'],
        }
        self.lepargs = {'tight':True}
        self.cutflow_labels = []
        #self.alternateIds, self.alternateIdMap = self.defineAlternateIds(period)
        self.doVBF = (period=='13')
        super(AnalyzerWZ_QCD, self).__init__(sample_name, file_list, out_file, period)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        highest pt lepton (eventually, match it to the trigger lepton)
        '''
        cands = []
        for l in permutations(self.objects):

            if l[0][0]=='j': continue
            lpt = getattr(rtrow,'%sPt' % l[0])
            jpt = getattr(rtrow,'%sPt' % l[1])
            cands.append((1./(lpt+jpt), l))

        if not len(cands): return 0

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        pt, leps = cands[0]

        return ([pt], leps)

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        #cuts.add(self.ID_tight)
        #cuts.add(self.ID_loose)
        cuts.add(self.ID_veto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
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
            triggers = ['muEPass', 'doubleMuPass', 'doubleEPass', 'eMuPass']

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

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))

    def mass3l(self,rtrow):
        return rtrow.Mass > 100.

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

    if args.analyzer == 'QCD': analyzer = AnalyzerWZ_QCD(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
