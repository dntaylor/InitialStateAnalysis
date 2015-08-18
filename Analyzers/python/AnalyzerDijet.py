#!/usr/bin/env python
'''
The WZ analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerDijet(AnalyzerBase):
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
        if not hasattr(self,'channel'): self.channel = 'Dijet'
        self.period = period
        self.final_states = ['ee','em','mm']
        self.initial_states = ['z1'] # in order of leptons returned in choose_objects
        self.other_states = [['w1', 'w2']]
        self.object_definitions = {
            'z1': ['em','em'],
            'w1': ['em','n'],
            'w2': ['em','n'],
        }
        self.lepargs = {'tight':True}
        self.cutflow_labels = ['Trigger','Fiducial','ID','Z Selection','W Selection']
        self.doVBF = (period=='13')
        super(AnalyzerWZ, self).__init__(sample_name, file_list, out_file, period,**kwargs)

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

            ordList = [l[1], l[0], l[2]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l

            if OS1 and l[0][0]==l[1][0]:
                cands.append((massdiff, list(ordList)))

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

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.metveto)
        cuts.add(self.wveto)
        cuts.add(self.zveto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.metveto)
        cuts.add(self.wveto)
        cuts.add(self.zveto)
        return cuts

    def getIdArgs(self,type):
        kwargs = {}
        if type=='Tight':
            kwargs['idDef'] = {
                'e':'Medium',
                'm':'Tight',
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
                'm':'Loose',
            }
            kwargs['isoCut'] = {
                'e':0.4,
                'm':0.4
            }
        if hasattr(self,'alternateIds'):
            if type in self.alternateIds:
                kwargs = self.alternateIdMap[type]
        return kwargs

    def trigger(self, rtrow):
        triggers = ["singleMuPass", "singleEPass"]

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

    def metveto(self,rtrow):
        if self.period=='8':
            if rtrow.type1_pfMetEt > 20.: return False
        else:
            if rtrow.pfMetEt > 20.: return False
        return True

    def wveto(self,rtrow):
        mtVar = 'PFMET' if self.period=='13' else 'PfMet_Ty1'
        for obj in self.objects:
            if obj[0] == 'm':
                if getattr(rtrow, "%sMtTo%s" % (obj, mtVar)) > 20.: return False
        return True

    def zveto(self,rtrow):
        leps = self.objects
        o = ordered(leps[0], leps[1])
        if o[0]!=o[1]: return True
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

    if args.analyzer == 'Dijet': analyzer = AnalyzerDijet(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
