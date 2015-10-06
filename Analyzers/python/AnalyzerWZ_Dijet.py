#!/usr/bin/env python
'''
The WZ analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerWZ_DijetFakeRate(AnalyzerBase):
    '''
    An implementation of the AnalyzerBase class for use in WZ fake rate analysis.
    '''

    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        if not hasattr(self,'channel'): self.channel = 'FakeRate'
        self.period = period
        self.final_states = ['e','m']
        self.initial_states = ['w1'] # in order of leptons returned in choose_objects
        self.object_definitions = {
            'w1': ['em','n'],
        }
        self.lepargs = {'tight':True}
        self.cutflow_labels = []
        self.doVBF = (period==13)
        super(AnalyzerWZ_DijetFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select lepton
        '''
        leps = self.objects
        return ([0.], leps)

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

    # overide good_to_store
    # will store via trigger selection
    #@staticmethod
    def good_to_store(self,rtrow, cand1, cand2):
        '''
        Select the trigger object
        '''
        #singleTrigMatch_leg1 = getattr(rtrow,'%sMatchesSingleE_leg1' %self.objCand[0]) if self.objCand[0][0]=='e' else getattr(rtrow,'%sMatchesSingleMu_leg1' %self.objCand[0])
        singleTrigMatch_leg2 = getattr(rtrow,'%sMatchesSingleE_leg2' %self.objCand[0]) if self.objCand[0][0]=='e' else getattr(rtrow,'%sMatchesSingleMu_leg2' %self.objCand[0])
        veto = (rtrow.eVeto + rtrow.muVeto == 0)
        return singleTrigMatch_leg2 and veto

    def getTriggerPrescale(self,rtrow):
        #singleTrigPrescale_leg1 = rtrow.singleE_leg1Prescale if self.objCand[0][0]=='e' else rtrow.singleMu_leg1Prescale
        singleTrigPrescale_leg2 = rtrow.singleE_leg2Prescale if self.objCand[0][0]=='e' else rtrow.singleMu_leg2Prescale
        return singleTrigPrescale_leg2

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        #cuts.add(self.ID_tight)
        #cuts.add(self.ID_loose)
        cuts.add(self.ID_veto)
        cuts.add(self.zVeto)
        cuts.add(self.jPsiVeto)
        cuts.add(self.wVeto)
        cuts.add(self.jetSelection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.zVeto)
        cuts.add(self.jPsiVeto)
        cuts.add(self.wVeto)
        cuts.add(self.jetSelection)
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
            if self.period==8:
                kwargs['idDef']['e'] = 'WZTight'
                kwargs['idDef']['m'] = 'WZTight'
            if self.period==13:
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
            if self.period==8:
                kwargs['idDef']['e'] = 'WZLoose'
                kwargs['idDef']['m'] = 'WZLoose'
            if self.period==13:
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
            if self.period==13:
                kwargs['isoCut']['e'] = 9999.
        if hasattr(self,'alternateIds'):
            if type in self.alternateIds:
                kwargs = self.alternateIdMap[type]
        return kwargs

    def trigger(self, rtrow):
        triggers = ['singleE_leg2Pass', 'singleMu_leg2Pass']

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

    def zVeto(self,rtrow):
        return getattr(rtrow,'%sNearestZMass' %self.objCand[0]) > 30.

    def jPsiVeto(self,rtrow):
        return getattr(rtrow,'%sLowestMll' %self.objCand[0]) > 20.

    def wVeto(self,rtrow):
        leps = self.objCand
        if self.period==8:
            if rtrow.type1_pfMetEt > 20.: return False
            if getattr(rtrow, "%sMtToPfMet_Ty1" % (leps[0])) > 25.: return False
        else:
            if rtrow.pfMetEt > 20.: return False
            if getattr(rtrow, "%sMtToPFMET" % (leps[0])) > 25.: return False
        return True

    def jetSelection(self,rtrow):
        leps = self.objCand
        if rtrow.jet1Pt<20: return False # moderate leading jet pt requirement
        lEta = getattr(rtrow,'%sEta' %leps[0])
        lPhi = getattr(rtrow,'%sPhi' %leps[0])
        jEta = rtrow.jet1Eta
        jPhi = rtrow.jet1Phi
        dr = deltaR(lEta,lPhi,jEta,jPhi)
        return dr>1. # jet far from lepton

class AnalyzerWZ_HZZDijetFakeRate(AnalyzerWZ_DijetFakeRate):
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerWZ_HZZDijetFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'HZZFakeRate'

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_loose)
        cuts.add(self.zVeto)
        cuts.add(self.jPsiVeto)
        cuts.add(self.wVeto)
        cuts.add(self.jetSelection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        cuts.add(self.zVeto)
        cuts.add(self.jPsiVeto)
        cuts.add(self.wVeto)
        cuts.add(self.jetSelection)
        return cuts

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
                'e':0.5,
                'm':0.4
            }
        if hasattr(self,'alternateIds'):
            if type in self.alternateIds:
                kwargs = self.alternateIdMap[type]
        return kwargs

    def good_to_store(self,rtrow, cand1, cand2):
        '''
        Select the trigger object
        '''
        #singleTrigMatch_leg1 = getattr(rtrow,'%sMatchesSingleE_leg1' %self.objCand[0]) if self.objCand[0][0]=='e' else getattr(rtrow,'%sMatchesSingleMu_leg1' %self.objCand[0])
        singleTrigMatch_leg2 = getattr(rtrow,'%sMatchesSingleE_leg2' %self.objCand[0]) if self.objCand[0][0]=='e' else getattr(rtrow,'%sMatchesSingleMu_leg2' %self.objCand[0])
        veto = (rtrow.eVetoHZZIso + rtrow.muVetoHZZIso == 0)
        return singleTrigMatch_leg2 and veto


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

    if args.analyzer == 'FakeRate': analyzer = AnalyzerWZ_DijetFakeRate(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer == 'HZZFakeRate': analyzer = AnalyzerWZ_HZZDijetFakeRate(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
