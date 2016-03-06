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
        self.tightW = False
        #self.doVBF = (period==13)
        super(AnalyzerWZ_DijetFakeRate, self).__init__(sample_name, file_list, out_file, period, **kwargs)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select lepton
        '''
        leps = self.objects
        pt = self.getObject(rtrow,'pt',leps[0])
        return ([-pt], leps)

    # overide good_to_store
    # will store via trigger selection
    #@staticmethod
    def good_to_store(self,rtrow, cand1, cand2):
        '''
        Select the trigger object
        '''
        obj = self.objCand[0]
        pt = self.getObject(rtrow,'pt',obj)
        #if obj=='e':
        #    if pt<20.: # use HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
        #        match = getattr(rtrow,'%sMatchesSingleE_leg2' %obj) > 0.5
        #    elif pt<30.: # use HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
        #        match = getattr(rtrow,'%sMatchesSingleE_leg1' %obj) > 0.5
        #    else:      # use HLT_Ele23_WPLoose_Gsf_v* TODO: this is actualy Ele27
        #        match = getattr(rtrow,'%sMatchesSingleE' %obj) > 0.5
        #elif obj=='m':
        #    if pt<20.: # use HLT_Mu8_TrkIsoVVL_v*
        #        match = getattr(rtrow,'%sMatchesSingleMu_leg2' %obj) > 0.5
        #    elif pt<25: # use HLT_Mu17_TrkIsoVVL_v*
        #        match = getattr(rtrow,'%sMatchesSingleMu_leg1' %obj) > 0.5
        #    else: # use HLT_IsoMu20_v*
        #        match = getattr(rtrow,'%sMatchesSingleMuIso20' %obj) > 0.5
        #else:
        #    match = 0.
        veto = self.veto(rtrow)
        #if obj=='m' and pt>20:
        #    print obj, pt, 'Match:', match, 'Num e:', rtrow.eVetoLoose, 'Num m:', rtrow.muVetoLoose
        #return match and veto
        return veto

    def getTriggerPrescale(self,rtrow):
        obj = self.objCand[0]
        pt = self.getObject(rtrow,'pt',obj)
        if obj=='e':
            #if pt<20: # use HLT_Ele12_CaloIdL_TrackIdL_IsoVL_v*
            prescale = 2263.552/4.174
            #elif pt<30: # use HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v*
            #else: # use HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v* not filled in MC
            #    prescale = 2263.552/45.941
            #else:      # use HLT_Ele23_WPLoose_Gsf_v* TODO this is actually Ele27
            #    prescale = 1.
        elif obj=='m':
            if pt<20: # use HLT_Mu8_TrkIsoVVL_v*
                prescale = 2263.552/1.330
            #elif pt<25: # use HLT_Mu17_TrkIsoVVL_v*
            else: # use HLT_Mu17_TrkIsoVVL_v*
                prescale = 2263.552/197.362
            #else: # use HLT_IsoMu20_v*
            #    prescale = 1.
        else:
            prescale = 1.
        return prescale

    ###########################
    ### Define preselection ###
    ###########################
    def cutTreeSelections(self):
        cutTree = CutTree()
        cutTree.add(self.returnTrue,'topology')
        cutTree.add(self.trigger,'trigger')
        cutTree.add(self.fiducial,'fiducial')
        cutTree.add(self.ID_loose,'looseID')
        cutTree.add(self.ID_medium,'mediumID')
        cutTree.add(self.ID_tight,'tightID')
        cutTree.add(self.zVeto,'zVeto')
        cutTree.add(self.jPsiVeto,'jPsiVeto')
        cutTree.add(self.wVeto,'wVeto')
        cutTree.add(self.jetSelection,'jetSelection')
        cutTree.add(self.veto,'veto2ndLepton')
        return cutTree

    def veto(self,rtrow):
        veto = (rtrow.eVetoLoose + rtrow.muVetoLoose == 0)
        if self.metShift and not self.isData:
            vetoMap = {
                'ees+' : (rtrow.eVetoLoose_eesUp   + rtrow.muVetoLoose == 0),
                'ees-' : (rtrow.eVetoLoose_eesDown + rtrow.muVetoLoose == 0),
                'mes+' : (rtrow.eVetoLoose         + rtrow.muVetoLoose_mesUp == 0),
                'mes-' : (rtrow.eVetoLoose         + rtrow.muVetoLoose_mesDown == 0),
            }
            if self.metShift in vetoMap:
                veto = vetoMap[self.metShift]
        return veto

    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        #cuts.add(self.ID_tight)
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
        #etriggers = ['singleE_leg2Pass', 'singleE_leg1Pass', 'singleEPass'] # TODO this is actually Ele27
        #etriggers = ['singleE_leg2Pass', 'singleE_leg1Pass'] # leg1 not filled in MC
        etriggers = ['singleE_leg2Pass','singleE_leg2Pass']
        #etriggers = ['singleE_leg2Pass']
        #mtriggers = ['singleMu_leg2Pass', 'singleMu_leg1Pass', 'singleIsoMu20Pass']
        mtriggers = ['singleMu_leg2Pass', 'singleMu_leg1Pass']
        #mtriggers = ['singleMu_leg2Pass']
        obj = self.objCand[0]

        triggers = etriggers if obj=='e' else mtriggers

        trigger = triggers[0]
        if getattr(rtrow, '%sPt' %self.objects[0])>=20: trigger = triggers[1]

        if getattr(rtrow,trigger)>0: return True
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

    def ID_loose(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Loose'))

    def ID_medium(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Medium'))

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))

    def zVeto(self,rtrow):
        return getattr(rtrow,'%sNearestZMass' %self.objCand[0]) > 30.

    def jPsiVeto(self,rtrow):
        return getattr(rtrow,'%sLowestMll' %self.objCand[0]) > 20.

    def wVeto(self,rtrow):
        leps = self.objCand
        #if self.getObject(rtrow,'pt','met') > 20.: return False
        #if self.getObject(rtrow, "mt", leps[0], 'met') > 20.: return False
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
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
