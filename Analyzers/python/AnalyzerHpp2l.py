#!/usr/bin/env python
'''
The doubly charged higgs associated production analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerHpp2l(AnalyzerBase):
    '''
    Dummy analyzer for 2 lepton final states.

    Objects:
        doubly charged higgs and signly charged higgs
        h1: H++ = 2 same sign leptons (not necessarily same flavor)
    Minimization function:
        By sT
    Selection:
        Trigger: double lepton
        Fiducial
        Trig threshold
        QCD suppression: M(ll) > 12
    '''

    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        runTau = kwargs.pop('runTau',False)
        #runTau=True
        self.channel = 'Hpp2l'
        self.final_states = ['ee','em','mm'] # no tau
        if runTau: self.final_states = ['ee','em','et','mm','mt','tt']
        self.initial_states = ['h1']
        self.other_states = []
        self.object_definitions = {
            'h1': ['em','em'],
        }
        if runTau:
            self.object_definitions['h1'] = ['emt', 'emt']
        self.lepargs = {'tight':True}
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','QCD Suppression']
        super(AnalyzerHpp2l, self).__init__(sample_name, file_list, out_file, period, **kwargs)

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

            if SS1:
                st = getattr(rtrow, "%sPt" % l[0]) + getattr(rtrow, "%sPt" % l[1])
                ordList = [l[1], l[0]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l
                cands.append([[1./st],list(ordList)]) # choose highest st

        if not len(cands): return 0

        cands.sort(key=lambda x: x[0])

        return (cands[0])

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.qcd_rejection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        #if self.isData: cuts.add(self.trigger)
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
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
            if self.period==8:
                kwargs['idDef']['e'] = 'WZTight'
                kwargs['idDef']['m'] = 'WZTight'
                #kwargs['idDef']['e'] = '4l'
                #kwargs['idDef']['m'] = '4l'
                #kwargs['isoCut']['e'] = 0.4
                #kwargs['isoCut']['m'] = 0.4
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
                #kwargs['idDef']['e'] = 'WZLoose'
                #kwargs['idDef']['m'] = 'WZLoose'
                kwargs['idDef']['e'] = '4l'
                kwargs['idDef']['m'] = '4l'
                kwargs['isoCut']['e'] = 0.4
                kwargs['isoCut']['m'] = 0.4
        return kwargs

    def trigger(self, rtrow):
        triggers = ["mu17ele8isoPass", "mu8ele17isoPass",
                    "doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == 13:
            triggers = ['muEPass', 'eMuPass', 'doubleMuPass',
                        'doubleEPass', 'tripleEPass']

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

    def overlap(self,rtrow):
        for l in permutations(self.objects):
            if lep_order(l[0],l[1]):
                continue
            dr = getattr(rtrow, '%s_%s_DR' % (l[0],l[1]))
            if dr < 0.1: return False
        return True


#######################
### Control regions ###
#######################
class AnalyzerHpp2l_Z(AnalyzerHpp2l):
    '''
    Z control region for Hpp2l
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerHpp2l_Z, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'Z'
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','Mass 3l','Z selection','W selection']
        self.initial_states = ['z1'] # in order of leptons returned in choose_objects
        self.other_states = []
        self.object_definitions = {
            'z1': ['em','em'],
        }

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select leptons that best fit the Z selection.
        The first two leptons are the Z.
        Z are then ordered in pt.
        We select combinatorics by closest to zmass.
        '''
        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            OS1 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) < 0.5 # select opposite sign
            SF1 = l[0][0]==l[1][0] # select same flavor
            mass = getattr(rtrow, "%s_%s_Mass" % (l[0], l[1]))
            massdiff = abs(ZMASS-mass)

            ordList = [l[1], l[0]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l

            if OS1 and SF1:
                cands.append((massdiff, list(ordList)))

        if not len(cands): return 0

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        massdiff, leps = cands[0]

        return ([massdiff], leps)

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        cuts.add(self.metVeto)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        #if self.isData: cuts.add(self.trigger)
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        cuts.add(self.metVeto)
        return cuts

    def zSelection(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o[0],o[1]))
        l0Pt = getattr(rtrow,'%sPt' %leps[0])
        return abs(m1-ZMASS)<20. and l0Pt>20.

    def metVeto(self,rtrow):
        if self.period==8:
            if rtrow.type1_pfMetEt > 30.: return False
        else:
            if rtrow.pfMetEt > 30.: return False
        return True

class AnalyzerHpp2l_Charge(AnalyzerHpp2l):
    '''
    Charge ID test control region for Hpp2l
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerHpp2l_Charge, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'Charge'
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','Mass 3l','Z selection','W selection']
        self.initial_states = ['z1'] # in order of leptons returned in choose_objects
        self.other_states = []
        self.object_definitions = {
            'z1': ['em','em'],
        }

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select leptons that best fit the Z selection.
        The first two leptons are the Z.
        Z are then ordered in pt.
        We select combinatorics by closest to zmass.
        '''
        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]):
                continue

            SF1 = l[0][0]==l[1][0] # select same flavor
            mass = getattr(rtrow, "%s_%s_Mass" % (l[0], l[1]))
            massdiff = abs(ZMASS-mass)

            ordList = [l[1], l[0]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l

            if SF1:
                cands.append((massdiff, list(ordList)))

        if not len(cands): return 0

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        massdiff, leps = cands[0]

        return ([massdiff], leps)

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        #if self.isData: cuts.add(self.trigger)
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        return cuts

    def zSelection(self,rtrow):
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o[0],o[1]))
        l0Pt = getattr(rtrow,'%sPt' %leps[0])
        return m1>60. and m1<120. and l0Pt>20.


class AnalyzerHpp2l_TT(AnalyzerHpp2l):
    '''
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerHpp2l_TT, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'TT'
        self.final_states = ['ee','em','mm'] 
        self.initial_states = ['z1']
        self.other_states = []
        self.object_definitions = {
            'z1': ['em','em'],
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

            OS1 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) < 0.5 # select opposite sign

            if OS1:
                st = getattr(rtrow, "%sPt" % l[0]) + getattr(rtrow, "%sPt" % l[1])
                ordList = [l[1], l[0]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l
                cands.append([[1./st],list(ordList)]) # choose highest st

        if not len(cands): return 0

        cands.sort(key=lambda x: x[0])

        return cands[0]

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.qcd_rejection)
        cuts.add(self.z_veto)
        cuts.add(self.metCut)
        cuts.add(self.jetCut)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        #if self.isData: cuts.add(self.trigger)
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.qcd_rejection)
        cuts.add(self.z_veto)
        cuts.add(self.metCut)
        cuts.add(self.jetCut)
        return cuts

    def z_veto(self,rtrow):
        '''Select Z candidate'''
        leps = self.objCand
        o = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o[0], o[1]))
        return abs(m1-ZMASS)>20.

    def metCut(self,rtrow):
        leps = self.objCand
        if self.period==8:
            if rtrow.type1_pfMetEt < 30.: return False
        else:
            if rtrow.pfMetEt < 30.: return False
        return True

    def jetCut(self,rtrow):
        return rtrow.jetVeto30 > 1


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

    if args.analyzer=='Hpp2l': analyzer = AnalyzerHpp2l(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer=='Z': analyzer = AnalyzerHpp2l_Z(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer=='Charge': analyzer = AnalyzerHpp2l_Charge(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer=='TT': analyzer = AnalyzerHpp2l_TT(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
