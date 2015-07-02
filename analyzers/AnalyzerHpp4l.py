#!/usr/bin/env python
'''
The doubly charged higgs associated production analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerHpp4l(AnalyzerBase):
    '''
    The primary analyzer for the doubly charged higgs associated production channel.

    Objects:
        two doubly charged higgs
        h1: H++ = 2 same sign leptons (not necessarily same flavor)
        h2: H-- = "
    Minimization function:
        abs(Mass(++)-Mass(--))
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
        #runTau = True
        self.channel = 'Hpp4l'
        self.final_states = ['eeee','eeem','eemm','emmm','mmmm'] # no tau
        if runTau: self.final_states = ['eeee','eeem','eeet','eemm','eemt','eett','emmm','emmt','emtt','ettt',\
                                        'mmmm','mmmt','mmtt','mttt','tttt']
        self.initial_states = ['h1','h2']
        self.other_states = [['z1','z2']]
        self.object_definitions = {
            'h1': ['em','em'],
            'h2': ['em','em'],
            'z1': ['em','em'],
            'z2': ['em','em'],
        }
        if runTau:
            self.object_definitions['h1'] = ['emt', 'emt']
            self.object_definitions['h2'] = ['emt', 'emt']
            self.object_definitions['z1'] = ['emt', 'emt']
            self.object_definitions['z2'] = ['emt', 'emt']
        self.lepargs = {'loose':True}
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','QCD Suppression']
        super(AnalyzerHpp4l, self).__init__(sample_name, file_list, out_file, period)

    ###############################
    ### Define Object selection ###
    ###############################
    def choose_objects(self, rtrow):
        '''
        Select candidate objects
        '''
        cands = []
        for l in permutations(self.objects):
            if lep_order(l[0], l[1]) or lep_order(l[2], l[3]):
                continue

            SS1 = getattr(rtrow, "%s_%s_SS" % (l[0], l[1])) > 0 # select same sign
            mass1 = getattr(rtrow, "%s_%s_Mass" % (l[0], l[1])) # select mass
            SS2 = getattr(rtrow, "%s_%s_SS" % (l[2], l[3])) > 0 # select same sign
            mass2 = getattr(rtrow, "%s_%s_Mass" % (l[2], l[3])) # select mass
            OS = getattr(rtrow, "%sCharge" % l[0]) != getattr(rtrow, "%sCharge" % l[2]) # select opposite sign
            massdiff = abs(mass1-mass2)

            if SS1 and SS2 and OS:
                cands.append([massdiff,list(l)]) # minimization is by mass diff

        if not len(cands): return 0

        cands.sort(key=lambda x: x[0])
        massdiff, leps = cands[0]

        return ([massdiff], leps)

    # override choose_alternative_objects
    def choose_alternative_objects(self, rtrow, state):
        '''
        Select alternative candidate.
        '''
        # ZZ
        if state == ['z1', 'z2']:
            bestZDiff = float('inf')
            bestSt = 0
            bestLeptons = []

            for l in permutations(self.objects):
                if lep_order(l[0],l[1]) or lep_order(l[2],l[3]):
                    continue

                os1 = getattr(rtrow,'%s_%s_SS' % (l[0], l[1])) < 0.5
                m1 = getattr(rtrow,'%s_%s_Mass' % (l[0], l[1]))
                os2 = getattr(rtrow,'%s_%s_SS' % (l[2], l[3])) < 0.5
                st2 = getattr(rtrow,'%sPt' %l[2]) + getattr(rtrow,'%sPt' %l[3])

                if l[0][0] == l[1][0] and os1 and l[2][0]==l[3][0] and os2:
                    if abs(m1-ZMASS) < bestZDiff:
                        bestZDiff = abs(m1-ZMASS)
                        bestSt = st2
                        bestLeptons = l
                    elif abs(m1-ZMASS)==bestZDiff and st2>bestSt:
                        bestZDiff = abs(m1-ZMASS)
                        bestSt = st2
                        bestLeptons = l

            return bestLeptons


    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_loose)
        cuts.add(self.qcd_rejection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
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
            if self.period=='8':
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
            if self.period=='8':
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
        return pts[0] > 25.0 and pts[1] > 15.0

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
class AnalyzerHpp4l_ZZ(AnalyzerHpp4l):
    '''
    ZZ control region for Hpp4l
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        super(AnalyzerHpp4l_ZZ, self).__init__(sample_name, file_list, out_file, period, **kwargs)
        self.channel = 'ZZ'
        self.cutflow_labels = ['Trigger','Fiducial','Trigger Threshold','ID','Mass 3l','Z selection','W selection']

    def preselection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_loose)
        cuts.add(self.zSelection)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        if self.isData: cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.overlap)
        cuts.add(self.trigger_threshold)
        cuts.add(self.ID_tight)
        cuts.add(self.zSelection)
        return cuts

    def zSelection(self,rtrow):
        leps = self.choose_alternative_objects(rtrow, ['z1','z2'])
        if not leps: return False
        o1 = ordered(leps[0], leps[1])
        m1 = getattr(rtrow,'%s_%s_Mass' % (o1[0],o1[1]))
        o2 = ordered(leps[2], leps[3])
        m2 = getattr(rtrow,'%s_%s_Mass' % (o2[0],o2[1]))
        l0Pt = getattr(rtrow,'%sPt' %leps[0])
        return abs(m1-ZMASS)<30. and abs(m2-ZMASS)<30. and l0Pt>20.


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

    if args.analyzer == 'Hpp4l': analyzer = AnalyzerHpp4l(args.sample_name,args.file_list,args.out_file,args.period)
    if args.analyzer == 'ZZ': analyzer = AnalyzerHpp4l_ZZ(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
