#!/usr/bin/env python
'''
The Z analyzer.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

class AnalyzerZ(AnalyzerBase):
    '''
    '''

    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        if not hasattr(self,'channel'): self.channel = 'Z'
        self.period = period
        self.final_states = ['ee','em','mm']
        self.initial_states = ['z1'] # in order of leptons returned in choose_objects
        self.object_definitions = {
            'z1': ['em','em'],
        }
        self.cutflow_labels = ['Trigger','Fiducial','ID','Z Selection','W Selection']
        self.doVBF = (period==13)
        super(AnalyzerZ, self).__init__(sample_name, file_list, out_file, period,**kwargs)

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
            mass = getattr(rtrow, "%s_%s_Mass" % (l[0], l[1]))
            massdiff = abs(ZMASS-mass)

            ordList = [l[1], l[0]] if getattr(rtrow,'%sPt' % l[0]) < getattr(rtrow,'%sPt' % l[1]) else l

            if OS1 and l[0][0]==l[1][0]:
                cands.append((massdiff, list(ordList)))

        if not len(cands): return 0

        # Sort by mass difference
        cands.sort(key=lambda x: x[0])
        massdiff, leps = cands[0]

        return ([massdiff], leps)


    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
        cuts.add(self.fiducial)
        cuts.add(self.ID_tight)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        cuts.add(self.trigger)
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
            if self.period==8:
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
            if self.period==8:
                kwargs['idDef']['e'] = 'WZLoose'
                kwargs['idDef']['m'] = 'WZLoose'
        if type=='Veto':
            kwargs['idDef'] = {
                'e':'Veto',
                'm':'Loose',
                't':'Loose'
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
        if self.period == 8:
            triggers = ["mu17ele8isoPass", "mu8ele17isoPass",
                        "doubleETightPass", "doubleMuPass", "doubleMuTrkPass"]

        if self.period == 13:
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

    def ID_loose(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Loose'))

    def ID_tight(self, rtrow):
        return self.ID(rtrow,*self.objects,**self.getIdArgs('Tight'))


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

    if args.analyzer == 'Z': analyzer = AnalyzerZ(args.sample_name,args.file_list,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
