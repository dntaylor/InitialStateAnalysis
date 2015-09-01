#!/usr/bin/env python
'''
An analyzer template.

Author: Devin N. Taylor, UW-Madison
'''

from AnalyzerBase import *

# TODO: give unique name for analyzer (must match file name). just find replace LABEL
class AnalyzerLABEL(AnalyzerBase):
    '''
    A template for an analysis in ISA.

    Objects:
        define in object_definitions
    Minimization function:
        define in choose_objects
    Loose Selection:
        define in preselection
    Tight selection:
        define in selection
    Tight ID:
        define in getIdArgs
    '''

    def __init__(self, sample_location, out_file, period, **kwargs):
        # TODO: fill out parameters for analysis
        self.channel = 'LABEL'
        self.final_states = [] # FSA final states
        self.initial_states = [] # custom names for initial states, in order of objects returned in choose_objects
        #self.other_states = [[]] # optional additional states to build in ntuple, requires definition os choose_alternative_objects
        self.object_definitions = { # definition of initial_states above, can use 'em', 'emt', 'g' (photon), 'j', 'n' (met)
        }
        self.cutflow_labels = [] # optional, length must match preselection cuts
        super(AnalyzerLABEL, self).__init__(sample_location, out_file, period, **kwargs)

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

            # TODO: define selection here

            # cands.append(([minimization variables], list(l)))

        if not len(cands): return 0

        # Sort candidates
        cands.sort(key=lambda x: x[0])

        # return selected candidate in form ([minimizing variables list],list(l))

    # optional: override choose_alternative_objects
    # def choose_alternative_objects(self, rtrow, state):
    #     if state == [ some alternative state ]

    ###########################
    ### Define preselection ###
    ###########################
    def preselection(self,rtrow):
        cuts = CutSequence()
        # TODO: define cuts
        #cuts.add(function)
        return cuts

    def selection(self,rtrow):
        cuts = CutSequence()
        # TODO: define cuts
        #cuts.add(function)
        return cuts

    def getIdArgs(self,type):
        kwargs = {}
        if type=='Tight':
           # TODO: define tight ID (with isolation)
           pass
        return kwargs

##########################
###### Command line ######
##########################
def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('in_sample', type=str)
    parser.add_argument('out_file', type=str)
    parser.add_argument('period', type=str)

    args = parser.parse_args(argv)
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    analyzer = AnalyzerLABEL(args.in_sample,args.out_file,args.period)
    with analyzer as thisAnalyzer:
        thisAnalyzer.analyze()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
