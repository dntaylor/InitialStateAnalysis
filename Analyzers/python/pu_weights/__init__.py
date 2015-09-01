import sys
import os
import json

from math import floor

sys.argv.append('-b')
import ROOT as rt
sys.argv.pop()


class PileupWeights(object):

    def __init__(self):
        path = os.path.join(os.path.dirname(__file__), 'pu_weights.json')
        with open(path, 'r') as pu_file:
            self.pu_weights = json.load(pu_file)

    def weight(self, rtrow):
        if rtrow.nTruePU < 0:
            return 1
        else:
            return self.pu_weights[str(int(floor(rtrow.nTruePU)))]
