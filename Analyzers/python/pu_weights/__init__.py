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
        path = os.path.join(os.path.dirname(__file__), 'pileup_RunIISpring2015_Run2015D_13TeV.root')
        self.scale_13tev = {}
        self.scale_13tev_up = {}
        self.scale_13tev_down = {}
        rootfile = rt.TFile(path)
        hist_scale = rootfile.Get('pileup_scale')
        for b in range(hist_scale.GetNbinsX()):
            self.scale_13tev[b] = hist_scale.GetBinContent(b+1)
        hist_scale = rootfile.Get('pileup_scale_up')
        for b in range(hist_scale.GetNbinsX()):
            self.scale_13tev_up[b] = hist_scale.GetBinContent(b+1)
        hist_scale = rootfile.Get('pileup_scale_down')
        for b in range(hist_scale.GetNbinsX()):
            self.scale_13tev_down[b] = hist_scale.GetBinContent(b+1)
        rootfile.Close()


    def weight(self, rtrow, **kwargs):
        period = kwargs.pop('period',8)
        if rtrow.nTruePU < 0:
            return [1,1,1]
        else:
            if period==8:
                val = self.pu_weights[str(int(floor(rtrow.nTruePU)))]
                return [val,val,val]
            else:
                val = self.scale_13tev[int(floor(rtrow.nTruePU))]
                up = self.scale_13tev_up[int(floor(rtrow.nTruePU))]
                down = self.scale_13tev_down[int(floor(rtrow.nTruePU))]
                return [val,up,down]
