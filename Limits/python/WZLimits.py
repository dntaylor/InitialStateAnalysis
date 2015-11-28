import sys
import os
import glob
import logging
import numpy as np
import ROOT

from .datacard import Datacard
import InitialStateAnalysis.Plotters.xsec as xsec
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, getChannels, getSigMap, getIntLumiMap, getMergeDict, getChannelBackgrounds
from InitialStateAnalysis.Utilities.utilities import python_mkdir

ROOT.gROOT.SetBatch(ROOT.kTRUE)

class WZLimits(object):

    def __init__(self, analysis, region, period, selection, ntuple_dir, out_dir, **kwargs):
        scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
        self.analysis = analysis
        self.region = region
        self.period = period
        self.selection = selection
        self.ntuple_dir = ntuple_dir
        python_mkdir(out_dir)
        self.out_dir = out_dir
        self.scalefactor = scalefactor
        self.datacard = Datacard(self.analysis)
        self.sample_groups = {}
        self.log = logging.getLogger(__name__)


    def add_systematics(self, syst_name, syst_type, **kwargs):
        self.log.debug('Adding systematic %s type %s' % (syst_name, syst_type))
        self.datacard.add_syst(syst_name, syst_type, **kwargs)

    def gen_card(self, file_name, **kwargs):
        doDataDriven = kwargs.pop('doDataDriven',True)

        # get the plotter
        nl = 3 
        sigMap = getSigMap(nl)
        intLumiMap = getIntLumiMap()
        mergeDict = getMergeDict(self.period)
        channelBackground = getChannelBackgrounds(self.period)
        channels, leptons = getChannels(nl)
        saves = '%s_%s_%iTeV' % (self.analysis,self.region,self.period)

        plotter = Plotter(self.region,ntupleDir=self.ntuple_dir,saveDir=saves,period=self.period,rootName='plots_limits_wz',mergeDict=mergeDict,scaleFactor=self.scalefactor,datadriven=True)
        plotter.initializeBackgroundSamples([sigMap[self.period][x] for x in channelBackground[self.region+'datadriven']])
        plotter.initializeDataSamples([sigMap[self.period]['data']])
        plotter.setIntLumi(intLumiMap[self.period])

        # set expected and observed yields
        sources = [x for x in channelBackground[self.region+'datadriven'] if x not in ['TT','T','DY','Z','Zfiltered']] + ['datadriven'] if doDataDriven else channelBackground[self.region]
        for bg in sources:
            self.log.info('Processing {0}'.format(bg))
            val = plotter.getNumEntries(self.selection,sigMap[self.period][bg])
            if bg in ['WZ']:
                self.datacard.add_sig(bg,val)
            else:
                self.datacard.add_bkg(bg,val)
        self.log.info('Processing observed')
        self.datacard.set_observed(plotter.getDataEntries(self.selection))

        # write out the datacard
        self.log.info("Saving card to file")

        with open("%s/%s" % (self.out_dir, file_name), 'w') as outfile:
            outfile.write(self.datacard.dump())

