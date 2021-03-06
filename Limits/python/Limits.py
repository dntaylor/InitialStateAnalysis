import sys
import os
import glob
import logging
import numpy as np
import ROOT

from .datacard import Datacard
import InitialStateAnalysis.Plotters.xsec as xsec
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import _3L_MASSES, _4L_MASSES, ZMASS, getChannels, getSigMap, getIntLumiMap, getMergeDict, getChannelBackgrounds, getNtupleDirectory

ROOT.gROOT.SetBatch(ROOT.kTRUE)

class Limits(object):

    def __init__(self, analysis, region, period, base_selections, ntuple_dir, out_dir,
                 channels=[], lumi=25.0, blinded=True, bgMode='mc', scalefactor='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',
                 sbcut='1', srcut='1'):
        self.base_selections = base_selections
        self.analysis = analysis
        self.region = region
        self.period = period
        self.out_dir = out_dir
        self.ntuple_dir = ntuple_dir
        self.sample_groups = {}
        self.lumi = lumi
        self.blinded = blinded
        self.sytematic_list = []
        self.datacard = Datacard(self.analysis)
        self.channels = channels
        self.xsecs = xsec.xsecs[period]
        self.bgMode = bgMode
        self.scalefactor = scalefactor
        self.sbcut = sbcut
        self.srcut = srcut

        self.log = logging.getLogger(__name__)

        os.system("mkdir -p %s" % self.out_dir)

    def getPlotter(self,analysis,region,runPeriod,mass,runTau,plotName,doFakes):
        nl = 3 if analysis=='Hpp3l' or analysis=='WZ' else 4
        ntuples = getNtupleDirectory(analysis,region,runPeriod)
        saves = '%s_%s_%iTeV' % (analysis,region,runPeriod)
        sigMap = getSigMap(nl,mass)
        intLumiMap = getIntLumiMap()
        mergeDict = getMergeDict(runPeriod)
        regionBackground = getChannelBackgrounds(runPeriod)
        channels, leptons = getChannels(nl,runTau=runTau)
    
        plotter = Plotter(region,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName=plotName,mergeDict=mergeDict,scaleFactor=self.scalefactor)
        if not doFakes: plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in regionBackground[analysis]])
        if runPeriod==8: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
    
        return plotter


    def add_group(self, group_name, *sample_names, **kwargs):
        samples = []
        for name in sample_names:
            samples += glob.glob("%s/%s.root" % (self.ntuple_dir, name))

        is_data = kwargs.get('isData', False)
        is_sig = kwargs.get('isSignal', False)
        scale = kwargs.get('scale', 1.0)

        self.sample_groups[group_name] = {
            'sample_names': [os.path.splitext(os.path.basename(x))[0]
                             for x in samples],
            'scale': scale,
            'isSig': is_sig,
            'isData': is_data}
        message = 'Added %s with scales %s' % (group_name, str(scale))
        if is_sig: message += ' (signal)'
        if is_data: message += ' (data)'
        self.log.debug(message)
        self.log.debug('Samples: %s' % str(self.sample_groups[group_name]['sample_names']))

    def add_systematics(self, syst_name, syst_type, **kwargs):
        self.log.debug('Adding systematic %s type %s' % (syst_name, syst_type))
        self.datacard.add_syst(syst_name, syst_type, **kwargs)

    def get_var_weights(self, sample, var, cut, scale):
        file = self.ntuple_dir+'/%s.root' % sample
        tfile = ROOT.TFile(file)
        cutflowHist = tfile.Get('cutflow')
        n_evts = cutflowHist.GetBinContent(1)
        sample_xsec = self.xsecs[sample]
        samplelumi = float(n_evts)/sample_xsec
        lumiscale = self.lumi/samplelumi
        val = 0
        tree = tfile.Get(self.region)
        #print tree.GetEntries()
        tree.Draw('1>>h%s(1,0,2)'%sample,'%s*(%s)' %(self.scalefactor,cut),'goff')
        if not ROOT.gDirectory.Get("h%s" %sample): return [0,0]
        hist = ROOT.gDirectory.Get("h%s" %sample).Clone("hnew%s" %sample)
        hist.Sumw2()
        val = hist.Integral()
        err = 0
        for i in range(hist.GetNbinsX()+1):
            err += hist.GetBinError(i)
        val = val * scale * lumiscale
        err = err * scale * lumiscale
        #print "val: %f" % val
        #print "err: %f" % err
        #if val < err: return err
        return [val, err]

    def gen_card(self, file_name, **kwargs):
        mass = kwargs.pop('mass',500)
        cuts = kwargs.pop('cuts','1')
        period = kwargs.pop('period',8)
        doAlphaTest = kwargs.pop('doAlphaTest',False)
        values = []
        weights = []

        chans = cuts

        var = 'mass'

        if isinstance(cuts, str):
            cut = self.base_selections + ' && ' + cuts
            cutMC_data = cut
            cutSig = cut
        else:
            cutMC_data = self.base_selections + ' && (' + ' || '.join(cuts) + ')' # full selection for bg and data
            cutSig = [self.base_selections + ' && ' + x for x in cuts]          # different selections or signal so it can be scaled
            cuts = '(' + ' || '.join(cuts) + ')'
        if self.region=='WZ': cuts += ' && select.PassTight'

        myCut = '1'

        plotter = self.getPlotter(self.analysis,self.region,self.period,mass,False,'plots_limits_temp',False)

        nSBDict = {}
        nSRDict = {}
        for background in plotter.backgrounds:
            nSBDict[background] = plotter.getNumEntries(self.sbcut,background,doError=True)
            nSRDict[background] = plotter.getNumEntries(self.srcut,background,doError=True)
        nSB = sum([x[0] for x in nSBDict.itervalues()])
        eSB = sum([x[1]*x[1] for x in nSBDict.itervalues()]) ** 0.5
        nSR = sum([x[0] for x in nSRDict.itervalues()])
        eSR = sum([x[1]*x[1] for x in nSRDict.itervalues()]) ** 0.5
        alpha = nSR/nSB if nSB else nSR
        # TODO: Check this
        if period == 13:
            nSBData, eSBData = (0., 0.)
            nSRData, eSRData = (0., 0.)
        else:
            nSBData, eSBData = plotter.getNumEntries(self.sbcut,*plotter.data,doError=True)
            nSRData, eSRData = plotter.getNumEntries(self.srcut,*plotter.data,doError=True)
        if self.blinded:
            nSRData = 0
            eSRData = 1
        nBGSR = alpha*(nSBData)
        eBGSR = alpha*((nSBData) ** 0.5)
        self.alpha = alpha
        self.nSBData = nSBData
        self.nBGSR = nBGSR

        # store values for values card
        sbVal = nBGSR
        sbStatErr = eBGSR
        dataVal = nSRData
        dataStatErr = eSRData
        dataSystErr = 0.0 # TODO
        mcVal = 0.0 # TODO
        mcStatErr = 0.0 # TODO
        pairVal = 0.0 # TODO
        pairStatErr = 0.0 # TODO
        assocVal = 0.0 # TODO
        assocStatErr = 0.0 # TODO

        # if nBGSR < stat uncertainty of MC in SR, nBGSR = stat uncertainty
        if nBGSR < eBGSR: nBGSR = eBGSR
        # print this out just to check
        #print "nSBMC: %0.4f" % nSB
        #print "eSBMC: %0.4f" % eSB
        #print "nSRMC: %0.4f" % nSR
        #print "eSRMC: %0.4f" % eSR
        #print "alpha: %0.4f" % alpha
        #print "nSBDa: %0.4f" % nSBData
        #print "eSBDa: %0.4f" % eSBData
        #print "nBGSR: %0.4f" % nBGSR
        #print "eBGSR: %0.4f" % eBGSR
        #if not self.blinded: print "nSRDa: %0.4f" % nSRData
        strToSave = ":".join(["%0.10f" % x for x in [nSB,eSB,nSR,eSR,alpha,nSBData,eSBData,nBGSR,eBGSR,nSRData,eSRData]])

        with open('%s/alphavalues_%s' % (self.out_dir, file_name), 'w') as file:
            file.write(strToSave)

        # calculate values
        bgMap = {}
        for key in self.sample_groups:
           is_data = self.sample_groups[key]['isData']
           scale = self.sample_groups[key]['scale']
           wgts = []
           if not is_data:
               self.log.info("Processing MC: %s" % key)
               for sample_name in self.sample_groups[key]["sample_names"]:
                   if type(scale) is list:
                       for s,c in zip(scale,cutSig):
                           self.log.debug("Scale applied: %f" % s)
                           self.log.debug("Cut applied: %s" % c)
                           thisWeight = self.get_var_weights(sample_name,var,c,s)
                           self.log.debug("Weight: %f +/- %f" % (thisWeight[0],thisWeight[1]))
                           wgts.append(thisWeight)
                   else:
                       self.log.debug("Scale applied: %f" % scale)
                       self.log.debug("Cut applied: %s" % cutMC_data)
                       thisWeight = self.get_var_weights(sample_name,var,cutMC_data,scale)
                       self.log.debug("Weight: %f +/- %f" % (thisWeight[0],thisWeight[1]))
                       wgts.append(thisWeight)
               bgMap[key] = [sum([x[0] for x in wgts]), sum([x[1]**2 for x in wgts])**0.5]

        for key in self.sample_groups:
            if self.sample_groups[key]['isSig']:
                val = bgMap[key]
                if 'AP' == key[-2:]:
                    assocVal = val[0]
                    assocStatErr = val[1]
                if 'PP' in key[-2:]:
                    pairVal = val[0]
                    pairStatErr = val[1]

        mcVal, mcStatErr = plotter.getBackgroundEntries(cutMC_data,doError=True)

        # here we decide what datacard format we want to output
        if self.bgMode=='sideband':
            # add the systematics for the alpha
            alphaSys = {'bg':self.alpha}
            self.add_systematics("alpha_%s" % file_name.split('.')[0],"gmN %i" %(self.nSBData), **alphaSys)
            for key in self.sample_groups:
                is_data = self.sample_groups[key]['isData']
                if self.sample_groups[key]['isSig']:
                    self.datacard.add_sig(key, bgMap[key][0])
                elif key=='bg': # don't add the mc bg
                    self.datacard.add_bkg(key, nBGSR)
                elif is_data:
                    self.datacard.set_observed(nSRData)

        if self.bgMode=='mc':
            for key in self.sample_groups:
                is_data = self.sample_groups[key]['isData']
                if not is_data:
                    if self.sample_groups[key]['isSig']:
                        self.datacard.add_sig(key, bgMap[key][0])
                    elif key!='bg': # dont add the sideband bg
                        self.datacard.add_bkg(key, bgMap[key][0])
                elif is_data:
                    self.datacard.set_observed(nSRData)

        self.log.info("Saving card to file")

        if doAlphaTest: file_name = 'alpha_' + file_name
        with open("%s/%s" % (self.out_dir, file_name), 'w') as outfile:
            outfile.write(self.datacard.dump())

        # output this to text file to be read later
        # format mcVal:mcStatErr:mcSystErr:sbVal:sbStatErr:sbSystErr:dataVal:dataStatErr:dataSystErr:\
        #        pairVal:pairStatErr:pairSystErr:assocVal:assocStatErr:assocSystErr
        with open('%s/values_%s' % (self.out_dir, file_name), 'w') as file:
            outString = '%f:%f' % (mcVal, mcStatErr)
            outString += ':%f:%f' % (sbVal, sbStatErr)
            if not self.blinded:
                outString += ':%f:%f:%f' % (dataVal, dataStatErr, dataSystErr)
            else:
                outString += ':%f:%f:%f' % (0, 0, 0)
            outString += ':%f:%f' % (pairVal, pairStatErr)
            outString += ':%f:%f' % (assocVal, assocStatErr)
            file.write(outString)
