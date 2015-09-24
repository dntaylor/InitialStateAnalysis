'''
Base plotting class using ntuples produced in the DblHAnalysis framework.
This class requires the definition of luminosity files claculated with lumiCalcy2.py, a cross section file,
and a data styles file for each mc sample. The xsecs should be defined in a dictionary found in xsec.py. 
The data styles should be defined in a disctionary in dataStyles.py. 

Author: Devin N. Taylor, UW-Madison
'''

import sys
import os
import errno
import glob
import ROOT
import json
import math
import logging
from array import array

from xsec import xsecs
from dataStyles import dataStyles
import CMS_lumi, tdrstyle
from plotUtils import *

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 2001;")
tdrstyle.setTDRStyle()

class PlotterBase(object):
    '''A Base plotting class to be used with flat histograms.'''
    def __init__(self,analysis,**kwargs):
        '''Initialize the plotter (optionally make the plots blinded).'''
        # get kwargs
        loglevel = kwargs.pop('loglevel','INFO')
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(getattr(logging,loglevel))
        saveDir = kwargs.pop('saveDir','')
        ntupleDir = kwargs.pop('ntupleDir','ntuples')
        period = kwargs.pop('period',13)
        blind = kwargs.pop('blind',False)
        rootName = kwargs.pop('rootName','plots')
        mergeDict = kwargs.pop('mergeDict',{})
        scaleFactor = kwargs.pop('scaleFactor','event.pu_weight*event.lep_scale*event.trig_scale')
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '%s' = %s" %(key,str(value)))

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        # first, setup our canvas
        self.W = 800
        self.H = 600
        self.T = 0.08
        self.B = 0.12
        self.L = 0.12
        self.R = 0.04
        self.canvas = ROOT.TCanvas("c1","c1",50,50,self.W,self.H)
        self.canvas.SetFillColor(0)
        self.canvas.SetBorderMode(0)
        self.canvas.SetFrameFillStyle(0)
        self.canvas.SetFrameBorderMode(0)
        self.canvas.SetLeftMargin( self.L )
        self.canvas.SetRightMargin( self.R )
        self.canvas.SetTopMargin( self.T )
        self.canvas.SetBottomMargin( self.B )
        self.canvas.SetTickx(0)
        self.canvas.SetTicky(0)

        # now, setup plotter conditions (some to be initalized later)
        self.j = 0 # global variable to prevent resusing histograms
        self.backgroundInitialized = False
        self.background = []
        self.dataInitialized = False
        self.data = []
        self.signalInitialized = False
        self.signal = []
        self.analysis = analysis
        self.blind = blind
        self.sqrts=period
        # worry about this later
        self.plot7TeV = self.sqrts==7
        self.plot8TeV = self.sqrts==8
        self.plot13TeV = self.sqrts==13
        self.xsecs = xsecs[self.sqrts]
        self.dataStyles = dataStyles
        self.ntupleDir = ntupleDir
        if saveDir=='': saveDir = self.analysis
        self.plotDir = 'plots/'+saveDir
        python_mkdir(self.plotDir)
        python_mkdir(self.plotDir+'/png')
        self.savefile = ROOT.TFile(self.plotDir+"/"+rootName+".root","recreate")
        self.samples = {}
        self.intLumi = 25000. # just a default 25 fb-1 for plotting without data
        self.sampleMergeDict = mergeDict
        self.scaleFactor = scaleFactor

    def reset(self):
        '''Reset the plotter class'''
        ROOT.gDirectory.Delete('h*') # clear histogram memory
        self.logger.info("Resetting the PlotterBase class")
        self.backgroundInitialized = False
        self.backgrounds = []
        self.dataInitialized = False
        self.data = []
        self.signalInitialized = False
        self.signal = []
        self.samples = {}
        self.intLumi = 25000.
        self.j = 0
        self.resetCanvas()

    def resetCanvas(self):
        '''Reset canvas after changes'''
        self.canvas.SetCanvasSize(self.W,self.H)
        self.canvas.SetLeftMargin( self.L )
        self.canvas.SetRightMargin( self.R )
        self.canvas.SetTopMargin( self.T )
        self.canvas.SetBottomMargin( self.B )
        self.canvas.SetTickx(0)
        self.canvas.SetTicky(0)
        

    def initializeBackgroundSamples(self,sampleList):
        '''Initialize the background samples.'''
        self.backgrounds = sampleList
        self.initializeSamples(sampleList)
        self.backgroundInitialized = True

    def initializeDataSamples(self,sampleList):
        '''Initialize the data samples.'''
        self.data = sampleList
        self.initializeSamples(sampleList)
        self.dataInitialized = True
        self.calculateIntLumi()

    def initializeSignalSamples(self,sampleList):
        '''Initialize the signal samples.'''
        self.signal = sampleList
        self.initializeSamples(sampleList)
        self.signalInitialized = True

    def initializeSamplesHelper(self,sample):
        '''initialize single sample'''
        self.samples[sample] = {} 
        file = self.ntupleDir+'/%s.root' % sample
        self.samples[sample]['file'] = ROOT.TFile(file)
        if 'data' in sample:
            lumifile = self.ntupleDir+'/%s.lumicalc.sum' % sample
        else:
            cutflowHist = self.samples[sample]['file'].Get('cutflow')
            n_evts = cutflowHist.GetBinContent(1)
            sample_xsec = self.xsecs[sample]
            self.samples[sample]['lumi'] = float(n_evts)/sample_xsec
            self.logger.debug('Initializing MC sample %s with %i events and xsec %f to lumi %f.'\
                  % (sample, n_evts, sample_xsec, self.samples[sample]['lumi']))

    def initializeSamples(self,sampleList):
        '''Initialize a list of samples to the sample dictionary.'''
        for sample in sampleList:
            if sample in self.sampleMergeDict:
                for s in self.sampleMergeDict[sample]:
                    self.initializeSamplesHelper(s)
            else:
                self.initializeSamplesHelper(sample)

    def calculateIntLumi(self):
        '''Calculate the integrated luminosity to scale the Monte Carlo'''
        if not self.dataInitialized: 
            self.logger.info("No data initialized, default to 25 fb-1")
            self.intLumi = 25000.
            return
        self.intLumi = 25000.

    def setIntLumi(self,intLumi):
        '''Set the integrated luminosity to scale MC to'''
        self.intLumi = intLumi

    def printInfo(self):
        if self.backgroundInitialized: self.logger.info('Backgrounds: ' + ' '.join(self.backgrounds))
        if self.dataInitialized: self.logger.info('Data: ' + ' '.join(self.data))
        if self.signalInitialized: self.logger.info('Signal: ' + ' '.join(self.signal))

    def setScaleFactor(self,scalefactor):
        '''Set Scale factor'''
        self.scaleFactor = scalefactor

    def getSignalEntries(self,selection,**kwargs):
        signal = kwargs.pop('signal','')
        doError = kwargs.pop('doError',False)
        result = 0
        err2 = 0
        if signal: # plot individual signal
             result,err = self.getNumEntries(selection,signal,doError=True,**kwargs)
        else:
            for sig in self.signal:
                temp = self.getNumEntries(selection,sig,doError=True,**kwargs)
                result += temp[0]
                err2 += temp[1]**2
            err = err2 ** 0.5
        if doError: return result, err
        return result

    def getBackgroundEntries(self,selection,**kwargs):
        doError = kwargs.pop('doError',False)
        result = 0
        err2 = 0
        for bg in self.backgrounds:
            temp = self.getNumEntries(selection,bg,doError=True,**kwargs)
            result += temp[0]
            err2 += temp[1]**2
        err = err2**0.5
        if doError: return result, err
        return result

    def getDataEntries(self,selection,**kwargs):
        doError = kwargs.pop('doError',False)
        result = 0
        err2 = 0
        for d in self.data:
            temp = self.getNumEntries(selection,d,doError=True,**kwargs)
            result += temp[0]
            err2 += temp[1]**2
        err = err2**0.5
        if doError: return result, err
        return result

    def getIndividualSampleEntries(self,sample,selection,scalefactor,**kwargs):
        tree = self.samples[sample]['file'].Get(self.analysis)
        tree.Draw('1>>h%s%i(1,0,2)'%(sample,self.j),'%s*(%s)' %(scalefactor,selection),'goff')
        if not ROOT.gDirectory.Get("h%s%i" %(sample,self.j)):
            self.logger.debug('%s: No entries'%sample)
            val = 0
        else:
            hist = ROOT.gDirectory.Get("h%s%i" %(sample,self.j)).Clone("hnew%s%i" %(sample,self.j))
            hist.Sumw2()
            val = hist.Integral()
            self.logger.debug('%s: Entries: %f'%(sample,val))
        err = val ** 0.5
        if 'data' in sample: return val, err
        lumi = self.samples[sample]['lumi']
        val = val * self.intLumi/lumi
        err = err * self.intLumi/lumi
        self.logger.debug('%s: Lumi scaled: %f'%(sample,val))
        return val, err


    def getNumEntries(self,selection,sample,**kwargs):
        '''Return the lumi scaled number of entries passing a given cut.'''
        doError = kwargs.pop('doError',False)
        scaleup = kwargs.pop('scaleup',False)
        unweighted = kwargs.pop('doUnweighted',False)
        doDataDriven = kwargs.pop('doDataDriven',False) # a custom weight (for datadriven background)
        customScale = kwargs.pop('customScale','')
        totalVal = 0
        totalErr2 = 0
        scalefactor = "event.lep_scale_up*event.trig_scale*event.pu_weight" if scaleup else self.scaleFactor
        if 'data' in sample:
            scalefactor = '1'
            if doDataDriven: scalefactor = 'event.datadriven_weight'
        if customScale: scalefactor = customScale
        self.j += 1
        self.logger.debug('Cut: %s'%selection)
        if sample in self.sampleMergeDict:
            for s in self.sampleMergeDict[sample]:
                thisCut = selection + ' & ' + self.sampleMergeDict[sample][s] if 'data' not in s else selection
                val, err = self.getIndividualSampleEntries(s,thisCut,scalefactor)
                totalVal += val
                totalErr2 += err*err
        else:
            val, err = self.getIndividualSampleEntries(sample,selection,scalefactor)
            totalVal += val
            totalErr2 += err*err
        totalErr = totalErr2 ** 0.5
        self.logger.debug('Total Integrated: %f'%(totalVal))
        if doError: return totalVal, totalErr
        return totalVal

    def getOverflowUnderflow(self,hist,**kwargs):
        '''Get the plot with overflow and underflow bins'''
        under = kwargs.pop('underflow',False)
        over = kwargs.pop('overflow',False)
        if not under and not over: return hist
        nx = hist.GetNbinsX()
        if under: nx += 1
        if over: nx += 1
        xbins = [0]*(nx+1)
        for i in range(nx):
            xbins[i]=hist.GetBinLowEdge(i+1)
        xbins[nx]=xbins[nx-1]+hist.GetBinWidth(nx)
        tempName = hist.GetName()+'OU%i' % self.j
        htmp = ROOT.TH1F(tempName, hist.GetTitle(), nx, array('d',xbins))
        htmp.Sumw2()
        for i in range(nx):
            htmp.SetBinContent(i+1, hist.GetBinContent(i+1))
            htmp.SetBinError(i+1, hist.GetBinError(i+1))
        htmp.SetBinContent(0, hist.GetBinContent(0))
        htmp.SetBinError(0, hist.GetBinError(0))
        htmp.SetEntries(hist.GetEntries())
        return htmp

    def getSingleVarHist2D(self,sample,var1,var2,bin1,bin2,cut,**kwargs):
        '''Plot a single sample hist with two variables'''
        tree = self.samples[sample]['file'].Get(self.analysis)
        zbin = kwargs.pop('zbin',[10,0,10])
        drawString = "%s:%s>>h%s%s%s(%s)" % (var2,var1,sample,var1,var2,', '.join(str(x) for x in bin1+bin2))
        if not cut: cut = '1'
        if 'data' not in sample:
            tree.Draw(drawString,'%s*(%s)' % (self.scaleFactor,cut),'goff')
        else:
            tree.Draw(drawString,cut,'goff')
        if not ROOT.gDirectory.Get("h%s%s%s" %(sample, var1,var2)):
            return 0
        hist = ROOT.gDirectory.Get("h%s%s%s" %(sample, var1,var2)).Clone("hmod%s%s%s"%(sample,var1,var2))
        hist.Sumw2()
        if 'data' not in sample: # if it is mc, scale to intLumi
            lumi = self.samples[sample]['lumi']
            theScale = float(self.intLumi)/lumi
            hist.Scale(theScale)
            hist.SetMarkerColor(4)
        hist.GetXaxis().SetLimits(bin1[1],bin1[2])
        hist.GetYaxis().SetLimits(bin2[1],bin2[2])
        #hist.GetZaxis().SetRangeUser(zbin[1],zbin[2])
        return hist

    def getHist2D(self, sample, var1, var2, bin1, bin2, cut, **kwargs):
        '''Return a histogram of a given variable from the given dataset with a cut'''
        normalize = kwargs.pop('normalize',False)
        hists = ROOT.TList()
        for v in range(len(var1)):
            if sample in self.sampleMergeDict:
                for s in self.sampleMergeDict[sample]:
                    if len(var1) != len(cut):
                        thisCut = cut + ' & ' + self.sampleMergeDict[sample][s]
                        hist = self.getSingleVarHist2D(s,var1[v],var2[v],bin1,bin2,thisCut,**kwargs)
                    else:
                        thisCut = cut[v] + ' & ' + self.sampleMergeDict[sample][s]
                        hist = self.getSingleVarHist2D(s,var1[v],var2[v],bin1,bin2,thisCut,**kwargs)
                    if hist:
                        hists.Add(hist)
            else:
                if len(var1) != len(cut):
                    hist = self.getSingleVarHist2D(sample,var1[v],var2[v],bin1,bin2,cut,**kwargs)
                else:
                    hist = self.getSingleVarHist2D(sample,var1[v],var2[v],bin1,bin2,cut[v],**kwargs)
                if hist:
                    hists.Add(hist)
        if hists.IsEmpty():
            return 0
        hist = hists[0].Clone("hmerged%s%s%s" % (sample, var1[0], var2[0]))
        hist.Reset()
        hist.Merge(hists)
        hist.SetTitle(self.dataStyles[sample]['name'])
        return hist


    def getSingleVarHist(self,sample,variable,binning,cut):
        '''Single variable, single sample hist'''
        tree = self.samples[sample]['file'].Get(self.analysis)
        self.j += 1
        histname = 'h%s%s' % (sample, variable.replace('(','_').replace(')','_'))
        if len(binning) == 3: # standard drawing
            drawString = "%s>>%s(%s)" % (variable, histname, ", ".join(str(x) for x in binning))
        else: # we will need to rebin TODO: this might cause issue in root6, double check later
            drawString = "%s>>%s()" % (variable, histname)
        if not cut: cut = '1'
        if 'data' not in sample:
            tree.Draw(drawString,'%s*(%s)' % (self.scaleFactor,cut),'goff')
        else:
            tree.Draw(drawString,cut,'goff')
        if not ROOT.gDirectory.Get(histname):
            return 0
        hist = ROOT.gDirectory.Get(histname).Clone(histname+'_mod')
        if len(binning) != 3: # variable binning (list of bin edges
            hist.Rebin(len(binning)-1,histname+'_new',array('d',binning))
            hist = ROOT.gDirectory.Get(histname+'_new').Clone(histname+'_newmod')
        hist.Sumw2()
        if 'data' not in sample: # if it is mc, scale to intLumi
            lumi = self.samples[sample]['lumi']
            theScale = float(self.intLumi)/lumi
            self.logger.debug('Scaling sample %s to %f with sample lumi of %f and scale factor %f.'\
                  % (sample, self.intLumi, lumi, theScale))
            self.logger.debug('The old integral was %f' % hist.Integral())
            hist.Scale(theScale)
            self.logger.debug('The new integral is %s' % hist.Integral())
        else:
            self.logger.debug('For data sample %s the integral is %f.' % (sample, hist.Integral()))
            
        return hist

    def getHist(self, sample, variables, binning, cut, noFormat=False, **kwargs):
        '''Return a histogram of a given variable from the given dataset with a cut'''
        normalize = kwargs.pop('normalize',False)
        hists = ROOT.TList()
        for v in range(len(variables)):
            if sample in self.sampleMergeDict:
                for s in self.sampleMergeDict[sample]:
                    if len(variables) != len(cut):
                        thisCut = cut + ' & ' + self.sampleMergeDict[sample][s]
                        hist = self.getSingleVarHist(s,variables[v],binning,thisCut)
                    else:
                        thisCut = cut[v] + ' & ' + self.sampleMergeDict[sample][s]
                        hist = self.getSingleVarHist(s,variables[v],binning,thisCut)
                    if hist:
                        hists.Add(hist)
            else:
                if len(variables) != len(cut):
                    hist = self.getSingleVarHist(sample,variables[v],binning,cut)
                else:
                    hist = self.getSingleVarHist(sample,variables[v],binning,cut[v])
                if hist:
                    hists.Add(hist)
        if hists.IsEmpty():
            return 0
        histname = 'h%s%s_merged' % (sample, variables[0].replace('(','_').replace(')','_'))
        hist = hists[0].Clone(histname)
        hist.Reset()
        hist.Merge(hists)
        self.logger.debug('The total integral for %s after merging is %f.' % (sample, hist.Integral()))
        hist = self.getOverflowUnderflow(hist,**kwargs)
        self.logger.debug('After overflow it is %f.' % (hist.Integral()))
        hist.Sumw2()
        if normalize:
            integral = hist.Integral()
            if integral: hist.Scale(1.0/integral)
        # set styles
        if not noFormat: hist.SetTitle(self.dataStyles[sample]['name'])
        if sample in self.data: return hist
        if not noFormat:
            hist.SetFillColor(self.dataStyles[sample]['fillcolor'])
            hist.SetLineColor(self.dataStyles[sample]['linecolor'])
            hist.SetFillStyle(self.dataStyles[sample]['fillstyle'])
        return hist

    def getData2D(self, var1, var2, bin1, bin2, cut, **kwargs):
        '''Return a histogram of data for the given variable'''
        hists = ROOT.TList()
        for sample in self.data:
            hist = self.getHist2D(sample, var1, var2, bin1, bin2, cut, **kwargs)
            hists.Add(hist)
        hist = hists[0].Clone("hdata%s%s" % (var1[0], var2[0]))
        hist.Reset()
        hist.Merge(hists)
        return hist

    def getData(self, variables, binning, cut, noFormat=False, **kwargs):
        '''Return a histogram of data for the given variable'''
        hists = ROOT.TList()
        for sample in self.data:
            hist = self.getHist(sample, variables, binning, cut, noFormat, **kwargs)
            hists.Add(hist)
        histname = 'h%s_data' % variables[0].replace('(','_').replace(')','_')
        hist = hists[0].Clone(histname)
        hist.Reset()
        hist.Merge(hists)
        return hist

    def getDataDrivenHist(self, variable, binning, cut, noFormat=False, **kwargs):
        '''Return a histogram corresponding to the data driven estimation of a background'''
        histName = 'hdataDriven%s' % variable
        binEdges = binning
        if len(binning)==3:
            hist = ROOT.TH1F(histName, histName, *binning)
            binEdges = range(binning[1],binning[2],(binning[2]-binning[1])/binning[0])
        else:
            hist = ROOT.TH1F(histName+'temp', histName+'temp', len(binning)-1, binning[0], binning[-1])
            hist.Rebin(len(binning)-1,histName+'new',array('d',binning))
            hist = ROOT.gDirectory.Get(histName+'new').Clone(histName)
        for b in range(len(binning)-1):
            bin = b+1
            cutString = '%s & %s >= %f & %s < %f' % (cut, variable, binEdges[b], variable, binEdges[b+1])
            vals = {}
            errs = {}
            for comb in itertools.product('PF', repeat=3): # hard code 3l for now
                combCut = cutString
                for l in range(3):
                    if comb[l] == 'P':
                        combCut += ' & l%i.passTight==1'
                    else:
                        combCut += ' & l%i.passTight==0'
                val = 0
                err2 = 0
                for sample in self.data:
                    # need to embed the weights, pt dependent, calculate in z->ll sample and QCD enriched smaple (prompt and fake rate, respectively)
                    tempval, temperr = self.getNumEntries(combCut,sample,doError=True,doDataDriven=True,**kwargs)
                    val += tempval
                    err2 += temperr**2
                err = err2**0.5
                vals[comb] = val
                errs[comp] = err
            binContent = vals['PPP'] - vals['PPF'] - vals['PFP'] + vals['PFF'] - vals['FPP'] + vals['FPF'] + vals['FFP'] - vals['FFF']
            binError = sum([x**2 for x in errs.itervalues()])**0.5
            hist.SetBinContent(bin,binContent)
            hist.SetBinError(bin,binError)

    def getMCStack2D(self, var1, var2, bin1, bin2, cut, **kwargs):
        '''Return a stack of MC histograms'''
        hists = ROOT.TList()
        for sample in self.backgrounds:
            hist = self.getHist2D(sample, var1, var2, bin1, bin2, cut, **kwargs)
            hists.Add(hist)
        hist = hists[0].Clone("hdata%s%s" % (var1[0], var2[0]))
        hist.Reset()
        hist.Merge(hists)
        return hist


    def getMCStack(self, variables, binning, cut, **kwargs):
        '''Return a stack of MC histograms'''
        nostack = kwargs.pop('nostack',False)
        mcstack = ROOT.THStack('hs%s' % variables[0],'mc stack')
        for sample in self.backgrounds:
            hist = self.getHist(sample, variables, binning, cut, **kwargs)
            if not hist: continue
            if nostack:
                hist.SetFillStyle(0)
                hist.SetLineWidth(2)
            mcstack.Add(hist)
        self.logger.debug('And the full stack integral is %f.' % mcstack.GetStack().Last().Integral())
        return mcstack

    def get_stat_err(self, hist):
        '''Create statistical errorbars froma histogram'''
        staterr = hist.Clone("staterr")
        staterr.Sumw2()
        staterr.SetFillColor(ROOT.kGray+3)
        staterr.SetLineColor(ROOT.kGray+3)
        staterr.SetLineWidth(0)
        staterr.SetMarkerSize(0)
        staterr.SetFillStyle(3013)
        return staterr

    def get_ratio(self, num, denom, label):
        '''Return a ratio histogram'''
        ratio = num.Clone(label)
        ratio.Sumw2()
        ratio.SetMarkerSize(0.8)
        ratio.Divide(num, denom, 1., 1., "")
        return ratio

    def get_ratio_stat_err(self, hist):
        '''Return a statistical error bars for a ratio plot'''
        ratiostaterr = hist.Clone("ratiostaterr")
        ratiostaterr.Sumw2()
        ratiostaterr.SetStats(0)
        ratiostaterr.SetTitle("")
        ratiostaterr.GetYaxis().SetTitle("Data/MC")
        ratiostaterr.SetMaximum(1.5)
        ratiostaterr.SetMinimum(0.5)
        ratiostaterr.SetMarkerSize(0)
        ratiostaterr.SetFillColor(ROOT.kGray+3)
        ratiostaterr.SetFillStyle(3013)
        ratiostaterr.GetXaxis().SetLabelSize(0.19)
        ratiostaterr.GetXaxis().SetTitleSize(0.21)
        ratiostaterr.GetXaxis().SetTitleOffset(1.0)
        ratiostaterr.GetYaxis().SetLabelSize(0.19)
        ratiostaterr.GetYaxis().SetTitleSize(0.21)
        ratiostaterr.GetYaxis().SetTitleOffset(0.27)
        ratiostaterr.GetYaxis().SetNdivisions(503)

        # bin by bin errors
        for i in range(hist.GetNbinsX()+2):
            ratiostaterr.SetBinContent(i, 1.0)
            if hist.GetBinContent(i)>1e-6:  # not empty
                binerror = hist.GetBinError(i) / hist.GetBinContent(i)
                ratiostaterr.SetBinError(i, binerror)
            else:
                ratiostaterr.SetBinError(i, 999.)

        return ratiostaterr

    def setStyle(self,position=11,plotdata=True,plotratio=False,preliminary=True):
        '''Set style for plots based on the CMS TDR style guidelines.
           https://twiki.cern.ch/twiki/bin/view/CMS/Internal/PubGuidelines#Figures_and_tables
           https://ghm.web.cern.ch/ghm/plots/'''
        # set period (used in CMS_lumi)
        # period : sqrts
        # 1 : 7, 2 : 8, 3 : 7+8, 4 : 13, ... 7 : 7+8+13
        self.period = 1*self.plot7TeV + 2*self.plot8TeV + 4*self.plot13TeV
        CMS_lumi.wrtieExtraText = preliminary
        CMS_lumi.extraText = "Preliminary" if plotdata else "Simulation Preliminary"
        CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        if self.intLumi < 1000:
            CMS_lumi.lumi_7TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
            CMS_lumi.lumi_8TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
            CMS_lumi.lumi_13TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
        CMS_lumi.CMS_lumi(self.plotpad if plotratio else self.canvas,self.period,position)

    def getLegend(self,plotdata,plotsig,plotratio,legendpos,mchist,datahist,sighists):

        numEntries = 0
        if mchist: numEntries += len(self.backgrounds)
        if plotsig: numEntries += len(sighists)
        if plotdata: numEntries += 1
        # setup legend position
        if legendpos % 10 == 1:   # on the left
            xend = 0.45
        elif legendpos % 10 == 2: # in the middle
            xend = 0.70
        else: # default (on right)
            xend = 0.95
        if legendpos//10 == 1:    # bottom
            yend = 0.35
        elif legendpos//10 == 2:  # middle
            yend = 0.59
        elif legendpos//10 == 4:  # very top (in line with CMS label)
            yend = 0.91
        else:                     # default, top, just below CMS label
            yend = 0.77
        xstart = xend-0.3
        ystart = yend-numEntries*0.045
        if plotratio: yend *= 0.95
        # create and draw legend
        leg = ROOT.TLegend(xstart,ystart,xend,yend,'','NDC')
        leg.SetTextFont(42)
        leg.SetTextSize(0.25/numEntries)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        if plotdata: leg.AddEntry(datahist,'Data','ep')
        if mchist:
            for hist in mchist.GetHists():
                leg.AddEntry(hist,hist.GetTitle(),'f')
        if plotsig:
            for s in self.signal:
                leg.AddEntry(sighists[s],sighists[s].GetTitle(),'f')
        return leg

    def save(self, savename):
        '''Save the canvas in multiple formats.'''
        for type in ['png','root']:
            name = "%s/%s/%s.%s" % (self.plotDir, type, savename, type)
            python_mkdir(os.path.dirname(name))
            self.canvas.Print(name)
        self.canvas.SetName(savename)
        self.savefile.WriteTObject(self.canvas)
        self.canvas.Clear()

