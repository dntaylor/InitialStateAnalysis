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
from InitialStateAnalysis.Utilities.utilities import *
#from InitialStateAnalysis.Limits.limitUtils import 
from systematicUncertainties import *

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
        region = kwargs.pop('region',analysis)
        saveDir = kwargs.pop('saveDir','')
        ntupleDir = kwargs.pop('ntupleDir','ntuples')
        period = kwargs.pop('period',13)
        blind = kwargs.pop('blind',False)
        rootName = kwargs.pop('rootName','plots')
        mergeDict = kwargs.pop('mergeDict',{})
        scaleFactor = kwargs.pop('scaleFactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
        dataScaleFactor = kwargs.pop('dataScaleFactor','1')
        datadriven = kwargs.pop('datadriven',False)
        baseSelection = kwargs.pop('baseSelection','')
        self.dontSave = kwargs.pop('dontSave',False)
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
        self.setupCanvas()

        # now, setup plotter conditions (some to be initalized later)
        self.j = 0 # global variable to prevent resusing histograms
        self.backgroundInitialized = False
        self.background = []
        self.dataInitialized = False
        self.data = []
        self.signalInitialized = False
        self.signal = []
        self.analysis = analysis
        self.region = region
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
        self.dataScaleFactor = dataScaleFactor
        self.datadriven = datadriven
        #self.baseSelection = baseSelection
        #if self.baseSelection:
        #    self.tempNtupleFileName = self.plotDir+"/"+rootName+"_temp.root"
        #    self.tempNtupleFile = ROOT.TFile(self.tempNtupleFileName,"recreate")
        self.period=period
        #self.proof = ROOT.TProof.Open("workers=8")
        # adapt this to turn off proof output
        #void PrintEmptyProgress(Long64_t, Long64_t, Float_t, Long64_t)
        #{
        #   return;
        #}
        #gProof->SetPrintProgress(&PrintEmptyProgress);

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        #if self.baseSelection:
        #    self.logger.debug('Cleaning up plotter')
        #    self.tempNtupleFile.Delete("all")
        #    os.remove(self.tempNtupleFileName)
        pass

    def __del__(self):
        #if self.baseSelection:
        #    self.logger.debug('Cleaning up plotter')
        #    self.tempNtupleFile.Delete("all")
        #    os.remove(self.tempNtupleFileName)
        pass

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

    def setupCanvas(self):
        '''Setup the intial canvas'''
        self.canvas.SetFillColor(0)
        self.canvas.SetBorderMode(0)
        self.canvas.SetFrameFillStyle(0)
        self.canvas.SetFrameBorderMode(0)
        self.canvas.SetLeftMargin( self.L )
        self.canvas.SetRightMargin( self.R )
        self.canvas.SetTopMargin( self.T )
        self.canvas.SetBottomMargin( self.B )
        self.canvas.SetTickx(1)
        self.canvas.SetTicky(1)

    def resetCanvas(self):
        '''Reset canvas after changes'''
        self.canvas.SetCanvasSize(self.W,self.H)
        self.canvas.SetLeftMargin( self.L )
        self.canvas.SetRightMargin( self.R )
        self.canvas.SetTopMargin( self.T )
        self.canvas.SetBottomMargin( self.B )
        self.canvas.SetTickx(1)
        self.canvas.SetTicky(1)
        

    def initializeBackgroundSamples(self,sampleList):
        '''Initialize the background samples.'''
        self.backgrounds = sampleList
        if self.datadriven:
            self.backgrounds = ['datadriven']
            self.backgrounds.extend([x for x in sampleList if x not in ['SingleTop', 'TTJets', 'ZJets', 'ZJetsFiltered']])
            sampleList = self.backgrounds
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
        if self.datadriven and sample == 'datadriven': return
        filename = self.ntupleDir+'/%s.root' % sample
        filedir = self.ntupleDir+'/%s' % sample
        if os.path.isfile(filename):
            filenames = [filename]
            self.logger.debug('Adding {0}'.format(filename))
        elif os.path.exists(filedir):
            filenames = glob.glob('{0}/*.root'.format(filedir))
            self.logger.debug('Adding {0} files from {1}'.format(len(filenames),filedir))
        else:
            filenames = []
            self.logger.warning('File {0} does not exist and directory {1} does not exists.'.format(filename,filedir))
        tchain = ROOT.TChain(self.analysis)
        for fname in filenames:
            tchain.Add(fname)
        #tchain.SetProof()
        self.samples[sample]['files'] = filenames
        #self.samples[sample]['file'] = ROOT.TFile(filename)
        #tree = self.samples[sample]['file'].GetTree(self.analysis)
        #if self.baseSelection:
        #    self.tempNtupleFile.cd()
        #    self.samples[sample]['tree'] = tree.CopyTree(self.baseSelection)
        #else:
        self.samples[sample]['tree'] = tchain
        if 'data' in sample:
            lumifile = self.ntupleDir+'/%s.lumicalc.sum' % sample
        else:
            n_evts = 0
            for fn in self.samples[sample]['files']:
                tfile = ROOT.TFile(fn)
                cutflowHist = tfile.Get('cutflow')
                n_evts += cutflowHist.GetBinContent(1)
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

    def getScaleFactor(self):
        '''Set Scale factor'''
        return self.scaleFactor

    def setScaleFactor(self,scalefactor):
        '''Set Scale factor'''
        self.scaleFactor = scalefactor

    def getAnalysis(self):
        return self.analysis

    def getPeriod(self):
        return self.period

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
        bp = kwargs.pop('bp','')
        if bp and 'HPlusPlus' in sample: # need to get scaled entries
            do4l = kwargs.pop('do4l',False)
            runTau = kwargs.pop('runTau',False)
            s = getScales(bp)
            genLeps = 4 if self.analysis=='Hpp4l' or do4l else 3
            recoLeps = 4 if self.analysis=='Hpp4l' else 3
            channelMap = getChannelMap(bp, genLeps, recoLeps, runTau=runTau)
            higgsChannels = channelMap['names']
            genChannelsMap = channelMap['genmap']
            recoChannelsMap = channelMap['recomap']
            allRecoChannels = channelMap['allreco']
            sf = getattr(s,'scale_%s'%self.analysis)
            if do4l: sf = getattr(s,'scale_Hpp4l')
            # setup the scaling
            genChannels = []
            genCuts = []
            genScales = []
            recoChannels = []
            for r in allRecoChannels:
                logging.debug('Adding reco channel %s' %r)
                recoChannels += [r]
                recoCut = 'channel=="%s"' % r
                for h in genChannelsMap:
                    if r not in recoChannelsMap[h]: continue
                    for g in genChannelsMap[h]:
                        scale = sf(g[:2],g[2:])
                        if not scale: continue
                        logging.debug('Adding gen channel %s with scale %f' %(g,scale))
                        genChannels += [g]
                        genCuts += ['(%s && genChannel=="%s")' % (recoCut,g)]
                        genScales += [scale]
            # now get the summed scaled entries
            totVal = 0.
            totErr2 = 0.
            for gencut, scale in zip(genCuts, genScales):
                thisCut = '{0} && {1}'.format(selection,gencut)
                val, err = self.getScaledIndividualSampleEntries(sample,thisCut,scalefactor,**kwargs)
                totVal += val*scale
                totErr2 += err*err*scale*scale
            return totVal, abs(totErr2)**0.5
        else:
            return self.getScaledIndividualSampleEntries(sample,selection,scalefactor,**kwargs)

    def getScaledIndividualSampleEntries(self,sample,selection,scalefactor,**kwargs):
        if 'datadriven' == sample:
            hist = self.getDataDrivenHist('1', [1,-10,10], selection, **kwargs)
            if hist:
                hist.Sumw2()
                val = hist.Integral()
                err = hist.GetBinError(1)
                self.logger.debug('%s: Entries: %f'%(sample,val))
            else:
                val = 0
                err = 0
        else:
            tree = self.samples[sample]['tree']
            histname = 'h%s%i' %(sample,self.j)
            tree.Draw('1>>%s(1,-10,10)'%(histname),'%s*(%s)' %(scalefactor,selection),'goff')
            #hist = self.proof.GetOutputList().FindObject(histname)
            if not ROOT.gDirectory.Get(histname):
                self.logger.debug('%s: No entries'%sample)
                val = 0
                err = 0
            else:
                hist = ROOT.gDirectory.Get(histname).Clone("hnew%s%i" %(sample,self.j))
                hist.Sumw2()
                val = hist.Integral()
                err = hist.GetBinError(1)
                self.logger.debug('%s: Entries: %f'%(sample,val))
        #err = abs(val) ** 0.5
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
        customScale = kwargs.pop('customScale','')
        totalVal = 0
        totalErr2 = 0
        scalefactor = "event.gen_weight*event.lep_scale_up*event.trig_scale*event.pu_weight" if scaleup else self.scaleFactor
        if 'data' in sample and sample!='datadriven':
            scalefactor = self.dataScaleFactor
        if customScale: scalefactor = customScale
        self.j += 1
        self.logger.debug('Cut: %s'%selection)
        if sample in self.sampleMergeDict:
            for s in self.sampleMergeDict[sample]:
                thisCut = selection + ' & ' + self.sampleMergeDict[sample][s] if 'data' not in s else selection
                val, err = self.getIndividualSampleEntries(s,thisCut,scalefactor,**kwargs)
                totalVal += val
                totalErr2 += err*err
        else:
            val, err = self.getIndividualSampleEntries(sample,selection,scalefactor,**kwargs)
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
        return hist # TODO fix
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
        #tree = self.samples[sample]['file'].Get(self.analysis)
        tree = self.samples[sample]['tree']
        zbin = kwargs.pop('zbin',[10,0,10])
        histname = 'h%s%s%s' % (sample,var1,var2)
        drawString = "%s:%s>>%s(%s)" % (var2,var1,histname,', '.join(str(x) for x in bin1+bin2))
        if not cut: cut = '1'
        if 'data' not in sample:
            tree.Draw(drawString,'%s*(%s)' % (self.scaleFactor,cut),'goff')
        else:
            tree.Draw(drawString,'%s*(%s)' % (self.dataScaleFactor,cut),'goff')
        #hist = self.proof.GetOutputList().FindObject(histname)
        if not ROOT.gDirectory.Get(histname):
            return 0
        hist = ROOT.gDirectory.Get(histname).Clone("hmod%s%s%s"%(sample,var1,var2))
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


    def getSingleVarHist(self,sample,variable,binning,cut,**kwargs):
        '''Single variable, single sample hist'''
        customScale = kwargs.pop('customScale','')
        #tree = self.samples[sample]['file'].Get(self.analysis)
        tree = self.samples[sample]['tree']
        self.j += 1
        histname = 'h%s%s' % (sample, variable.replace('(','_').replace(')','_'))
        if len(binning) == 3: # standard drawing
            drawString = "%s>>%s(%s)" % (variable, histname, ", ".join(str(x) for x in binning))
        else: # we will need to rebin TODO: this might cause issue in root6, double check later
            drawString = "%s>>%s()" % (variable, histname)
        if not cut: cut = '1'
        if customScale:
            tree.Draw(drawString,'%s*(%s)' % (customScale,cut),'goff')
        else:
            if 'data' not in sample:
                tree.Draw(drawString,'%s*(%s)' % (self.scaleFactor,cut),'goff')
            else:
                tree.Draw(drawString,'%s*(%s)' % (self.dataScaleFactor,cut),'goff')
        #hist = self.proof.GetOutputList().FindObject(histname)
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
        if sample=='datadriven': return self.getDataDrivenHist(variables, binning, cut, noFormat=False, normalize=normalize, **kwargs)
        hists = ROOT.TList()
        for v in range(len(variables)):
            if sample in self.sampleMergeDict:
                for s in self.sampleMergeDict[sample]:
                    if len(variables) != len(cut):
                        thisCut = cut + ' & ' + self.sampleMergeDict[sample][s]
                        hist = self.getSingleVarHist(s,variables[v],binning,thisCut,**kwargs)
                    else:
                        thisCut = cut[v] + ' & ' + self.sampleMergeDict[sample][s]
                        hist = self.getSingleVarHist(s,variables[v],binning,thisCut,**kwargs)
                    if hist:
                        hists.Add(hist)
            else:
                if len(variables) != len(cut):
                    hist = self.getSingleVarHist(sample,variables[v],binning,cut,**kwargs)
                else:
                    hist = self.getSingleVarHist(sample,variables[v],binning,cut[v],**kwargs)
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
        #hist = self.getPoissonErrors(hist)
        return hist

    def getPoissonErrors(self,hist):
        #return hist
        # adapted from rootpy to get asymmetric poisson errors
        graph = ROOT.TGraphAsymmErrors(hist.GetNbinsX())
        #graph.SetLineWidth(self.GetLineWidth())
        #graph.SetMarkerSize(self.GetMarkerSize())
        chisqr = ROOT.TMath.ChisquareQuantile
        npoints = 0
        for bin in range(hist.GetNbinsX()):
            entries = hist.GetBinContent(bin+1)
            if entries <= 0:
                continue
            ey_low = entries - 0.5 * chisqr(0.1586555, 2. * entries)
            ey_high = 0.5 * chisqr(
                1. - 0.1586555, 2. * (entries + 1)) - entries
            ex = hist.GetBinWidth(bin+1) / 2.
            graph.SetPoint(npoints, hist.GetBinCenter(bin+1), hist.GetBinContent(bin+1))
            #graph.SetPointEXlow(npoints, ex)
            graph.SetPointEXlow(npoints, 0)
            #graph.SetPointEXhigh(npoints, ex)
            graph.SetPointEXhigh(npoints, 0)
            graph.SetPointEYlow(npoints, ey_low)
            graph.SetPointEYhigh(npoints, ey_high)
            npoints += 1
        graph.Set(npoints)
        return graph

    def getDataDrivenHist(self, variables, binning, cut, noFormat=False, **kwargs):
        '''Return a histogram corresponding to the data driven estimation of a background'''
        getPrompt = kwargs.pop('getPrompt',False)
        doSimple = kwargs.pop('doSimple',True)
        nameMap = {
            0: 'z1.PassTight1',
            1: 'z1.PassTight2',
            2: 'w1.PassTight1',
        }
        effMap = {
            0: 'z1.LepEffTight1',
            1: 'z1.LepEffTight2',
            2: 'w1.LepEffTight1',
        }
        fakeMap = {
            0: 'z1.LepFake1',
            1: 'z1.LepFake2',
            2: 'w1.LepFake1',
        }
        if self.analysis in ['Hpp3l']:
            nameMap = {
                0: 'h1.PassTight1',
                1: 'h1.PassTight2',
                2: 'h2.PassTight1',
            }
            effMap = {
                0: 'h1.LepEffTight1',
                1: 'h1.LepEffTight2',
                2: 'h2.LepEffTight1',
            }
            fakeMap = {
                0: 'h1.LepFake1',
                1: 'h1.LepFake2',
                2: 'h2.LepFake1',
            }

        # remove any passTight cuts, assumes these are there, since it only returns the all tight stuff
        allCuts = getPassTightDefinition(self.analysis,self.region,self.period)
        if self.analysis in ['Hpp3l']:
            allCuts = '1' # no special cuts at select.passTight
        if type(cut) is list:
            for c in range(len(cut)):
                for l in range(3):
                    cut[c] = cut[c].replace(nameMap[l],'1')
                    cut[c] = cut[c].replace('l{0}.PassTight'.format(l),'1')
                cut[c] = cut[c].replace('select.passTight',allCuts)
        else:
            for l in range(3):
                cut = cut.replace(nameMap[l],'1')
                cut = cut.replace('l{0}.PassTight'.format(l),'1')
            cut = cut.replace('select.passTight',allCuts)
        hists = ROOT.TList()

        # the new way, dont do it a million times
        if type(cut) is list:
            custCut = ['{0} && fakeChannel!="PPP"'.format(x) for x in cut]
        else:
            custCut = '{0} && fakeChannel!="PPP"'.format(cut)
        scalefactor = 'event.fakerate'
        # get contribution from data
        for sample in self.data:
            hist = self.getHist(sample,variables,binning,custCut,customScale=scalefactor,**kwargs)
            hists.Add(hist)
        # subtract contribution from signal MC
        mcScalefactor = self.scaleFactor + '*(-1*' + scalefactor + ')'
        for sample in self.backgrounds:
            if sample=='datadriven': continue
            if 'HPlusPlus' in sample: continue # dont include signal in subtraction
            hist = self.getHist(sample,variables,binning,custCut,customScale=mcScalefactor,**kwargs)
            hists.Add(hist)


        #for comb in itertools.product('PF', repeat=3): # hard code 3l for now
        #    if comb.count('P')==3: continue # signal region
        #    if type(cut) is list:
        #        custCut = ['{0} && fakeChannel=="{1}"'.format(x,''.join(comb)) for x in cut]
        #    else:
        #        custCut = '{0} && fakeChannel=="{1}"'.format(cut,''.join(comb))
        #    scalefactor = 'event.fakerate'
        #    #for l in range(3):
        #    #    if type(cut) is list:
        #    #        for c in range(len(cut)):
        #    #            if doSimple:
        #    #                custCut[c] += ' && {0}==1'.format(nameMap[l]) if comb[l] == 'P' else ' && {0}==0'.format(nameMap[l])
        #    #            else:
        #    #                custCut[c] += ' && {0}==1'.format(nameMap[l]) if comb[l] == 'P' else ' && {0}==0'.format(nameMap[l])
        #    #    else:
        #    #        if doSimple:
        #    #            custCut += ' && {0}==1'.format(nameMap[l]) if comb[l] == 'P' else ' && {0}==0'.format(nameMap[l])
        #    #        else:
        #    #            custCut += ' && {0}==1'.format(nameMap[l]) if comb[l] == 'P' else ' && {0}==0'.format(nameMap[l])
        #    #denom = '1./(({0}-{1})*({2}-{3})*({4}-{5}))'.format(effMap[0],fakeMap[0],effMap[1],fakeMap[1],effMap[2],fakeMap[2])
        #    #num = '1' if comb.count('P') in [1,3] else '-1'
        #    #for l in range(3):
        #    #    if comb[l] == 'P':
        #    #        num += '*(1-{0})'.format(fakeMap[l])
        #    #    else:
        #    #        num += '*({0})'.format(fakeMap[l])
        #    #scalefactor = num + '*' + denom
        #    #if doSimple: # sascha's way
        #    #    num = '1' if comb.count('F') in [1,3] else '-1'
        #    #    for l in range(3):
        #    #        if comb[l] == 'F': num += '*{0}'.format(fakeMap[l])
        #    #    scalefactor = num
        #    #    if comb.count('P')==3: continue
        #    # get contribution from data
        #    for sample in self.data:
        #        hist = self.getHist(sample,variables,binning,custCut,customScale=scalefactor,**kwargs)
        #        hists.Add(hist)
        #    # subtract contribution from signal MC
        #    mcScalefactor = self.scaleFactor + '*(-1*' + scalefactor + ')'
        #    for sample in self.backgrounds:
        #        if sample=='datadriven': continue
        #        if 'HPlusPlus' in sample: continue # dont include signal in subtraction
        #        hist = self.getHist(sample,variables,binning,custCut,customScale=mcScalefactor,**kwargs)
        #        hists.Add(hist)

        histname = 'h%s_datadriven' % variables[0].replace('(','_').replace(')','_')
        if hists.IsEmpty(): return 0
        hist = hists[0].Clone(histname)
        hist.Reset()
        hist.Merge(hists)
        if getPrompt:
            result = self.getData(variables,binning,cut,noFormat=True,**kwargs)
            result.Add(hist,-1)
        else:
            result = hist
        if not noFormat:
            result.SetTitle(self.dataStyles['datadriven']['name'])
            result.SetFillColor(self.dataStyles['datadriven']['fillcolor'])
            result.SetLineColor(self.dataStyles['datadriven']['linecolor'])
            result.SetFillStyle(self.dataStyles['datadriven']['fillstyle'])

        return result

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

    def addSystematicUncertainty(self,hist,sample):
        unc = getSystUncertaintyMap(self.analysis,self.region,self.period,sample)
        totsyst2 = 0.
        for u in unc:
            totsyst2 += unc[u]**2
        totsyst = totsyst2**0.5
        nbins = hist.GetNbinsX()
        for n in range(nbins):
            val = hist.GetBinContent(n+1)
            err = hist.GetBinError(n+1)
            toterr = ((val*totsyst)**2 + err**2)**0.5
            hist.SetBinError(n+1,toterr)
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
            histsyst = self.addSystematicUncertainty(hist,sample)
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
        #ratio.Divide(num.GetHistogram(), denom, "pois")
        return ratio

    def getPoissonRatio(self,num,denom,label):
        # get ratio between two hists with poisson errors
        graph = ROOT.TGraphAsymmErrors(num.GetNbinsX())
        chisqr = ROOT.TMath.ChisquareQuantile
        npoints = 0
        for bin in range(num.GetNbinsX()):
            entries = num.GetBinContent(bin+1)
            denomentries = denom.GetBinContent(bin+1)
            if entries <= 0 or denomentries <= 0:
                continue
            ey_low = entries - 0.5 * chisqr(0.1586555, 2. * entries)
            ey_high = 0.5 * chisqr(
                1. - 0.1586555, 2. * (entries + 1)) - entries
            ex = num.GetBinWidth(bin+1) / 2.
            graph.SetPoint(npoints, num.GetBinCenter(bin+1), num.GetBinContent(bin+1)/denomentries)
            graph.SetPointEXlow(npoints, 0)
            graph.SetPointEXhigh(npoints, 0)
            graph.SetPointEYlow(npoints, ey_low/denomentries)
            graph.SetPointEYhigh(npoints, ey_high/denomentries)
            npoints += 1
        graph.Set(npoints)
        return graph


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
        self.period_int = 1*self.plot7TeV + 2*self.plot8TeV + 4*self.plot13TeV
        CMS_lumi.wrtieExtraText = preliminary
        CMS_lumi.extraText = "Preliminary" if plotdata else "Simulation Preliminary"
        CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        if self.intLumi < 1000:
            CMS_lumi.lumi_7TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
            CMS_lumi.lumi_8TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
            CMS_lumi.lumi_13TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
        CMS_lumi.CMS_lumi(self.plotpad if plotratio else self.canvas,self.period_int,position)

    def getLegend(self,mchist,datahist,sighists,**kwargs):
        legendpos = kwargs.pop('legendpos',33)
        plotdata = kwargs.pop('plotdata',True)
        plotsig = kwargs.pop('plotsig',True)
        plotmc = kwargs.pop('plotmc',True)
        plotratio = kwargs.pop('plotratio',True)
        numcol = kwargs.pop('numcol',1)

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
            xend = 0.92
        if legendpos//10 == 1:    # bottom
            yend = 0.35
        elif legendpos//10 == 2:  # middle
            yend = 0.59
        elif legendpos//10 == 4:  # very top (in line with CMS label)
            yend = 0.91
        else:                     # default, top, just below CMS label
            yend = 0.77
        xstart = xend-0.15*numcol-0.15
        ystart = yend-math.ceil(float(numEntries)/numcol)*0.06
        if plotratio: yend *= 0.95
        # create and draw legend
        leg = ROOT.TLegend(xstart,ystart,xend,yend,'','NDC')
        if numcol>1: leg.SetNColumns(int(numcol))
        leg.SetTextFont(42)
        #leg.SetTextSize(0.25/numEntries)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        if plotdata: leg.AddEntry(datahist,'Data','ep')
        if mchist:
            for hist in reversed(mchist.GetHists()):
                leg.AddEntry(hist,hist.GetTitle(),'f')
        if plotsig:
            for s in self.signal:
                leg.AddEntry(sighists[s],sighists[s].GetTitle(),'f')
        return leg

    def save(self, savename):
        '''Save the canvas in multiple formats.'''
        self.canvas.SetName(savename)
        if self.dontSave: return
        for type in ['png','root','pdf']:
            name = "%s/%s/%s.%s" % (self.plotDir, type, savename, type)
            python_mkdir(os.path.dirname(name))
            self.canvas.Print(name)
        self.savefile.WriteTObject(self.canvas)
        #self.canvas.Clear()

