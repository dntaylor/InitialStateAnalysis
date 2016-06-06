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
from multiprocessing import Pool
#import copy_reg
#import types
import time

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

# Non class functions for parallelization
# python can't pickle class methods, which is required by Pool
def getHist(args):
    histname,tree,variable,binning,scalefactor,cut = args
    drawString = "%s>>%s(%s)" % (variable, histname, ", ".join(str(x) for x in binning))
    tree.Draw(drawString,'%s*(%s)' % (scalefactor,cut),'goff')
    if not ROOT.gDirectory.Get(histname):
        return histname, 0
    hist = ROOT.gDirectory.Get(histname)
    hist.Sumw2()
    return histname, hist

def getNumEntries(args):
    histname,tree,scalefactor,cut = args
    variable = '1'
    binning = [1,-10,10]
    _, hist = getHist([histname,tree,variable,binning,scalefactor,cut])
    if not hist: return histname, [0.,0.]
    val = hist.Integral()
    err = hist.GetBinError(1)
    return histname, [val,err]

def getSummedNumEntries(histname,jobs):
    try:
        results = pool.map_async(getNumEntries, jobs).get()
    except KeyboardInterrupt:
        pool.terminate()
        print 'getSummedEntries cancelled'
        sys.exit(1)
    #print histname
    #for h in results:
    #    print '    {0:20.20}: {1:8.4f} {2:8.4f}'.format(h[0],h[1][0],h[1][1])
    val = sum([h[1][0] for h in results])
    err = sum([h[1][1]**2 for h in results])**0.5
    #print '    {0:20.20}: {1:8.4f} {2:8.4f}'.format('',val,err)
    return histname, [val,err]

def getMergedHist(histname,jobs):
    hists = ROOT.TList()
    results = [getHist(j) for j in jobs]
    #try:
    #    results = pool.map_async(getHist,jobs).get()
    #except KeyboardInterrupt:
    #    pool.terminate()
    #    print 'getMergedHist cancelled'
    #    sys.exit(1)
    for h in results:
        if not h[1]: continue
        #print '    ', h[0], h[1].Integral()
        hists.Add(h[1])
    if hists.IsEmpty():
        return histname, 0
    hist = hists[0].Clone(histname)
    hist.Reset()
    hist.Merge(hists)
    return histname, hist

# basic plotter class
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
        self.tightW = kwargs.pop('tightW',True)
        self.allMedium = kwargs.pop('allMedium',False)
        self.fakeMode = kwargs.pop('fakeMode','fakerate')
        self.numLeptons = kwargs.pop('numLeptons',3)
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
        canvas = ROOT.TCanvas("c1","c1",50,50,self.W,self.H)
        canvas = self.setupCanvas(canvas)

        # setup the scale factors for alternate scaling
        # hard coded for now
        self.scales = {
            # name : [filename, varNames],
            'generatorWeight' : ['',['event.gen_weight']], # no filename, use tree directly
        }

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

        # status 1 cut
        self.mcCut = ' && '.join(['fabs(l{0}.GenPatPdgId)<100'.format(l+1) for l in range(self.numLeptons)]) if self.datadriven else '1'

        # setup datadriven stuff
        self.nameMap = {
            0: 'z1.PassTight1',
            1: 'z1.PassTight2',
            2: 'w1.PassTight1',
        }
        self.scaleMap = {
            0: 'z1.LepScaleTight1',
            1: 'z1.LepScaleTight2',
            2: 'w1.LepScaleTight1',
        }
        self.scaleLooseMap = {
            0: 'z1.LepScaleLoose1',
            1: 'z1.LepScaleLoose2',
            2: 'w1.LepScaleLoose1',
        }
        self.effMap = {
            0: 'z1.LepEffTight1',
            1: 'z1.LepEffTight2',
            2: 'w1.LepEffTight1',
        }
        self.fakeMap = {
            0: 'z1.LepFakeTight1',
            1: 'z1.LepFakeTight2',
            2: 'w1.LepFakeTight1',
        }
        self.fakeMapUp = {
            0: 'z1.LepFakeTight_up1',
            1: 'z1.LepFakeTight_up2',
            2: 'w1.LepFakeTight_up1',
        }
        self.fakeMapDown = {
            0: 'z1.LepFakeTight_down1',
            1: 'z1.LepFakeTight_down2',
            2: 'w1.LepFakeTight_down1',
        }
        self.fakeMapMC = {
            0: 'z1.LepFakeMCTight1',
            1: 'z1.LepFakeMCTight2',
            2: 'w1.LepFakeMCTight1',
        }
        self.fakeMapMCUp = {
            0: 'z1.LepFakeMCTight_up1',
            1: 'z1.LepFakeMCTight_up2',
            2: 'w1.LepFakeMCTight_up1',
        }
        self.fakeMapMCDown = {
            0: 'z1.LepFakeMCTight_down1',
            1: 'z1.LepFakeMCTight_down2',
            2: 'w1.LepFakeMCTight_down1',
        }
        self.promptMap = {
            0: 'z1.GenIsPrompt1',
            1: 'z1.GenIsPrompt2',
            2: 'w1.GenIsPrompt1',
        }
        if self.analysis in ['Hpp3l']:
            self.nameMap = {
                0: 'h1.PassTight1',
                1: 'h1.PassTight2',
                2: 'h2.PassTight1',
            }
            self.effMap = {
                0: 'h1.LepEffTight1',
                1: 'h1.LepEffTight2',
                2: 'h2.LepEffTight1',
            }
            self.fakeMap = {
                0: 'h1.LepFake1',
                1: 'h1.LepFake2',
                2: 'h2.LepFake1',
            }
            self.fakeMapMC = {
                0: 'h1.LepFakeMC1',
                1: 'h1.LepFakeMC2',
                2: 'h2.LepFakeMC1',
            }
        if self.tightW:
            self.nameMap = {
                0: 'z1.PassMedium1',
                1: 'z1.PassMedium2',
                2: 'w1.PassTight1',
            }
            self.scaleMap = {
                0: 'z1.LepScaleMedium1',
                1: 'z1.LepScaleMedium2',
                2: 'w1.LepScaleTight1',
            }
            self.scaleLooseMap = {
                0: 'z1.LepScaleLoose1',
                1: 'z1.LepScaleLoose2',
                2: 'w1.LepScaleLoose1',
            }
            self.effMap = {
                0: 'z1.LepEffMedium1',
                1: 'z1.LepEffMedium2',
                2: 'w1.LepEffTight1',
            }
            self.fakeMap = {
                0: 'z1.LepFakeMedium1',
                1: 'z1.LepFakeMedium2',
                2: 'w1.LepFakeTight1',
            }
            self.fakeMapUp = {
                0: 'z1.LepFakeMedium_up1',
                1: 'z1.LepFakeMedium_up2',
                2: 'w1.LepFakeTight_up1',
            }
            self.fakeMapDown = {
                0: 'z1.LepFakeMedium_down1',
                1: 'z1.LepFakeMedium_down2',
                2: 'w1.LepFakeTight_down1',
            }
            self.fakeMapMC = {
                0: 'z1.LepFakeMCMedium1',
                1: 'z1.LepFakeMCMedium2',
                2: 'w1.LepFakeMCTight1',
            }
            self.fakeMapMCUp = {
                0: 'z1.LepFakeMCMedium_up1',
                1: 'z1.LepFakeMCMedium_up2',
                2: 'w1.LepFakeMCTight_up1',
            }
            self.fakeMapMCDown = {
                0: 'z1.LepFakeMCMedium_down1',
                1: 'z1.LepFakeMCMedium_down2',
                2: 'w1.LepFakeMCTight_down1',
            }
        if self.allMedium:
            self.nameMap = {
                0: 'z1.PassMedium1',
                1: 'z1.PassMedium2',
                2: 'w1.PassMedium1',
            }
            self.scaleMap = {
                0: 'z1.LepScaleMedium1',
                1: 'z1.LepScaleMedium2',
                2: 'w1.LepScaleMedium1',
            }
            self.scaleLooseMap = {
                0: 'z1.LepScaleLoose1',
                1: 'z1.LepScaleLoose2',
                2: 'w1.LepScaleLoose1',
            }
            self.effMap = {
                0: 'z1.LepEffMedium1',
                1: 'z1.LepEffMedium2',
                2: 'w1.LepEffMedium1',
            }
            self.fakeMap = {
                0: 'z1.LepFakeMedium1',
                1: 'z1.LepFakeMedium2',
                2: 'w1.LepFakeMedium1',
            }
            self.fakeMapMC = {
                0: 'z1.LepFakeMCMedium1',
                1: 'z1.LepFakeMCMedium2',
                2: 'w1.LepFakeMCMedium1',
            }


    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.__cleanup()

    def __del__(self):
        self.__cleanup()

    def __cleanup(self):
        self.savefile.Close()

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

    def setupCanvas(self,canvas):
        '''Setup the intial canvas'''
        canvas.SetFillColor(0)
        canvas.SetBorderMode(0)
        canvas.SetFrameFillStyle(0)
        canvas.SetFrameBorderMode(0)
        canvas.SetCanvasSize(self.W,self.H)
        canvas.SetLeftMargin( self.L )
        canvas.SetRightMargin( self.R )
        canvas.SetTopMargin( self.T )
        canvas.SetBottomMargin( self.B )
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        return canvas

    def resetCanvas(self,canvas):
        '''Reset canvas after changes'''
        return setupCanvas(canvas)
        

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
        # find files
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
        self.samples[sample]['files'] = filenames
        # get total lumi of sample
        if 'data' in sample:
            lumifile = self.ntupleDir+'/%s.lumicalc.sum' % sample
        else:
            n_evts = 0
            tfiles = {}
            cutflowHists = {}
            for fn in self.samples[sample]['files']:
                #tfile = ROOT.TFile(fn)
                tfiles[fn] = ROOT.TFile.Open(fn)
                cutflowHists[fn] = tfiles[fn].Get('cutflow')
                n_evts += cutflowHists[fn].GetBinContent(1)
                #tfile.Close()
            sample_xsec = self.xsecs[sample]
            self.samples[sample]['lumi'] = float(n_evts)/sample_xsec
            self.logger.debug('Initializing MC sample %s with %i events and xsec %f to lumi %f.'\
                  % (sample, n_evts, sample_xsec, self.samples[sample]['lumi']))
        # create tchain
        tchain = ROOT.TChain(self.analysis)
        for fname in filenames:
            tchain.Add(fname)
        #tchain.SetProof()
        #self.samples[sample]['file'] = ROOT.TFile(filename)
        #tree = self.samples[sample]['file'].GetTree(self.analysis)
        #if self.baseSelection:
        #    self.tempNtupleFile.cd()
        #    self.samples[sample]['tree'] = tree.CopyTree(self.baseSelection)
        #else:
        self.samples[sample]['tree'] = tchain

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

    def __buildAsync_getDataDrivenNumEntriesJobs(self,cut,**kwargs):
        getPrompt = kwargs.pop('getPrompt',False)
        singleBackground = kwargs.pop('singleBackground','')
        jobs = []

        (cutLoose,promptCut,mcLooseCut,fakeScaleFactor,fullMCScale,nonPromptScale) = self.__getDataDrivenHistParameters(cut,getPrompt=getPrompt,**kwargs)

        # get contribution from data
        if not singleBackground or singleBackground=='data':
            for sample in self.data:
                jobs += self.__buildAsync_getNumEntriesJobs(cutLoose,sample,customScale=fakeScaleFactor,**kwargs)
                # if not prompt, subtract from signal region
                #if not getPrompt:
                #    jobs += self.__buildAsync_getNumEntriesJobs(cut,sample,**kwargs)

        # subtract contribution from signal MC if fakerate method used
        if self.fakeMode == 'fakerate':
            for sample in self.backgrounds:
                if singleBackground and singleBackground!=sample: continue
                if 'datadriven' in sample: continue
                if 'HPlusPlus' in sample: continue # dont include signal in subtraction
                jobs += self.__buildAsync_getNumEntriesJobs(mcLooseCut,sample,customScale=fullMCScale,**kwargs)
                #if not getPrompt:
                #    jobs += self.__buildAsync_getNumEntriesJobs(promptCut,sample,customScale=nonPromptScale,**kwargs)

        return jobs

    def __buildAsync_getDataDrivenNumEntriesJobs_mcClosure(self,cut,**kwargs):
        getPrompt = kwargs.pop('getPrompt',False)
        singleBackground = kwargs.pop('singleBackground','')
        jobs = []

        (cutLoose,promptCut,mcLooseCut,fakeScaleFactor,fullMCScale,nonPromptScale) = self.__getDataDrivenHistParameters(cut,getPrompt=getPrompt,mcFake=True,**kwargs)

        # get contribution from control regions
        for sample in self.backgrounds:
            jobs += self.__buildAsync_getNumEntriesJobs(cutLoose,sample,customScale=fullMCScale+'*(-1.)',**kwargs)

        return jobs


    def __buildAsync_getNumEntriesJobs(self,selection,sample,**kwargs):
        doError = kwargs.pop('doError',False)
        scaleup = kwargs.pop('scaleup',False)
        unweighted = kwargs.pop('doUnweighted',False)
        customScale = kwargs.pop('customScale','')
        self.j += 1
        jobs = []
        scalefactor = "event.gen_weight*event.lep_scale_up*event.pu_weight*event.trig_scale" if scaleup else self.scaleFactor
        if 'data' in sample and 'datadriven' not in sample:
            scalefactor = self.dataScaleFactor
        if customScale: scalefactor = customScale
        if sample in self.sampleMergeDict:
            for s in self.sampleMergeDict[sample]:
                thisCut = selection + ' && ' + self.sampleMergeDict[sample][s] + ' && ' + self.mcCut if 'data' not in s else selection
                jobs += [self.__buildAsync_getSingleNumEntriesJob(s,thisCut,scalefactor,**kwargs)]
        else:
            thisCut = selection + ' && ' + self.mcCut if 'data' not in sample else selection
            if 'datadriven' in sample:
                jobs += self.__buildAsync_getDataDrivenNumEntriesJobs(thisCut,**kwargs)
            elif 'mcClosure' in sample:
                jobs += self.__buildAsync_getDataDrivenNumEntriesJobs_mcClosure(thisCut,**kwargs)
            else:
                jobs += [self.__buildAsync_getSingleNumEntriesJob(sample,thisCut,scalefactor,**kwargs)]
        return jobs

    def __buildAsync_getSingleNumEntriesJob(self,sample,selection,scalefactor,**kwargs):
        tree = self.samples[sample]['tree']
        self.j += 1
        histname = 'h_{0}_numEntries_{1}'.format(sample,self.j)
        if 'data' not in sample:
            lumiscale = self.intLumi/self.samples[sample]['lumi']
            scalefactor = '{0}*{1}'.format(scalefactor,lumiscale)
        job = (histname,tree,scalefactor,selection)
        return job

    def __getEntries_async(self,name,samples,selection,**kwargs):
        jobs = []
        for s in samples:
            jobs += self.__buildAsync_getNumEntriesJobs(selection,s,**kwargs)
        self.j += 1
        histname = 'h_{0}_numEntries_{1}'.format(name,self.j)
        result = getSummedNumEntries(histname, jobs)
        val,err = result[1]
        return result

    def getNumEntries(self,selection,sample,**kwargs):
        doError = kwargs.pop('doError',False)
        doSyst = kwargs.pop('doSyst',False)
        kwargs['doError'] = True
        _, result = self.__getEntries_async(sample,[sample],selection,**kwargs)
        if doSyst: # also include syst
            syst = self.getSystematicUncertainty(result[0],sample,**kwargs)
            result += [syst]
        return result if doError else results[0]

    def getSignalEntries(self,selection,**kwargs):
        signal = kwargs.pop('signal','')
        doError = kwargs.pop('doError',False)
        kwargs['doError'] = True
        signals = [signal] if signal else self.signal
        _, result = self.__getEntries_async('signal',signals,selection,**kwargs)
        return result if doError else results[0]

    def getBackgroundEntries(self,selection,**kwargs):
        doError = kwargs.pop('doError',False)
        kwargs['doError'] = True
        _, result = self.__getEntries_async('background',self.backgrounds,selection,**kwargs)
        return result if doError else results[0]

    def getDataEntries(self,selection,**kwargs):
        doError = kwargs.pop('doError',False)
        kwargs['doError'] = True
        _, result = self.__getEntries_async('data',self.data,selection,**kwargs)
        return result if doError else result[0]


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



    def __getDataDrivenHistParameters(self,cut,**kwargs):
        getPrompt = kwargs.pop('getPrompt',False)
        baseScaleFactor = kwargs.pop('baseScaleFactor','event.pu_weight*event.gen_weight*event.trig_scale')
        singleComponent = kwargs.pop('singleComponent','')
        electronOnly = kwargs.pop('electronOnly',False)
        muonOnly = kwargs.pop('muonOnly',False)
        singleChannel = kwargs.pop('channel','')
        mcFake = kwargs.pop('mcFake',False)
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)

        # remove any passTight cuts, assumes these are there, since it only returns the all tight stuff
        allCuts = getPassTightDefinition(self.analysis,self.region,self.period)
        if self.analysis in ['Hpp3l']:
            allCuts = '1' # no special cuts at select.passTight
        cutLoose = cut
        if type(cut) is list:
            for c in range(len(cut)):
                for l in range(self.numLeptons):
                    cutLoose[c] = cutLoose[c].replace(self.nameMap[l],'1')
                    cutLoose[c] = cutLoose[c].replace('l{0}.PassTight'.format(l),'1')
                    cutLoose[c] = cutLoose[c].replace('l{0}.PassMedium'.format(l),'1')
                cutLoose[c] = cutLoose[c].replace('select.passTight',allCuts)
        else:
            for l in range(self.numLeptons):
                cutLoose = cutLoose.replace(self.nameMap[l],'1')
                cutLoose = cutLoose.replace('l{0}.PassTight'.format(l),'1')
                cutLoose = cutLoose.replace('l{0}.PassMedium'.format(l),'1')
            cutLoose = cutLoose.replace('select.passTight',allCuts)

        fakechan = 'fakeChannel'
        if self.tightW: fakechan = 'fakeChannel_tightW'
        if self.allMedium: fakechan = 'fakeChannel_allMedium'

        # exclude signal region
        if self.numLeptons==3:
            cutLoose += ' && {0}!="PPP"'.format(fakechan)
            modes = ['PPP','PPF','PFP','FPP','PFF','FPF','FFP','FFF']
            modeMap = {
                'l': {
                    'eee': {'PPP': 1., 'PPF': 1., 'PFP': 1., 'FPP': 1., 'PFF': 1.,  'FPF': 1.,  'FFP': 1., 'FFF': 1.,},
                    'eem': {'PPP': 1., 'PPF': 1., 'PFP': 1., 'FPP': 1., 'PFF': 1.,  'FPF': 1.,  'FFP': 1., 'FFF': 1.,},
                    'mme': {'PPP': 1., 'PPF': 1., 'PFP': 1., 'FPP': 1., 'PFF': 1.,  'FPF': 1.,  'FFP': 1., 'FFF': 1.,},
                    'mmm': {'PPP': 1., 'PPF': 1., 'PFP': 1., 'FPP': 1., 'PFF': 1.,  'FPF': 0.,  'FFP': 1., 'FFF': 1.,},
                },
                'e': {
                    'eee': {'PPP': 1., 'PPF': 1., 'PFP': 1., 'FPP': 1., 'PFF': 1.,  'FPF': 1.,  'FFP': 1., 'FFF': 1.,},
                    'eem': {'PPP': 1., 'PPF': 0., 'PFP': 1., 'FPP': 1., 'PFF': 0.5, 'FPF': 0.5, 'FFP': 1., 'FFF': 2./3,},
                    'mme': {'PPP': 1., 'PPF': 1., 'PFP': 0., 'FPP': 0., 'PFF': 0.5, 'FPF': 0.5, 'FFP': 0., 'FFF': 1./3,},
                    'mmm': {'PPP': 1., 'PPF': 0., 'PFP': 0., 'FPP': 0., 'PFF': 0.,  'FPF': 0.,  'FFP': 0., 'FFF': 0.,},
                },
                'm':{
                    'eee': {'PPP': 1., 'PPF': 0., 'PFP': 0., 'FPP': 0., 'PFF': 0.,  'FPF': 0.,  'FFP': 0., 'FFF': 0.,},
                    'eem': {'PPP': 1., 'PPF': 1., 'PFP': 0., 'FPP': 0., 'PFF': 0.5, 'FPF': 0.5, 'FFP': 0., 'FFF': 1./3,},
                    'mme': {'PPP': 1., 'PPF': 0., 'PFP': 1., 'FPP': 1., 'PFF': 0.5, 'FPF': 0.5, 'FFP': 1., 'FFF': 2./3,},
                    'mmm': {'PPP': 1., 'PPF': 1., 'PFP': 1., 'FPP': 1., 'PFF': 1.,  'FPF': 1.,  'FFP': 1., 'FFF': 1.,},
                },
            }
        elif self.numLeptons==2:
            cutLoose += ' && {0}!="PP"'.format(fakechan)
            modes = ['PP','PF','FP','FF']
            modeMap = {
                'l': {
                    'ee': {'PP': 1., 'PF': 1., 'FP': 1., 'FF': 1.,  },
                    'em': {'PP': 1., 'PF': 1., 'FP': 1., 'FF': 1.,  },
                    'mm': {'PP': 1., 'PF': 1., 'FP': 1., 'FF': 1.,  },
                },
                'e': {
                    'ee': {'PP': 1., 'PF': 1., 'FP': 1., 'FF': 1.,  },
                    'em': {'PP': 1., 'PF': 0., 'FP': 1., 'FF': 0.5, },
                    'mm': {'PP': 1., 'PF': 0., 'FP': 0., 'FF': 0., },
                },
                'm':{
                    'ee': {'PP': 1., 'PF': 0., 'FP': 0., 'FF': 0.,  },
                    'em': {'PP': 1., 'PF': 1., 'FP': 0., 'FF': 0.5, },
                    'mm': {'PP': 1., 'PF': 1., 'FP': 1., 'FF': 1., },
                },
            }

        modesToRun = modes[1:]
        if singleComponent: modesToRun = [singleComponent]

        # get mc scale factors
        mcScaleFactorBase = baseScaleFactor
        mcFactors = {}
        for mode in modes:
            factors = [mcScaleFactorBase] if mode in modesToRun else ['0']
            for i,f in enumerate(mode):
                factors += [self.scaleMap[i] if f=='P' else self.scaleLooseMap[i]]
            factor = '*'.join(factors)
            mcFactors[mode] = factor
        mcScaleFactor = '*'.join(['({0}=="{1}" ? ({2}) : (1.))'.format(fakechan,c,mcFactors[c]) for c in modes[1:]])

        fakeMap = self.fakeMapMC if mcFake else self.fakeMap
        if shiftUp: fakeMap = self.fakeMapMCUp if mcFake else self.fakeMapUp
        if shiftDown: fakeMap = self.fakeMapMCDown if mcFake else self.fakeMapDown

        # matrix method
        #coeff = '1./(({0}-{1})*({2}-{3})*({4}-{5}))'.format(self.effMap[0],fakeMap[0],self.effMap[1],fakeMap[1],self.effMap[2],fakeMap[2])
        coeff = '1./({0})'.format('*'.join(['({0}-{1})'.format(self.effMap[i],fakeMap[i]) for i in range(self.numLeptons)]))
        fakefactors = {}
        for mode in modes:
            factors = [coeff]
            for i,f in enumerate(mode):
                factors += ['(1-{0})'.format(fakeMap[i]) if f=='P' else fakeMap[i]]
            factor = '*'.join(factors)
            if mode.count('F') in [1,3]: factor += '*(-1.)'
            fakefactors[mode] = factor
        # prompt lepton scale factor
        matrixMethodScaleFactor = '*'.join(['({0}=="{1}" ? ({2}) : (1.))'.format(fakechan,c,fakefactors[c]) for c in modes])

        # fake rate method
        fakefactors = {}
        for mode in modes:
            mapMode = 'l'
            if electronOnly: mapMode='e'
            if muonOnly: mapMode='m'
            factors = ['{0}'.format(modeMap[mapMode][singleChannel][mode] if singleChannel else 1.)]
            if mode not in modesToRun: factors = ['0']
            for i,f in enumerate(mode):
                if f=='F': factors += ['({0}/(1-{0}))'.format(fakeMap[i])]
                # my single bin fakerates
                #if f=='F':
                #    if i in [0,1]: # z1 medium
                #        factors += ['(z1Flv=="ee" ? 0.14617 : 0.45367)']
                #    else: # w1 tight
                #        factors += ['(w1Flv=="e" ? 0.05382 : 0.45367)']
            factor = '*'.join(factors)
            if mode.count('F') in [2]: factor += '*(-1.)'
            fakefactors[mode] = factor
        # prompt lepton scale factor
        fakerateMethodScaleFactor = '*'.join(['({0}=="{1}" ? ({2}) : (1.))'.format(fakechan,c,fakefactors[c]) for c in modes[1:]])

        ## try it a different way:
        #mapMode = 'l'
        #if electronOnly: mapMode='e'
        #if muonOnly: mapMode='m'
        #for l in range(self.numLeptons):
        #    factors = 

        fakeScaleFactor = '1'
        if self.fakeMode == 'matrix': fakeScaleFactor = matrixMethodScaleFactor
        if self.fakeMode == 'fakerate': fakeScaleFactor = fakerateMethodScaleFactor
        #if not getPrompt: fakeScaleFactor += '*(-1.)'

        mcLooseCut = '{0} && {1}'.format(self.mcCut,cutLoose)
        promptCut = '{0} && {1}'.format(self.mcCut,cut)
        fullMCScale = mcScaleFactor + '*' + fakeScaleFactor + '*(-1.)'
        nonPromptScale = mcScaleFactor+'*(-1.)'

        return (cutLoose,promptCut,mcLooseCut,fakeScaleFactor,fullMCScale,nonPromptScale)

    def __buildAsync_getDataDrivenHistJobs(self,variables,binning,cut,**kwargs):
        getPrompt = kwargs.pop('getPrompt',False)
        jobs = []

        (cutLoose,promptCut,mcLooseCut,fakeScaleFactor,fullMCScale,nonPromptScale) = self.__getDataDrivenHistParameters(cut,getPrompt=getPrompt,**kwargs)

        # get contribution from data
        for sample in self.data:
            jobs += self.__buildAsync_getMergedHistJobs(sample,variables,binning,cutLoose,customScale=fakeScaleFactor,**kwargs)
            # if not prompt, subtract from signal region
            #if not getPrompt:
            #    jobs += self.__buildAsync_getMergedHistJobs(sample,variables,binning,cut,**kwargs)

        # subtract contribution from signal MC if fakerate method used
        if self.fakeMode == 'fakerate':
            for sample in self.backgrounds:
                if 'datadriven' in sample: continue
                if 'HPlusPlus' in sample: continue # dont include signal in subtraction
                jobs += self.__buildAsync_getMergedHistJobs(sample,variables,binning,mcLooseCut,customScale=fullMCScale,**kwargs)
                #if not getPrompt:
                #    jobs += self.__buildAsync_getMergedHistJobs(sample,variables,binning,promptCut,customScale=nonPromptScale,**kwargs)

        return jobs

    def __buildAsync_getDataDrivenHistJobs_mcClosure(self,variables,binning,cut,**kwargs):
        getPrompt = kwargs.pop('getPrompt',False)
        jobs = []

        (cutLoose,promptCut,mcLooseCut,fakeScaleFactor,fullMCScale,nonPromptScale) = self.__getDataDrivenHistParameters(cut,getPrompt=getPrompt,mcFake=True,**kwargs)

        # get contribution from control regions
        for sample in self.backgrounds:
            jobs += self.__buildAsync_getMergedHistJobs(sample,variables,binning,cutLoose,customScale=fullMCScale+'*(-1.)',**kwargs)

        return jobs

    def __buildAsync_getMergedHistJobs(self,sample,variables,binning,cut,**kwargs):
        jobs = []
        if isinstance(variables, basestring): variables = [variables]
        if isinstance(cut, basestring): cut = [cut]
        for v in range(len(variables)):
            if sample in self.sampleMergeDict:
                for s in self.sampleMergeDict[sample]:
                    if len(variables) != len(cut):
                        thisCut = '{0} && {1}'.format(cut[0],self.sampleMergeDict[sample][s])
                    else:
                        thisCut = '{0} && {1}'.format(cut[v],self.sampleMergeDict[sample][s])
                    if 'data' not in s: thisCut += ' && {0}'.format(self.mcCut)
                    jobs += [self.__buildAsync_getHistJob(s,variables[v],binning,thisCut,**kwargs)]
            else:
                if len(variables) != len(cut):
                    thisCut = '{0} && {1}'.format(cut[0],self.mcCut) if 'data' not in sample else cut[0]
                else:
                    thisCut = '{0} && {1}'.format(cut[v],self.mcCut) if 'data' not in sample else cut[v]
                if 'datadriven' in sample: 
                    jobs += self.__buildAsync_getDataDrivenHistJobs(variables[v],binning,thisCut,**kwargs) # need to do some special stuff
                elif 'mcClosure' in sample: 
                    jobs += self.__buildAsync_getDataDrivenHistJobs_mcClosure(variables[v],binning,thisCut,**kwargs) # for the closure test
                else:
                    jobs += [self.__buildAsync_getHistJob(sample,variables[v],binning,thisCut,**kwargs)]
        return jobs

    def __buildAsync_getHistJob(self,sample,variable,binning,cut,**kwargs):
        customScale = kwargs.pop('customScale','')
        tree = self.samples[sample]['tree']
        self.j += 1
        histname = 'h_{0}_{1}'.format(sample,self.j)
        scalefactor = self.scaleFactor
        if 'data' in sample: scalefactor = self.dataScaleFactor
        if customScale: scalefactor = customScale
        if 'data' not in sample:
            scalefactor = '{0}*{1}'.format(scalefactor, float(self.intLumi)/self.samples[sample]['lumi'])
        job = (histname,tree,variable,binning,scalefactor,cut)
        return job

    def __getHist_async(self,sample,variables,binning,cut,noFormat=False,**kwargs):
        normalize = kwargs.pop('normalize',False)
        #print sample
        jobs = self.__buildAsync_getMergedHistJobs(sample,variables,binning,cut,**kwargs)
        self.j += 1
        histname = 'h_{0}_getHist_{1}'.format(sample,self.j)
        result = getMergedHist(histname, jobs)
        hist = result[1]
        if not hist: return 0
        #print '    ', histname, hist.Integral()
        hist.SetTitle(self.dataStyles[sample]['name'])
        if sample in self.data: return hist
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
            #hist = self.getHist(sample, variables, binning, cut, noFormat, **kwargs)
            hist = self.__getHist_async(sample,variables,binning,cut,**kwargs)
            hists.Add(hist)
        histname = 'h%s_data' % variables[0].replace('(','_').replace(')','_')
        hist = hists[0].Clone(histname)
        hist.Reset()
        hist.Merge(hists)
        #hist = self.getPoissonErrors(hist)
        return hist

    def getHist(self, sample, variables, binning, cut, noFormat=False, **kwargs):
        return self.__getHist_async(sample,variables,binning,cut,**kwargs)

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
                #continue
                entries = 0
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

    def getSystematicUncertainty(self,val,sample,**kwargs):
        chan = kwargs.pop('chan','wz')
        unc = getSystUncertaintyMap(self.analysis,self.region,self.period,sample)
        totsyst2 = 0.
        for u in unc:
            totsyst2 += unc[u]**2
        totsyst = totsyst2**0.5
        return val*totsyst

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
        histname = kwargs.pop('histname','mcstack')
        mcstack = ROOT.THStack('hs%s' % variables[0],histname)
        plotsig = kwargs.pop('plotsig',False)
        samples = self.backgrounds
        if plotsig: samples = self.backgrounds + self.signal
        for sample in samples:
            hist = self.__getHist_async(sample,variables,binning,cut,**kwargs)
            if not hist: continue
            if nostack:
                hist.SetFillStyle(0)
                hist.SetLineWidth(2)
            #histsyst = self.addSystematicUncertainty(hist,sample)
            #mcstack.Add(histsyst)
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
            if entries <= 0:
                entries = 0
            #if denomentries <= 0:
            #    continue
            ey_low = entries - 0.5 * chisqr(0.1586555, 2. * entries)
            ey_high = 0.5 * chisqr(
                1. - 0.1586555, 2. * (entries + 1)) - entries
            ex = num.GetBinWidth(bin+1) / 2.
            if denomentries > 0:
                graph.SetPoint(npoints, num.GetBinCenter(bin+1), num.GetBinContent(bin+1)/denomentries)
                graph.SetPointEXlow(npoints, 0)
                graph.SetPointEXhigh(npoints, 0)
                graph.SetPointEYlow(npoints, ey_low/denomentries)
                graph.SetPointEYhigh(npoints, ey_high/denomentries)
            else:
                graph.SetPoint(npoints, num.GetBinCenter(bin+1), 1)
                graph.SetPointEXlow(npoints, 0)
                graph.SetPointEXhigh(npoints, 0)
                graph.SetPointEYlow(npoints, 0)
                graph.SetPointEYhigh(npoints, 0)
            npoints += 1
        graph.Set(npoints)
        return graph


    def get_ratio_stat_err(self, hist, **kwargs):
        '''Return a statistical error bars for a ratio plot'''
        ratiomin = kwargs.pop('ratiomin',0.5)
        ratiomax = kwargs.pop('ratiomax',1.5)
        ratiostaterr = hist.Clone("ratiostaterr")
        ratiostaterr.Sumw2()
        ratiostaterr.SetStats(0)
        ratiostaterr.SetTitle("")
        ratiostaterr.GetYaxis().SetTitle("Data / MC")
        ratiostaterr.SetMaximum(ratiomax)
        ratiostaterr.SetMinimum(ratiomin)
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
                ratiostaterr.SetBinError(i, 0.)

        return ratiostaterr

    def setStyle(self,pad,position=11,plotdata=True,preliminary=True):
        '''Set style for plots based on the CMS TDR style guidelines.
           https://twiki.cern.ch/twiki/bin/view/CMS/Internal/PubGuidelines#Figures_and_tables
           https://ghm.web.cern.ch/ghm/plots/'''
        # set period (used in CMS_lumi)
        # period : sqrts
        # 1 : 7, 2 : 8, 3 : 7+8, 4 : 13, ... 7 : 7+8+13
        self.period_int = 1*self.plot7TeV + 2*self.plot8TeV + 4*self.plot13TeV
        CMS_lumi.writeExtraText = preliminary
        CMS_lumi.extraText = "Preliminary" if plotdata else "Simulation Preliminary"
        CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (float(self.intLumi)/1000.)
        if self.intLumi < 1000:
            CMS_lumi.lumi_7TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
            CMS_lumi.lumi_8TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
            CMS_lumi.lumi_13TeV = "%0.1f pb^{-1}" % (float(self.intLumi))
        CMS_lumi.CMS_lumi(pad,self.period_int,position)

    def getLegend(self,mchist,datahist,sighists,**kwargs):
        legendpos = kwargs.pop('legendpos',43)
        plotdata = kwargs.pop('plotdata',True)
        plotsig = kwargs.pop('plotsig',True)
        plotmc = kwargs.pop('plotmc',True)
        plotratio = kwargs.pop('plotratio',True)
        numcol = kwargs.pop('numcol',1)

        numEntries = 0
        if mchist: numEntries += len(self.backgrounds) + len(self.signal)
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
            yend = 0.85
        else:                     # default, top, just below CMS label
            yend = 0.77
        xstart = xend-0.17*numcol-0.1
        ystart = yend-math.ceil(float(numEntries)/numcol)*0.07
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
            for s in sighists:
                leg.AddEntry(s,s.GetTitle(),'f')
        return leg

    def save(self, canvas, savename):
        '''Save the canvas in multiple formats.'''
        canvas.SetName(savename)
        if self.dontSave: return
        #for type in ['png','root','pdf']:
        for type in ['pdf','root','png']:
            name = "%s/%s/%s.%s" % (self.plotDir, type, savename, type)
            python_mkdir(os.path.dirname(name))
            canvas.Print(name)
        self.savefile.WriteTObject(canvas)

pool = Pool(16)
