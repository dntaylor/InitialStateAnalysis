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
from array import array

from xsec import xsecs
from dataStyles import dataStyles
import CMS_lumi, tdrstyle
from plotUtils import *

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()

class PlotterBase(object):
    '''A Base plotting class to be used with flat histograms.'''
    def __init__(self,analysis,**kwargs):
        '''Initialize the plotter (optionally make the plots blinded).'''
        # get kwargs
        saveDir = kwargs.pop('saveDir','')
        ntupleDir = kwargs.pop('ntupleDir','ntuples')
        period = kwargs.pop('period',13)
        blind = kwargs.pop('blind',False)
        rootName = kwargs.pop('rootName','plots')
        mergeDict = kwargs.pop('mergeDict',{})
        for key, value in kwargs.iteritems():
            print "Unrecognized parameter '" + key + "' = " + str(value)

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

    def reset(self):
        '''Reset the plotter class'''
        print "Resetting the PlotterBase class"
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
            print "No data initialized, default to 25 fb-1"
            self.intLumi = 25000.
            return
        self.intLumi = 25000.

    def setIntLumi(self,intLumi):
        '''Set the integrated luminosity to scale MC to'''
        self.intLumi = intLumi

    def getNumEntries(self,selection,sample,**kwargs):
        '''Return the lumi scaled number of entries passing a given cut.'''
        doError = kwargs.pop('doError',False)
        scaleup = kwargs.pop('scaleup',False)
        unweighted = kwargs.pop('doUnweighted',False)
        totalVal = 0
        totalErr2 = 0
        if sample in self.sampleMergeDict:
            for s in self.sampleMergeDict[sample]:
                tree = self.samples[s]['file'].Get(self.analysis)
                if 'data' not in s and not unweighted:
                    #if scaleup: tree.Draw('event.pu_weight>>h%s()'%s,'event.lep_scale_up*event.trig_scale*(%s)' %selection,'goff')
                    #if not scaleup: tree.Draw('event.pu_weight>>h%s()'%s,'event.lep_scale*event.trig_scale*(%s)' %selection,'goff')
                    tree.Draw('event.pu_weight>>h%s()'%s,'event.lep_scale*(%s)' %selection,'goff')
                    if not ROOT.gDirectory.Get("h%s" %s):
                        val = 0
                    else:
                        hist = ROOT.gDirectory.Get("h%s" %s).Clone("hnew%s" %s)
                        hist.Sumw2()
                        val = hist.Integral()
                    err = val ** 0.5
                    lumi = self.samples[s]['lumi']
                    val = val * self.intLumi/lumi
                    err = err * self.intLumi/lumi
                else:
                    val = tree.GetEntries(selection)
                    err = val ** 0.5
                totalVal += val
                totalErr2 += err*err
        else:
            tree = self.samples[sample]['file'].Get(self.analysis)
            if 'data' not in sample and not unweighted:
                #if scaleup: tree.Draw('event.pu_weight>>h%s()'%sample,'event.lep_scale_up*event.trig_scale*(%s)' %selection,'goff')
                #if not scaleup: tree.Draw('event.pu_weight>>h%s()'%sample,'event.lep_scale*event.trig_scale*(%s)' %selection,'goff')
                tree.Draw('event.pu_weight>>h%s()'%sample,'event.lep_scale*(%s)' %selection,'goff')
                if not ROOT.gDirectory.Get("h%s" %sample):
                    val = 0
                else:
                    hist = ROOT.gDirectory.Get("h%s" %sample).Clone("hnew%s" %sample)
                    hist.Sumw2()
                    val = hist.Integral()
                err = val ** 0.5
                lumi = self.samples[sample]['lumi']
                val = val * self.intLumi/lumi
                err = err * self.intLumi/lumi
            else:
                val = tree.GetEntries(selection)
                err = val ** 0.5
            totalVal += val
            totalErr2 += err*err
        totalErr = totalErr2 ** 0.5
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
        for i in range(nx):
            htmp.Fill(htmp.GetBinCenter(i+1), hist.GetBinContent(i+1))
        htmp.Fill(hist.GetBinLowEdge(1)-1, hist.GetBinContent(0))
        htmp.SetEntries(hist.GetEntries())
        return htmp

    def getSingleVarHist(self,tree,sample,variable,binning,cut):
        '''Single variable, single sample hist'''
        self.j += 1
        if 'data' not in sample: lumi = self.samples[sample]['lumi']
        if len(binning) == 3: # standard drawing
            drawString = "%s>>h%s%s(%s)" % (variable, sample, variable, ", ".join(str(x) for x in binning))
        else: # we will need to rebin
            drawString = "%s>>h%s%s()" % (variable, sample, variable)
        if not cut: cut = '1'
        if 'data' not in sample and self.sqrts != 13: # TODO: dont forget to remove when we have data!
            tree.Draw(drawString,'(event.pu_weight*event.lep_scale*event.trig_scale)*('+cut+')','goff')
        else:
            tree.Draw(drawString,cut,'goff')
        if not ROOT.gDirectory.Get("h%s%s" %(sample, variable)):
            return 0
        hist = ROOT.gDirectory.Get("h%s%s" %(sample, variable)).Clone("hmod%s%s"%(sample,variable))
        if len(binning) != 3: # variable binning (list of bin edges
            hist.Rebin(len(binning)-1,"hnew%s%s" %(sample,variable),array('d',binning))
            hist = ROOT.gDirectory.Get("hnew%s%s" %(sample,variable)).Clone("hnewmod%s%s"%(sample,variable))
        if 'data' not in sample: # if it is mc, scale to intLumi
            hist.Scale(self.intLumi/lumi)
        return hist

    def getHist(self, sample, variables, binning, cut, noFormat=False, **kwargs):
        '''Return a histogram of a given variable from the given dataset with a cut'''
        normalize = kwargs.pop('normalize',False)
        hists = ROOT.TList()
        for v in range(len(variables)):
            if sample in self.sampleMergeDict:
                for s in self.sampleMergeDict[sample]:
                    tree = self.samples[s]['file'].Get(self.analysis)
                    if len(variables) != len(cut):
                        hist = self.getSingleVarHist(tree,s,variables[v],binning,cut)
                    else:
                        hist = self.getSingleVarHist(tree,s,variables[v],binning,cut[v])
                    if hist:
                        hists.Add(hist)
            else:
                tree = self.samples[sample]['file'].Get(self.analysis)
                if len(variables) != len(cut):
                    hist = self.getSingleVarHist(tree,sample,variables[v],binning,cut)
                else:
                    hist = self.getSingleVarHist(tree,sample,variables[v],binning,cut[v])
                if hist:
                    hists.Add(hist)
        if hists.IsEmpty():
            return 0
        hist = hists[0].Clone("hmerged%s%s" % (sample, variables[0]))
        hist.Reset()
        hist.Merge(hists)
        hist = self.getOverflowUnderflow(hist,**kwargs)
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

    def getData(self, variables, binning, cut, noFormat=False, **kwargs):
        '''Return a histogram of data for the given variable'''
        hists = ROOT.TList()
        for sample in self.data:
            hist = self.getHist(sample, variables, binning, cut, noFormat, **kwargs)
            hists.Add(hist)
        hist = hists[0].Clone("hdata%s" % variables[0])
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
        return mcstack

    def get_stat_err(self, hist):
        '''Create statistical errorbars froma histogram'''
        staterr = hist.Clone("staterr")
        staterr.Sumw2()
        staterr.SetFillColor(ROOT.EColor.kGray+3)
        staterr.SetLineColor(ROOT.EColor.kGray+3)
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
        ratiostaterr.SetMaximum(2.)
        ratiostaterr.SetMinimum(0)
        ratiostaterr.SetMarkerSize(0)
        ratiostaterr.SetFillColor(ROOT.EColor.kGray+3)
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
        CMS_lumi.CMS_lumi(self.plotpad if plotratio else self.canvas,self.period,position)

    def drawLegend(self,plotdata,plotsignal,plotratio,legendpos):
        '''Plot a legend using the new PUB COMM recommendations.'''
        # initialize latex for labels and legend
        latex = ROOT.TLatex()
        # setup number of things to include in legend
        n_ = len(self.backgrounds)
        if plotdata: n_ += 1
        if plotsignal: n_ += len(self.signal)
        # create legend position
        if legendpos % 10 == 1: # on the left
            x1_l = 0.45
        elif legendpos % 10 == 2: # in the middle
            x1_l = 0.70
        else: # default (on right)
            x1_l = 0.95
        dx_l = 0.3
        x0_l = x1_l-dx_l
        if legendpos//10 == 1: # bottom
            y1_l = 0.35
        elif legendpos//10 == 2: # middle
            y1_l = 0.59
        elif legendpos//10 == 4: # very top (in line with CMS label)
            y1_l = 0.91
        else: # default, top, just below CMS label
            y1_l = 0.77
        if plotratio: y1_l *= 0.95
        dy_l = 0.045*n_ # number of things in legend
        #if plotratio: dy_l *= 0.85
        y0_l = y1_l-dy_l
        legend = ROOT.TPad("legend_0","legend_0",x0_l,y0_l,x1_l,y1_l)
        legend.Draw()
        # now setup the position of the boxes in the legend
        legend.cd()
        ar_l = dy_l/dx_l
        gap_ = 1./(n_+1)
        bwx_ = 0.12
        bwy_ = gap_/1.5
        x_l = [1.2*bwx_]
        y_l = [1-gap_]
        ex_l = [0]
        ey_l = [0.04/ar_l]
        x_l = array("d",x_l)
        ex_l = array("d",ex_l)
        y_l = array("d",y_l)
        ey_l = array("d",ey_l)
        latex.SetTextFont(42)
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.EColor.kBlack)    
        latex.SetTextSize(0.7/n_) 
        #latex.SetTextSize(0.25) 
        latex.SetTextAlign(12) 
        # plot data label
        if plotdata:
            gr_l =  ROOT.TGraphErrors(1, x_l, y_l, ex_l, ey_l)
            ROOT.gStyle.SetEndErrorSize(0)
            gr_l.SetMarkerStyle(20)
            gr_l.SetMarkerSize(0.9)
            gr_l.Draw("0P")
            xx_ = x_l[0]
            yy_ = y_l[0]
            latex.DrawLatex(xx_+1.*bwx_,yy_,"Data")
        # plot mc labels
        allMC = self.backgrounds + self.signal if plotsignal else self.backgrounds
        box_ = ROOT.TBox()
        for s in range(len(allMC)):
            xx_ = x_l[0]
            yy_ = y_l[0]
            yy_ -= gap_*(s+1) if plotdata else gap_*s
            box_.SetFillColor( self.dataStyles[allMC[s]]['fillcolor'] )
            if allMC[s] in self.signal:
                box_.SetFillStyle( 0 )
            else:
                box_.SetFillStyle( self.dataStyles[allMC[s]]['fillstyle'] )
            box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
            box_.SetLineStyle( ROOT.kSolid )
            if allMC[s] in self.signal:
                box_.SetLineWidth( 2 )
            else:
                box_.SetLineWidth( 1 )
            box_.SetLineColor( self.dataStyles[allMC[s]]['linecolor'] )
            box_.SetFillStyle( 0 )
            box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
            latex.DrawLatex(xx_+1.*bwx_,yy_,self.dataStyles[allMC[s]]['name'])
        self.canvas.Update()

    def save(self, savename):
        '''Save the canvas in multiple formats.'''
        for type in ['png']:
            name = "%s/%s/%s.%s" % (self.plotDir, type, savename, type)
            python_mkdir(os.path.dirname(name))
            self.canvas.Print(name)
        self.canvas.SetName(savename)
        self.savefile.WriteTObject(self.canvas)
        self.canvas.Clear()

