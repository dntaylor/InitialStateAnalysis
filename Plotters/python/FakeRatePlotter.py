'''
A cut flow plotter class.

Author: Devin N. Taylor, UW-Madison
'''

from PlotterBase import PlotterBase
from array import array
import ROOT
import math

ROOT.gStyle.SetPalette(1)

class FakeRatePlotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    def getFakeRateProjection(self, numString, denomString, bins, var, savename, **kwargs):
        '''Get 1d histogram of fakerates'''
        dataDriven = kwargs.pop('dataDriven',True)
        subtractSamples = kwargs.pop('subtractSamples',[])
        fakeHist = ROOT.TH1F(savename,'',len(bins)-1,array('d',bins))
        for b in range(len(bins)-1):
            #print '{0}: [{1},{2}]'.format(var,bins[b],bins[b+1])
            kinCut = 'abs({var})>={lowbin} && abs({var})<{highbin}'.format(var=var,lowbin=bins[b],highbin=bins[b+1])
            numCut = '{0} && {1}'.format(numString,kinCut)
            denomCut = '{0} && {1}'.format(denomString,kinCut)
            #print numCut
            #print denomCut
            num = 0
            denom = 0
            numErr2 = 0
            denomErr2 = 0
            samples = self.data if dataDriven else self.backgrounds
            for sample in samples:
                n, nErr = self.getNumEntries(numCut, sample, doError=True)
                d, dErr = self.getNumEntries(denomCut, sample, doError=True)
                #print 'Sample {0}: num: {1}; denom: {2}'.format(sample,n,d)
                num += n
                numErr2 += nErr ** 2
                denom += d
                denomErr2 += dErr ** 2
            for sample in subtractSamples:
                n, nErr = self.getNumEntries(numCut, sample, doError=True)
                d, dErr = self.getNumEntries(denomCut, sample, doError=True)
                #print 'Subtract Sample {0}: num: {1}; denom: {2}'.format(sample,n,d)
                num -= n
                numErr2 += nErr ** 2
                denom -= d
                denomErr2 += dErr ** 2
            if num < 0: num = 0
            if denom and num:
                fakerate = float(num)/denom
                err = fakerate * (numErr2/(num**2) + denomErr2/(denom**2)) ** 0.5
            else:
                fakerate = 0
                err = 0
            fakeHist.SetBinContent(b+1,fakerate)
            fakeHist.SetBinError(b+1,err)
        return fakeHist

    def getFakeRate(self,numString, denomString, ptBins, etaBins, ptVar, etaVar, savename, **kwargs):
        '''Get 2d histogram of fakerates'''
        dataDriven = kwargs.pop('dataDriven',True)
        subtractSamples = kwargs.pop('subtractSamples',[])
        fakeHist = ROOT.TH2F(savename,'',len(ptBins)-1,array('d',ptBins),len(etaBins)-1,array('d',etaBins))
        for p in range(len(ptBins)-1):
            for e in range(len(etaBins)-1):
                #print '{0}: [{1},{2}]; {3}: [{4},{5}]'.format(ptVar,ptBins[p],ptBins[p+1],etaVar,etaBins[e],etaBins[e+1])
                kinCut = '%s>=%f & %s<%f & abs(%s)>=%f & abs(%s)<%f' %\
                         (ptVar, ptBins[p], ptVar, ptBins[p+1], etaVar, etaBins[e], etaVar, etaBins[e+1])
                numCut = '%s && %s' % (kinCut, numString)
                denomCut = '%s && %s' % (kinCut, denomString)
                #print numCut
                #print denomCut
                num = 0
                denom = 0
                numErr2 = 0
                denomErr2 = 0
                samples = self.data if dataDriven else self.backgrounds
                for sample in samples:
                    n, nErr = self.getNumEntries(numCut, sample, doError=True)
                    d, dErr = self.getNumEntries(denomCut, sample, doError=True)
                    #print 'Sample {0}: num: {1}; denom: {2}'.format(sample,n,d)
                    num += n
                    numErr2 += nErr ** 2
                    denom += d
                    denomErr2 += dErr ** 2
                for sample in subtractSamples:
                    n, nErr = self.getNumEntries(numCut, sample, doError=True)
                    d, dErr = self.getNumEntries(denomCut, sample, doError=True)
                    #print 'Subtract Sample {0}: num: {1}; denom: {2}'.format(sample,n,d)
                    num -= n
                    numErr2 += nErr ** 2
                    denom -= d
                    denomErr2 += dErr ** 2
                if num < 0: num = 0
                if denom and num:
                    fakerate = float(num)/denom
                    err = fakerate * (numErr2/(num**2) + denomErr2/(denom**2)) ** 0.5
                else:
                    fakerate = 0
                    err = 0
                fakeHist.SetBinContent(p+1,e+1,fakerate)
                fakeHist.SetBinError(p+1,e+1,err)
        return fakeHist

    def plotFakeRate(self, passSelection, failSelection, savename, **kwargs):
        '''A function to calculate and plot the fake rate for a given selection.
           kwargs accepts:
               cut         string           applied with all selections
               ptBins      list (float)     list of pt bin edges for fakerate
               etaBins     list (float)     list of eta bin edges for fakerate
               ptVar       string           probe pt variable
               etaVar      string           probe eta variable
               logy        bool             set logy plot
               logx        bool             set logx plot
               lumitext    int              location of lumitext (from CMS_lumi)
               isprelim    bool             The plot is CMS preliminary'''
        cut = kwargs.pop('cut', '1')
        logy = kwargs.pop('logy', 0)
        logx = kwargs.pop('logx', 0)
        lumitext = kwargs.pop('lumitext', 11)
        isprelim = kwargs.pop('isprelim', 1)
        ptBins = kwargs.pop('ptBins', [10,15,20,25,30,35,40,50])
        etaBins = kwargs.pop('etaBins', [0,1,1.479,2,2.5])
        ptVar = kwargs.pop('ptVar','w1.Pt1')
        etaVar = kwargs.pop('etaVar','w1.Eta1')
        xaxis = kwargs.pop('xaxis','p_{T} (GeV)')
        yaxis = kwargs.pop('yaxis','#eta')
        dataDriven = kwargs.pop('dataDriven',True)
        subtractSamples = kwargs.pop('subtractSamples',[])
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        self.canvas.SetLogy(logy)
        self.canvas.SetLogx(logx)
        self.canvas.SetRightMargin(0.14)

        # calculate fake rate
        fakeRateHist = self.getFakeRate(passSelection,failSelection,ptBins,etaBins,ptVar,etaVar,savename,dataDriven=dataDriven,subtractSamples=subtractSamples)
        fakeRateHist.GetXaxis().SetTitle(xaxis)
        fakeRateHist.GetYaxis().SetTitle(yaxis)
        fakeRateHist.GetYaxis().SetTitleOffset(1.)
        fakeRateHist.SetTitle('')
        fakeRateHist.GetZaxis().SetRangeUser(0,1)
        self.savefile.WriteTObject(fakeRateHist)

        # plot fakerate
        fakeRateHist.Draw('colz text')

        # draw cms lumi
        #self.setStyle(lumitext,plotdata,plotratio,isprelim)
        #self.canvas.cd()
        #self.canvas.Update()
        #self.canvas.RedrawAxis()
        #frame = self.canvas.GetFrame()
        #frame.Draw()

        # save everything
        self.canvas.cd()
        self.save(savename)
        self.resetCanvas()

    def plotFakeRateProjection(self, passSelection, failSelection, savename, projection, **kwargs):
        '''A function to calculate and plot the fake rate for a given selection.
           kwargs accepts:
               cut         string           applied with all selections
               ptBins      list (float)     list of pt bin edges for fakerate
               etaBins     list (float)     list of eta bin edges for fakerate
               ptVar       string           probe pt variable
               etaVar      string           probe eta variable
               logy        bool             set logy plot
               logx        bool             set logx plot
               lumitext    int              location of lumitext (from CMS_lumi)
               isprelim    bool             The plot is CMS preliminary'''
        cut = kwargs.pop('cut', '1')
        logy = kwargs.pop('logy', 0)
        logx = kwargs.pop('logx', 0)
        lumitext = kwargs.pop('lumitext', 11)
        isprelim = kwargs.pop('isprelim', 1)
        ptBins = kwargs.pop('ptBins', [10,15,20,25,30,35,40,50])
        etaBins = kwargs.pop('etaBins', [0,1,1.479,2,2.5])
        ptVar = kwargs.pop('ptVar','w1.Pt1')
        etaVar = kwargs.pop('etaVar','w1.Eta1')
        xaxis = kwargs.pop('xaxis','')
        yaxis = kwargs.pop('yaxis','Fake Rate')
        dataDriven = kwargs.pop('dataDriven',True)
        subtractSamples = kwargs.pop('subtractSamples',[])
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        if projection not in ['pt','eta']:
            self.logger.warning("Must project in either pt or eta")

        notproj = 'eta' if projection=='pt' else 'pt'

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        self.canvas.SetLogy(logy)
        self.canvas.SetLogx(logx)

        # calculate fake rate
        fakeBins, iterBins = (ptBins, etaBins) if projection=='pt' else (etaBins, ptBins)
        fakeVar, iterVar = (ptVar, etaVar) if projection=='pt' else (etaVar, ptVar)
        for b in range(len(iterBins)):
            if b < len(iterBins)-1:
                binCut = 'abs({var}) >= {lowbin} && abs({var}) < {highbin}'.format(var=iterVar,lowbin=iterBins[b],highbin=iterBins[b+1])
                name = '{0}_{1}_{2}'.format(savename,notproj,b)
                thisNum = '{0} && {1}'.format(passSelection,binCut)
                thisDenom = '{0} && {1}'.format(failSelection,binCut)
            else:
                name = savename
                thisNum = passSelection
                thisDenom = failSelection
            fakeRateHist = self.getFakeRateProjection(thisNum,thisDenom,fakeBins,fakeVar,name,dataDriven=dataDriven,subtractSamples=subtractSamples)
            fakeRateHist.GetXaxis().SetTitle(xaxis)
            fakeRateHist.GetYaxis().SetTitle(yaxis)
            fakeRateHist.GetYaxis().SetTitleOffset(1.)
            fakeRateHist.SetTitle('')
            fakeRateHist.GetYaxis().SetRangeUser(0,1)
            #fakeRateHist.SetLineColor(ROOT.kBlue)
            self.savefile.WriteTObject(fakeRateHist)

            # plot fakerate
            fakeRateHist.Draw()

            # plot legend
            leg = ROOT.TLegend(0.65,0.72,0.95,0.77,'','NDC')
            leg.SetTextFont(42)
            leg.SetBorderSize(0)
            leg.SetFillColor(0)
            leg.AddEntry(fakeRateHist,'Fake Rate','ep')
            leg.Draw()

            # draw cms lumi
            self.setStyle(lumitext,True,False,True)
            self.canvas.cd()
            self.canvas.Update()
            self.canvas.RedrawAxis()
            frame = self.canvas.GetFrame()
            frame.Draw()

            # save everything
            self.canvas.cd()
            self.save(name)

