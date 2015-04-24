'''
A cut flow efficiency plotter class.

Author: Devin N. Taylor, UW-Madison
'''

from PlotterBase import PlotterBase
import ROOT
import math

class EfficiencyPlotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    def getSummedCutSampleYield(self,samples,cut,**kwargs):
        '''Get efficiency for many samples.'''
        count = 0
        err2 = 0
        if not isinstance(samples,list): samples = [samples]
        for sample in samples:
            n, nErr = self.getNumEntries(cut,sample,doError=True)
            count += n
            err2 += nErr**2
        return count, err2 ** 0.5

    def getSummedCutSampleEfficiency(self,samples,numCut,denomCut,**kwargs):
        '''Get efficiency for many samples.'''
        numCount = 0
        denomCount = 0
        numErr2 = 0
        denErr2 = 0
        if not isinstance(samples,list): samples = [samples]
        for sample in samples:
            num, nErr = self.getNumEntries(numCut,sample,doError=True)
            den, dErr = self.getNumEntries(denomCut,sample,doError=True)
            numCount += num
            denomCount += den
            numErr2 += nErr**2
            denErr2 += dErr**2
        if denomCount:
            eff = float(numCount)/denomCount
            if numCount:
                err = eff * (numErr2/(numCount**2) + denErr2/(denomCount**2)) ** 0.5
            else:
                err = 0.
            return eff, err
        return 0, 0

    def getEfficiencyCutflow(self,samples,cutflow,**kwargs):
        '''Return histogram of efficiency v cutflow.'''
        effMode = kwargs.pop('effMode','default') # default, relative, pre, post
        effs = []
        errors = []
        for c in range(len(cutflow)):
            if effMode == 'default':
                numCut = ' & '.join(cutflow[:c+1])
                denomCut = cutflow[0]
            elif effMode == 'relative':
                numCut = ' & '.join(cutflow[:c+1]) if c else cutflow[0]
                denomCut = ' & '.join(cutflow[:c]) if c else cutflow[0]
            elif effMode == 'pre':
                numCut = '%s & %s' % (cutflow[0], cutflow[c]) if c else cutflow[0]
                denomCut = cutflow[0]
            elif effMode == 'post':
                numCut = ' & '.join(cutflow)
                numCut = ' & '.join(cutflow[:c] + cutflow[c+1:])
            else:
                print 'Unrecognized Efficieny Mode'
                numCut = ''
                denomCut = ''
            eff, err = self.getSummedCutSampleEfficiency(samples,numCut,denomCut)
            effs.append(eff)
            errors.append(err)
        hist = ROOT.TH1F('hEffCutFlow', 'EfficiencyCutFlow', len(effs),0,len(effs))
        for eff in range(len(effs)):
            hist.SetBinContent(eff+1,effs[eff])
            hist.SetBinError(eff+1,errors[eff])
        return hist
        

    # several aliases for fast plotting configuration
    def plotCutFlowMC(self, selections, savename, **kwargs):
        '''Plot Monte Carlo'''
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        self.plotCutFlowMCDataSignalRatio(selections, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotCutFlowMCSignal(self, selections, savename, **kwargs):
        '''Plot Monte Carlo with signal overlay'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        self.plotCutFlowMCDataSignalRatio(selections, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotCutFlowMCData(self, selections, savename, **kwargs):
        '''Plot Monte Carlo with data'''
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 0)
        self.plotCutFlowMCDataSignalRatio(selections, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotCutFlowMCDataSignal(self, selections, savename, **kwargs):
        '''Plot Monte Carlo with data and signal overlay'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 0)
        self.plotCutFlowMCDataSignalRatio(selections, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotCutFlowMCDataRatio(self, selections, savename, **kwargs):
        '''Plot Monte Carlo with data and a ratio plot'''
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 1)
        self.plotCutFlowMCDataSignalRatio(selections, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotCutFlowMCSignalRatio(self, selections, savename, **kwargs):
        '''Plot Monte Carlo with signal and a ratio plot'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 1)
        self.plotCutFlowMCDataSignalRatio(selections, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotCutFlowMCDataSignalRatio(self, selections, savename, **kwargs):
        '''A function to plot cut flow efficiency of selections. Selections is a list of cuts to apply in order.
           Each subsequent cut is the logical and of the previous.
           kwargs accepts:
               labels      list (string)    bin labels for each cut
               cut         string           applied with all selections
               doYields    bool             do yields rather than efficiencies
               blinder     list (double)    range to blind (2 elements)
               logy        bool             set logy plot
               plotsig     bool             plot signal
               plotdata    bool             plot data
               plotratio   bool             make ratio plot
               lumitext    int              location of lumitext (from CMS_lumi)
               legendpos   int              location of legendtext AB (A=012=LCR, B=012=TMB)
               signalscale int              factor to scale signal by
               nosum       bool             Don't sum cut selections
               isprecf     bool             Do the preselection cutflow
               isprelim    bool             The plot is CMS preliminary'''
        labels = kwargs.pop('labels',[])
        cut = kwargs.pop('cut', '1')
        blinder = kwargs.pop('blinder', [])
        logy = kwargs.pop('logy', 0)
        plotsig = kwargs.pop('plotsig', 1)
        plotdata = kwargs.pop('plotdata', 1)
        plotratio = kwargs.pop('plotratio', 1)
        lumitext = kwargs.pop('lumitext', 11)
        legendpos = kwargs.pop('legendpos', 33)
        signalscale = kwargs.pop('signalscale',1)
        isprecf = kwargs.pop('isprecf',False)
        isprelim = kwargs.pop('isprelim', 1)
        for key, value in kwargs.iteritems():
            print "Unrecognized parameter '" + key + "' = " + str(value)

        self.canvas.SetLogy(logy)

        numSelections = len(selections)
        bgEff = self.getEfficiencyCutflow(self.backgrounds,selections)
        bgEff.SetLineColor(ROOT.EColor.kBlue)
        bgEff.SetMarkerStyle(0)
        bgEff.Draw('hist')
        bgEff.GetYaxis().SetTitle('Efficiency')
        bgEff.GetYaxis().SetTitleOffset(1)
        bgEff.SetMaximum(1.25)
        bgEff.SetMinimum(0)
        if labels:
            for bin in range(numSelections):
                bgEff.GetXaxis().SetBinLabel(bin+1,labels[bin])
        bgEffRel = self.getEfficiencyCutflow(self.backgrounds,selections,effMode='relative')
        bgEffRel.SetLineColor(ROOT.EColor.kBlue)
        bgEffRel.SetLineStyle(2)
        bgEffRel.SetMarkerStyle(0)
        bgEffRel.Draw('same')

        if plotsig:
            sigEff = self.getEfficiencyCutflow(self.signal,selections)
            sigEff.SetLineColor(ROOT.EColor.kRed)
            sigEff.SetMarkerStyle(0)
            sigEff.Draw('hist same')
            sigEffRel = self.getEfficiencyCutflow(self.signal,selections,effMode='relative')
            sigEffRel.SetLineColor(ROOT.EColor.kRed)
            sigEffRel.SetLineStyle(2)
            sigEffRel.SetMarkerStyle(0)
            sigEffRel.Draw('same')

        if plotdata: 
            dataEff = self.getEfficiencyCutflow(self.data,selections)
            dataEff.SetLineColor(ROOT.EColor.kBlack)
            dataEff.Draw('esamex0')

        # draw cms lumi
        self.setStyle(lumitext,plotdata,plotratio,isprelim)
        self.canvas.cd()
        self.canvas.Update()
        self.canvas.RedrawAxis()
        frame = self.canvas.GetFrame()
        frame.Draw()

        # legend
        #self.drawLegend(plotdata,plotsig,plotratio,legendpos)
        legend = ROOT.TLegend(0.3,0.77,0.6,0.90)
        legend.SetFillColor(0)
        legend.AddEntry(bgEff,'Background Efficiency')
        legend.AddEntry(bgEffRel,'Relative Background Efficiency')
        if plotsig:
            legend.AddEntry(sigEff,'Signal Efficiency')
            legend.AddEntry(sigEffRel,'Relative Signal Efficiency')
        if plotdata: legend.AddEntry(dataEff,'Data Efficiency')
        legend.Draw('same')

        # save everything
        self.canvas.cd()
        self.save(savename)


