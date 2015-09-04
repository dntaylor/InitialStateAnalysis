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

    def getSummedCutSampleEfficiency(self,sampleCategory,numCut,denomCut,**kwargs):
        '''Get efficiency for many samples.'''
        if sampleCategory=='bg':
            method = getattr(self,'getBackgroundEntries')
        elif sampleCategory=='sig':
            method = getattr(self,'getSignalEntries')
        elif sampleCategory=='data':
            method = getattr(self,'getDataEntries')
        else:
            return 0,0
        num, nErr = method(numCut,doError=True)
        den, dErr = method(denomCut,doError=True)
        if den:
            eff = float(num)/den
            if num:
                err = eff * (nErr**2/(num**2) + dErr**2/(den**2)) ** 0.5
            else:
                err = 0.
            return eff, err
        return 0, 0

    def getEfficiencyCutflow(self,sampleCategory,cutflow,**kwargs):
        '''Return histogram of efficiency v cutflow.'''
        effMode = kwargs.pop('effMode','default') # default, relative, pre, post
        effs = []
        errors = []
        for c in range(len(cutflow)):
            if effMode == 'default':
                numCut = ' && '.join(cutflow[:c+1])
                denomCut = cutflow[0]
            elif effMode == 'relative':
                numCut = ' && '.join(cutflow[:c+1]) if c else cutflow[0]
                denomCut = ' && '.join(cutflow[:c]) if c else cutflow[0]
            elif effMode == 'pre':
                numCut = '%s && %s' % (cutflow[0], cutflow[c]) if c else cutflow[0]
                denomCut = cutflow[0]
            elif effMode == 'post':
                numCut = ' && '.join(cutflow)
                denomCut = ' && '.join(cutflow[:c] + cutflow[c+1:])
            else:
                print 'Unrecognized Efficieny Mode'
                numCut = ''
                denomCut = ''
            eff, err = self.getSummedCutSampleEfficiency(sampleCategory,numCut,denomCut)
            effs.append(eff)
            errors.append(err)
        hist = ROOT.TH1F('hEffCutFlow', 'EfficiencyCutFlow', len(effs),0,len(effs))
        for eff in range(len(effs)):
            hist.SetBinContent(eff+1,effs[eff])
            #hist.SetBinError(eff+1,errors[eff])
            hist.SetBinError(eff+1,0.00001)
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
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        self.canvas.SetLogy(logy)

        numSelections = len(selections)
        bgEff = self.getEfficiencyCutflow('bg',selections)
        bgEff.SetLineColor(ROOT.kBlue)
        bgEff.SetMarkerStyle(0)
        bgEff.Draw('hist')
        bgEff.GetYaxis().SetTitle('Efficiency')
        bgEff.GetYaxis().SetTitleOffset(1)
        bgEff.SetMaximum(1.25)
        bgEff.SetMinimum(0)
        if labels:
            for bin in range(numSelections):
                bgEff.GetXaxis().SetBinLabel(bin+1,labels[bin])
        bgEffPre = self.getEfficiencyCutflow('bg',selections,effMode='pre')
        bgEffPre.SetLineColor(ROOT.kBlue)
        bgEffPre.SetLineStyle(2)
        bgEffPre.SetMarkerStyle(0)
        bgEffPre.Draw('same')
        bgEffPost = self.getEfficiencyCutflow('bg',selections,effMode='post')
        bgEffPost.SetLineColor(ROOT.kBlue)
        bgEffPost.SetLineStyle(3)
        bgEffPost.SetMarkerStyle(0)
        bgEffPost.Draw('same')

        if plotsig:
            sigEff = self.getEfficiencyCutflow('sig',selections)
            sigEff.SetLineColor(ROOT.kRed)
            sigEff.SetMarkerStyle(0)
            sigEff.Draw('hist same')
            sigEffPre = self.getEfficiencyCutflow('sig',selections,effMode='pre')
            sigEffPre.SetLineColor(ROOT.kRed)
            sigEffPre.SetLineStyle(2)
            sigEffPre.SetMarkerStyle(0)
            sigEffPre.Draw('same')
            sigEffPost = self.getEfficiencyCutflow('sig',selections,effMode='post')
            sigEffPost.SetLineColor(ROOT.kRed)
            sigEffPost.SetLineStyle(3)
            sigEffPost.SetMarkerStyle(0)
            sigEffPost.Draw('same')

        if plotdata: 
            dataEff = self.getEfficiencyCutflow('data',selections)
            dataEff.SetLineColor(ROOT.kBlack)
            dataEff.Draw('esamex0')

        # draw cms lumi
        self.setStyle(lumitext,plotdata,plotratio,isprelim)
        self.canvas.cd()
        self.canvas.Update()
        self.canvas.RedrawAxis()
        frame = self.canvas.GetFrame()
        frame.Draw()

        # legend
        legend = ROOT.TLegend(0.3,0.77,0.6,0.90)
        legend.SetFillColor(0)
        legend.AddEntry(bgEff,'Background Efficiency')
        legend.AddEntry(bgEffPre,'Background Efficiency (first cut)')
        legend.AddEntry(bgEffPost,'Background Efficiency (last cut)')
        if plotsig:
            legend.AddEntry(sigEff,'Signal Efficiency')
            legend.AddEntry(sigEffPre,'Signal Efficiency (first cut)')
            legend.AddEntry(sigEffPost,'Signal Efficiency (last cut)')
        if plotdata: legend.AddEntry(dataEff,'Data Efficiency')
        legend.Draw('same')

        # save everything
        self.canvas.cd()
        self.save(savename)


