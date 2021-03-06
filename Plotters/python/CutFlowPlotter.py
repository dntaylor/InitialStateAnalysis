'''
A cut flow plotter class.

Author: Devin N. Taylor, UW-Madison
'''

from PlotterBase import PlotterBase
import ROOT
import math

class CutFlowPlotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    def writeCutString(self,sample,cutList):
        '''A function to write out the number of events for a sample after a series of selections'''
        cutString = '{:20s}'.format(sample if len(sample)<20 else sample[0:18])
        for cut in cutList:
            cutString += '{:20.3f}'.format(cut)
        cutString += '\n'
        with open(self.cutFlowFile,'a') as txtfile:
            txtfile.write(cutString)

    def getSampleCutFlow(self,selections,cut,sample,**kwargs):
        '''Return a cut flow histogram with style'''
        sumEntries = kwargs.pop('sumEntries',True)
        printString = kwargs.pop('printString',True)
        hist = ROOT.TH1F('h%sCutFlow' % sample, 'CutFlow', len(selections),0,len(selections))
        valList = []
        errList = []
        for bin in range(len(selections)):
            cutString = cut + ' && ' + selections[bin] if not sumEntries else\
                        cut + ' && ' + ' && '.join(selections[:bin+1])
            val,err = self.getNumEntries(cutString,sample,doError=True,**kwargs)
            #val = self.getNumEntries(selections[:bin+1],sample)
            valList += [val]
            errList += [err]
            hist.SetBinContent(bin+1,val)
            hist.SetBinError(bin+1,err)
        if printString: self.writeCutString(sample,valList)
        hist.SetTitle(self.dataStyles[sample]['name'])
        if sample in self.data: return hist
        hist.SetFillColor(self.dataStyles[sample]['fillcolor'])
        hist.SetLineColor(self.dataStyles[sample]['linecolor'])
        hist.SetFillStyle(self.dataStyles[sample]['fillstyle'])
        return hist

    def getDataCutFlow(self,selections,cut,**kwargs):
        hists = ROOT.TList()
        for sample in self.data:
            hist = self.getSampleCutFlow(selections,cut,sample,**kwargs)
            hists.Add(hist)
        hist = hists[0].Clone("hdataCutFlow")
        hist.Reset()
        hist.Merge(hists)
        return hist

    def getMCStackCutFlow(self,selections,cut,**kwargs):
        '''Return a cut flow stack'''
        mcstack = ROOT.THStack('hsCutFlow','cutFlowStack')
        for sample in self.backgrounds:
            hist = self.getSampleCutFlow(selections,cut,sample,**kwargs)
            if not hist: continue
            mcstack.Add(hist)
        return mcstack

    def getSingleSampleCutFlow_Preselection(self,sample):
        entries = []
        eventsfile = self.ntupleDir+'/%s.num.txt' % sample
        eventsdata = open(eventsfile)
        n_evts = eventsdata.readline()
        eventsdata.close()
        if 'data' not in sample:
            lumi = self.samples[sample]['lumi']
            n_evts = float(n_evts) * self.intLumi/lumi
        entries += [float(n_evts)]
        cutflowfile = self.ntupleDir+'/%s.cutflow.txt' % sample
        with open(cutflowfile) as f:
            cf = f.readlines()
            for val in cf:
                if 'data' not in sample:
                    lumi = self.samples[sample]['lumi']
                    val = float(val) * self.intLumi/lumi
                entries += [float(val)]
        return entries

    def getSampleCutFlow_Preselection(self,sample,**kwargs):
        '''Return preselection cutflow (extracted from sampleName.cutflow.txt and sampleName.num.txt'''
        printString = kwargs.pop('printString',True)
        entries = []
        if sample in self.sampleMergeDict:
            for s in self.sampleMergeDict[sample]:
                newentries = self.getSingleSampleCutFlow_Preselection(s,**kwargs)
                if entries:
                    for e in range(max(len(entries),len(newentries))):
                        if e in range(len(entries)) and e in range(len(newentries)):
                            entries[e] += newentries[e]
                        elif e in range(len(newentries)):
                            entries += [newentries[e]]
                        else:
                            continue
                else:
                    entries = newentries
        else:
            entries = self.getSingleSampleCutFlow_Preselection(sample,**kwargs)
        if printString: self.writeCutString(sample,entries)
        numEntries = len(entries)
        hist = ROOT.TH1F('h%sCutFlow_Preselection' % sample, 'CutFlow', numEntries,0,numEntries)
        for bin in range(numEntries):
            hist.SetBinContent(bin+1,entries[bin])
        hist.SetTitle(self.dataStyles[sample]['name'])
        if sample in self.data: return hist
        hist.SetFillColor(self.dataStyles[sample]['fillcolor'])
        hist.SetLineColor(self.dataStyles[sample]['linecolor'])
        hist.SetFillStyle(self.dataStyles[sample]['fillstyle'])
        return hist

    def getDataCutFlow_Preselection(self,**kwargs):
        hists = ROOT.TList()
        for sample in self.data:
            hist = self.getSampleCutFlow_Preselection(sample,**kwargs)
            hists.Add(hist)
        hist = hists[0].Clone("hdataCutFlow_Preselection")
        hist.Reset()
        hist.Merge(hists)
        return hist

    def getMCStackCutFlow_Preselection(self,**kwargs):
        '''Return a cut flow stack'''
        mcstack = ROOT.THStack('hsCutFlow_Preselection','cutFlowStack')
        for sample in self.backgrounds:
            hist = self.getSampleCutFlow_Preselection(sample,**kwargs)
            if not hist: continue
            mcstack.Add(hist)
        return mcstack

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
        '''A function to plot cut flow of selections. Selections is a list of cuts to apply in order.
           Each subsequent cut is the logical and of the previous.
           kwargs accepts:
               labels      list (string)    bin labels for each cut
               cut         string           applied with all selections
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
        numcol = kwargs.pop('numcol',1)
        signalscale = kwargs.pop('signalscale',1)
        nosum = kwargs.pop('nosum',False)
        isprecf = kwargs.pop('isprecf',False)
        isprelim = kwargs.pop('isprelim', 1)
        yscale = kwargs.pop('yscale',1.2)
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        self.cutFlowFile = self.plotDir+'/'+savename.replace('/','_')+'.txt'
        cutString = '{0: <20}'.format(self.analysis)
        for label in labels:
            cutString += '{0: <20}'.format(label if len(label)<20 else label[0:18])
        cutString += '\n'
        with open(self.cutFlowFile,'w') as txtfile:
            txtfile.write(cutString)

        #print savename, 'Canvas'
        canvas = ROOT.TCanvas(savename,savename,50,50,self.W,self.H)
        canvas = self.setupCanvas(canvas)

        canvas.SetLogy(logy)

        # hack to show both mc and data on same plot
        if plotdata:
            dataHist = self.getDataCutFlow_Preselection(printString=False) if isprecf else self.getDataCutFlow(selections,cut,sumEntries=not nosum,printString=False)
            datamax = dataHist.GetMaximum()
            datamin = dataHist.GetMinimum()
        if plotsig:
            sigHists = {}
            sigmax = 0.
            sigmin = 999999.
            for signal in self.signal:
                sigHists[signal] = self.getSampleCutFlow_Preselection(signal) if isprecf else self.getSampleCutFlow(selections,cut,signal,sumEntries=not nosum)
                sigHists[signal].Scale(signalscale)
                sigmax = max(sigmax,sigHists[signal].GetMaximum())
                sigmin = min(sigmin,sigHists[signal].GetMinimum())

        numSelections = len(selections)
        mc = self.getMCStackCutFlow_Preselection() if isprecf else self.getMCStackCutFlow(selections,cut,sumEntries=not nosum)
        mc.Draw('hist')
        mc.GetYaxis().SetTitle('Events')
        mc.GetYaxis().SetTitleOffset(1)
        newymax = max(datamax,mc.GetMaximum()) if plotdata else mc.GetMaximum()
        newymax = max(sigmax,newymax) if plotsig else newymax
        newymin = mc.GetMinimum()
        if plotdata:
            newymin = min(datamin,newymin) if datamin>0 else newymin
        if plotsig:
            newymin = min(sigmin,newymin) if sigmin>0 else newymin
        mc.SetMaximum(yscale*newymax)
        if isprecf: mc.SetMinimum(1)
        mc.SetMinimum(0.1) if logy else mc.SetMinimum(0)
        if newymin>0 and logy: mc.SetMinimum(newymin)
        if labels:
            for bin in range(numSelections):
                mc.GetHistogram().GetXaxis().SetBinLabel(bin+1,labels[bin])

        staterr = self.get_stat_err(mc.GetStack().Last())
        staterr.Draw("e2 same")

        if plotsig:
            #sigHists = {}
            for signal in self.signal:
                #sigHists[signal] = self.getSampleCutFlow_Preselection(signal) if isprecf else self.getSampleCutFlow(selections,cut,signal,sumEntries=not nosum)
                #sigHists[signal].Scale(signalscale)
                sigHists[signal].SetFillStyle(0)
                sigHists[signal].SetLineWidth(2)
                if signalscale != 1:
                    sigHists[signal].SetTitle(self.dataStyles[signal]['name'] + ' (x%i)' % signalscale)
                sigHists[signal].Draw('hist same')

        if plotdata: 
            #dataHist = self.getDataCutFlow_Preselection() if isprecf else self.getDataCutFlow(selections,cut,sumEntries=not nosum)
            datapois = self.getPoissonErrors(dataHist)
            dataHist.SetMarkerStyle(20)
            dataHist.SetMarkerSize(1.0)
            dataHist.SetLineColor(ROOT.kBlack)
            #dataHist.Draw('esamex0')
            datapois.SetMarkerStyle(20)
            datapois.SetMarkerSize(1.0)
            datapois.SetLineColor(ROOT.kBlack)
            datapois.Draw('0P')

        # draw cms lumi
        self.setStyle(canvas,lumitext,plotdata,isprelim)
        #canvas.cd()
        canvas.Update()
        canvas.RedrawAxis()
        frame = canvas.GetFrame()
        frame.Draw()

        # legend
        if not plotdata: dataHist = 0
        if not plotsig: sigHists = 0
        leg = self.getLegend(mc,dataHist,sigHists,plotdata=plotdata,plotsig=plotsig,plotratio=plotratio,legendpos=legendpos,numcol=numcol)
        leg.Draw()

        # save everything
        #canvas.cd()
        self.save(canvas,savename)

    def plotMCClosure(self, selections, savename, **kwargs):
        '''A function to plot cut flow of selections. Selections is a list of cuts to apply in order.
           Each subsequent cut is the logical and of the previous.
           kwargs accepts:
               labels      list (string)    bin labels for each cut
               cut         string           applied with all selections
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
        numcol = kwargs.pop('numcol',1)
        signalscale = kwargs.pop('signalscale',1)
        nosum = kwargs.pop('nosum',False)
        isprecf = kwargs.pop('isprecf',False)
        isprelim = kwargs.pop('isprelim', 1)
        yscale = kwargs.pop('yscale',1.2)
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        self.cutFlowFile = self.plotDir+'/'+savename.replace('/','_')+'.txt'
        cutString = '{0: <20}'.format(self.analysis)
        for label in labels:
            cutString += '{0: <20}'.format(label if len(label)<20 else label[0:18])
        cutString += '\n'
        with open(self.cutFlowFile,'w') as txtfile:
            txtfile.write(cutString)

        #print savename, 'Canvas'
        canvas = ROOT.TCanvas(savename,savename,50,50,self.W,self.H)
        canvas = self.setupCanvas(canvas)

        canvas.SetLogy(logy)


        numSelections = len(selections)
        mcClosure = self.getSampleCutFlow(selections,cut,'mcClosure',sumEntries=False)
        mcClosureUp = self.getSampleCutFlow(selections,cut,'mcClosure',sumEntries=False,shiftUp=True)
        mcClosureDown = self.getSampleCutFlow(selections,cut,'mcClosure',sumEntries=False,shiftDown=True)
        mc = self.getMCStackCutFlow(selections,cut,sumEntries=False)
        mc.Draw('hist')
        mc.GetYaxis().SetTitle('Events')
        mc.GetYaxis().SetTitleOffset(1)
        newymax = max(mc.GetMaximum(),mcClosure.GetMaximum(),mcClosureUp.GetMaximum(),mcClosureDown.GetMaximum())
        mc.SetMaximum(1.2*newymax)
        mc.SetMinimum(0.1) if logy else mc.SetMinimum(0)
        if labels:
            for bin in range(numSelections):
                mc.GetHistogram().GetXaxis().SetBinLabel(bin+1,labels[bin])

        staterr = self.get_stat_err(mc.GetStack().Last())
        staterr.Draw("e2 same")

        mcClosure.SetFillStyle(0)
        mcClosure.SetLineWidth(2)
        mcClosure.Draw('hist same')

        mcClosureUp.SetFillStyle(0)
        mcClosureUp.SetLineWidth(2)
        mcClosureUp.SetLineStyle(2)
        mcClosureUp.Draw('hist same')

        mcClosureDown.SetFillStyle(0)
        mcClosureDown.SetLineWidth(2)
        mcClosureDown.SetLineStyle(3)
        mcClosureDown.Draw('hist same')

        # draw cms lumi
        self.setStyle(canvas,lumitext,plotdata,isprelim)
        #canvas.cd()
        canvas.Update()
        canvas.RedrawAxis()
        frame = canvas.GetFrame()
        frame.Draw()

        # legend
        dataHist = 0
        sigHists = [mcClosure]
        leg = self.getLegend(mc,dataHist,sigHists,plotdata=False,plotsig=plotsig,plotratio=plotratio,legendpos=legendpos,numcol=numcol)
        leg.Draw()

        # save everything
        #canvas.cd()
        self.save(canvas,savename)
