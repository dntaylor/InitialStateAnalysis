'''
A plotter class.
'''
from PlotterBase import *

class Plotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    # several aliases for fast plotting configuration
    def plotMC(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo'''
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotMCSignal(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo with signal overlay'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotMCData(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo with data'''
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotMCDataSignal(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo with data and signal overlay'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotMCDataRatio(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo with data and a ratio plot'''
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 1)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotMCSignalRatio(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo with signal and a ratio plot'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 1)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    # our primary plotting script
    def plotMCDataSignalRatio(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo with data and signal overlay and a ratio plot
           variables is a list of variables to combine plot (will combine all into same hist)
           binning is either a 3 element list [numBins, binLow, binHigh] or a list of bin edges.
           kwargs accepts:
               cut         string           string for plotting
               xaxis       string           label on xaxis
               yaxis       string           label on yaxis
               xrange      list (double)    range of xaxis (2 elements)
               xmin        double           minimum of xaxis
               xmax        double           maximum of xaxis
               ymin        double           minimum of yaxis
               ymax        double           maximum of yaxis
               overflow    bool             plot overflow bin
               underflow   bool             plot underflow bin
               nostack     bool             do not stack mc histograms
               normalize   bool             normalize histograms to 1
               blinder     list (double)    range to blind (2 elements)
               logy        bool             set logy plot
               logx        bool             set logx plot
               plotsig     bool             plot signal
               plotdata    bool             plot data
               plotratio   bool             make ratio plot
               lumitext    int              location of lumitext (from CMS_lumi)
               legendpos   int              location of legendtext AB (A=012=LCR, B=012=TMB)
               signalscale int              factor to scale signal by
               isprelim    bool             The plot is CMS preliminary'''
        cut = kwargs.pop('cut', '')
        xaxis = kwargs.pop('xaxis', '')
        yaxis = kwargs.pop('yaxis', '')
        xrange = kwargs.pop('xrange', [])
        xmin = kwargs.pop('xmin', None)
        xmax = kwargs.pop('xmax', None)
        ymin = kwargs.pop('xmin', None)
        ymax = kwargs.pop('xmax', None)
        overflow = kwargs.pop('overflow',False)
        underflow = kwargs.pop('underflow',False)
        nostack = kwargs.pop('nostack',False)
        normalize = kwargs.pop('normalize',False)
        blinder = kwargs.pop('blinder', [])
        logy = kwargs.pop('logy', 1)
        logx = kwargs.pop('logx', 0)
        plotsig = kwargs.pop('plotsig', 1)
        plotdata = kwargs.pop('plotdata', 1)
        plotratio = kwargs.pop('plotratio', 1)
        lumitext = kwargs.pop('lumitext', 11)
        legendpos = kwargs.pop('legendpos', 33)
        signalscale = kwargs.pop('signalscale',1)
        isprelim = kwargs.pop('isprelim', 1)
        if not xmin and len(xrange)==2: xmin = xrange[0]
        if not xmax and len(xrange)==2: xmax = xrange[1]
        if xmin or xmax: xrange = [xmin, xmax]
        for key, value in kwargs.iteritems():
            print "Unrecognized parameter '" + key + "' = " + str(value)

        if type(variables) is not list: variables = [variables]

        if plotratio:
            self.canvas.SetCanvasSize(796,666)
            #self.canvas.SetLeftMargin(0)
            #self.canvas.SetRightMargin(0)
            #self.canvas.SetTopMargin(0)
            #self.canvas.SetBottomMargin(0)
            #self.canvas.SetBorderMode(0)
            plotpad = ROOT.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
            #plotpad.SetFillColor(0)
            #plotpad.SetBorderMode(0)
            #plotpad.SetFrameFillStyle(0)
            #plotpad.SetFrameBorderMode(0)
            plotpad.SetLeftMargin(self.L)
            plotpad.SetRightMargin(self.R)
            plotpad.SetTopMargin(0.0875)
            plotpad.SetBottomMargin(28./666.)
            plotpad.SetTickx(0)
            plotpad.SetTicky(0)
            plotpad.Draw()
            self.plotpad = plotpad
            ratiopad = ROOT.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
            ratiopad.SetTopMargin(0.)
            ratiopad.SetBottomMargin(0.5)
            ratiopad.SetLeftMargin(self.L)
            ratiopad.SetRightMargin(self.R)
            ratiopad.SetFillColor(0)
            #ratiopad.SetBorderMode(0)
            #ratiopad.SetFrameFillStyle(0)
            #ratiopad.SetFrameBorderMode(0)
            ratiopad.SetTickx(0)
            ratiopad.SetTicky(0)
            ratiopad.Draw()
            plotpad.cd()
            plotpad.SetLogy(logy)
            plotpad.SetLogx(logx)
            ratiopad.SetLogx(logx)
        else:
            self.canvas.SetLogy(logy)
            self.canvas.SetLogx(logx)

        # hack to show both mc and data on same plot
        if plotdata:
            data = self.getData(variables, binning, cut, overflow=overflow, underflow=underflow, normalize=normalize)
            datamax = data.GetMaximum()
        

        # plot monte carlo
        stack = self.getMCStack(variables,binning,cut,overflow=overflow,underflow=underflow,nostack=nostack,normalize=normalize)
        stack.SetTitle("")
        stack.Draw("hist nostack") if nostack else stack.Draw("hist")
        stack.GetXaxis().SetTitle(xaxis)
        stack.GetYaxis().SetTitle(yaxis)
        stack.GetYaxis().SetTitleOffset(1)
        if len(xrange)==2:
            stack.GetXaxis().SetRangeUser(xrange[0],xrange[1])
        if ymin: stack.SetMinimum(ymin)
        if ymax: stack.SetMaximum(ymax)
        else:
            newymax = max(datamax,stack.GetMaximum()) if plotdata else stack.GetMaximum()
            if not nostack: stack.SetMaximum(1.25*newymax)
        if plotratio:
            stack.GetHistogram().GetXaxis().SetLabelOffset(999)

        # add errors
        if not nostack:
            staterr = self.get_stat_err(stack.GetStack().Last())
            staterr.Draw("e2 same")

        # plot signal
        if plotsig:
            sigLabels = {}
            for signal in self.signal:
                sighist = self.getHist(signal,variables,binning,cut,overflow=overflow,underflow=underflow,normalzie=normalize)
                sighist.Scale(signalscale)
                sighist.SetFillStyle(0)
                sighist.SetLineWidth(2)
                sighist.Draw('hist same')
                if signalscale != 1:
                    sigLabels[signal] = dataStyles[signal]['name']
                    dataStyles[signal]['name'] += ' (x%i)' % signalscale

        # plot data
        if plotdata:
            data = self.getData(variables, binning, cut, overflow=overflow, underflow=underflow, normalize=normalize)
            data.SetMarkerStyle(20)
            data.SetMarkerSize(1.0)
            data.SetLineColor(ROOT.EColor.kBlack)
            if len(blinder)==2:
                datablind = data.Clone("datablind")
                start = datablind.FindBin(blinder[0])
                end = datablind.FindBin(blinder[1])
                for i in range(start,end+1):
                    datablind.SetBinContent(i,0)
                    datablind.SetBinError(i,0)
                datablind.Draw("esamex0")
            else:
                data.Draw("esamex0")
            self.dataHist = data

        if plotratio:
            ratiopad.cd()
            mchist = stack.GetStack().Last().Clone("mchist%s" % savename)
            if plotdata:
                ratio = self.get_ratio(data,mchist,"ratio%s" % savename)
            if plotsig:
                sig = sighist.Clone("sig%s" % savename)
                sig.Add(mchist)
                ratiosig = self.get_ratio(sig,mchist,"ratiosig%s" % savename)
                ratiosig.SetLineWidth(1)

            ratiostaterr = self.get_ratio_stat_err(mchist)
            ratiostaterr.GetXaxis().SetTitle(xaxis)
            if len(xrange)==2:
                ratiostaterr.GetXaxis().SetRangeUser(xrange[0],xrange[1])

            if len(xrange)==2:
                ratiounity = ROOT.TLine(xrange[0],1,xrange[1],1)
            else:
                ratiounity = ROOT.TLine(stack.GetXaxis().GetXmin(),1,stack.GetXaxis().GetXmax(),1)
            ratiounity.SetLineStyle(2)

            # draw ratios
            ratiopad.cd()
            ratiopad.SetGridy(0)
            ratiostaterr.Draw("e2")
            ratiostaterr.Draw("e2 same")
            ratiounity.Draw("same")
            if plotdata: ratio.Draw("e0 same")
            if plotsig: ratiosig.Draw("hist same")

        # draw cms lumi
        self.setStyle(lumitext,plotdata,plotratio,isprelim)

        # legend
        self.drawLegend(plotdata,plotsig,plotratio,legendpos)

        # save everything
        self.canvas.cd()
        self.save(savename)

        # reset signal names
        if plotsig:
            if signalscale != 1:
                for signal in self.signal:
                    dataStyles[signal]['name'] = sigLabels[signal]

        if plotratio:
            self.resetCanvas()

