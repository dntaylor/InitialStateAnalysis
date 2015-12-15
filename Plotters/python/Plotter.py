'''
A plotter class.
'''
from PlotterBase import *

class Plotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    # several aliases for fast plotting configuration
    def plotSignal(self, variables, binning, savename, **kwargs):
        '''Plot Monte Carlo'''
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        nb = kwargs.pop('nobg',1)
        self.plotMCDataSignalRatio(variables, binning, savename, plotsig=ps, plotdata=pd, plotratio=pr, nobg=nb, **kwargs)

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
        yaxis = kwargs.pop('yaxis', 'Events')
        xrange = kwargs.pop('xrange', [])
        xmin = kwargs.pop('xmin', None)
        xmax = kwargs.pop('xmax', None)
        ymin = kwargs.pop('ymin', None)
        ymax = kwargs.pop('ymax', None)
        overflow = kwargs.pop('overflow',False)
        underflow = kwargs.pop('underflow',False)
        nostack = kwargs.pop('nostack',False)
        normalize = kwargs.pop('normalize',False)
        blinder = kwargs.pop('blinder', [])
        boxes = kwargs.pop('boxes',[])
        logy = kwargs.pop('logy', 0)
        logx = kwargs.pop('logx', 0)
        nobg = kwargs.pop('nobg',0)
        plotsig = kwargs.pop('plotsig', 1)
        plotdata = kwargs.pop('plotdata', 1)
        plotratio = kwargs.pop('plotratio', 1)
        lumitext = kwargs.pop('lumitext', 11)
        legendpos = kwargs.pop('legendpos', 33)
        numcol = kwargs.pop('numcol',1)
        signalscale = kwargs.pop('signalscale',1)
        isprelim = kwargs.pop('isprelim', 1)
        scalefactor = kwargs.pop('scalefactor','')
        yscale = kwargs.pop('yscale',1.25)
        if not xmin and len(xrange)==2: xmin = xrange[0]
        if not xmax and len(xrange)==2: xmax = xrange[1]
        if xmin or xmax: xrange = [xmin, xmax]
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        if type(variables) is not list: variables = [variables]

        if scalefactor:
            oldscalefactor = self.getScaleFactor()
            self.setScaleFactor(scalefactor)

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        if plotratio:
            self.canvas.SetCanvasSize(796,666)
            plotpad = ROOT.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
            plotpad.SetLeftMargin(self.L)
            plotpad.SetRightMargin(self.R)
            plotpad.SetTopMargin(0.0875)
            plotpad.SetBottomMargin(28./666.)
            plotpad.SetTickx(1)
            plotpad.SetTicky(1)
            plotpad.Draw()
            self.plotpad = plotpad
            ratiopad = ROOT.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
            ratiopad.SetTopMargin(0.)
            ratiopad.SetBottomMargin(0.5)
            ratiopad.SetLeftMargin(self.L)
            ratiopad.SetRightMargin(self.R)
            ratiopad.SetFillColor(0)
            ratiopad.SetTickx(1)
            ratiopad.SetTicky(1)
            ratiopad.Draw()
            plotpad.cd()
            plotpad.SetLogy(logy)
            plotpad.SetLogx(logx)
            ratiopad.SetLogx(logx)
            curPad = plotpad
        else:
            self.canvas.SetLogy(logy)
            self.canvas.SetLogx(logx)
            curPad = self.canvas

        # hack to show both mc and data on same plot
        if plotdata:
            data = self.getData(variables, binning, cut, overflow=overflow, underflow=underflow, normalize=normalize)
            datamax = data.GetMaximum()
        

        # plot monte carlo
        if not nobg:
            stack = self.getMCStack(variables,binning,cut,overflow=overflow,underflow=underflow,nostack=nostack,normalize=normalize)
            stack.SetTitle("")
            stack.Draw("hist nostack") if nostack else stack.Draw("hist")
            stack.GetXaxis().SetTitle(xaxis)
            stack.GetYaxis().SetTitle(yaxis)
            stack.GetYaxis().SetTitleOffset(1)
            if len(xrange)==2:
                stack.GetXaxis().SetRangeUser(xrange[0],xrange[1])
            if ymin: stack.SetMinimum(ymin)
            if ymax: 
                y2 = ymax
                stack.SetMaximum(ymax)
            else:
                newymax = max(datamax,stack.GetMaximum()) if plotdata else stack.GetMaximum()
                y2 = yscale*newymax
                if not nostack: stack.SetMaximum(yscale*newymax)
            if plotratio:
                stack.GetHistogram().GetXaxis().SetLabelOffset(999)

            # add errors
            if not nostack:
                staterr = self.get_stat_err(stack.GetStack().Last())
                staterr.Draw("e2 same")

        # plot signal
        if plotsig:
            sighists = {}
            sigcolors = [
                ROOT.TColor.GetColor('#000000'),
                #ROOT.TColor.GetColor('#1A0000'),
                ROOT.TColor.GetColor('#330000'),
                #ROOT.TColor.GetColor('#4C0000'),
                ROOT.TColor.GetColor('#660000'),
                ROOT.TColor.GetColor('#800000'),
                ROOT.TColor.GetColor('#990000'),
                ROOT.TColor.GetColor('#B20000'),
                ROOT.TColor.GetColor('#CC0000'),
                #ROOT.TColor.GetColor('#E60000'),
                ROOT.TColor.GetColor('#FF0000'),
                #ROOT.TColor.GetColor('#FF1919'),
                ROOT.TColor.GetColor('#FF3333'),
                #ROOT.TColor.GetColor('#FF4D4D'),
                ROOT.TColor.GetColor('#FF6666'),
                ROOT.TColor.GetColor('#FF8080'),
                ROOT.TColor.GetColor('#FF9999'),
                ROOT.TColor.GetColor('#FFB2B2'),
                ROOT.TColor.GetColor('#FFCCCC'),
            ]
            colorCount = 0
            for signal in self.signal:
                sighists[signal] = self.getHist(signal,variables,binning,cut,normalize=normalize,overflow=overflow,underflow=underflow)
                sighists[signal].Scale(signalscale)
                sighists[signal].SetFillStyle(0)
                sighists[signal].SetLineWidth(3)
                if len(self.signal) > 1:
                    sighists[signal].SetLineColor(sigcolors[colorCount])
                    dataStyles[signal]['linecolor'] = sigcolors[colorCount]
                    colorCount += 1
                if nobg and colorCount < 2:
                    sighists[signal].Draw('hist')
                    newymax = max(datamax,sighists[signal].GetMaximum()) if plotdata else sighists[signal].GetMaximum()
                    sighists[signal].SetMaximum(1.25*newymax)
                    y2 = 1.25*newymax
                else:
                    sighists[signal].Draw('hist same')
                if nobg and colorCount < 2:
                    sighists[signal].GetXaxis().SetTitle(xaxis)
                    sighists[signal].GetYaxis().SetTitle(yaxis)
                    sighists[signal].GetYaxis().SetTitleOffset(1)
                if signalscale != 1:
                    sighists[signal].SetTitle(dataStyles[signal]['name']+' (x%i)' % signalscale)

        # plot data
        if plotdata:
            data = self.getData(variables, binning, cut, overflow=overflow, underflow=underflow, normalize=normalize)
            datapois = self.getPoissonErrors(data)
            data.SetMarkerStyle(20)
            data.SetMarkerSize(1.0)
            data.SetLineColor(ROOT.kBlack)
            datapois.SetMarkerStyle(20)
            datapois.SetMarkerSize(1.0)
            datapois.SetLineColor(ROOT.kBlack)
            if len(blinder)==2:
                datablind = data.Clone("datablind")
                start = datablind.FindBin(blinder[0])
                end = datablind.FindBin(blinder[1])
                for i in range(start,end+1):
                    datablind.SetBinContent(i,0)
                    datablind.SetBinError(i,0)
                datablind.Draw("esamex0")
            else:
                #non poisson
                #data.Draw("esamex0")
                # poisson
                datapois.Draw("0P")
            self.dataHist = data

        if boxes:
            box = ROOT.TBox()
            y1 = curPad.GetY1()
            for b in boxes:
                box.SetFillColor(b[2])
                box.SetFillStyle(3002)
                box.SetLineColor(b[2])
                if logy:
                    box.DrawBox(b[0],0,b[1],1.55*467)
                else:
                    box.DrawBox(b[0],0,b[1],467)



        if plotratio:
            ratiopad.cd()
            mchist = stack.GetStack().Last().Clone("mchist%s" % savename)
            if plotdata:
                ratio = self.get_ratio(data,mchist,"ratio%s" % savename)
                ratiopois = self.getPoissonRatio(data,mchist,"ratiopois%s" % savename)
                if len(blinder)==2:
                    ratioblind = ratio.Clone("ratioblind")
                    start = ratioblind.FindBin(blinder[0])
                    end = ratioblind.FindBin(blinder[1])
                    for i in range(start,end+1):
                        ratioblind.SetBinContent(i,999)
                        ratioblind.SetBinError(i,0)
            if plotsig:
                sig = sighists[self.signal[0]].Clone("sig%s" % savename)
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
            if plotdata:
                if len(blinder)==2:
                    ratioblind.Draw("e0 same")
                else:
                    #ratio.Draw("e0 same")
                    ratiopois.Draw("0P same")
            if plotsig: ratiosig.Draw("hist same")

            if boxes:
                for b in boxes:
                    box.SetFillColor(b[2])
                    box.SetFillStyle(3002)
                    box.SetLineColor(b[2])
                    box.DrawBox(b[0],0.5,b[1],1.5)


        if plotratio:
            plotpad.cd()
        else:
            self.canvas.cd()

        # legend
        if not plotdata: data = 0
        if not plotsig: sighists = 0
        if nobg: stack = 0
        leg = self.getLegend(stack,data,sighists,plotdata=plotdata,plotsig=plotsig,plotratio=plotratio,legendpos=legendpos,numcol=numcol)
        leg.Draw()

        # draw cms lumi
        self.setStyle(lumitext,plotdata,plotratio,isprelim)

        # save everything
        self.canvas.cd()
        self.save(savename)

        if plotratio:
            self.resetCanvas()

        if scalefactor:
            self.setScaleFactor(oldscalefactor)

    def plotMCDataSignalRatio2D(self, var1, var2, bin1, bin2, savename, **kwargs):
        cut = kwargs.pop('cut', '')
        zbin = kwargs.pop('zbin',[10,0,10])
        xaxis = kwargs.pop('xaxis', '')
        yaxis = kwargs.pop('yaxis', '')
        xrange = kwargs.pop('xrange', [])
        xmin = kwargs.pop('xmin', None)
        xmax = kwargs.pop('xmax', None)
        ymin = kwargs.pop('ymin', None)
        ymax = kwargs.pop('ymax', None)
        overflow = kwargs.pop('overflow',False)
        underflow = kwargs.pop('underflow',False)
        nostack = kwargs.pop('nostack',False)
        normalize = kwargs.pop('normalize',False)
        blinder = kwargs.pop('blinder', [])
        logy = kwargs.pop('logy', 1)
        logx = kwargs.pop('logx', 0)
        plotmc = kwargs.pop('plotmc', 1)
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
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        if type(var1) is not list: var1 = [var1]
        if type(var2) is not list: var2 = [var2]

        ROOT.gDirectory.Delete('h*') # clear histogram memory

        self.canvas.SetRightMargin(0.14)

        # plot monte carlo
        if plotmc:
            stack = self.getMCStack2D(var1,var2,bin1,bin2,cut,zbin=zbin)
            stack.SetTitle("")
            stack.Draw("colz")
            stack.GetXaxis().SetTitle(xaxis)
            stack.GetYaxis().SetTitle(yaxis)
            stack.GetYaxis().SetTitleOffset(1)

        # plot data
        if plotdata:
            data = self.getData2D(var1, var2, bin1, bin2, cut,zbin=zbin)
            data.SetTitle("")
            data.Draw("colz")
            data.GetXaxis().SetTitle(xaxis)
            data.GetYaxis().SetTitle(yaxis)
            data.GetYaxis().SetTitleOffset(1)

        if plotsig:
            for signal in self.signal:
                sig = self.getHist2D(signal, var1, var2, bin1, bin2, cut, zbin=zbin)
                sig.SetTitle("")
                sig.Draw("colz")
                sig.GetXaxis().SetTitle(xaxis)
                sig.GetYaxis().SetTitle(yaxis)
                sig.GetYaxis().SetTitleOffset(1)

        # save everything
        self.canvas.cd()
        self.save(savename)

