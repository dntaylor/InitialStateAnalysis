'''
A plotter class.
'''
from PlotterBase import *

class MCClosurePlotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    # our primary plotting script
    def plotRatio(self, variables, binning, savename, **kwargs):
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
        ratiomin = kwargs.pop('ratiomin',0.5)
        ratiomax = kwargs.pop('ratiomax',1.5)
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

        #print savename, 'Delete'
        ROOT.gDirectory.Delete('h*') # clear histogram memory

        #print savename, 'Canvas'
        canvas = ROOT.TCanvas(savename,savename,50,50,self.W,self.H)
        canvas = self.setupCanvas(canvas)

        if plotratio:
            #print savename, 'Plotpad'
            canvas.SetCanvasSize(796,666)
            plotpad = ROOT.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
            plotpad.SetLeftMargin(self.L)
            plotpad.SetRightMargin(self.R)
            plotpad.SetTopMargin(0.0875)
            plotpad.SetBottomMargin(0.04)
            plotpad.SetTickx(1)
            plotpad.SetTicky(1)
            plotpad.Draw()
            plotpad = plotpad
            #print savename, 'Ratiopad'
            ratiopad = ROOT.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
            ratiopad.SetTopMargin(0.06)
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
            canvas.SetLogy(logy)
            canvas.SetLogx(logx)
            curPad = canvas

        # plot monte carlo
        #print savename, 'getMCStack'
        mcClosure = self.getHist('mcClosure',variables,binning,cut,normalize=normalize,overflow=overflow,underflow=underflow)
        mcClosure.SetFillStyle(0)
        mcClosure.SetLineWidth(3)
        highest = mcClosure.GetMaximum()

        stack = self.getMCStack(variables,binning,cut,overflow=overflow,underflow=underflow,nostack=nostack,normalize=normalize,plotsig=not plotsig)
        stack.SetTitle("")
        stack.Draw("hist nostack") if nostack else stack.Draw("hist")
        stack.GetXaxis().SetTitle(xaxis)
        stack.GetYaxis().SetTitle(yaxis)
        stack.GetYaxis().SetTitleOffset(1)
        highest = max(highest, stack.GetMaximum())
        stack.SetMaximum(yscale*highest)
        if plotratio:
            stack.GetHistogram().GetXaxis().SetLabelOffset(999)

        # add errors
        staterr = self.get_stat_err(stack.GetStack().Last())
        staterr.Draw("e2 same")

        mcClosure.Draw('hist same')

        data = 0
        sighists = [mcClosure]
        #print savename, 'getLegend'
        leg = self.getLegend(stack,data,sighists,plotdata=False,plotsig=True,plotratio=plotratio,legendpos=legendpos,numcol=numcol)
        leg.Draw()

        # draw cms lumi
        pad = plotpad if plotratio else canvas
        #print savename, 'setStyle'
        self.setStyle(pad,lumitext,plotdata,isprelim)


        if plotratio:
            #print savename, 'plotratio'
            ratiopad.cd()
            #if plotsig:
            #    mchist = stack.GetStack().Last().Clone("mchist%s" % savename)
            #else:
            #    newstack = self.getMCStack(variables,binning,cut,overflow=overflow,underflow=underflow,nostack=nostack,normalize=normalize,plotsig=False,histname='newstack')
            mchist = stack.GetStack().Last().Clone("mchist%s" % savename)
            ratioClosure = self.get_ratio(mcClosure,mchist,"ratiosig%s" % savename)
            ratioClosure.SetLineWidth(1)

            ratiostaterr = self.get_ratio_stat_err(mchist,ratiomin=ratiomin,ratiomax=ratiomax)
            ratiostaterr.SetXTitle(xaxis)
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
            ratiounity.Draw("same")
            ratioClosure.Draw("hist same")


        if plotratio:
            plotpad.cd()
        else:
            canvas.cd()


        #print savename, 'save'
        self.save(canvas,savename)


        if scalefactor:
            self.setScaleFactor(oldscalefactor)

