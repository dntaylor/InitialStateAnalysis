'''
A plotter class.
'''
from PlotterBase import *

class ShapePlotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    def getShape(self,variables,binning,thisCut,plotType,**kwargs):
        signal = kwargs.pop('signal',0)
        plotTypeMap = {
            'MC' : 'getMCStack',
            'Data' : 'getData',
            'Signal' : 'getHist',
        }
        histMethod = getattr(self,plotTypeMap[plotType])
        if plotType == 'MC':
            stack = histMethod(variables,binning,thisCut)
            shape = stack.GetStack().Last()
        elif plotType == 'Data':
            shape = histMethod(variables,binning,thisCut)
        elif plotType == 'Signal' and signal:
            shape = histMethod(signal,variables,binning,thisCut)
        else:
            print 'Error, invalid plot type.'
            return 0
        integral = shape.Integral()
        if integral: shape.Scale(1./integral)
        return shape.Clone()

    # several aliases for fast plotting configuration
    def plotMC(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo'''
        pm = kwargs.pop('plotmc', 1)
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, plotCut, binning, savename, plotmc=pm, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotSignal(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo with signal overlay'''
        pm = kwargs.pop('plotmc', 0)
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, plotCut, binning, savename, plotmc=pm, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotData(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo with data'''
        pm = kwargs.pop('plotmc', 0)
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 0)
        self.plotMCDataSignalRatio(variables, plotCut, binning, savename, plotmc=pm, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotMCRatio(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo with data and a ratio plot'''
        pm = kwargs.pop('plotmc', 1)
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 1)
        self.plotMCDataSignalRatio(variables, plotCut, binning, savename, plotmc=pm, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotSignalRatio(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo with signal and a ratio plot'''
        pm = kwargs.pop('plotmc', 0)
        ps = kwargs.pop('plotsig', 1)
        pd = kwargs.pop('plotdata', 0)
        pr = kwargs.pop('plotratio', 1)
        self.plotMCDataSignalRatio(variables, plotCut, binning, savename, plotmc=pm, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    def plotDataRatio(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo with signal and a ratio plot'''
        pm = kwargs.pop('plotmc', 0)
        ps = kwargs.pop('plotsig', 0)
        pd = kwargs.pop('plotdata', 1)
        pr = kwargs.pop('plotratio', 1)
        self.plotMCDataSignalRatio(variables, plotCut, binning, savename, plotmc=pm, plotsig=ps, plotdata=pd, plotratio=pr, **kwargs)

    # our primary plotting script
    def plotMCDataSignalRatio(self, variables, plotCut, binning, savename, **kwargs):
        '''Plot Monte Carlo with data and signal overlay and a ratio plot
           variables is a list of variables to combine plot (will combine all into same hist)
           binning is either a 3 element list [numBins, binLow, binHigh] or a list of bin edges.
           kwargs accepts:
               cutNames    list (string)    string for labelling histogram for each cut
               cut         string           string for plotting (applied to all histograms)
               xaxis       string           label on xaxis
               yaxis       string           label on yaxis
               xrange      list (double)    range of xaxis (2 elements)
               xmin        double           minimum of xaxis
               xmax        double           maximum of xaxis
               ymin        double           minimum of yaxis
               ymax        double           maximum of yaxis
               blinder     list (double)    range to blind (2 elements)
               logy        bool             set logy plot
               logx        bool             set logx plot
               plotmc      bool             plot mc
               plotsig     bool             plot signal
               plotdata    bool             plot data
               plotratio   bool             make ratio plot
               lumitext    int              location of lumitext (from CMS_lumi)
               legendpos   int              location of legendtext AB (A=012=LCR, B=012=TMB)
               signalscale int              factor to scale signal by
               isprelim    bool             The plot is CMS preliminary'''
        cutNames = kwargs.pop('cutNames',[])
        cut = kwargs.pop('cut', '')
        xaxis = kwargs.pop('xaxis', '')
        yaxis = kwargs.pop('yaxis', '')
        xrange = kwargs.pop('xrange', [])
        xmin = kwargs.pop('xmin', None)
        xmax = kwargs.pop('xmax', None)
        ymin = kwargs.pop('xmin', None)
        ymax = kwargs.pop('xmax', None)
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
            print "Unrecognized parameter '" + key + "' = " + str(value)

        if type(variables) is not list: variables = [variables]

        #if plotratio:
        #    self.canvas.SetCanvasSize(796,666)
        #    plotpad = ROOT.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
        #    plotpad.SetLeftMargin(self.L)
        #    plotpad.SetRightMargin(self.R)
        #    plotpad.SetTopMargin(0.0875)
        #    plotpad.SetBottomMargin(28./666.)
        #    plotpad.SetTickx(0)
        #    plotpad.SetTicky(0)
        #    plotpad.Draw()
        #    self.plotpad = plotpad
        #    ratiopad = ROOT.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
        #    ratiopad.SetTopMargin(0.)
        #    ratiopad.SetBottomMargin(0.5)
        #    ratiopad.SetLeftMargin(self.L)
        #    ratiopad.SetRightMargin(self.R)
        #    ratiopad.SetFillColor(0)
        #    ratiopad.SetTickx(0)
        #    ratiopad.SetTicky(0)
        #    ratiopad.Draw()
        #    plotpad.cd()
        #    plotpad.SetLogy(logy)
        #    plotpad.SetLogx(logx)
        #    ratiopad.SetLogx(logx)
        #else:
        self.canvas.SetLogy(logy)
        self.canvas.SetLogx(logx)

        colors = [ROOT.EColor.kRed, ROOT.EColor.kBlue, ROOT.EColor.kGreen, ROOT.EColor.kBlack]

        legend = ROOT.TLegend(0.65, 0.6, 0.85, 0.88)
        legend.SetFillColor(0)
        legend.SetLineColor(0)
        legend.SetShadowColor(0)
        legend.SetTextFont(62)
        legend.SetTextSize(0.03)
        legend.SetBorderSize(1)
        
        shapes = []
        for c in range(len(plotCut)):
            thisCut = '%s & %s' % (cut, plotCut[c])
            name = cutNames[c] if cutNames else ''
            color = colors[c]
            plotType = 'error'
            lineStyle = 0
            if plotmc: 
                plotType = 'MC'
                lineStyle = 2
            if plotdata: 
                plotType = 'Data'
                lineStyle = 1
            if plotsig:
                plotType = 'Signal'
                lineStyle = 8
                print 'Signal plotting not supported'
                return 0
            shape = self.getShape(variables,binning,thisCut,plotType)
            shape.SetFillColor(color)
            shape.SetLineColor(color)
            shape.SetLineStyle(lineStyle)
            shape.SetFillStyle(0)
            shape.SetTitle(name + ' ' + plotType)
            shape.SetName(shape.GetName()+name+plotType)
            legend.AddEntry(shape)
            shapes += [shape]

        shape = shapes[0]
        shape.Draw()
        shape.GetXaxis().SetTitle(xaxis)
        shape.GetYaxis().SetTitle(yaxis)
        shape.GetYaxis().SetTitleOffset(1)
        if len(xrange)==2:
            shape.GetXaxis().SetRangeUser(xrange[0],xrange[1])
        if ymin: shape.SetMinimum(ymin)
        if ymax: shape.SetMaximum(ymax)
        else: shape.SetMaximum(1.25*shape.GetMaximum())
        #if plotratio:
        #    shape.GetHistogram().GetXaxis().SetLabelOffset(999)
        for shape in shapes[1:]:
            shape.Draw("same")

        legend.Draw('same')

        #if plotratio:
        #    ratiopad.cd()
        #    if len(xrange)==2:
        #        ratiounity = ROOT.TLine(xrange[0],1,xrange[1],1)
        #    else:
        #        ratiounity = ROOT.TLine(stack.GetXaxis().GetXmin(),1,stack.GetXaxis().GetXmax(),1)
        #    ratiounity.SetLineStyle(2)

        #    ratiounity.Draw("same")

        # draw cms lumi
        self.setStyle(lumitext,plotdata,plotratio,isprelim)
        frame = self.canvas.GetFrame()
        frame.Draw()

        # save everything
        self.canvas.cd()
        self.save(savename)

        if plotratio:
            self.resetCanvas()

