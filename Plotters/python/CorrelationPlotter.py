'''
A correlation plotter class.

Author: Devin N. Taylor, UW-Madison
'''

from PlotterBase import PlotterBase
import ROOT

class CorrelationPlotter(PlotterBase):
    def __init__(self,analysis,**kwargs):
        PlotterBase.__init__(self,analysis,**kwargs)

    def getCorrelation(self, samples, cuts):
        '''Get correlation of cuts'''
        hist = ROOT.TH2F('correlation','correlation',len(cuts),0,len(cuts),len(cuts),0,len(cuts))
        for x,xcut in enumerate(cuts):
            for y,ycut in enumerate(cuts):
                num = 0.
                numErr2 = 0.
                denom = 0.
                denomErr2 = 0.
                for sample in samples:
                    val, err = self.getNumEntries('%s && %s' % (xcut,ycut), sample, doError=True)
                    num += val
                    numErr2 += err**2
                    val, err = self.getNumEntries(xcut, sample, doError=True)
                    denom += val
                    denomErr2 += err**2
                val = num / denom
                err = val * (numErr2/num**2 + denomErr2/denom**2)**0.5
                hist.SetBinContent(x+1,y+1,val)
                hist.SetBinError(x+1,y+1,err)
        return hist

    # aliases for plotting
    def plotMCCorrelation(self,selections, savename, **kwargs):
        plottype = kwargs.pop('plottype','mc')
        self.plotCorrelation(selections, savename, plottype=plottype, **kwargs)

    def plotDataCorrelation(self,selections, savename, **kwargs):
        plottype = kwargs.pop('plottype','data')
        self.plotCorrelation(selections, savename, plottype=plottype, **kwargs)

    def plotSignalCorrelation(self,selections, savename, **kwargs):
        plottype = kwargs.pop('plottype','sig')
        self.plotCorrelation(selections, savename, plottype=plottype, **kwargs)

    def plotCorrelation(self, selections, savenames, **kwargs):
        '''Plot the correlation between a series of cuts. Correlation is defined as:
               x bin = base selection
               y bin = correlated selection
               z value = events passing y selection and x selection / events passing x selection
           kwargs accepts:
               labels      list (string)    bin labels for each cut
               cut         string           applied with all selections
               plottype    [mc,data,sig]    plot type
               lumitext    int              location of lumitext (from CMS_lumi)
               legendpos   int              location of legendtext AB (A=012=LCR, B=012=TMB)
               signalscale int              factor to scale signal by
               isprelim    bool             The plot is CMS preliminary'''
        labels = kwargs.pop('labels',[])
        cut = kwargs.pop('cut', '1')
        plottype = kwargs.pop('plottype', 'mc')
        sample = kwargs.pop('sample','')
        samples = kwargs.pop('samples',[])
        lumitext = kwargs.pop('lumitext', 11)
        legendpos = kwargs.pop('legendpos', 33)
        signalscale = kwargs.pop('signalscale',1)
        isprelim = kwargs.pop('isprelim', 1)
        for key, value in kwargs.iteritems():
            self.logger.warning("Unrecognized parameter '" + key + "' = " + str(value))

        if plottype not in ['mc','data','sig']:
            self.logger.error('Allowed plot types are: \'mc\', \'data\', \'sig\'.')

        if type(savenames) is not list: savenames = [savenames]

        self.canvas.SetRightMargin(0.14)

        cuts = ['%s && %s' % (cut, x) for x in selections]
        hists = []
        if sample or samples: # a specific sample was passed
            if type(samples) is not list: samples = [samples]
            if sample: samples += [sample]
            for sample in samples:
                hists += [self.getCorrelation(sample, cuts)]
        else:
            if plottype=='mc':
                hists += [self.getCorrelation(self.backgrounds, cuts)]
            elif plottype=='data':
                hists += [self.getCorrelation(self.data, cuts)]
            elif plottype=='sig':
                for signal in self.signal:
                    hists += [self.getCorrelation([signal], cuts)]
            else:
                pass

        if len(hists) != len(savenames):
            self.logger.error('Number of outfiles does not match number of histograms')
            return

        for hist,savename in zip(hists,savenames):
            # set bin labels
            if labels:
                for bin in range(len(cuts)):
                    hist.GetXaxis().SetBinLabel(bin+1,labels[bin])
                    hist.GetYaxis().SetBinLabel(bin+1,labels[bin])
            hist.Draw('colz')
            # set styles
            # save everything
            self.canvas.cd()
            self.save(savename)

