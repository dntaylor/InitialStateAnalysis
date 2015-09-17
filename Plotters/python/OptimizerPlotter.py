import ROOT
import pickle
from array import array

import CMS_lumi, tdrstyle
from plotUtils import *
from plotUtils import _3L_MASSES, _4L_MASSES

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPalette(1)

def plotOptimization(analysis,period,variable,numTaus):
    nl = 3 if analysis in ['Hpp3l'] else 4
    masses = _3L_MASSES if nl==3 else _4L_MASSES
    filenames = ['pickles/%s_%iTeV_%s/optimize_%iTau_%s_%i.pkl' %(analysis,period,analysis,numTaus,variable,mass) for mass in masses]
    optVals = {}
    for mass,filename in zip(masses,filenames):
        with open(filename,'rb') as f:
            optVals[mass] = pickle.load(f)

    functions = {
        'Hpp3l': {
            'hmassUnder' : [['x-0.9*x'],['x-x/2'],['x-x/2']],
            'hmassOver'  : [['1.1*x-x'],['1.1*x-x'],['1.1*x-x']],
            'zmass'      : [['80'],['80'],['50']],
            'dPhi'       : [['x/600+1.95'],['x/200+1.15'],['2.1']],
            'dR'         : [['x/380+2.06','x/1200+2.77'],['x/380+1.96','x/1000+2.6'],['x/380+1.86','x/1000+2.37']],
            'met'        : [['0'],['20'],['40']],
            'st'         : [['1.07*x+45'],['0.72*x+50'],['0.44*x+65']],
        },
        'Hpp4l': {
            'hmassUnder' : [['x-0.9*x'],['x-x/2'],['x-x/2']],
            'hmassOver'  : [['1.1*x-x'],['1.1*x-x'],['1.1*x-x']],
            'zmass'      : [['0'],['0'],['0']],
            'dPhi'       : [['0'],['0'],['0']],
            'dR'         : [['0'],['0'],['0']],
            'met'        : [['0'],['0'],['0']],
            'st'         : [['0.6*x+130.'],['0'],['0']],
        }
    }

    canvas = ROOT.TCanvas('c','c',50,50,800,600)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.14)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)

    for plotType in ['sigEff','bgEff','significance1','significance2','significance3']:
        # draw the values
        numMasses = len(masses)
        extMasses = [150] + masses + [800]
        mbinning = [float(extMasses[x]+extMasses[x+1])/2. for x in range(len(masses)+1)]
        massBins = array('d',mbinning)
        numBins = len(optVals[masses[0]][plotType])
        binning = optVals[masses[0]]['cuts']
        extBinning = [binning[0] - (binning[1]-binning[0])] + binning
        bins = array('d',extBinning)
        hist = ROOT.TH2F(plotType,plotType,numMasses,massBins,numBins,bins)
        for m,mass in enumerate(masses):
            for b,val in enumerate(optVals[mass][plotType]):
                hist.SetBinContent(m+1,b+1,val[0])
                hist.SetBinError(m+1,b+1,val[1])
        hist.GetXaxis().SetTitle('M(\\ell^{\\pm}\\ell^{\\pm}) (GeV/c^2)')
        hist.Draw('colz goff')
        # draw a fit line
        if 'significance' in plotType:
            maxVals = {}
            maxIndex = {}
            downTen = {}
            upTen = {}
            for mass in masses:
                valList = [x[0] for x in optVals[mass][plotType]]
                maxVals[mass] = max(valList)
                maxIndex[mass] = valList.index(maxVals[mass])
                downTen[mass] = 0
                canContinue = False
                for v,val in enumerate(valList):
                    if val>0.9*maxVals[mass] and not downTen[mass] and v<maxIndex[mass]: # find lowest within 10%
                        downTen[mass] = v
                    if v>maxIndex[mass] and val<0.9*maxVals[mass]: # find highest within 10%
                        upTen[mass] = v
                        break
                if mass not in upTen: upTen[mass] = len(valList)-1
            # create TGraph for points
            maxVal = ROOT.TGraph(numMasses)
            down = ROOT.TGraph(numMasses)
            up = ROOT.TGraph(numMasses)
            for m,mass in enumerate(masses):
                maxVal.SetPoint(m,mass,binning[maxIndex[mass]])
                down.SetPoint(m,mass,binning[downTen[mass]])
                up.SetPoint(m,mass,binning[upTen[mass]])
            maxVal.SetMarkerStyle(0)
            maxVal.SetFillStyle(0)
            down.SetMarkerStyle(0)
            down.SetFillStyle(0)
            down.SetLineStyle(7)
            up.SetMarkerStyle(0)
            up.SetFillStyle(0)
            up.SetLineStyle(7)
            maxVal.Draw('same')
            down.Draw('same')
            up.Draw('same')

        if variable in functions[analysis]:
            funcs = {}
            for func in functions[analysis][variable][numTaus]:
                funcs[func] = ROOT.TF1("7TeV_%s"%(plotType),func,extMasses[0],extMasses[-1])
                funcs[func].SetLineColor(ROOT.kBlack)
                funcs[func].SetLineWidth(3)
                funcs[func].Draw("lsame")


        savedir = 'plots/%s_%s_%iTeV/png/optimization' % (analysis,analysis,period)
        python_mkdir(savedir)
        savename = '%s/%s_%iTau_%s.png' % (savedir,variable,numTaus,plotType)
        canvas.Print(savename)
