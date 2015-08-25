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

def plotOptimization(variable,numTaus,nl):
    masses = _3L_MASSES if nl==3 else _4L_MASSES
    filenames = ['pickles/optimize_%iTau_%s_%i.pkl' %(numTaus,variable,mass) for mass in masses]
    optVals = {}
    for mass,filename in zip(masses,filenames):
        with open(filename,'rb') as f:
            optVals[mass] = pickle.load(f)

    functions = {
        'hmassUnder' : [['x-0.9*x'],['x-x/2'],['x-(x/2-20)']],
        'hmassOver'  : [['1.1*x-x'],['1.1*x-x'],['1.1*x-x']],
        'zmass'      : [['80'],['80'],['50']],
        'dPhi'       : [['x/600+1.95'],['x/200+1.15'],['2.1']],
        'dR'         : [['x/1400.+2.43'],['x/1400.+2.43'],['x/1400.+2.43']],
        'met'        : [['0'],['20'],['40']],
        'st'         : [['1.1*x+60'],['0.85*x+125'],['x-10','200']]
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

        if variable in functions:
            funcs = {}
            for func in functions[variable][numTaus]:
                funcs[func] = ROOT.TF1("7TeV_%s"%(plotType),func,extMasses[0],extMasses[-1])
                funcs[func].SetLineColor(ROOT.kBlack)
                funcs[func].SetLineWidth(3)
                funcs[func].Draw("lsame")


        savedir = 'plots/Hpp3l_Hpp3l_8TeV/png/optimization'
        python_mkdir(savedir)
        savename = '%s/%s_%iTau_%s.png' % (savedir,variable,numTaus,plotType)
        canvas.Print(savename)
