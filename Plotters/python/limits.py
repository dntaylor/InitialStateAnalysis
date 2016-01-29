'''
Class for plotting limits.
'''

import sys
import os
import errno
import numpy as np
import CMS_lumi, tdrstyle
from plotUtils import _3L_MASSES, _4L_MASSES, python_mkdir
from xsec import xsecs

sys.argv.append('-b')
import ROOT
sys.argv.pop()

ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()

# bp labels
bpLabels = {
    'ee100' : '100% #Phi^{#pm#pm} #rightarrow ee',
    'em100' : '100% #Phi^{#pm#pm} #rightarrow e#mu',
    'et100' : '100% #Phi^{#pm#pm} #rightarrow e#tau',
    'mm100' : '100% #Phi^{#pm#pm} #rightarrow #mu#mu',
    'mt100' : '100% #Phi^{#pm#pm} #rightarrow #mu#tau',
    'tt100' : '100% #Phi^{#pm#pm} #rightarrow #tau#tau',
    'BP1' : 'Benchmark 1',
    'BP2' : 'Benchmark 2',
    'BP3' : 'Benchmark 3',
    'BP4' : 'Benchmark 4',
}

def save(savename,saveDir,canvas):
    '''Save the limits into root file and images.'''
    #for type in ['png', 'pdf', 'eps']:
    for type in ['png','pdf','root']:
        name = "%s/%s/%s.%s" % (saveDir, type, savename, type)
        python_mkdir(os.path.dirname(name))
        canvas.Print(name)
    #canvas.SetName(savename)
    #savefile.WriteTObject(self.canvas)
    #canvas.Clear()

def plot_limits(analysis, region, period, savename, **kwargs):
    '''Plot limits and get exclusion limits'''
    datacardBaseDir = kwargs.pop('datacardBaseDir','datacards')
    limitDataBaseDir = kwargs.pop('limitDataBaseDir','limitData')
    saveDir = kwargs.pop('saveDir','plots/limits')
    blind = not kwargs.pop('unblind',False)
    bp = kwargs.pop('branchingPoint','')
    bgMode = kwargs.pop('bgMode','sideband')
    do4l = kwargs.pop('do4l',False)
    limitMode = kwargs.pop('limitMode','fullCLs')

    saveDir += '/' + limitMode

    # directories    
    datacardDir = '%s/%s_%itev_%s' % (datacardBaseDir, analysis, period, region)
    if bp: datacardDir += '/%s' % bp
    limitDataDir = '%s/%s_%itev_%s' % (limitDataBaseDir, analysis, period, region)
    if bp: limitDataDir += '/%s' % bp
    datacardString = '' if bgMode == "sideband" else "_{0}".format(bgMode)
    if do4l: datacardString += '_4l'

    if analysis == 'HppComb': analysisName = 'Combined'
    if analysis == 'HppAP': analysisName = 'AP'
    if analysis == 'HppPP': analysisName = 'PP'
    if analysis == 'Hpp3l': analysisName = '3lAP'
    if analysis == 'Hpp3l' and do4l: analysisName = '3lPP'
    if analysis == 'Hpp4l': analysisName = '4lPP'

    # masses to include
    masses = _3L_MASSES if analysis in ['Hpp3l','HppAP'] else _4L_MASSES
    if do4l: masses = _4L_MASSES
    if period==13: masses = [500]
    n = len(masses)

    # get cross sections
    xsecMap = {}
    xsecGraph = ROOT.TGraph(n)
    for i,mass in enumerate(masses):
        sample = 'HPlusPlusHMinusMinusHTo4L_M-{0}_8TeV-pythia6' if analysis in ['Hpp4l','HppComb','HppPP'] or do4l else 'HPlusPlusHMinusHTo3L_M-{0}_8TeV-calchep-pythia6'
        xsecMap[mass] = xsecs[period][sample.format(mass)]
        if analysis in ['HppComb'] and mass >= 170:
            xsecMap[mass] += xsecs[period]['HPlusPlusHMinusHTo3L_M-{0}_8TeV-calchep-pythia6'.format(mass)]
        xsecGraph.SetPoint(i,mass,xsecMap[mass])
    xsecGraph.SetMarkerStyle(0)
    xsecGraph.SetFillStyle(0)
    xsecGraph.SetLineColor(ROOT.kBlue)

    # get limit values
    quartiles = np.empty((6, len(masses)), dtype=float)
    for j, mass in enumerate(masses):
        #fname = os.path.join(limitDataDir, "higgsCombineTest.Asymptotic.mH%i%s.root" % (mass,datacardString))
        #file = ROOT.TFile(fname,"READ")
        #tree = file.Get("limit")
        #if not tree: continue
        #for i, row in enumerate(tree):
        #    quartiles[i,j] = row.limit
        limitFilename = '{0}/{1}/{2}/{3}/limits.txt'.format(limitMode,analysisName,bp,mass)
        with open(limitFilename,'r') as f:
            limString = f.readline()
            limList = limString.split(' ')
            for i in range(6):
                quartiles[i,j] = limList[i]

    twoSigma = ROOT.TGraph(2*n)
    oneSigma = ROOT.TGraph(2*n)
    expected  = ROOT.TGraph(n)
    oneSigma_low = ROOT.TGraph(n)
    oneSigma_high = ROOT.TGraph(n)
    if not blind: observed  = ROOT.TGraph(n)
    xtwoSigma = ROOT.TGraph(2*n)
    xoneSigma = ROOT.TGraph(2*n)
    xexpected  = ROOT.TGraph(n)
    xoneSigma_low = ROOT.TGraph(n)
    xoneSigma_high = ROOT.TGraph(n)
    if not blind: xobserved  = ROOT.TGraph(n)
    for i, mass in enumerate(masses):
        twoSigma.SetPoint(i,masses[i],quartiles[4][i])
        twoSigma.SetPoint(n+i,masses[n-i-1],quartiles[0][n-i-1])
        oneSigma.SetPoint(i,masses[i],quartiles[3][i])
        oneSigma.SetPoint(n+i,masses[n-i-1],quartiles[1][n-i-1])
        oneSigma_low.SetPoint(i,masses[i],quartiles[3][i])
        oneSigma_high.SetPoint(i,masses[i],quartiles[1][i])
        expected.SetPoint(i,masses[i],quartiles[2][i])
        if not blind: observed.SetPoint(i,masses[i],quartiles[5][i])
        xtwoSigma.SetPoint(i,masses[i],quartiles[4][i]*xsecMap[masses[i]])
        xtwoSigma.SetPoint(n+i,masses[n-i-1],quartiles[0][n-i-1]*xsecMap[masses[n-i-1]])
        xoneSigma.SetPoint(i,masses[i],quartiles[3][i]*xsecMap[masses[i]])
        xoneSigma.SetPoint(n+i,masses[n-i-1],quartiles[1][n-i-1]*xsecMap[masses[n-i-1]])
        xoneSigma_low.SetPoint(i,masses[i],quartiles[3][i]*xsecMap[masses[i]])
        xoneSigma_high.SetPoint(i,masses[i],quartiles[1][i]*xsecMap[masses[i]])
        xexpected.SetPoint(i,masses[i],quartiles[2][i]*xsecMap[masses[i]])
        if not blind: xobserved.SetPoint(i,masses[i],quartiles[5][i]*xsecMap[masses[i]])
    twoSigma.SetFillColor(ROOT.kYellow)
    twoSigma.SetLineColor(ROOT.kYellow)
    twoSigma.SetMarkerStyle(0)
    oneSigma.SetFillColor(ROOT.kSpring)
    oneSigma.SetLineColor(ROOT.kSpring)
    oneSigma.SetMarkerStyle(0)
    expected.SetLineStyle(7)
    expected.SetMarkerStyle(0)
    expected.SetFillStyle(0)
    if not blind:
        observed.SetMarkerStyle(0)
        observed.SetFillStyle(0)
    xtwoSigma.SetFillColor(ROOT.kYellow)
    xtwoSigma.SetLineColor(ROOT.kYellow)
    xtwoSigma.SetMarkerStyle(0)
    xoneSigma.SetFillColor(ROOT.kSpring)
    xoneSigma.SetLineColor(ROOT.kSpring)
    xoneSigma.SetMarkerStyle(0)
    xexpected.SetLineStyle(7)
    xexpected.SetMarkerStyle(0)
    xexpected.SetFillStyle(0)
    if not blind:
        xobserved.SetMarkerStyle(0)
        xobserved.SetFillStyle(0)
    
    # create canvas
    canvas = ROOT.TCanvas('c%s'%bp,'c%s'%bp,50,50,800,600)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.SetLogy(1)

    # create CLs plot
    expected.GetXaxis().SetLimits(170,masses[-1])
    expected.GetXaxis().SetTitle('#Phi^{++} Mass (GeV)')
    expected.GetYaxis().SetTitle('95% CLs Upper Limit on #sigma/#sigma_{model}')
    expected.GetYaxis().SetTitleOffset(1.)
    expected.GetYaxis().SetTitleSize(0.05)

    expected.Draw()
    twoSigma.Draw('f')
    oneSigma.Draw('f')
    expected.Draw('same')
    ROOT.gPad.RedrawAxis()
    if not blind: observed.Draw('same')

    ratiounity = ROOT.TLine(expected.GetXaxis().GetXmin(),1,expected.GetXaxis().GetXmax(),1)
    ratiounity.Draw()

    legend = ROOT.TLegend(0.60,0.15,0.90,0.45)
    legend.SetFillColor(0)
    if not blind: legend.AddEntry(observed, 'Observed','l')
    legend.AddEntry(expected, 'Expected','l')
    legend.AddEntry(twoSigma, 'Expected 2#sigma', 'F')
    legend.AddEntry(oneSigma, 'Expected 1#sigma', 'F')
    legend.AddEntry(None,bpLabels[bp],'')

    legend.Draw('same')

    lumiperiod = 2 if period == 8 else 4
    CMS_lumi.wrtieExtraText = True
    CMS_lumi.extraText = "Preliminary" if not blind else "Simulation Preliminary"
    CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (4.9)
    CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (19.7)
    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (25.0)
    CMS_lumi.CMS_lumi(canvas,lumiperiod,11)

    # plot label
    #label = ROOT.TLatex(0.15,0.7,bpLabels[bp])
    #label.SetNDC(True)
    #label.SetTextFont(52)
    #label.Draw()

    if do4l: savename += '_4l'
    save(savename,saveDir,canvas)

    canvas.Clear()

    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.SetLogy(1)

    # create the limit on cross section plot
    xsecGraph.GetXaxis().SetLimits(masses[0],masses[-1])
    xsecGraph.GetXaxis().SetTitle('#Phi^{++} Mass (GeV)')
    xsecGraph.GetYaxis().SetTitle('#sigma*BR (pb)')
    if analysis in ['Hpp4l'] or do4l: xsecGraph.GetYaxis().SetTitle('#sigma*BR^{2} (pb)')
    xsecGraph.GetYaxis().SetTitleOffset(1.)
    xsecGraph.GetYaxis().SetTitleSize(0.05)

    xsecGraph.Draw('')
    xsecGraph.SetMaximum(1)
    xtwoSigma.Draw('f')
    xoneSigma.Draw('f')
    xexpected.Draw('same')
    xsecGraph.Draw('same')
    ROOT.gPad.RedrawAxis()
    if not blind: xobserved.Draw('same')

    legend = ROOT.TLegend(0.5,0.6,0.90,0.85)
    legend.SetFillColor(0)
    if not blind: legend.AddEntry(xobserved, 'Observed', 'l')
    legend.AddEntry(xexpected, 'Expected', 'l')
    legend.AddEntry(xtwoSigma, 'Expected 2#sigma', 'F')
    legend.AddEntry(xoneSigma, 'Expected 1#sigma', 'F')
    name = 'Pair Production Cross Section' if analysis in ['Hpp4l','HppPP'] or do4l else 'Associated Production Cross Section'
    if analysis in ['HppComb']: name = 'Cross Section'
    legend.AddEntry(xsecGraph, name, 'l')
    legend.AddEntry(None,bpLabels[bp],'')

    legend.Draw('same')

    lumiperiod = 2 if period == 8 else 4
    CMS_lumi.wrtieExtraText = True
    CMS_lumi.extraText = "Preliminary" if not blind else "Simulation Preliminary"
    CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (4.9)
    CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (19.7)
    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (25.0)
    CMS_lumi.CMS_lumi(canvas,lumiperiod,11)

    # plot label
    #label = ROOT.TLatex(0.15,0.7,bpLabels[bp])
    #label.SetNDC(True)
    #label.SetTextFont(52)
    #label.Draw()

    savename += '_crossSection'
    save(savename,saveDir,canvas)

    # get expected limit
    y = 0
    for x in range(1,1000):
        y = expected.Eval(x)
        if y > 1: break
    y = 0
    for l in range(1,1000):
        y = oneSigma_low.Eval(l)
        if y > 1: break
    y = 0
    for h in range(1,1000):
        y = oneSigma_high.Eval(h)
        if y > 1: break

    # get observed limit
    if not blind:
        dy = 0
        for dx in range(1,1000):
            dy = observed.Eval(dx)
            if dy > 1: break
    else:
        dx = 0


    print "Expected Limit: %i GeV (+%i, -%i)" % (x, h-x, x-l)
    if not blind: print "Observed Limit: %i GeV" % (dx)
    return [x,h-x,x-l,dx]

def plot_combined_limits(period, savename, **kwargs):
    '''Plot limits and get exclusion limits'''
    datacardBaseDir = kwargs.pop('datacardBaseDir','datacards')
    limitDataBaseDir = kwargs.pop('limitDataBaseDir','limitData')
    saveDir = kwargs.pop('saveDir','plots/limits')
    blind = not kwargs.pop('unblind',False)
    bp = kwargs.pop('branchingPoint','')
    bgMode = kwargs.pop('bgMode','sideband')
    limitMode = kwargs.pop('limitMode','fullCLs')

    saveDir += '/' + limitMode

    # directories
    datacardDir = {}
    limitDataDir = {}
    for analysis in ['AP','PP','Comb']:
        datacardDir[analysis] = '%s/Hpp%s_%itev_Hpp%s/%s' % (datacardBaseDir, analysis, period, analysis, bp)
        limitDataDir[analysis] = '%s/Hpp%s_%itev_Hpp%s/%s' % (limitDataBaseDir, analysis, period, analysis, bp)

    datacardString = '' if bgMode == "sideband" else "_{0}".format(bgMode)

    # masses to include
    masses = _4L_MASSES
    n = len(masses)

    # get cross sections
    xsecMap = {'AP':{},'PP':{}}
    xsecGraph = {'AP':ROOT.TGraph(n),'PP':ROOT.TGraph(n)}
    for i,mass in enumerate(masses):
        sample_4l = 'HPlusPlusHMinusMinusHTo4L_M-{0}_8TeV-pythia6'
        sample_3l = 'HPlusPlusHMinusHTo3L_M-{0}_8TeV-calchep-pythia6'
        xsecMap['PP'][mass] = xsecs[period][sample_4l.format(mass)]
        xsecMap['AP'][mass] = xsecs[period][sample_3l.format(mass)] if mass in _3L_MASSES else 0.
        xsecGraph['AP'].SetPoint(i,mass,xsecMap['AP'][mass])
        xsecGraph['PP'].SetPoint(i,mass,xsecMap['PP'][mass])
    for analysis in ['AP','PP']:
        xsecGraph[analysis].SetMarkerStyle(0)
        xsecGraph[analysis].SetFillStyle(0)
        xsecGraph[analysis].SetLineColor(ROOT.kBlue)

    # get limit values
    quartiles = {}
    for analysis in ['AP','PP','Comb']:
        if analysis == 'Comb': analysisName = 'Combined'
        if analysis == 'AP': analysisName = 'AP'
        if analysis == 'PP': analysisName = 'PP'
        quartiles[analysis] = np.empty((6, len(masses)), dtype=float)
        for j, mass in enumerate(masses):
            if mass not in _3L_MASSES and analysis == 'AP': continue
            #fname = os.path.join(limitDataDir[analysis], "higgsCombineTest.Asymptotic.mH%i%s.root" % (mass,datacardString))
            #file = ROOT.TFile(fname,"READ")
            #tree = file.Get("limit")
            #if not tree: continue
            #for i, row in enumerate(tree):
            #    quartiles[analysis][i,j] = row.limit
            limitFilename = '{0}/{1}/{2}/{3}/limits.txt'.format(limitMode,analysisName,bp,mass)
            with open(limitFilename,'r') as f:
                limString = f.readline()
                limList = limString.split(' ')
                for i in range(6):
                    quartiles[analysis][i,j] = limList[i]


    twoSigma = {}
    oneSigma = {}
    expected = {}
    oneSigma_low = {}
    oneSigma_high = {}
    observed = {}
    xtwoSigma = {}
    xoneSigma = {}
    xexpected = {}
    xoneSigma_low = {}
    xoneSigma_high = {}
    xobserved = {}
    for analysis in ['AP','PP','Comb']:
        twoSigma[analysis] = ROOT.TGraph(2*n)
        oneSigma[analysis] = ROOT.TGraph(2*n)
        expected[analysis]  = ROOT.TGraph(n)
        oneSigma_low[analysis] = ROOT.TGraph(n)
        oneSigma_high[analysis] = ROOT.TGraph(n)
        if not blind: observed[analysis]  = ROOT.TGraph(n)
        if analysis in ['AP','PP']:
            xtwoSigma[analysis] = ROOT.TGraph(2*n)
            xoneSigma[analysis] = ROOT.TGraph(2*n)
            xexpected[analysis]  = ROOT.TGraph(n)
            xoneSigma_low[analysis] = ROOT.TGraph(n)
            xoneSigma_high[analysis] = ROOT.TGraph(n)
            if not blind: xobserved[analysis]  = ROOT.TGraph(n)
        for i, mass in enumerate(masses):
            twoSigma[analysis].SetPoint(i,masses[i],quartiles[analysis][4][i])
            twoSigma[analysis].SetPoint(n+i,masses[n-i-1],quartiles[analysis][0][n-i-1])
            oneSigma[analysis].SetPoint(i,masses[i],quartiles[analysis][3][i])
            oneSigma[analysis].SetPoint(n+i,masses[n-i-1],quartiles[analysis][1][n-i-1])
            oneSigma_low[analysis].SetPoint(i,masses[i],quartiles[analysis][3][i])
            oneSigma_high[analysis].SetPoint(i,masses[i],quartiles[analysis][1][i])
            expected[analysis].SetPoint(i,masses[i],quartiles[analysis][2][i])
            if not blind: observed[analysis].SetPoint(i,masses[i],quartiles[analysis][5][i])
            if analysis in ['AP','PP']:
                xtwoSigma[analysis].SetPoint(i,masses[i],quartiles[analysis][4][i]*xsecMap[analysis][masses[i]])
                xtwoSigma[analysis].SetPoint(n+i,masses[n-i-1],quartiles[analysis][0][n-i-1]*xsecMap[analysis][masses[n-i-1]])
                xoneSigma[analysis].SetPoint(i,masses[i],quartiles[analysis][3][i]*xsecMap[analysis][masses[i]])
                xoneSigma[analysis].SetPoint(n+i,masses[n-i-1],quartiles[analysis][1][n-i-1]*xsecMap[analysis][masses[n-i-1]])
                xoneSigma_low[analysis].SetPoint(i,masses[i],quartiles[analysis][3][i]*xsecMap[analysis][masses[i]])
                xoneSigma_high[analysis].SetPoint(i,masses[i],quartiles[analysis][1][i]*xsecMap[analysis][masses[i]])
                xexpected[analysis].SetPoint(i,masses[i],quartiles[analysis][2][i]*xsecMap[analysis][masses[i]])
                if not blind: xobserved[analysis].SetPoint(i,masses[i],quartiles[analysis][5][i]*xsecMap[analysis][masses[i]])
        twoSigma[analysis].SetFillColor(ROOT.kYellow)
        twoSigma[analysis].SetLineColor(ROOT.kYellow)
        twoSigma[analysis].SetMarkerStyle(0)
        oneSigma[analysis].SetFillColor(ROOT.kSpring)
        oneSigma[analysis].SetLineColor(ROOT.kSpring)
        oneSigma[analysis].SetMarkerStyle(0)
        expected[analysis].SetLineStyle(7)
        expected[analysis].SetMarkerStyle(0)
        expected[analysis].SetFillStyle(0)
        if not blind:
            observed[analysis].SetMarkerStyle(0)
            observed[analysis].SetFillStyle(0)
        if analysis in ['AP','PP']:
            xtwoSigma[analysis].SetFillColor(ROOT.kYellow)
            xtwoSigma[analysis].SetLineColor(ROOT.kYellow)
            xtwoSigma[analysis].SetMarkerStyle(0)
            xoneSigma[analysis].SetFillColor(ROOT.kSpring)
            xoneSigma[analysis].SetLineColor(ROOT.kSpring)
            xoneSigma[analysis].SetMarkerStyle(0)
            xexpected[analysis].SetLineStyle(7)
            xexpected[analysis].SetMarkerStyle(0)
            xexpected[analysis].SetFillStyle(0)
            if not blind:
                xobserved[analysis].SetMarkerStyle(0)
                xobserved[analysis].SetFillStyle(0)

    # create canvas
    canvas = ROOT.TCanvas('cxs%s'%bp,'cxs%s'%bp,50,50,800,600)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    #canvas.SetLeftMargin(0.12)
    #canvas.SetRightMargin(0.04)
    #canvas.SetTopMargin(0.08)
    #canvas.SetBottomMargin(0.12)
    #canvas.SetLogy(1)

    appad = ROOT.TPad("appad", "top pad", 0.0, 0.5, 1.0, 1.0)
    appad.SetLeftMargin(0.12)
    appad.SetRightMargin(0.04)
    appad.SetTopMargin(0.08)
    appad.SetBottomMargin(0.12)
    appad.SetTickx(1)
    appad.SetTicky(1)
    appad.SetLogy(1)
    appad.Draw()
    pppad = ROOT.TPad("pppad", "bottom pad", 0.0, 0.0, 1.0, 0.5)
    pppad.SetLeftMargin(0.12)
    pppad.SetRightMargin(0.04)
    pppad.SetTopMargin(0.08)
    pppad.SetBottomMargin(0.12)
    pppad.SetFillColor(0)
    pppad.SetTickx(1)
    pppad.SetTicky(1)
    pppad.SetLogy(1)
    pppad.Draw()

    # create the limit on cross section plot
    legend = {}
    for analysis in ['AP','PP']:
        if analysis in ['AP']: appad.cd()
        if analysis in ['PP']: pppad.cd()
        xsecGraph[analysis].GetXaxis().SetLimits(_3L_MASSES[0],_3L_MASSES[-1])
        xsecGraph[analysis].GetXaxis().SetTitle('#Phi^{++} Mass (GeV)')
        if analysis in ['AP']: xsecGraph[analysis].GetYaxis().SetTitle('#sigma*BR (pb)')
        if analysis in ['PP']: xsecGraph[analysis].GetYaxis().SetTitle('#sigma*BR^{2} (pb)')
        xsecGraph[analysis].GetYaxis().SetTitleOffset(1.)
        xsecGraph[analysis].GetYaxis().SetTitleSize(0.06)

        xsecGraph[analysis].Draw('')
        xsecGraph[analysis].SetMaximum(1)
        xtwoSigma[analysis].Draw('f')
        xoneSigma[analysis].Draw('f')
        xexpected[analysis].Draw('same')
        xsecGraph[analysis].Draw('same')
        ROOT.gPad.RedrawAxis()
        if not blind: xobserved[analysis].Draw('same')

        legend[analysis] = ROOT.TLegend(0.5,0.52,0.90,0.88)
        legend[analysis].SetFillColor(0)
        if not blind: legend[analysis].AddEntry(xobserved[analysis], 'Observed', 'l')
        legend[analysis].AddEntry(xexpected[analysis], 'Expected', 'l')
        legend[analysis].AddEntry(xtwoSigma[analysis], 'Expected 2#sigma', 'F')
        legend[analysis].AddEntry(xoneSigma[analysis], 'Expected 1#sigma', 'F')
        name = 'Pair Production Cross Section' if analysis in ['PP'] else 'Associated Production Cross Section'
        legend[analysis].AddEntry(xsecGraph[analysis], name, 'l')
        legend[analysis].AddEntry(None,bpLabels[bp],'')
        legend[analysis].SetTextSize(0.05)

        legend[analysis].Draw('same')

        lumiperiod = 2 if period == 8 else 4
        CMS_lumi.wrtieExtraText = True
        CMS_lumi.extraText = "Preliminary" if not blind else "Simulation Preliminary"
        CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (4.9)
        CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (19.7)
        CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (25.0)
        if analysis in ['AP']: CMS_lumi.CMS_lumi(appad,lumiperiod,11)
        if analysis in ['PP']: CMS_lumi.CMS_lumi(pppad,lumiperiod,11)

    canvas.cd()

    # plot label
    #label = ROOT.TLatex(0.15,0.7,bpLabels[bp])
    #label.SetNDC(True)
    #label.SetTextFont(52)
    #label.Draw()

    save(savename,saveDir,canvas)

