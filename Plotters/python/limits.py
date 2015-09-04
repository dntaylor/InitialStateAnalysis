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

def save(savename,saveDir,canvas):
    '''Save the limits into root file and images.'''
    #for type in ['png', 'pdf', 'eps']:
    for type in ['png']:
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
    blind = kwargs.pop('blind',True)
    bp = kwargs.pop('branchingPoint','')
    bgMode = kwargs.pop('bgMode','sideband')
    do4l = kwargs.pop('do4l',False)

    # directories    
    datacardDir = '%s/%s_%itev_%s' % (datacardBaseDir, analysis, period, region)
    if bp: datacardDir += '/%s' % bp
    limitDataDir = '%s/%s_%itev_%s' % (limitDataBaseDir, analysis, period, region)
    if bp: limitDataDir += '/%s' % bp
    datacardString = '' if bgMode == "sideband" else "_{0}".format(bgMode)
    if do4l: datacardString += '_4l'

    # masses to include
    masses = _3L_MASSES if analysis == 'Hpp3l' else _4L_MASSES
    if do4l: masses = _4L_MASSES
    if period==13: masses = [500]
    n = len(masses)

    # get cross sections
    xsecMap = {}
    xsecGraph = ROOT.TGraph(n)
    for i,mass in enumerate(masses):
        sample = 'HPlusPlusHMinusMinusHTo4L_M-{0}_8TeV-pythia6' if analysis in ['Hpp4l','HppComb'] or do4l else 'HPlusPlusHMinusHTo3L_M-{0}_8TeV-calchep-pythia6'
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
        fname = os.path.join(limitDataDir, "higgsCombineTest.Asymptotic.mH%i%s.root" % (mass,datacardString))
        file = ROOT.TFile(fname,"READ")
        tree = file.Get("limit")
        if not tree: continue
        for i, row in enumerate(tree):
            quartiles[i,j] = row.limit

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
    expected.GetXaxis().SetLimits(masses[0],masses[-1])
    expected.GetXaxis().SetTitle('#Phi^{++} Mass (GeV)')
    expected.GetYaxis().SetTitle('95% CLs Upper Limit on #sigma/#sigma_{model}')
    expected.GetYaxis().SetTitleOffset(1.)
    expected.GetYaxis().SetTitleSize(0.05)

    expected.Draw()
    twoSigma.Draw('f')
    oneSigma.Draw('f')
    expected.Draw('same')
    ROOT.gPad.RedrawAxis()
    if not blind: observed.Draw()

    ratiounity = ROOT.TLine(expected.GetXaxis().GetXmin(),1,expected.GetXaxis().GetXmax(),1)
    ratiounity.Draw()

    legend = ROOT.TLegend(0.65,0.2,0.90,0.4)
    legend.SetFillColor(0)
    if not blind: legend.AddEntry(observed, 'Observed','l')
    legend.AddEntry(expected, 'Expected','l')
    legend.AddEntry(twoSigma, 'Expected 2#sigma', 'F')
    legend.AddEntry(oneSigma, 'Expected 1#sigma', 'F')

    legend.Draw('same')

    lumiperiod = 2 if period == 8 else 4
    CMS_lumi.wrtieExtraText = True
    CMS_lumi.extraText = "Preliminary" if not blind else "Simulation Preliminary"
    CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (4.9)
    CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (19.7)
    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (25.0)
    CMS_lumi.CMS_lumi(canvas,lumiperiod,11)

    if do4l: savename += '_4l'
    save(savename,saveDir,canvas)

    canvas.Clear()

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
    if not blind: xobserved.Draw()

    legend = ROOT.TLegend(0.5,0.6,0.90,0.85)
    legend.SetFillColor(0)
    if not blind: legend.AddEntry(xobserved, 'Observed', 'l')
    legend.AddEntry(xexpected, 'Expected', 'l')
    legend.AddEntry(xtwoSigma, 'Expected 2#sigma', 'F')
    legend.AddEntry(xoneSigma, 'Expected 1#sigma', 'F')
    name = 'Pair Production Cross Section' if analysis in ['Hpp4l'] or do4l else 'Associated Production Cross Section'
    if analysis in ['HppComb']: name = 'Cross Section'
    legend.AddEntry(xsecGraph, name, 'l')

    legend.Draw('same')

    lumiperiod = 2 if period == 8 else 4
    CMS_lumi.wrtieExtraText = True
    CMS_lumi.extraText = "Preliminary" if not blind else "Simulation Preliminary"
    CMS_lumi.lumi_7TeV = "%0.1f fb^{-1}" % (4.9)
    CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (19.7)
    CMS_lumi.lumi_13TeV = "%0.1f fb^{-1}" % (25.0)
    CMS_lumi.CMS_lumi(canvas,lumiperiod,11)

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

    print "Expected Limit: %i GeV (+%i, -%i)" % (x, h-x, x-l)

