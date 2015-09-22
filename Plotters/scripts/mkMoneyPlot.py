#!/usr/bin/env python

from InitialStateAnalysis.Plotters import CMS_lumi, tdrstyle
import ROOT
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()


# sub categories
limitNames = ['ee100','em100','mm100','et100','mt100','BP1','BP2','BP3','BP4']

## without charge id and with tau
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#}
#
## with charge id and with tau
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
#}
#
## without charge id and without tau
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : 533, 'AP3l' : 509, 'PP3l' : 507, 'AP': 509, 'PP': 592, 'Comb' : 620,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : 545, 'AP3l' : 510, 'PP3l' : 513, 'AP': 510, 'PP': 604, 'Comb' : 624,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 515, 'PP3l' : 517, 'AP': 515, 'PP': 609, 'Comb' : 633,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : 406, 'AP3l' : 430, 'PP3l' : 351, 'AP': 430, 'PP': 434, 'Comb' : 494,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : 416, 'AP3l' : 456, 'PP3l' : 386, 'AP': 456, 'PP': 462, 'Comb' : 512,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : 482, 'AP3l' : 476, 'PP3l' : 444, 'AP': 476, 'PP': 519, 'Comb' : 557,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : 430, 'AP3l' : 445, 'PP3l' : 384, 'AP': 445, 'PP': 464, 'Comb' : 510,},
#}

# with charge id and without tau
limits = {
    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'BP1'   : {'name' : 'BP1',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'BP2'   : {'name' : 'BP2',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'BP3'   : {'name' : 'BP3',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
    'BP4'   : {'name' : 'BP4',                                'PP4l' : , 'AP3l' : , 'PP3l' : , 'AP': , 'PP': , 'Comb' : ,},
}


# now make the plot
canvas = ROOT.TCanvas('cMoneyPlot','cMoneyPlot',50,50,800,600)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin(0.22)
canvas.SetRightMargin(0.04)
canvas.SetTopMargin(0.08)
canvas.SetBottomMargin(0.12)


nl = len(limitNames)
h =  ROOT.TH2F("h", "h; Excluded Masses (GeV/c^{2}); ", 1,0,1000,nl+2,0.5,nl+2.5)
h.GetYaxis().SetRangeUser(2,nl+1)
h.Draw()

binPos = nl+2
limitBars = {}
for l,lim in enumerate(limitNames):
    binPos -= 1
    h.GetYaxis().SetBinLabel(binPos,limits[lim]['name'])
    # draw combined
    limitBars['combined%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['combined%s' %(lim)]
    bar.SetPoint(0,0,binPos-.45)
    bar.SetPoint(1,limits[lim]['Comb'],binPos-.45)
    bar.SetPoint(2,limits[lim]['Comb'],binPos+.45)
    bar.SetPoint(3,0,binPos+.45)
    bar.SetFillColor(ROOT.kOrange+8)
    bar.SetLineColor(ROOT.kOrange+8)
    bar.Draw('f')
    # draw pair
    limitBars['pair%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['pair%s' %(lim)]
    bar.SetPoint(0,0,binPos+.15)
    bar.SetPoint(1,limits[lim]['PP4l'],binPos+.15)
    bar.SetPoint(2,limits[lim]['PP4l'],binPos+.45)
    bar.SetPoint(3,0,binPos+.45)
    bar.SetFillColor(ROOT.kOrange-1)
    bar.SetLineColor(ROOT.kOrange-1)
    bar.Draw('f')
    # draw pair 3l
    limitBars['pair3l%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['pair3l%s' %(lim)]
    bar.SetPoint(0,0,binPos-.15)
    bar.SetPoint(1,limits[lim]['PP3l'],binPos-.15)
    bar.SetPoint(2,limits[lim]['PP3l'],binPos+.15)
    bar.SetPoint(3,0,binPos+.15)
    bar.SetFillColor(ROOT.kOrange+2)
    bar.SetLineColor(ROOT.kOrange+2)
    bar.Draw('f')
    # draw associated
    limitBars['associated%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['associated%s' %(lim)]
    bar.SetPoint(0,0,binPos-.45)
    bar.SetPoint(1,limits[lim]['AP3l'],binPos-.45)
    bar.SetPoint(2,limits[lim]['AP3l'],binPos-.15)
    bar.SetPoint(3,0,binPos-.15)
    bar.SetFillColor(ROOT.kOrange)
    bar.SetLineColor(ROOT.kOrange)
    bar.Draw('f')

#set the colors and size for the legend
markerSize  = 1.0

latex = ROOT.TLatex()
n_ = 3

x1_l = 0.92
y1_l = 0.80

dx_l = 0.20
dy_l = 0.04+0.06*n_
x0_l = x1_l-dx_l
y0_l = y1_l-dy_l

legend =  ROOT.TPad("legend_0","legend_0",x0_l,y0_l,x1_l, y1_l )
#legend.SetFillColor( rt.kGray )
legend.Draw()
legend.cd()

ar_l = dy_l/dx_l
#gap_ = 0.09/ar_l
gap_ = 1./(n_+1)
bwx_ = 0.12
bwy_ = gap_/1.5
    
x_l = [1.2*bwx_]
#y_l = [1-(1-0.10)/ar_l]
y_l = [1-gap_]
ex_l = [0]
ey_l = [0.04/ar_l]

#array must be converted 
x_l = array("f",x_l)
ex_l = array("f",ex_l)
y_l = array("f",y_l)
ey_l = array("f",ey_l)

latex.SetTextFont(42)
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)
latex.SetTextSize(0.1)
latex.SetTextAlign(12)

box_ = ROOT.TBox()
xx_ = x_l[0]
yy_ = y_l[0]

dataStyles = {
    'Comb': {'name': 'Combined',                   'color' : ROOT.kOrange+8,},
    'PP4l': {'name': 'Pair Production (4l)',       'color' : ROOT.kOrange-1,},
    'PP3l': {'name': 'Pair Production (3l)',       'color' : ROOT.kOrange+2,},
    'AP3l': {'name': 'Associated Production (3l)', 'color' : ROOT.kOrange,},
}

for n in ['Comb','PP4l','PP3l','AP3l']:
    box_.SetLineStyle( ROOT.kSolid )
    box_.SetLineWidth( 1 )
    # box_.SetLineColor( kBlack )
    box_.SetLineColor( dataStyles[n]['color'])
    box_.SetFillColor( dataStyles[n]['color'])
    box_.SetFillStyle(1001)
    box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
    box_.SetFillStyle(0)
    box_.DrawBox( xx_-bwx_/2, yy_-bwy_/2, xx_+bwx_/2, yy_+bwy_/2 )
    #Draw Z->ee text
    latex.DrawLatex(xx_+1.*bwx_,yy_,dataStyles[n]['name'])
    yy_ -= gap_


lumiperiod = 2
CMS_lumi.wrtieExtraText = True
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_8TeV = "%0.1f fb^{-1}" % (19.7)
CMS_lumi.CMS_lumi(canvas,lumiperiod,33)

canvas.RedrawAxis()

canvas.Print('plots/limits/png/moneyPlot.png')
