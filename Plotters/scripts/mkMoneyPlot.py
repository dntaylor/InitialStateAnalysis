#!/usr/bin/env python

from InitialStateAnalysis.Plotters import CMS_lumi, tdrstyle
import ROOT
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()


# sub categories
limitNames = ['ee100','em100','mm100','et100','mt100','BP1','BP2','BP3','BP4']

# without charge id and with tau
limits = {
    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : 533, 'AP3l' : 509, 'PP3l' : 507, 'AP': 509, 'PP': 592, 'Comb' : 620,},
    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : 545, 'AP3l' : 510, 'PP3l' : 513, 'AP': 510, 'PP': 604, 'Comb' : 626,},
    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 515, 'PP3l' : 517, 'AP': 515, 'PP': 609, 'Comb' : 635,},
    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' : 261, 'AP3l' : 313, 'PP3l' : 383, 'AP': 313, 'PP': 394, 'Comb' : 422,},
    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : 278, 'AP3l' : 326, 'PP3l' : 406, 'AP': 326, 'PP': 417, 'Comb' : 443,},
    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' :   1, 'AP3l' : 141, 'PP3l' :  54, 'AP': 141, 'PP':   1, 'Comb' :   1,},
    'BP1'   : {'name' : 'BP1',                                'PP4l' : 372, 'AP3l' : 417, 'PP3l' : 449, 'AP': 417, 'PP': 478, 'Comb' : 511,},
    'BP2'   : {'name' : 'BP2',                                'PP4l' : 462, 'AP3l' : 469, 'PP3l' : 498, 'AP': 469, 'PP': 535, 'Comb' : 574,},
    'BP3'   : {'name' : 'BP3',                                'PP4l' : 482, 'AP3l' : 476, 'PP3l' : 477, 'AP': 476, 'PP': 533, 'Comb' : 575,},
    'BP4'   : {'name' : 'BP4',                                'PP4l' : 432, 'AP3l' : 447, 'PP3l' : 477, 'AP': 447, 'PP': 513, 'Comb' : 544,},
}

## with charge id and with tau
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : 511, 'AP3l' : 498, 'PP3l' : 492, 'AP': 498, 'PP': 557, 'Comb' : 603,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : 529, 'AP3l' : 504, 'PP3l' : 507, 'AP': 504, 'PP': 593, 'Comb' : 615,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 513, 'PP3l' : 517, 'AP': 513, 'PP': 609, 'Comb' : 632,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' : 257, 'AP3l' : 309, 'PP3l' : 378, 'AP': 309, 'PP': 388, 'Comb' : 406,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : 278, 'AP3l' : 323, 'PP3l' : 405, 'AP': 323, 'PP': 415, 'Comb' : 428,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' :   1, 'AP3l' : 145, 'PP3l' :  42, 'AP': 145, 'PP': 60, 'Comb' : 165,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : 407, 'AP3l' : 428, 'PP3l' : 482, 'AP': 428, 'PP': 510, 'Comb' : 530,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : 402, 'AP3l' : 443, 'PP3l' : 458, 'AP': 443, 'PP': 489, 'Comb' : 519,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : 471, 'AP3l' : 468, 'PP3l' : 469, 'AP': 468, 'PP': 524, 'Comb' : 557,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : 422, 'AP3l' : 438, 'PP3l' : 468, 'AP': 438, 'PP': 506, 'Comb' : 529,},
#}

## without charge id and without tau
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : 533, 'AP3l' : 509, 'PP3l' : 507, 'AP': 509, 'PP': 592, 'Comb' : 620,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : 545, 'AP3l' : 510, 'PP3l' : 513, 'AP': 510, 'PP': 604, 'Comb' : 624,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 515, 'PP3l' : 517, 'AP': 515, 'PP': 609, 'Comb' : 633,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : 370, 'AP3l' : 416, 'PP3l' : 307, 'AP': 416, 'PP': 401, 'Comb' : 463,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : 461, 'AP3l' : 469, 'PP3l' : 328, 'AP': 469, 'PP': 504, 'Comb' : 542,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : 482, 'AP3l' : 476, 'PP3l' : 444, 'AP': 476, 'PP': 519, 'Comb' : 557,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : 430, 'AP3l' : 445, 'PP3l' : 384, 'AP': 445, 'PP': 464, 'Comb' : 510,},
#}

## with charge id and without tau
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : 511, 'AP3l' : 498, 'PP3l' : 492, 'AP': 498, 'PP': 557, 'Comb' : 603,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : 529, 'AP3l' : 504, 'PP3l' : 507, 'AP': 504, 'PP': 593, 'Comb' : 615,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 513, 'PP3l' : 517, 'AP': 513, 'PP': 609, 'Comb' : 632,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' :   0, 'AP3l' :   0, 'PP3l' :   0, 'AP':   0, 'PP':   0, 'Comb' :   0,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : 405, 'AP3l' : 427, 'PP3l' : 350, 'AP': 427, 'PP': 432, 'Comb' : 490,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : 401, 'AP3l' : 443, 'PP3l' : 368, 'AP': 443, 'PP': 443, 'Comb' : 498,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : 471, 'AP3l' : 468, 'PP3l' : 435, 'AP': 468, 'PP': 511, 'Comb' : 545,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : 420, 'AP3l' : 437, 'PP3l' : 376, 'AP': 437, 'PP': 457, 'Comb' : 503,},
#}


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

pp4l = ROOT.TColor.GetColor('#FFAE03')
pp3l = ROOT.TColor.GetColor('#E67F0D')
ap3l = ROOT.TColor.GetColor('#FE4E00')
comb = ROOT.TColor.GetColor('#E9190F')

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
    bar.SetFillColor(comb)
    bar.SetLineColor(comb)
    bar.Draw('f')
    # draw pair
    limitBars['pair%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['pair%s' %(lim)]
    bar.SetPoint(0,0,binPos+.15)
    bar.SetPoint(1,limits[lim]['PP4l'],binPos+.15)
    bar.SetPoint(2,limits[lim]['PP4l'],binPos+.45)
    bar.SetPoint(3,0,binPos+.45)
    bar.SetFillColor(pp4l)
    bar.SetLineColor(pp4l)
    bar.Draw('f')
    # draw pair 3l
    limitBars['pair3l%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['pair3l%s' %(lim)]
    bar.SetPoint(0,0,binPos-.15)
    bar.SetPoint(1,limits[lim]['PP3l'],binPos-.15)
    bar.SetPoint(2,limits[lim]['PP3l'],binPos+.15)
    bar.SetPoint(3,0,binPos+.15)
    bar.SetFillColor(pp3l)
    bar.SetLineColor(pp3l)
    bar.Draw('f')
    # draw associated
    limitBars['associated%s' %(lim)] = ROOT.TGraph(4)
    bar = limitBars['associated%s' %(lim)]
    bar.SetPoint(0,0,binPos-.45)
    bar.SetPoint(1,limits[lim]['AP3l'],binPos-.45)
    bar.SetPoint(2,limits[lim]['AP3l'],binPos-.15)
    bar.SetPoint(3,0,binPos-.15)
    bar.SetFillColor(ap3l)
    bar.SetLineColor(ap3l)
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
    'Comb': {'name': 'Combined',                   'color' : comb},
    'PP4l': {'name': 'Pair Production (4l)',       'color' : pp4l},
    'PP3l': {'name': 'Pair Production (3l)',       'color' : pp3l},
    'AP3l': {'name': 'Associated Production (3l)', 'color' : ap3l},
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
