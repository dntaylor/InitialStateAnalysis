#!/usr/bin/env python

from InitialStateAnalysis.Plotters import CMS_lumi, tdrstyle
import ROOT
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()


# sub categories
limitNames = ['ee100','em100','mm100','et100','mt100','BP1','BP2','BP3','BP4']

# expected limites
#limits = {
#    'ee100' : {'name' : '100% H^{#pm#pm}#rightarrow ee',      'PP4l' : 533, 'AP3l' : 544, 'PP3l' : 507, 'AP': 544, 'PP': 592, 'Comb' : 636,},
#    'em100' : {'name' : '100% H^{#pm#pm}#rightarrow e#mu',    'PP4l' : 544, 'AP3l' : 544, 'PP3l' : 512, 'AP': 544, 'PP': 604, 'Comb' : 642,},
#    'mm100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 553, 'PP3l' : 516, 'AP': 553, 'PP': 609, 'Comb' : 653,},
#    'et100' : {'name' : '100% H^{#pm#pm}#rightarrow e#tau',   'PP4l' : 262, 'AP3l' : 339, 'PP3l' : 382, 'AP': 339, 'PP': 398, 'Comb' : 437,},
#    'mt100' : {'name' : '100% H^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : 285, 'AP3l' : 354, 'PP3l' : 407, 'AP': 354, 'PP': 421, 'Comb' : 457,},
#    'tt100' : {'name' : '100% H^{#pm#pm}#rightarrow #tau#tau','PP4l' :   1, 'AP3l' : 118, 'PP3l' :   6, 'AP': 118, 'PP':   1, 'Comb' :   1,},
#    'BP1'   : {'name' : 'BP1',                                'PP4l' : 372, 'AP3l' : 448, 'PP3l' : 450, 'AP': 448, 'PP': 477, 'Comb' : 522,},
#    'BP2'   : {'name' : 'BP2',                                'PP4l' : 462, 'AP3l' : 505, 'PP3l' : 499, 'AP': 505, 'PP': 535, 'Comb' : 597,},
#    'BP3'   : {'name' : 'BP3',                                'PP4l' : 482, 'AP3l' : 512, 'PP3l' : 477, 'AP': 512, 'PP': 533, 'Comb' : 601,},
#    'BP4'   : {'name' : 'BP4',                                'PP4l' : 432, 'AP3l' : 481, 'PP3l' : 477, 'AP': 481, 'PP': 513, 'Comb' : 563,},
#}
# observed limits
limits = {
    'ee100' : {'name' : '100% #Phi^{#pm#pm}#rightarrow ee',      'PP4l' : 533, 'AP3l' : 546, 'PP3l' : 508, 'AP': 546, 'PP': 593, 'Comb' : 636,},
    'em100' : {'name' : '100% #Phi^{#pm#pm}#rightarrow e#mu',    'PP4l' : 545, 'AP3l' : 550, 'PP3l' : 517, 'AP': 550, 'PP': 607, 'Comb' : 645,},
    'mm100' : {'name' : '100% #Phi^{#pm#pm}#rightarrow #mu#mu',  'PP4l' : 562, 'AP3l' : 557, 'PP3l' : 518, 'AP': 557, 'PP': 611, 'Comb' : 654,},
    'et100' : {'name' : '100% #Phi^{#pm#pm}#rightarrow e#tau',   'PP4l' : 266, 'AP3l' : 313, 'PP3l' : 344, 'AP': 313, 'PP': 355, 'Comb' : 371,},
    'mt100' : {'name' : '100% #Phi^{#pm#pm}#rightarrow #mu#tau', 'PP4l' : 289, 'AP3l' : 317, 'PP3l' : 350, 'AP': 317, 'PP': 392, 'Comb' : 419,},
    'tt100' : {'name' : '100% #Phi^{#pm#pm}#rightarrow #tau#tau','PP4l' :   1, 'AP3l' : 123, 'PP3l' :   1, 'AP': 123, 'PP':   1, 'Comb' :   1,},
    'BP1'   : {'name' : 'Benchmark 1',                           'PP4l' : 373, 'AP3l' : 458, 'PP3l' : 457, 'AP': 458, 'PP': 483, 'Comb' : 529,},
    'BP2'   : {'name' : 'Benchmark 2',                           'PP4l' : 463, 'AP3l' : 509, 'PP3l' : 503, 'AP': 509, 'PP': 539, 'Comb' : 601,},
    'BP3'   : {'name' : 'Benchmark 3',                           'PP4l' : 483, 'AP3l' : 515, 'PP3l' : 484, 'AP': 515, 'PP': 537, 'Comb' : 603,},
    'BP4'   : {'name' : 'Benchmark 4',                           'PP4l' : 434, 'AP3l' : 494, 'PP3l' : 489, 'AP': 494, 'PP': 520, 'Comb' : 572,},
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
h =  ROOT.TH2F("h", "h; Excluded Masses (GeV}); ", 1,0,1000,nl+2,0.5,nl+2.5)
h.GetYaxis().SetRangeUser(2,nl+1)
h.Draw()

pp3l = ROOT.TColor.GetColor('#00BD39')
pp4l = ROOT.TColor.GetColor('#760BAA')
ap3l = ROOT.TColor.GetColor('#0B5FA5')
comb = ROOT.TColor.GetColor('#FF9400')

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
    'Comb': {'name': '#splitline{Associated Production}{ + Pair Production}',                   'color' : comb},
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
