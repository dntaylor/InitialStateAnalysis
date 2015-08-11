#!/usr/bin/env python
from plotters.plotUtils import *
from plotters.plotUtils import _3L_MASSES, _4L_MASSES, ZMASS
from plotters.Plotter import Plotter
from plotters import tdrstyle
from math import sqrt
import ROOT

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()
ROOT.gStyle.SetPalette(1)

def getVals(plotter,passCut,denomCut):
    bgPass = 0
    bgTotal = 0
    for bg in plotter.backgrounds:
        bgPass += plotter.getNumEntries(passCut,bg,doError=True)[0]
        bgTotal += plotter.getNumEntries(denomCut,bg,doError=True)[0]
    sigPass = 0
    sigTotal = 0
    for sig in plotter.signal:
        sigPass += plotter.getNumEntries(passCut,sig,doError=True)[0]
        sigTotal += plotter.getNumEntries(denomCut,sig,doError=True)[0]
    p = sigPass/sigTotal if sigTotal else sigPass
    f = bgPass/bgTotal if bgTotal else bgPass
    val = sigPass/sqrt(bgPass) if bgPass else sigPass
    return (p,f,val)

def getVal(plotter,cut,genCut):
    bgPass = 0
    for bg in plotter.backgrounds:
        bgPass += plotter.getNumEntries(cut,bg,doError=True)[0]
    sigPass = 0
    for sig in plotter.signal:
        sigPass += plotter.getNumEntries(cut+' & '+genCut,sig,doError=True)[0]
    return (bgPass,sigPass)


def initializePlotter(analysis, period, mass, plotName, nl, runTau):
    ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,period,analysis)
    saves = '%s_%s_%iTeV' % (analysis,analysis,period)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    regionBackground = getChannelBackgrounds(period)
    channels, leptons = getChannels(nl,runTau=runTau)
    mergeDict = getMergeDict(period)
    scaleFactor = 'event.pu_weight*event.lep_scale*event.trig_scale'
    plotter = Plotter(analysis,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName=plotName)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in regionBackground[analysis]])
    plotter.initializeSignalSamples([sigMap[period]['Sig']])
    plotter.setIntLumi(intLumiMap[period])
    return plotter

def optimizeSelections(analysis, period):
    nl = 3 if analysis=='Hpp3l' else 4

    hNumTaus = range(3)
    
    masses = _3L_MASSES
    mbinning = [17,175,1025]
    allMasses = range(200,1050,50)
    allMasses = masses
    
    myCut = 'select.passTight'
    
    tauFlavor = {
        0: ['ee','em','mm'],
        1: ['et','mt'],
        2: ['tt']
    }

    genChannels = {
        'ee': ['eee','eem','eet'],
        'em': ['eme','emm','emt'],
        'et': ['ete','etm','ett'],
        'mm': ['mme','mmm','mmt'],
        'mt': ['mte','mtm','mtt'],
        'tt': ['tte','ttm','ttt'],
    }

    recoChannels = {
        'ee': ['eee','eem'],
        'em': ['eme','emm','mee','mem'],
        'mm': ['mme','mmm'],
        'et': ['eee','eme','eem','emm','mee','mem'],
        'mt': ['mee','mem','mme','mmm','eme','emm'],
        'tt': ['eee','eem','eme','emm','mee','mem','mme','mmm'],
    }

    plotTypes = {
        'h1mass'     : [50,0,500],
        'h1massUnder': [50,0,500],
        'h1massOver' : [50,0,500],
        #'h2mass'     : [50,0,500],
        #'h2massUnder': [50,0,500],
        #'h2massOver' : [50,0,500],
        'z1mass'     : [20,0,100],
        #'z2mass'     : [20,0,100],
        'dPhi1'      : [32,0,3.2],
        #'dPhi2'      : [32,0,3.2],
        'dR1'        : [60,0,6.],
        #'dR2'        : [60,0,6.],
        'met'        : [20,0,100],
        'st'         : [150,0,1500],
    }

    axisLabels = {
        'h1mass'     : '|m(l^{#pm}l^{#pm})-Mass hyp.|< (GeV)',
        'h1massUnder': 'Mass hyp.-m(l^{#pm}l^{#pm})< (GeV)',
        'h1massOver' : 'm(l^{#pm}l^{#pm})-Mass hyp.< (GeV)',
        'h2mass'     : '|m(l^{#pm}l^{#pm})-Mass hyp.|<  (GeV)',
        'h2massUnder': 'Mass hyp.-m(l^{#pm}l^{#pm})< (GeV)',
        'h2massOver' : 'm(l^{#pm}l^{#pm})-Mass hyp.< (GeV)',
        'z1mass'     : '|m(l^{#pm}l^{#mp})-m_{Z}|>  (GeV)',
        'z2mass'     : '|m(l^{#pm}l^{#mp})-m_{Z}|>  (GeV)',
        'dPhi1'      : '|#Delta#phi(l^{#pm}l^{#pm})|<',
        'dPhi2'      : '|#Delta#phi(l^{#pm}l^{#pm})|<',
        'dR1'        : '#DeltaR(l^{#pm}l^{#pm})<',
        'dR2'        : '#DeltaR(l^{#pm}l^{#pm})<',
        'met'        : 'E_{T}^{miss}> (GeV)',
        'st'         : 's_{T}> (GeV)',
    }

    sevenTeVFunctions = {
        'h1massUnder': [['x-0.9*x'],['x-x/2'],['x-(x/2-20)']],
        'h1massOver' : [['1.1*x-x'],['1.1*x-x'],['1.1*x-x']],
        'z1mass'     : [['80'],['80'],['50']],
        'dPhi1'      : [['x/600+1.95'],['x/200+1.15'],['2.1']],
        'met'        : [['0'],['20'],['40']],
        'st'         : [['1.1*x+60'],['0.85*x+125'],['x-10','200']]
    }
    if nl==4:
        sevenTeVFunctions = {
            'h1massUnder': [['x-0.9*x'],['x-x/2'],['0']],
            'h1massOver' : [['1.1*x-x'],['1.1*x-x'],['0']],
            'z1mass'     : [['0'],['10'],['50']],
            'dPhi1'      : [['3.14159'],['3.14159'],['2.5']],
            'st'         : [['0.6*x+130'],['x+100','400'],['120']]
        }

    vals = {}
    valTypes = ['p','f','sOverSqrtB']
    
    for nTaus in hNumTaus:
        vals[nTaus] = {}
        genSelect = '(' + ' | '.join(['genChannel=="%s"'%x for y in tauFlavor[nTaus] for x in genChannels[y] + ['aaa']]) + ')'
        recoSelect = '(' + ' | '.join(['channel=="%s"'%x for y in tauFlavor[nTaus] for x in recoChannels[y]]) + ')'

        for plotType,binning in plotTypes.iteritems():
            if nl==3 and (plotType=='h2mass' or plotType=='z2mass' or plotType=='dR2'): continue
            vals[nTaus][plotType] = {}
            #for m in range(mbinning[0]):
            for m in range(len(allMasses)):
                mass = allMasses[m]
                vals[nTaus][plotType][m] = {}
                if mass not in masses: continue
                plotter = initializePlotter(analysis, period, mass, 'optimize%i'%mass, nl, True)
                print 'Optimize %s %i TeV %s %i Taus %i GeV' %(analysis,period,plotType,nTaus,mass)

                denomCut = myCut + '&&' + recoSelect
                if plotType=='h1massUnder': denomCut += '&&h1.mass<%f' % mass
                if plotType=='h1massOver' : denomCut += '&&h1.mass>%f' % mass
                if plotType=='h2massUnder': denomCut += '&&h2.mass<%f' % mass
                if plotType=='h2massOver' : denomCut += '&&h2.mass>%f' % mass
                bins = [x*binning[2]/binning[0] for x in range(binning[0])]
                bgDenom, sigDenom = getVal(plotter,denomCut,genSelect)
                for cut in range(binning[0]):
                    if plotType=='h1mass'     : passCut = 'fabs(h1.mass-%f)<%f&&%s' % (mass, bins[cut], denomCut)
                    if plotType=='h1massUnder': passCut = '%f-h1.mass<%f&&h1.mass<%f&&%s' % (mass, bins[cut], mass, denomCut)
                    if plotType=='h1massOver' : passCut = 'h1.mass-%f<%f&&h1.mass>%f&&%s' % (mass, bins[cut], mass, denomCut)
                    if plotType=='h2mass'     : passCut = 'fabs(h2.mass-%f)<%f&&%s' % (mass, bins[cut], denomCut)
                    if plotType=='h2massUnder': passCut = '%f-h2.mass<%f&&h2.mass<%f&&%s' % (mass, bins[cut], mass, denomCut)
                    if plotType=='h2massOver' : passCut = 'h2.mass-%f<%f&&h2.mass>%f&&%s' % (mass, bins[cut], mass, denomCut)
                    if plotType=='z1mass'     : passCut = 'fabs(z1.mass-%f)>%f&&%s' % (ZMASS, bins[cut], denomCut)
                    if plotType=='z2mass'     : passCut = 'fabs(z2.mass-%f)>%f&&%s' % (ZMASS, bins[cut], denomCut)
                    if plotType=='dPhi1'      : passCut = 'fabs(h1.dPhi)<%f&&%s' % (bins[cut], denomCut)
                    if plotType=='dPhi2'      : passCut = 'fabs(h2.dPhi)<%f&&%s' % (bins[cut], denomCut)
                    if plotType=='dR1'        : passCut = 'h1.dR<%f&&%s' % (bins[cut], denomCut)
                    if plotType=='dR2'        : passCut = 'h2.dR<%f&&%s' % (bins[cut], denomCut)
                    if plotType=='met'        : passCut = 'finalstate.met>%f&&%s' % (bins[cut], denomCut)
                    if plotType=='st'         : passCut = 'finalstate.sT>%f&&%s' % (bins[cut], denomCut)
                    bgNum, sigNum = getVal(plotter,passCut,genSelect)
                    p = sigNum/sigDenom
                    f = bgNum/bgDenom
                    if not bgNum: bgNum = 0.0001
                    val = sigNum/sqrt(bgNum)
                    vals[nTaus][plotType][m][cut] = (p,f,val)

            canvas = ROOT.TCanvas('%s_%i'%(plotType,nTaus),'%s_%i'%(plotType,nTaus),50,50,800,600)
            canvas.SetFillColor(0)
            canvas.SetBorderMode(0)
            canvas.SetFrameFillStyle(0)
            canvas.SetFrameBorderMode(0)
            canvas.SetLeftMargin(0.12)
            canvas.SetRightMargin(0.14)
            canvas.SetTopMargin(0.08)
            canvas.SetBottomMargin(0.12)

            for v in range(2):
                valType = valTypes[v]
                plot = ROOT.TH2F('%s_%i_%s' %(plotType,nTaus,valType),\
                                 '%s_%i_%s; H^{#pm#pm} Mass (GeV); %s'%(plotType,nTaus,valType,axisLabels[plotType]),\
                                 *mbinning+binning)
                savedir = 'plots/%s_optimization_%iTeV/png' % (analysis, period)
                python_mkdir(savedir)
                savename = '%s/%s_%s_%i_%s.png' % (savedir, analysis, plotType, nTaus, valType)
                #for m in range(mbinning[0]):
                for m in range(len(allMasses)):
                    for c,cutVals in vals[nTaus][plotType][m].iteritems():
                        plot.SetBinContent(m+1,c+1,cutVals[v])
                plot.GetYaxis().SetTitleOffset(1.)
                if v<2: plot.GetZaxis().SetRangeUser(0.,1.)
                plot.Draw('colz')
                if plotType in sevenTeVFunctions:
                    funcs = {}
                    for func in sevenTeVFunctions[plotType][nTaus]:
                        funcs[func] = ROOT.TF1("7TeV_%i_%s"%(nTaus,plotType),func,mbinning[1],mbinning[2])
                        funcs[func].SetLineColor(ROOT.kBlack)
                        funcs[func].SetLineWidth(3)
                        funcs[func].Draw("lsame")
                canvas.Print(savename)
        

optimizeSelections('Hpp3l',8)
#optimizeSelections('Hpp4l',13)
