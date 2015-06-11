#!/usr/bin/env python

from plotters.Plotter import Plotter
from plotters.ShapePlotter import ShapePlotter
from plotters.CutFlowPlotter import CutFlowPlotter
from plotters.EfficiencyPlotter import EfficiencyPlotter
from plotters.FakeRatePlotter import FakeRatePlotter
from plotters.plotUtils import *
import argparse
import itertools
import sys

ZMASS = 91.1876

def plotKinematicsMethod(plotMethod,variables,savename,cuts,**kwargs):
    savedir = kwargs.pop('savedir','')
    if savedir: savedir += '/'
    if type(variables) is not list: variables = [variables]
    for var, c in zip(variables, cuts):
        plotMethod(var+'.Pt',[50,0,250],savedir+var+savename+'Pt',cut=c,yaxis='Events/5.0 GeV/c^{2}',xaxis='p_{T} (GeV/c)',legendpos=43,**kwargs)
        plotMethod(var+'.Eta',[30,-3.0,3.0],savedir+var+savename+'Eta',cut=c,yaxis='Events',xaxis='#eta',legendpos=43,**kwargs)
        plotMethod(var+'.Phi',[30,-3.14159,3.14159],savedir+var+savename+'Phi',cut=c,yaxis='Events',xaxis='#phi',legendpos=43,**kwargs)


def plotDistributions(plotMethod,myCut,nl,isControl,**kwargs):
    savedir = kwargs.pop('savedir','')
    analysis = kwargs.pop('analysis','')
    region = kwargs.pop('region','')
    if savedir: savedir += '/'
    plotMethod('finalstate.sT',[40,0,1000],savedir+'sT',yaxis='Events/25.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.sT',[50,0,500],savedir+'sT_zoom',yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    #plotMethod('finalstate.jetVeto20',[8,0,8],savedir+'numJets20',yaxis='Events',xaxis='Number of Jets (p_{T}>20 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.jetVeto30',[8,0,8],savedir+'numJets30',yaxis='Events',xaxis='Number of Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.elecVetoLoose',[8,0,8],savedir+'numElectronsLoose',yaxis='Events',xaxis='Number of Electrons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.muonVetoLoose',[8,0,8],savedir+'numMuonsLoose',yaxis='Events',xaxis='Number of Muons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.elecVetoTight',[8,0,8],savedir+'numElectronsTight',yaxis='Events',xaxis='Number of Electrons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.muonVetoTight',[8,0,8],savedir+'numMuonsTight',yaxis='Events',xaxis='Number of Muons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.bjetVeto30Loose',[8,0,8],savedir+'numBJets30Loose',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.bjetVeto30Medium',[8,0,8],savedir+'numBJets30Medium',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.bjetVeto30Tight',[8,0,8],savedir+'numBJets30Tight',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.jetVeto40',[8,0,8],savedir+'numJets40',yaxis='Events',xaxis='Number of Jets (p_{T}>40 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.muonVeto5',[8,0,8],savedir+'muonVeto5',yaxis='Events',xaxis='Muon Veto (p_{T}>5 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.muonVeto10Loose',[8,0,8],savedir+'muonVeto10',yaxis='Events',xaxis='Muon Veto (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.muonVeto15',[8,0,8],savedir+'muonVeto15',yaxis='Events',xaxis='Muon Veto (p_{T}>15 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.elecVeto10',[8,0,8],savedir+'elecVeto10',yaxis='Events',xaxis='Electron Veto (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.met',[40,0,200],savedir+'met',yaxis='Events/5.0 GeV',xaxis='E_{T}^{miss} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.mass',[40,0,400],savedir+'mass',yaxis='Events/10.0 GeV',xaxis='Mass (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.mass',[150,0,300],savedir+'mass_zoom',yaxis='Events/2.0 GeV',xaxis='Mass (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('event.nvtx',[50,0,50],savedir+'puVertices',yaxis='Events',xaxis='Number PU Vertices',legendpos=43,logy=0,cut=myCut,**kwargs)
    if analysis in ['Hpp3l','Hpp4l'] or region in ['Hpp2l']:
        plotMethod('h1.mass',[24,0,600],savedir+'hppMass',yaxis='Events/25.0 GeV/c^{2}',xaxis='M(l^{+}l^{+}) (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.mass',[24,0,600],savedir+'hppMass_mod',yaxis='Events/25.0 GeV/c^{2}',xaxis='M(l^{+}l^{+}) (GeV/c^{2})',lumitext=33,legendpos=41,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.dPhi',[32,0,3.2],savedir+'hppDphi',yaxis='Events/0.1 rad',xaxis='#Delta#phi(l^{+}l^{+}) (rad)',legendpos=41,lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('h1.Pt',[40,0,400],savedir+'hppPt',yaxis='Events/10.0 GeV',xaxis='p_{T}^{#Phi^{++}} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Pt1',[40,0,200],savedir+'hppLeadingLeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{#Phi^{++} Leading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Pt2',[40,0,200],savedir+'hppSubleadingLeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{#Phi^{++} Subleading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Iso1',[50,0,0.5],savedir+'hppLeadingIso',yaxis='Events',xaxis='Iso/p_{T} (#Phi^{++} Leading Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Iso2',[50,0,0.5],savedir+'hppSubleadingIso',yaxis='Events',xaxis='Iso/p_{T} (#Phi^{++} Subleading Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.dR',[60,0,6],savedir+'hppdR',yaxis='Events',xaxis='#DeltaR(l^{+}l^{+})',legendpos=43,logy=0,cut=myCut,**kwargs)
    if analysis in ['Z', 'Hpp3l', 'Hpp4l', 'WZ'] or region in ['Z', 'TT']:
        plotMethod('z1.mass',[42,70,112],savedir+'z1Mass',yaxis='Events/1.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass',[7,80.5,101.5],savedir+'z1Mass_wideBin',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass',[80,0,240],savedir+'z1Mass_fullWindow',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass',[80,0,240],savedir+'z1Mass_fullWindow_log',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=1,cut=myCut,**kwargs)
        plotMethod('z1.Pt',[40,0,400],savedir+'z1Pt',yaxis='Events/10.0 GeV',xaxis='p_{T}^{Z1} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Pt1',[40,0,200],savedir+'z1LeadingLeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{Z1 Leading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Pt2',[40,0,200],savedir+'z1SubleadingLeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{Z1 Subleading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Iso1',[50,0,0.5],savedir+'z1LeadingIso',yaxis='Events',xaxis='Iso/p_{T} (Z1 Leading Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Iso2',[50,0,0.5],savedir+'z1SubleadingIso',yaxis='Events',xaxis='Iso/p_{T} (Z1 Subleading Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.dR',[60,0,6],savedir+'z1dR',yaxis='Events',xaxis='#DeltaR',legendpos=43,logy=0,cut=myCut,**kwargs)
    if analysis in ['Hpp4l']:
        plotMethod('z2.mass',[42,70,112],savedir+'z2Mass',yaxis='Events/1.0 GeV',xaxis='M(l^{+}l^{-}) (Z2) (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[7,80.5,101.5],savedir+'z2Mass_wideBin',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z2) (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[80,0,240],savedir+'z2Mass_fullWindow',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z2) (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[80,0,240],savedir+'z2Mass_fullWindow_log',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z2) (GeV)',legendpos=43,logy=1,cut=myCut,**kwargs)
        plotMethod('z2.Pt',[40,0,400],savedir+'z2Pt',yaxis='Events/10.0 GeV',xaxis='p_{T}^{Z2} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Pt1',[40,0,200],savedir+'z2LeadingLeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{Z2 Leading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Pt2',[40,0,200],savedir+'z2SubleadingLeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{Z2 Subleading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Iso1',[50,0,0.5],savedir+'z2LeadingIso',yaxis='Events',xaxis='Iso/p_{T} (Z2 Leading Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Iso2',[50,0,0.5],savedir+'z2SubleadingIso',yaxis='Events',xaxis='Iso/p_{T} (Z2 Subleading Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
    if analysis in ['Hpp3l', 'WZ']:
        plotMethod('w1.Pt',[40,0,400],savedir+'w1Pt',yaxis='Events/10.0 GeV',xaxis='p_{T}^{W} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.Pt1',[40,0,200],savedir+'w1LeptonPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{W Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.Iso1',[50,0,0.5],savedir+'w1Iso',yaxis='Events',xaxis='Iso/p_{T} (W Lepton)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.mass',[40,0,200],savedir+'w1Mass',yaxis='Events/5.0 GeV',xaxis='M_{T}^{W} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.dPhi',[32,0,3.2],savedir+'w1dPhi',yaxis='Events/0.1 rad',xaxis='#Delta#phi(W lepton, E_{T}^{miss}) (rad)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('w1.dR1_z1_1',[60,0,6],savedir+'w1dR1_1',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{leading lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('w1.dR1_z1_2',[60,0,6],savedir+'w1dR1_2',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{subleading lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.mll_z1_1',[80,0,240],savedir+'dilepton_mass_1_ss',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#pm}) (Z_{l1},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1==z1.Chg1 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod('w1.mll_z1_2',[80,0,240],savedir+'dilepton_mass_2_ss',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#pm}) (Z_{l2},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1==z1.Chg2 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod('w1.mll_z1_1',[80,0,240],savedir+'dilepton_mass_1_os',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#mp}) (Z_{l1},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1!=z1.Chg1 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod('w1.mll_z1_2',[80,0,240],savedir+'dilepton_mass_2_os',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#mp}) (Z_{l2},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1!=z1.Chg2 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod(['w1.mll_z1_1','w1.mll_z1_2'],[80,0,240],savedir+'dilepton_mass',yaxis='Events/3.0 GeV',xaxis='M(ll) (Z_{l},W_{l}) (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)


def plotRegion(analysis,channel,runPeriod,**kwargs):
    '''A function to simplify plotting multiple channels and run periods.'''
    blind = kwargs.pop('blind',True)
    mass = kwargs.pop('mass',500)
    runTau = kwargs.pop('runTau',False)
    myCut = kwargs.pop('myCut','1')
    plotFinalStates = kwargs.pop('plotFinalStates',False)
    plotJetBins = kwargs.pop('plotJetBins',False)
    plotOverlay = kwargs.pop('plotOverlay',False)
    plotShapes = kwargs.pop('plotShapes',False)
    plotCutFlow = kwargs.pop('plotCutFlow',False)
    finalStatesToPlot = kwargs.pop('finalStates','all')
    nostack = kwargs.pop('nostack',False)
    normalize = kwargs.pop('normalize',False)
    scaleFactor = kwargs.pop('scaleFactor','event.pu_weight*event.lep_scale*event.trig_scale')
    useSignal = analysis in ['Hpp3l','Hpp4l']
    for key, value in kwargs.iteritems():
        print "Unrecognized parameter '" + key + "' = " + str(value)
        return 0

    if useSignal: print "MKPLOTS:%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass)
    isControl = analysis != channel
    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ'   : 3,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[analysis]
    ntuples = 'ntuples%s_%stev_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = {
        #'Hpp2l' : ['T', 'TT', 'TTV', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Hpp2l' : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        #'Z'     : ['T', 'TT', 'TTV', 'W', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Z'     : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        #'WZ'    : ['T', 'TT', 'TTV', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'WZ'    : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        #'TT'    : ['T', 'TT', 'TTV', 'W', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'TT'    : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        #'Hpp3l' : ['T', 'TT', 'TTV', 'Z', 'ZG', 'VVV', 'DB'],
        'Hpp3l' : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Hpp4l' : ['TT', 'Z', 'DB']
    }
    if runPeriod==13:
        channelBackground = {
            'WZ'    : ['T', 'TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'TT'    : ['T', 'TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'Hpp3l' : ['T', 'TT', 'TTV', 'Z', 'DB'],
            'Hpp4l' : ['T', 'TT', 'Z', 'TTV', 'DB']
        }

    finalStates, leptons = getChannels(nl,runTau=runTau)
    if finalStatesToPlot=='all':
        fsToPlot = finalStates
    else:
        fsToPlot = finalStatesToPlot.split(',')
    print 'MKPLOTS:%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut)
    dataplot = (isControl or not blind) and runPeriod in [7,8]
    mergeDict = getMergeDict(runPeriod)

    # Plotting discriminating variables
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotMCDataRatio' if dataplot else 'plotMC'
    plotMethod = getattr(plotter,plotMode)
    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables" % (analysis,channel, runPeriod)
    plotDistributions(plotMethod,myCut,nl,isControl,analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
    #if channel in ['Z', 'TT']: # plot the 3rd lepton veto selection
    #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables loose e veto" % (analysis,channel, runPeriod)
    #    plotDistributions(plotMethod,myCut+' & finalstate.elecVetoLoose==1 & finalstate.muonVetoLoose==0',nl,isControl,savedir='eVetoLoose',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
    #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables loose m veto" % (analysis,channel, runPeriod)
    #    plotDistributions(plotMethod,myCut+' & finalstate.muonVetoLoose==1 & finalstate.elecVetoLoose==0',nl,isControl,savedir='mVetoLoose',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
    #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables tight e veto" % (analysis,channel, runPeriod)
    #    plotDistributions(plotMethod,myCut+' & finalstate.elecVetoTight==1 & finalstate.muonVetoLoose==0',nl,isControl,savedir='eVetoTight',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
    #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables tight m veto" % (analysis,channel, runPeriod)
    #    plotDistributions(plotMethod,myCut+' & finalstate.muonVetoTight==1 & finalstate.elecVetoLoose==0',nl,isControl,savedir='mVetoTight',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)

    # plotting kinematic plots
    #print "MKPLOTS:%s:%s:%iTeV: Plotting kinematics" % (analysis, channel, runPeriod)
    #plotKinematicsMethod(plotMethod,leptons,'Lepton',[myCut]*nl,logy=0)
    #plotKinematicsMethod(plotMethod,leptons,'Electron',['%sFlv=="e"&&%s' %(x,myCut) for x in leptons],logy=0)
    #plotKinematicsMethod(plotMethod,leptons,'Muon',['%sFlv=="m"&&%s' %(x,myCut) for x in leptons],logy=0)
    #if runTau: plotKinematicsMethod(plotMethod,leptons,'Tau',['%sFlv=="t"&&%s' %(x,myCut) for x in leptons],logy=0)

    # each channel
    if plotFinalStates:
        print "MKPLOTS:%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod)
        for c in fsToPlot:
            print "MKPLOTS:%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c)
            plotDistributions(plotMethod,myCut+'&&channel=="%s"'%c,nl,isControl,savedir=c,analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
            #if channel in ['Z', 'TT']: # plot the 3rd lepton veto selection
            #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables loose e veto" % (analysis,channel, runPeriod)
            #    plotDistributions(plotMethod,myCut+' & finalstate.elecVetoLoose==1 & finalstate.muonVetoLoose==0 & channel=="%s"'%c,nl,isControl,savedir=c+'/eVetoLoose',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
            #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables loose m veto" % (analysis,channel, runPeriod)
            #    plotDistributions(plotMethod,myCut+' & finalstate.muonVetoLoose==1 & finalstate.elecVetoLoose==0 & channel=="%s"'%c,nl,isControl,savedir=c+'/mVetoLoose',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
            #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables tight e veto" % (analysis,channel, runPeriod)
            #    plotDistributions(plotMethod,myCut+' & finalstate.elecVetoTight==1 & finalstate.muonVetoLoose==0 & channel=="%s"'%c,nl,isControl,savedir=c+'/eVetoTight',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
            #    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables tight m veto" % (analysis,channel, runPeriod)
            #    plotDistributions(plotMethod,myCut+' & finalstate.muonVetoTight==1 & finalstate.elecVetoLoose==0 & channel=="%s"'%c,nl,isControl,savedir=c+'/mVetoTight',analysis=analysis,region=channel,nostack=nostack,normalize=normalize)

    # some partially blind plots for h++
    if runPeriod==8 and not dataplot and analysis in ['Hpp3l']:
        print "MKPLOTS:%s:%s:%iTeV: Plotting partially blinded variables" % (analysis, channel, runPeriod)
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotModeUnblind = 'plotMCDataRatio'
        plotMethodUnblind = getattr(plotter,plotModeUnblind)
        plotMethodUnblind('h1.mass',[24,0,600],'hppMass_unblind',yaxis='Events/25.0 GeV/c^{2}',xaxis='M(l^{+}l^{+}) (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,blinder=[150,99999])
        plotMethodUnblind('finalstate.sT',[50,0,500],'sT_unblind',yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,blinder=[200,99999])
        plotMethodUnblind('z1.mass',[42,70,112],'z1Mass_unblind',yaxis='Events/1.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut,blinder=[112,99999])
        plotMethodUnblind('z1.mass',[80,0,240],'z1Mass_fullWindow_unblind',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut,blinder=[112,99999])
        if plotFinalStates:
            print "MKPLOTS:%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod)
            for c in fsToPlot:
                print "MKPLOTS:%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c)
                plotMethodUnblind('h1.mass',[24,0,600],c+'/hppMass_unblind',yaxis='Events/25.0 GeV/c^{2}',xaxis='M(l^{+}l^{+}) (GeV/c^{2})',lumitext=33,logy=0,cut=myCut+'&&channel=="%s"'%c,overflow=True,blinder=[150,99999])
                plotMethodUnblind('finalstate.sT',[50,0,500],c+'/sT_unblind',yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut+'&&channel=="%s"'%c,overflow=True,blinder=[200,99999])
                plotMethodUnblind('z1.mass',[42,70,112],c+'/z1Mass_unblind',yaxis='Events/1.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut+'&&channel=="%s"'%c,blinder=[112,99999])
                plotMethodUnblind('z1.mass',[80,0,240],c+'/z1Mass_fullWindow_unblind',yaxis='Events/3.0 GeV',xaxis='M(l^{+}l^{-}) (Z1) (GeV)',legendpos=43,logy=0,cut=myCut+'&&channel=="%s"'%c,blinder=[112,99999])

    # setup signal overlay plots
    if plotOverlay and useSignal:
        plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_overlay',mergeDict=mergeDict,scaleFactor=scaleFactor)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotMCDataSignalRatio' if dataplot else 'plotMCSignalRatio'
        plotMethod = getattr(plotter,plotMode)
        # plot the signal overlay
        print "MKPLOTS:%s:%s:%iTeV: Plotting signal overlay discriminating variables" % (analysis, channel, runPeriod)
        plotDistributions(plotMethod,myCut,nl,isControl,savedir='overlay',signalscale=100,analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
        #print "MKPLOTS:%s:%s:%iTeV: Plotting signal overlay kinematics" % (analysis, channel, runPeriod)
        #plotKinematicsMethod(plotMethod,leptons,'Lepton',[myCut]*nl,signalscale=1000,savedir='overlay',logy=1)
        #plotKinematicsMethod(plotMethod,leptons,'Electron',['%sFlv=="e"&&%s' %(x,myCut) for x in leptons],signalscale=100,savedir='overlay',logy=1)
        #plotKinematicsMethod(plotMethod,leptons,'Muon',['%sFlv=="m"&&%s' %(x,myCut) for x in leptons],signalscale=100,savedir='overlay',logy=1)

    # plot shapes
    if plotShapes:
        print "MKPLOTS:%s:%s:%iTeV: Plotting shapes" % (analysis, channel, runPeriod)
        plotter = ShapePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_shapes',mergeDict=mergeDict,scaleFactor=scaleFactor)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotter.plotMC('z1.mass',['channel=="mmm"','channel=="emm"'],[42,70,112],'zMass_mc_mm',yaxis='Normalized',xaxis='M(l^{+}l^{-}) (GeV)',logy=0,cut=myCut,cutNames=['mmm','emm'])
        if dataplot: plotter.plotData('z1.mass',['channel=="mmm"','channel=="emm"'],[42,70,112],'zMass_data_mm',yaxis='Normalized',xaxis='M(l^{+}l^{-}) (GeV)',logy=0,cut=myCut,cutNames=['mmm','emm'])
        plotter.plotMC('z1.mass',['channel=="eee"','channel=="eme"'],[42,70,112],'zMass_mc_ee',yaxis='Normalized',xaxis='M(l^{+}l^{-}) (GeV)',logy=0,cut=myCut,cutNames=['eee','eme'])
        if dataplot: plotter.plotData('z1.mass',['channel=="eee"','channel=="eme"'],[42,70,112],'zMass_data_ee',yaxis='Normalized',xaxis='M(l^{+}l^{-}) (GeV)',logy=0,cut=myCut,cutNames=['eee','eme'])

    # plot cut flows (each cut)
    print "MKPLOTS:%s:%s:%iTeV: Plotting cut flow" % (analysis, channel, runPeriod)
    cutFlowMap = {}
    cutFlowMap[channel] = defineCutFlowMap(channel,finalStates,mass)
    #print cutFlowMap
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutFlowSelections',mergeDict=mergeDict,scaleFactor=scaleFactor)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotMCDataRatio' if dataplot else 'plotMC'
    plotMethod = getattr(plotter,plotMode)
    if plotCutFlow:
        for i in range(len(cutFlowMap[channel]['cuts'])):
            print 'MKPLOTS:%s:%s:%iTeV: Plotting cut flow selections %s' % (analysis, channel, runPeriod, cutFlowMap[channel]['labels_simple'][i])
            thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
            plotDistributions(plotMethod,'%s&%s'%(myCut,thisCut),nl,isControl,savedir='cutflow/%s'%cutFlowMap[channel]['labels_simple'][i],analysis=analysis,region=channel,nostack=nostack,normalize=normalize)
            if plotFinalStates:
                print "MKPLOTS:%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod)
                for c in fsToPlot:
                    print "MKPLOTS:%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c)
                    plotDistributions(plotMethod,'%s&channel=="%s"&%s'%(myCut,c,thisCut),nl,isControl,savedir='cutflow/%s/%s' %(cutFlowMap[channel]['labels_simple'][i],c),analysis=analysis,region=channel,nostack=nostack,normalize=normalize)


    # plot cut flows on same plot
    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflow',mergeDict=mergeDict,scaleFactor=scaleFactor)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    if useSignal: plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
    plotMethod = getattr(plotter,plotMode)
    plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'cutFlow',labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)
    if plotFinalStates:
        for c in fsToPlot:
            print "MKPLOTS:%s:%s:%iTeV: Plotting cut flow  %s" % (analysis, channel, runPeriod, c)
            plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)

    # setup individual channel cuts on same plot
    plotChannelStrings, plotChannelCuts = getChannelStringsCuts(channel,finalStates)
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    plotMethod = getattr(plotter,plotMode)
    plotMethod([myCut]+['%s&&%s' %(x,myCut) for x in plotChannelCuts],'individualChannels',labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0)
    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflowSelectionsChannels',mergeDict=mergeDict,scaleFactor=scaleFactor)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    plotMethod = getattr(plotter,plotMode)
    for i in range(len(cutFlowMap[channel]['cuts'])):
        thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
        plotMethod(['%s&&%s'%(myCut,thisCut)]+['%s&&%s&&%s' %(x,myCut,thisCut) for x in plotChannelCuts],'cutflow/%s/individualChannels'%cutFlowMap[channel]['labels_simple'][i],labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0)

    # plot efficiencies
    if analysis in ['WZ']:
        print "MKPLOTS:%s:%s:%iTeV: Plotting efficiency" % (analysis, channel, runPeriod)
        plotter = EfficiencyPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_efficiency',mergeDict=mergeDict,scaleFactor=scaleFactor)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel] if x != 'WZ'])
        plotter.initializeSignalSamples([sigMap[runPeriod]['WZ']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
        plotMethod = getattr(plotter,plotMode)
        plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'efficiency',labels=cutFlowMap[channel]['labels'],lumitext=33)
        if plotFinalStates:
            for c in fsToPlot:
                print "MKPLOTS:%s:%s:%iTeV: Plotting efficiency  %s" % (analysis, channel, runPeriod, c)
                plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/efficiency'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)

        # plot cut flows overlays
        plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflow_overlay',mergeDict=mergeDict,scaleFactor=scaleFactor)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel] if x not in ['WZ']])
        plotter.initializeSignalSamples([sigMap[runPeriod]['WZ']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
        plotMethod = getattr(plotter,plotMode)
        plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'cutFlow_overlay',labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)
        if plotFinalStates:
            for c in fsToPlot:
                print "MKPLOTS:%s:%s:%iTeV: Plotting cut flow overlay  %s" % (analysis, channel, runPeriod, c)
                plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow_overlay'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)


def plotFakeRate(analysis,channel,runPeriod,**kwargs):
    '''Plot fake rate for an analysis.'''
    ntuples = 'ntuples%s_%stev_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV_Fakes' % (analysis,channel,runPeriod)
    nl = 3 if analysis in ['Hpp3l','WZ'] else 4
    mass = kwargs.pop('mass',500)
    isControl = False
    sigMap = getSigMap(nl,mass)
    mergeDict = getMergeDict(runPeriod)
    intLumiMap = getIntLumiMap()
    channelBackground = {
        'WZ' : ['T','TT', 'TTV', 'Z', 'VVV', 'ZZ','WZ'],
        'Hpp3l' : ['T','TT', 'TTV','Z','VVV','DB'],
        'Hpp4l' : ['TT','Z','DB'],
        'FakeRate' : ['T','TT', 'TTV','Z','VVV','DB'],
    }
    if runPeriod==13:
        channelBackground = {
            'WZ' : ['T','TT', 'TTV', 'Z', 'ZZ','WZ'],
            'Hpp3l' : ['T', 'TT', 'TTV','Z','DB'],
            'Hpp4l' : ['T', 'TT', 'Z', 'TTV','DB'],
            'FakeRate' : ['T', 'TT', 'TTV','Z','DB'],
        }

    # setup selections
    denomSelection = 'select.pass_elecLoose0p20 & select.pass_muonLoose0p20 & finalstate.bjetVeto30==0 & (z1.Pt1>20.&z1.Pt2>10.) & fabs(z1.mass-91.1876)<20. & w1.Pt1>20. & w1.met<30.'
    eDenom = denomSelection+' & w1Flv=="e"'
    mDenom = denomSelection+' & w1Flv=="m"'
    ePtBins = [20,30,40,60,100,1000]
    mPtBins = [20,30,40,60,100,1000]
    eEtaBins = [0,1.479,2.5]
    mEtaBins = [0,1.2,2.4]

    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    #plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotMode = 'plotMC'
    plotMethod = getattr(plotter,plotMode)
    print "MKPLOTS:%s:%s:%iTeV: Plotting discriminating variables" % (analysis,channel, runPeriod)
    plotDistributions(plotMethod,denomSelection,nl,isControl,analysis=analysis)
    plotDistributions(plotMethod,eDenom,nl,isControl,analysis=analysis,savedir='elec')
    plotDistributions(plotMethod,mDenom,nl,isControl,analysis=analysis,savedir='muon')

    # now plot the fake rates
    print "MKPLOTS:%s:%s:%iTeV: Computing fake rates" % (analysis,channel, runPeriod)
    plotter = FakeRatePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='fakerates',mergeDict=mergeDict)
    # this should be done in data... using mc until we get some!
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    #plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])
    plotter.plotFakeRate('w1.pass_elecMedium0p15_1 & '+eDenom, eDenom, 'fakeRate_elec', ptBins=ePtBins, etaBins=eEtaBins, logx=1)
    plotter.plotFakeRate('w1.pass_muonTight0p12_1 & '+mDenom, mDenom, 'fakeRate_muon', ptBins=mPtBins, etaBins=mEtaBins, logx=1)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['Z','WZ','Hpp2l','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['Z','WZ','TT','Hpp2l','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-pf','--plotFinalStates',action='store_true',help='Plot individual final states')
    parser.add_argument('-pj','--plotJetBins',action='store_true',help='Plot jet bins')
    parser.add_argument('-po','--plotOverlay',action='store_true',help='Plot overlay')
    parser.add_argument('-ps','--plotShapes',action='store_true',help='Plot shapes')
    parser.add_argument('-pcf','--plotCutFlow',action='store_true',help='Plot cutflow distributions')
    parser.add_argument('-rt','--runTau',action='store_true',help='Run Tau finalStates (not implemented)')
    parser.add_argument('-ub','--unblind',action='store_false',help='Unblind signal channel')
    parser.add_argument('-ns','--nostack',action='store_true',help='Plot histograms unstacked')
    parser.add_argument('-no','--normalize',action='store_true',help='Plot histograms normalized to 1')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500)
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses')
    parser.add_argument('-ac','--allControls',action='store_true',help='Run over all controls for a given analysis (3l, 4l)')
    parser.add_argument('-fr','--doFakeRate',action='store_true',help='Make fake rate plots and output fake rate histograms')
    parser.add_argument('-fs','--finalStates',type=str,default='all',help='Only run given channels (ie: "eee,emm")')
    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to plots.')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for plots.')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    massLists = {
        13 : {
            'Hpp3l' : [500],
            'Hpp4l' : [500]
        }, 
        8 : {
            '3l' : [170, 200, 250, 300, 350, 400, 450, 500, 600, 700],
            '4l' : [110, 130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]
        }
    }

    controlList = {
        'Hpp3l' : ['Hpp3l', 'WZ'],
        'Hpp4l' : ['Hpp4l']
    }

    if args.period == 7:
        print "7 TeV not implemented"
    elif args.doFakeRate:
        plotFakeRate(args.analysis,args.channel,args.period,mass=args.mass)
    else:
        plotRegion(args.analysis,args.channel,args.period,plotFinalStates=args.plotFinalStates,runTau=args.runTau,blind=args.unblind,mass=args.mass,plotJetBins=args.plotJetBins,plotOverlay=args.plotOverlay,plotShapes=args.plotShapes,plotCutFlow=args.plotCutFlow,myCut=args.cut,finalStates=args.finalStates,nostack=args.nostack,normalize=args.normalize,scaleFactor=args.scaleFactor)

    return 0


if __name__ == "__main__":
    main()
