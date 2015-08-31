#!/usr/bin/env python

import argparse
import itertools
import sys
import pickle
import json
import logging
from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.ShapePlotter import ShapePlotter
from InitialStateAnalysis.Plotters.CutFlowPlotter import CutFlowPlotter
from InitialStateAnalysis.Plotters.EfficiencyPlotter import EfficiencyPlotter
from InitialStateAnalysis.Plotters.FakeRatePlotter import FakeRatePlotter
from InitialStateAnalysis.Plotters.CorrelationPlotter import CorrelationPlotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS, _3L_MASSES, _4L_MASSES

def plotDistributions(plotMethod,myCut,nl,isControl,**kwargs):
    savedir = kwargs.pop('savedir','')
    analysis = kwargs.pop('analysis','')
    region = kwargs.pop('region','')
    mass = kwargs.pop('mass',500)
    if savedir: savedir += '/'
    plotMethod('finalstate.sT',[40,0,1000],savedir+'sT',yaxis='Events/25.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',lumitext=33,logy=1,ymin=0.1,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.sT',[50,0,500],savedir+'sT_zoom',yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    #plotMethod('finalstate.jetVeto20',[8,0,8],savedir+'numJets20',yaxis='Events',xaxis='Number of Jets (p_{T}>20 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.jetVeto30',[8,0,8],savedir+'numJets30',yaxis='Events',xaxis='Number of Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.elecVetoLoose',[8,0,8],savedir+'numElectronsLoose',yaxis='Events',xaxis='Number of Electrons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.muonVetoLoose',[8,0,8],savedir+'numMuonsLoose',yaxis='Events',xaxis='Number of Muons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.elecVetoTight',[8,0,8],savedir+'numElectronsTight',yaxis='Events',xaxis='Number of Electrons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.muonVetoTight',[8,0,8],savedir+'numMuonsTight',yaxis='Events',xaxis='Number of Muons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.bjetVeto30Loose',[8,0,8],savedir+'numBJets30Loose',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.bjetVeto30Medium',[8,0,8],savedir+'numBJets30Medium',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.bjetVeto30Tight',[8,0,8],savedir+'numBJets30Tight',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.jetVeto40',[8,0,8],savedir+'numJets40',yaxis='Events',xaxis='Number of Jets (p_{T}>40 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.muonVeto5',[8,0,8],savedir+'muonVeto5',yaxis='Events',xaxis='Muon Veto (p_{T}>5 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.muonVeto10Loose',[8,0,8],savedir+'muonVeto10',yaxis='Events',xaxis='Muon Veto (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.muonVeto15',[8,0,8],savedir+'muonVeto15',yaxis='Events',xaxis='Muon Veto (p_{T}>15 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    #plotMethod('finalstate.elecVeto10',[8,0,8],savedir+'elecVeto10',yaxis='Events',xaxis='Electron Veto (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.met',[40,0,200],savedir+'met',yaxis='Events/5.0 GeV',xaxis='E_{T}^{miss} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.mass',[40,0,400],savedir+'mass',yaxis='Events/10.0 GeV',xaxis='M_{3\\ell} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.mass',[150,0,300],savedir+'mass_zoom',yaxis='Events/2.0 GeV',xaxis='Mass (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('event.nvtx',[50,0,50],savedir+'puVertices',yaxis='Events',xaxis='Number PU Vertices',legendpos=43,logy=0,cut=myCut,**kwargs)
    if analysis in ['WZ_FakeRate']:
        plotMethod('finalstate.leadJetPt',[40,0,200],savedir+'JetPt',yaxis='Events/5.0 GeV',xaxis='p_{T}^{jet} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('finalstate.leadJetEta',[50,-5.0,5.0],savedir+'JetEta',yaxis='Events',xaxis='\\eta^{jet}',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('finalstate.leadJetPhi',[30,-3.14159,3.14159],savedir+'JetPhi',yaxis='Events',xaxis='\\phi^{jet}',legendpos=43,logy=0,cut=myCut,**kwargs)
    # plot lepton kinematics
    for l in range(nl):
        name = 'l%i' % (l+1)
        plotMethod('%s.Pt' %name,[40,0,200],savedir+'%sPt' %name,yaxis='Events/5.0 GeV',xaxis='p_{T}^{\\ell%i} (GeV)' %(l+1),legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('%s.Eta' %name,[30,-3.0,3.0],savedir+'%sEta' %name,yaxis='Events',xaxis='\\eta^{\\ell%i}' %(l+1),legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('%s.Phi' %name,[30,-3.14159,3.14159],savedir+'%sPhi' %name,yaxis='Events',xaxis='\\phi^{\\ell%i}' %(l+1),legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('%s.Iso' %name,[50,0,.5],savedir+'%sIso' %name,yaxis='Events',xaxis='Relative Isolation (\\ell%i)' %(l+1),legendpos=43,logy=0,cut=myCut,**kwargs)    
        plotMethod('%s.ChargeConsistent' %name,[3,-1.5,1.5],savedir+'%sChargeId' %name,yaxis='Events',xaxis='Charge ID (\\ell%i)' %(l+1),legendpos=43,logy=0,cut=myCut,**kwargs)    
    names = {
        'e': 'Elec',
        'm': 'Mu',
        't': 'Tau',
    }
    tex = {
        'e': 'e',
        'm': '\\mu',
        't': '\\tau',
    }
    for l in ['e','m']:
        name = names[l]
        t = tex[l]
        cuts = ['%s & %s' %(myCut,'l%iFlv=="%s"' %((x+1),l)) for x in range(nl)]
        plotMethod(['l%i.Pt'  %(x+1) for x in range(nl)], [40,0,200],            savedir+'%sPt'  %name, yaxis='Events/5.0 GeV', xaxis='p_{T}^{%s} (GeV)' %t,        legendpos=43, logy=0, cut=cuts, overflow=True, **kwargs)
        plotMethod(['l%i.Eta' %(x+1) for x in range(nl)], [30,-3.0,3.0],         savedir+'%sEta' %name, yaxis='Events',         xaxis='\\eta^{%s}' %t,              legendpos=43, logy=0, cut=cuts, **kwargs)
        plotMethod(['l%i.Phi' %(x+1) for x in range(nl)], [30,-3.14159,3.14159], savedir+'%sPhi' %name, yaxis='Events',         xaxis='\\phi^{%s}' %t,              legendpos=43, logy=0, cut=cuts, **kwargs)
        plotMethod(['l%i.Iso' %(x+1) for x in range(nl)], [50,0,.5],             savedir+'%sIso' %name, yaxis='Events',         xaxis='Relative Isolation (%s)' %t, legendpos=43, logy=0, cut=cuts, **kwargs)    
        plotMethod(['l%i.ChargeConsistent' %(x+1) for x in range(nl)], [3,-1.5,1.5], savedir+'%sChargeId' %name, yaxis='Events', xaxis='Charge ID (%s)' %t,         legendpos=43, logy=0, cut=cuts, **kwargs)    

    # plot doubly charged higgs stuff
    if analysis in ['Hpp3l','Hpp4l'] or region in ['Hpp2l']:
        plotMethod('h1.mass', [24,0,600],savedir+'hppMass',              yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',        lumitext=33,logy=1,cut=myCut,overflow=True,**kwargs)
        #plotMethod('h1.mass', [24,0,600],savedir+'hppMass_mod',          yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{+}\\ell^{+}} (GeV/c^{2})',              lumitext=33,legendpos=41,logy=1,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.mass', [32,0,800],savedir+'hppMass_alpha',        yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})', lumitext=33,logy=1,cut=myCut,boxes=[[12,0.9*mass,1],[1.1*mass,800,1],[0.9*mass,1.1*mass,2]],**kwargs)
        plotMethod('h1.dPhi', [32,0,3.2],savedir+'hppDphi',              yaxis='Events/0.1 rad',       xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',   legendpos=41,lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('h1.Pt',   [40,0,400],savedir+'hppPt',                yaxis='Events/10.0 GeV',      xaxis='p_{T}^{\\Phi^{\\pm\\pm}} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Pt1',  [40,0,200],savedir+'hppLeadingLeptonPt',   yaxis='Events/5.0 GeV',       xaxis='p_{T}^{\\Phi^{\\pm\\pm} Leading Lepton} (GeV)',   legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Pt2',  [40,0,200],savedir+'hppSubleadingLeptonPt',yaxis='Events/5.0 GeV',       xaxis='p_{T}^{\\Phi^{\\pm\\pm} Subleading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Iso1', [50,0,0.5],savedir+'hppLeadingIso',        yaxis='Events',               xaxis='Iso/p_{T} (\\Phi^{\\pm\\pm} Leading Lepton)',     legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.Iso2', [50,0,0.5],savedir+'hppSubleadingIso',     yaxis='Events',               xaxis='Iso/p_{T} (\\Phi^{\\pm\\pm} Subleading Lepton)',  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.dR',   [60,0,6], savedir+'hppdR',                 yaxis='Events',               xaxis='\\Delta R(\\ell^{\\pm}\\ell^{\\pm})',             legendpos=43,logy=0,cut=myCut,**kwargs)
    # plot Z stuff
    if analysis in ['Z', 'Hpp3l', 'Hpp4l', 'WZ', 'WZ_W'] or region in ['Z', 'TT']:
        plotMethod('z1.mass', [42,70,112],   savedir+'z1Mass',               yaxis='Events/1.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass', [60,60,120],   savedir+'z1Mass_newWidth',      yaxis='Events/1.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass', [7,80.5,101.5],savedir+'z1Mass_wideBin',       yaxis='Events/3.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass', [80,0,240],    savedir+'z1Mass_fullWindow',    yaxis='Events/3.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.mass', [80,0,240],    savedir+'z1Mass_fullWindow_log',yaxis='Events/3.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=1,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Pt',   [40,0,400],    savedir+'z1Pt',                 yaxis='Events/10.0 GeV',xaxis='p_{T}^{Z} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Pt1',  [40,0,200],    savedir+'z1LeadingLeptonPt',    yaxis='Events/5.0 GeV', xaxis='p_{T}^{Z Leading Lepton} (GeV)',   legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Pt2',  [40,0,200],    savedir+'z1SubleadingLeptonPt', yaxis='Events/5.0 GeV', xaxis='p_{T}^{Z Subleading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Iso1', [50,0,0.5],    savedir+'z1LeadingIso',         yaxis='Events',         xaxis='Iso/p_{T} (Z Leading Lepton)',     legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.Iso2', [50,0,0.5],    savedir+'z1SubleadingIso',      yaxis='Events',         xaxis='Iso/p_{T} (Z Subleading Lepton)',  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z1.dR',   [60,0,6],      savedir+'z1dR',                 yaxis='Events',         xaxis='\\Delta R(\\ell^{+}\\ell^{-})',    legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.JetPt1',[40,0,400],savedir+'z1JetPtLead',yaxis='Events/10.0 GeV',xaxis='p_{T}^{Z Lead Lepton Jet} (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.JetBTag1',[40,-1,1],savedir+'z1JetBTagLead',yaxis='Events',xaxis='Z Lead Lepton Jet BTag',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.GenPdgId1',[80,-40,40],savedir+'z1GenPdgIdLead',yaxis='Events',xaxis='Z Lead Lepton PDG ID',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.MotherGenPdgId1',[80,-40,40],savedir+'z1MotherGenPdgIdLead',yaxis='Events',xaxis='Z Lead Lepton Mother PDG ID',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.Dxy1',[100,-.1,.1],savedir+'z1DXY1',yaxis='Events',xaxis='#Delta(lepton,pv)_{XY}',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.Dz1',[100,-.1,.1],savedir+'z1DZ1',yaxis='Events',xaxis='#Delta(lepton,pv)_{Z}',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.Dxy2',[100,-.1,.1],savedir+'z1DXY2',yaxis='Events',xaxis='#Delta(lepton,pv)_{XY}',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.Dz2',[100,-.1,.1],savedir+'z1DZ2',yaxis='Events',xaxis='#Delta(lepton,pv)_{Z}',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.JetPt2',[40,0,400],savedir+'z1JetPtSubLead',yaxis='Events/10.0 GeV',xaxis='p_{T}^{Z SubLead Lepton Jet} (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.JetBTag2',[40,-1,1],savedir+'z1JetBTagSubLead',yaxis='Events',xaxis='Z Lead SubLepton Jet BTag',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.GenPdgId2',[80,-40,40],savedir+'z1GenPdgIdSubLead',yaxis='Events',xaxis='Z SubLead Lepton PDG ID',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z1.MotherGenPdgId2',[80,-40,40],savedir+'z1MotherGenPdgIdSubLead',yaxis='Events',xaxis='Z SubLead Lepton Mother PDG ID',legendpos=43,logy=0,cut=myCut,**kwargs)
    # plot second z stuff
    if analysis in ['Hpp4l']:
        plotMethod('z2.mass',[42,70,112],   savedir+'z2Mass',               yaxis='Events/1.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[7,80.5,101.5],savedir+'z2Mass_wideBin',       yaxis='Events/3.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[80,0,240],    savedir+'z2Mass_fullWindow',    yaxis='Events/3.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[80,0,240],    savedir+'z2Mass_fullWindow_log',yaxis='Events/3.0 GeV', xaxis='M_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=1,cut=myCut,**kwargs)
        plotMethod('z2.Pt',  [40,0,400],    savedir+'z2Pt',                 yaxis='Events/10.0 GeV',xaxis='p_{T}^{Z2} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Pt1', [40,0,200],    savedir+'z2LeadingLeptonPt',    yaxis='Events/5.0 GeV', xaxis='p_{T}^{Z2 Leading Lepton} (GeV)',   legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Pt2', [40,0,200],    savedir+'z2SubleadingLeptonPt', yaxis='Events/5.0 GeV', xaxis='p_{T}^{Z2 Subleading Lepton} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Iso1',[50,0,0.5],    savedir+'z2LeadingIso',         yaxis='Events',         xaxis='Iso/p_{T} (Z2 Leading Lepton)',     legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('z2.Iso2',[50,0,0.5],    savedir+'z2SubleadingIso',      yaxis='Events',         xaxis='Iso/p_{T} (Z2 Subleading Lepton)',  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
    # plot W stuff
    if analysis in ['Hpp3l', 'WZ', 'WZ_W', 'WZ_FakeRate']:
        plotMethod('w1.Pt',  [40,0,400],savedir+'w1Pt',      yaxis='Events/10.0 GeV',xaxis='p_{T}^{W} (GeV)',                           legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.Pt1', [40,0,200],savedir+'w1LeptonPt',yaxis='Events/5.0 GeV', xaxis='p_{T}^{W Lepton} (GeV)',                    legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.Iso1',[50,0,0.5],savedir+'w1Iso',     yaxis='Events',         xaxis='Iso/p_{T} (W Lepton)',                      legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.mass',[40,0,200],savedir+'w1Mass',    yaxis='Events/5.0 GeV', xaxis='M_{T}^{W} (GeV)',                           legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.dPhi',[32,0,3.2],savedir+'w1dPhi',    yaxis='Events/0.1 rad', xaxis='\\Delta\\phi(W lepton, E_{T}^{miss}) (rad)',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.Dxy1',[100,-.1,.1],savedir+'w1DXY',yaxis='Events',xaxis='#Delta(lepton,pv)_{XY}',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.Dz1',[100,-.1,.1],savedir+'w1DZ',yaxis='Events',xaxis='#Delta(lepton,pv)_{Z}',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.JetPt1',[40,0,400],savedir+'w1JetPt',yaxis='Events/10.0 GeV',xaxis='p_{T}^{W Lepton Jet} (GeV)',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.JetBTag1',[40,-1,1],savedir+'w1JetBTag',yaxis='Events',xaxis='W Lepton Jet BTag',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.GenPdgId1',[80,-40,40],savedir+'w1GenPdgId',yaxis='Events',xaxis='W Lepton PDG ID',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.MotherGenPdgId1',[80,-40,40],savedir+'w1MotherGenPdgId',yaxis='Events',xaxis='W Lepton Mother PDG ID',legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('w1.mll_z1_1',[80,0,240],savedir+'dilepton_mass_1_ss',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#pm}) (Z_{l1},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1==z1.Chg1 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod('w1.mll_z1_2',[80,0,240],savedir+'dilepton_mass_2_ss',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#pm}) (Z_{l2},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1==z1.Chg2 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod('w1.mll_z1_1',[80,0,240],savedir+'dilepton_mass_1_os',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#mp}) (Z_{l1},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1!=z1.Chg1 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod('w1.mll_z1_2',[80,0,240],savedir+'dilepton_mass_2_os',yaxis='Events/3.0 GeV',xaxis='M(l^{#pm}l^{#mp}) (Z_{l2},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1!=z1.Chg2 & %s' %myCut,overflow=True,**kwargs)
        #plotMethod(['w1.mll_z1_1','w1.mll_z1_2'],[80,0,240],savedir+'dilepton_mass',yaxis='Events/3.0 GeV',xaxis='M(ll) (Z_{l},W_{l}) (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
    #if analysis in ['Hpp3l', 'WZ']:
    #    plotMethod('w1.dR1_z1_1',[60,0,6],savedir+'w1dR1_1',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{leading lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)
    #    plotMethod('w1.dR1_z1_2',[60,0,6],savedir+'w1dR1_2',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{subleading lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)


def plotRegion(analysis,channel,runPeriod,**kwargs):
    '''A function to simplify plotting multiple channels and run periods.'''
    logger = logging.getLogger(__name__)
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
    loglevel = kwargs.pop('loglevel','INFO')
    for key, value in kwargs.iteritems():
        logger.warning("Unrecognized parameter '" + key + "' = " + str(value))
        return 0

    if useSignal: logger.info("%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass))
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
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = getChannelBackgrounds(runPeriod)

    finalStates, leptons = getChannels(nl,runTau=runTau)
    if finalStatesToPlot=='all':
        fsToPlot = finalStates
    else:
        fsToPlot = finalStatesToPlot.split(',')
    logger.info('%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut))
    dataplot = (isControl or not blind)
    mergeDict = getMergeDict(runPeriod)
    cutFlowMap = {}
    cutFlowMap[channel] = defineCutFlowMap(channel,finalStates,mass)

    genChannels = {
        'ee': ['eee','eem','eet'],
        'em': ['eme','emm','emt'],
        'et': ['ete','etm','ett'],
        'mm': ['mme','mmm','mmt'],
        'mt': ['mte','mtm','mtt'],
        'tt': ['tte','ttm','ttt'],
    }

    hppChannels = ['ee','em','et','mm','mt','tt']
    hpChannels = ['e','m','t']

    recoChannels = {
        'ee': ['eee','eem'],
        'em': ['eme','emm','mee','mem'],
        'mm': ['mme','mmm'],
        'et': ['eee','eme','eem','emm','mee','mem'],
        'mt': ['mee','mem','mme','mmm','eme','emm'],
        'tt': ['eee','eem','eme','emm','mee','mem','mme','mmm'],
    }

    customFinalStates = {
        'Hpp3l' : {},
        'Hpp4l' : {},
    }

    numTauCuts = {
        'Hpp3l' : {},
        'Hpp4l' : {},
    }

    for nt in range(3):
        theCut = '(' + ' || '.join(['genChannel=="%s"' %gChan for hChan in genChannels for gChan in genChannels[hChan] if hChan.count('t')==nt]) + ')'
        numTauCuts['Hpp3l'][nt] = theCut

    for c in genChannels:
        genCut = '(' + ' | '.join(['genChannel=="%s"'%x for x in genChannels[c] + ['aaa']]) + ')'
        recoCut = '(' + ' | '.join(['channel=="%s"'%x for x in recoChannels[c]]) + ')'
        customFinalStates['Hpp3l'][c] = genCut + ' && ' + recoCut

    # plotting correlation
    plotter = CorrelationPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName='plots_correlation',loglevel=loglevel)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])

    logger.info('%s:%s:%iTeV: Plotting correlation' % (analysis, channel, runPeriod))
    plotMethod = getattr(plotter,'plotCorrelation')
    plotMethod(cutFlowMap[channel]['cuts'][1:], 'correlation/mc', cut=myCut, labels=cutFlowMap[channel]['labels'][1:], plottype='mc')
    if useSignal: plotMethod(cutFlowMap[channel]['cuts'][1:], 'correlation/sig', cut=myCut, labels=cutFlowMap[channel]['labels'][1:], plottype='sig')

    # do variables on same plot
    if useSignal:
        logger.info("%s:%s:%iTeV: Plotting signal" % (analysis, channel, runPeriod))
        plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_signal',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        masses = _3L_MASSES if nl==3 else _4L_MASSES
        plotter.initializeSignalSamples([sigMap[runPeriod][x] for x in masses])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotSignal'
        plotMethod = getattr(plotter,plotMode)
        plotMethod('h1.mass',       [1000,0,1000],'signal/hppMass',          yaxis='A.U.',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     lumitext=33,cut=myCut,logy=0, normalize=1)
        plotMethod('h1.dPhi',       [100,0,5],    'signal/hppDphi',          yaxis='A.U.',xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',lumitext=33,logy=0,cut=myCut,normalize=1)
        plotMethod('h1.dR',         [100,0,6.28], 'signal/hppDR',            yaxis='A.U.',xaxis='\\Delta R_{\\ell^{\\pm}\\ell^{\\pm}}',         lumitext=33,logy=0,cut=myCut,normalize=1)
        plotMethod('finalstate.sT', [500,0,2000], 'signal/sT',               yaxis='A.U.',xaxis='S_{T} (GeV/c^{2})',                            lumitext=33,logy=0,cut=myCut,overflow=True,normalize=1)
        plotMethod('z1.mass',       [250,0,1000], 'signal/z1Mass_fullWindow',yaxis='A.U.',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             legendpos=43,logy=0,cut=myCut,overflow=True,normalize=1)
        plotMethod('finalstate.met',[200,0,1000], 'signal/met',              yaxis='A.U.',xaxis='E_{T}^{miss} (GeV/c^{2})',                     lumitext=33,logy=0,cut=myCut,overflow=True,normalize=1)
        for nt in range(3):
            if analysis not in ['Hpp3l']: continue
            theCut = numTauCuts['Hpp3l'][nt] + ' && ' + myCut
            plotMethod('h1.mass',       [1000,0,1000],'signal/hppMass_%iTau'%nt,          yaxis='A.U.',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     lumitext=33,cut=theCut,logy=0, normalize=1)
            plotMethod('h1.dPhi',       [100,0,5],    'signal/hppDphi_%iTau'%nt,          yaxis='A.U.',xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',lumitext=33,logy=0,cut=theCut,normalize=1)
            plotMethod('h1.dR',         [100,0,6.28], 'signal/hppDR_%iTau'%nt,            yaxis='A.U.',xaxis='\\Delta R_{\\ell^{\\pm}\\ell^{\\pm}}',         lumitext=33,logy=0,cut=theCut,normalize=1)
            plotMethod('finalstate.sT', [500,0,2000], 'signal/sT_%iTau'%nt,               yaxis='A.U.',xaxis='S_{T} (GeV/c^{2})',                            lumitext=33,logy=0,cut=theCut,overflow=True,normalize=1)
            plotMethod('z1.mass',       [250,0,1000], 'signal/z1Mass_fullWindow_%iTau'%nt,yaxis='A.U.',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             legendpos=43,logy=0,cut=theCut,overflow=True,normalize=1)
            plotMethod('finalstate.met',[200,0,1000], 'signal/met_%iTau'%nt,              yaxis='A.U.',xaxis='E_{T}^{miss} (GeV/c^{2})',                     lumitext=33,logy=0,cut=theCut,overflow=True,normalize=1)
        #plotter.initializeSignalSamples([allSigMap[runPeriod][mass]])
        #plotter.setIntLumi(intLumiMap[runPeriod])
        #plotMode = 'plotSignal'
        #plotMethod = getattr(plotter,plotMode)
        #plotMethod('h1.mass',[400,300,700],'signal/hppMass_cut',yaxis='Events/0.1 GeV/c^{2}',xaxis='M_{\\ell^{+}\\ell^{+}} (GeV/c^{2})',lumitext=33,logy=0,cut=myCut,boxes=[[450,550,2]])
        names = {
            'e': 'Elec',
            'm': 'Mu',
            't': 'Tau',
        }
        tex = {
            'e': 'e',
            'm': '\\mu',
            't': '\\tau',
        }
        for l in ['e','m']:
            name = names[l]
            t = tex[l]
            cuts = ['%s & %s' %(myCut,'l%iFlv=="%s"' %((x+1),l)) for x in range(nl)]
            plotMethod(['l%i.Pt'  %(x+1) for x in range(nl)], [100,0,1000], 'signal/%sPt'%name, yaxis='A.U.', xaxis='p_{T}^{%s} (GeV)' %t, legendpos=43, logy=0, cut=cuts, overflow=True, normalize=1)
            for nt in range(3):
                if analysis not in ['Hpp3l']: continue
                theCuts = [numTauCuts['Hpp3l'][nt] + ' && ' + c for c in cuts]
                plotMethod(['l%i.Pt'  %(x+1) for x in range(nl)], [100,0,1000], 'signal/%sPt_%iTau'%(name,nt), yaxis='A.U.', xaxis='p_{T}^{%s} (GeV)' %t, legendpos=43, logy=0, cut=theCuts, overflow=True, normalize=1)




    # Plotting discriminating variables
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,rootName='plots_2d',loglevel=loglevel)
    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])


    if analysis in ['Hpp3l'] and not blind:
        plotter.plotMCDataSignalRatio2D('z1.mass',       'finalstate.mass',[70,0,140], [100,0,200],'m3l_v_z1_mc',     xaxis='M_{\\ell^{+}\\ell^{-}}',                 yaxis='M_{3\\ell}',        cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('z1.mass',       'finalstate.mass',[70,0,140], [100,0,200],'m3l_v_z1_data',   xaxis='M_{\\ell^{+}\\ell^{-}}',                 yaxis='M_{3\\ell}',        cut=myCut, plotdata=1, plotmc=0, plotsig=0)
    if analysis in ['Hpp3l']:
        #plotter.setScaleFactor('event.lep_scale*event.trig_scale')
        plotter.plotMCDataSignalRatio2D('h1.mass',       'event.nvtx',[24,0,600], [50,0,50],'hppMass_v_pu_sig',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('z1.mass',       'event.nvtx',[80,0,240], [50,0,50],'z1Mass_v_pu_sig', xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('h1.dPhi',       'event.nvtx',[32,0,3.2], [50,0,50],'hppDphi_v_pu_sig',xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('finalstate.sT', 'event.nvtx',[40,0,1000],[50,0,50],'sT_v_pu_sig',     xaxis='S_{T} (GeV/c^{2})',                            yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('finalstate.met','event.nvtx',[40,0,200], [50,0,50],'met_v_pu_sig',    xaxis='E_{T}^{miss} (GeV)',                           yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=0, plotsig=1)
        plotter.plotMCDataSignalRatio2D('h1.mass',       'event.nvtx',[24,0,600], [50,0,50],'hppMass_v_pu_mc', xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',     yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('z1.mass',       'event.nvtx',[80,0,240], [50,0,50],'z1Mass_v_pu_mc',  xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',             yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('h1.dPhi',       'event.nvtx',[32,0,3.2], [50,0,50],'hppDphi_v_pu_mc', xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('finalstate.sT', 'event.nvtx',[40,0,1000],[50,0,50],'sT_v_pu_mc',      xaxis='S_{T} (GeV/c^{2})',                            yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)
        plotter.plotMCDataSignalRatio2D('finalstate.met','event.nvtx',[40,0,200], [50,0,50],'met_v_pu_mc',     xaxis='E_{T}^{miss} (GeV)',                           yaxis='Number PU Vertices',cut=myCut, plotdata=0, plotmc=1, plotsig=0)

    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    if useSignal:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]+['Sig']])
    else:
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    plotter.setIntLumi(intLumiMap[runPeriod])

    plotMode = 'plotMCDataRatio' if dataplot else 'plotMC'
    plotMethod = getattr(plotter,plotMode)
    logger.info("%s:%s:%iTeV: Plotting discriminating variables" % (analysis,channel, runPeriod))
    plotDistributions(plotMethod,myCut,nl,isControl,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)

    # each channel
    if plotFinalStates:
        logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
        for c in fsToPlot:
            logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
            plotDistributions(plotMethod,myCut+'&&channel=="%s"'%c,nl,isControl,savedir=c,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
               sel = customFinalStates[analysis][c]
               logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
               plotDistributions(plotMethod,'%s & %s'%(myCut,sel),nl,isControl,savedir=c,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)

    # some partially blind plots for h++
    if runPeriod==8 and not dataplot and analysis in ['Hpp3l']:
        logger.info("%s:%s:%iTeV: Plotting partially blinded variables" % (analysis, channel, runPeriod))
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotModeUnblind = 'plotMCDataRatio'
        plotMethodUnblind = getattr(plotter,plotModeUnblind)
        plotMethodUnblind('h1.mass',      [24,0,600], 'hppMass_unblind',          yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',lumitext=33,logy=1,cut=myCut,overflow=True,blinder=[150,99999])
        plotMethodUnblind('finalstate.sT',[50,0,500], 'sT_unblind',               yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',                       lumitext=33,logy=0,cut=myCut,overflow=True,blinder=[200,99999])
        plotMethodUnblind('z1.mass',      [42,70,112],'z1Mass_unblind',           yaxis='Events/1.0 GeV',       xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',        legendpos=43,logy=0,cut=myCut,blinder=[112,99999])
        plotMethodUnblind('z1.mass',      [80,0,240], 'z1Mass_fullWindow_unblind',yaxis='Events/3.0 GeV',       xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',        legendpos=43,logy=0,cut=myCut,blinder=[112,99999])
        if plotFinalStates:
            logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
            for c in fsToPlot:
                logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                plotMethodUnblind('h1.mass',      [24,0,600], c+'/hppMass_unblind',          yaxis='Events/25.0 GeV/c^{2}',xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}} (GeV/c^{2})',lumitext=33,logy=1,cut=myCut+'&&channel=="%s"'%c,overflow=True,blinder=[150,99999])
                plotMethodUnblind('finalstate.sT',[50,0,500], c+'/sT_unblind',               yaxis='Events/10.0 GeV/c^{2}',xaxis='S_{T} (GeV/c^{2})',                       lumitext=33,logy=0,cut=myCut+'&&channel=="%s"'%c,overflow=True,blinder=[200,99999])
                plotMethodUnblind('z1.mass',      [42,70,112],c+'/z1Mass_unblind',           yaxis='Events/1.0 GeV',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',               legendpos=43,logy=0,cut=myCut+'&&channel=="%s"'%c,blinder=[112,99999])
                plotMethodUnblind('z1.mass',      [80,0,240], c+'/z1Mass_fullWindow_unblind',yaxis='Events/3.0 GeV',xaxis='M_{\\ell^{+}\\ell^{-}} (Z) (GeV)',               legendpos=43,logy=0,cut=myCut+'&&channel=="%s"'%c,blinder=[112,99999])

    # setup signal overlay plots
    if useSignal:
        plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_overlay',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotMCDataSignalRatio' if dataplot else 'plotMCSignalRatio'
        plotMethod = getattr(plotter,plotMode)
        if analysis in ['Hpp3l']:
            plotter.plotMCDataSignalRatio2D('h1.mass','h1.dPhi', [80,0,800], [32,0,3.2], 'h1mass_v_h1dphi_mc',  plotdata=0, plotmc=1, plotsig=0, cut=myCut, xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}}',yaxis='#Delta#phi_{\\ell^{\\pm}\\ell^{\\pm}}')
            plotter.plotMCDataSignalRatio2D('h1.mass','h1.dPhi', [80,0,800], [32,0,3.2], 'h1mass_v_h1dphi_sig', plotdata=0, plotmc=0, plotsig=1, cut=myCut, xaxis='M_{\\ell^{\\pm}\\ell^{\\pm}}',yaxis='#Delta#phi_{\\ell^{\\pm}\\ell^{\\pm}}')
    if plotOverlay and useSignal:
        # plot the signal overlay
        logger.info("%s:%s:%iTeV: Plotting signal overlay discriminating variables" % (analysis, channel, runPeriod))
        plotDistributions(plotMethod,myCut,nl,isControl,savedir='overlay',signalscale=100,analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)


    # plot shapes
    if plotShapes:
        logger.info("%s:%s:%iTeV: Plotting shapes" % (analysis, channel, runPeriod))
        plotter = ShapePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_shapes',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        if useSignal: plotter.initializeSignalSamples([sigMap[runPeriod]['Sig']])
        if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotter.plotMC('z1.mass',['channel=="mmm"','channel=="emm"'],[42,70,112],'zMass_mc_mm',yaxis='Normalized',xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,cutNames=['mmm','emm'])
        if dataplot: plotter.plotData('z1.mass',['channel=="mmm"','channel=="emm"'],[42,70,112],'zMass_data_mm',yaxis='Normalized',xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,cutNames=['mmm','emm'])
        plotter.plotMC('z1.mass',['channel=="eee"','channel=="eme"'],[42,70,112],'zMass_mc_ee',yaxis='Normalized',xaxis='M(l^{+}l^{-}) (GeV)',logy=0,cut=myCut,cutNames=['eee','eme'])
        if dataplot: plotter.plotData('z1.mass',['channel=="eee"','channel=="eme"'],[42,70,112],'zMass_data_ee',yaxis='Normalized',xaxis='M_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,cutNames=['eee','eme'])

    # plot cut flows (each cut)
    logger.info("%s:%s:%iTeV: Plotting cut flow" % (analysis, channel, runPeriod))
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutFlowSelections',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
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
            logger.info('%s:%s:%iTeV: Plotting cut flow selections %s' % (analysis, channel, runPeriod, cutFlowMap[channel]['labels_simple'][i]))
            thisCut = '&&'.join(cutFlowMap[channel]['cuts'][:i+1])
            plotDistributions(plotMethod,'%s&%s'%(myCut,thisCut),nl,isControl,savedir='cutflow/%s'%cutFlowMap[channel]['labels_simple'][i],analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)
            if plotFinalStates:
                logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
                for c in fsToPlot:
                    logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                    plotDistributions(plotMethod,'%s&channel=="%s"&%s'%(myCut,c,thisCut),nl,isControl,savedir='cutflow/%s/%s' %(cutFlowMap[channel]['labels_simple'][i],c),analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)
            thisCut = cutFlowMap[channel]['cuts'][i]
            plotDistributions(plotMethod,'%s&%s'%(myCut,thisCut),nl,isControl,savedir='cutflow/%s_only'%cutFlowMap[channel]['labels_simple'][i],analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)
            if plotFinalStates:
                logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
                for c in fsToPlot:
                    logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                    plotDistributions(plotMethod,'%s&channel=="%s"&%s'%(myCut,c,thisCut),nl,isControl,savedir='cutflow/%s_only/%s' %(cutFlowMap[channel]['labels_simple'][i],c),analysis=analysis,region=channel,nostack=nostack,normalize=normalize,mass=mass)

    # plot cut flows on same plot
    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflow',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
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
            logger.info("%s:%s:%iTeV: Plotting cut flow  %s" % (analysis, channel, runPeriod, c))
            plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
               sel = customFinalStates[analysis][c]
               logger.info("%s:%s:%iTeV: Plotting cut flow  %s" % (analysis, channel, runPeriod, c))
               plotMethod(['%s&& %s &&%s' %(x,sel,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=1)


    # setup individual channel cuts on same plot
    plotChannelStrings, plotChannelCuts = getChannelStringsCuts(channel,finalStates)
    plotMode = 'plotCutFlowMCData' if dataplot else 'plotCutFlowMC'
    if useSignal: plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
    plotMethod = getattr(plotter,plotMode)
    plotMethod([myCut]+['%s&&%s' %(x,myCut) for x in plotChannelCuts],'individualChannels',labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0,signalscale=1000)
    if plotFinalStates:
        logger.info("%s:%s:%iTeV: Plotting individual finalStates" % (analysis, channel, runPeriod))
        if analysis in customFinalStates:
            plotGenChannelStrings, plotGenChannelCuts = getGenChannelStringsCuts(channel,customFinalStates[analysis])
            for c in fsToPlot:
                logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
                plotMethod(['%s&&channel=="%s"'%(myCut,c)]+['%s&&%s&&channel=="%s"' %(x,myCut,c) for x in plotGenChannelCuts],'%s/individualGenChannels'%c,labels=['Total']+plotGenChannelStrings,nosum=True,lumitext=33,logy=0,signalscale=1000)
        if analysis in customFinalStates:
            for c in customFinalStates[analysis]:
               sel = customFinalStates[analysis][c]
               logger.info("%s:%s:%iTeV: Channel %s" % (analysis, channel, runPeriod, c))
               plotMethod(['%s&&%s' %(myCut,sel)]+['%s&&%s&&%s' %(x,myCut,sel) for x in plotChannelCuts],'%s/individualChannels'%c,labels=['Total']+plotChannelStrings,nosum=True,lumitext=33,logy=0,signalscale=1000)


    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflowSelectionsChannels',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
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
    #if analysis in ['WZ']:
    #    logger.info("%s:%s:%iTeV: Plotting efficiency" % (analysis, channel, runPeriod))
    #    plotter = EfficiencyPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_efficiency',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    #    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel] if x != 'WZ'])
    #    plotter.initializeSignalSamples([sigMap[runPeriod]['WZ']])
    #    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    #    plotter.setIntLumi(intLumiMap[runPeriod])
    #    plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
    #    plotMethod = getattr(plotter,plotMode)
    #    plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'efficiency',labels=cutFlowMap[channel]['labels'],lumitext=33)
    #    if plotFinalStates:
    #        for c in fsToPlot:
    #            logger.info("%s:%s:%iTeV: Plotting efficiency  %s" % (analysis, channel, runPeriod, c))
    #            plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/efficiency'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)

    #    # plot cut flows overlays
    #    plotter = CutFlowPlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='plots_cutflow_overlay',mergeDict=mergeDict,scaleFactor=scaleFactor,loglevel=loglevel)
    #    plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel] if x not in ['WZ']])
    #    plotter.initializeSignalSamples([sigMap[runPeriod]['WZ']])
    #    if dataplot: plotter.initializeDataSamples([sigMap[runPeriod]['data']])
    #    plotter.setIntLumi(intLumiMap[runPeriod])
    #    plotMode = 'plotCutFlowMCDataSignal' if dataplot else 'plotCutFlowMCSignal'
    #    plotMethod = getattr(plotter,plotMode)
    #    plotMethod([x+'&&'+myCut for x in cutFlowMap[channel]['cuts']],'cutFlow_overlay',labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)
    #    if plotFinalStates:
    #        for c in fsToPlot:
    #            logger.info("%s:%s:%iTeV: Plotting cut flow overlay  %s" % (analysis, channel, runPeriod, c))
    #            plotMethod(['%s&&channel=="%s"&&%s' %(x,c,myCut) for x in cutFlowMap[channel]['cuts']],'%s/cutFlow_overlay'%c,labels=cutFlowMap[channel]['labels'],lumitext=33,logy=0)


def plotFakeRate(analysis,channel,runPeriod,**kwargs):
    '''Plot fake rate for an analysis.'''
    logger = logging.getLogger(__name__)
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
        logger.warning("Unrecognized parameter '" + key + "' = " + str(value))
        return 0

    if useSignal: logger.info("%s:%s:%iTeV: Mass: %i" % (analysis,channel,runPeriod,mass))
    isControl = analysis != channel
    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ_W' : 2,
        'WZ_FakeRate' : 1,
        'WZ'   : 3,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[analysis]
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,runPeriod,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,runPeriod)
    sigMap = getSigMap(nl,mass)
    intLumiMap = getIntLumiMap()
    channelBackground = getChannelBackgrounds(runPeriod)

    finalStates, leptons = getChannels(nl,runTau=runTau)
    if finalStatesToPlot=='all':
        fsToPlot = finalStates
    else:
        fsToPlot = finalStatesToPlot.split(',')
    logger.info('%s:%s:%iTeV: Cuts to be applied: %s' % (analysis, channel, runPeriod, myCut))
    dataplot = (isControl or not blind)
    mergeDict = getMergeDict(runPeriod)

    # TODO: load from pickle file
    # define fake regions
    lepName = {'e': 'Elec', 'm': 'Muon', 't': 'Tau'}
    fakeRegions = {}
    fakeRegions['WZ'] = {}
    for f in ['e', 'm']:
        for p in ['Loose', 'Tight']:
            # select leading Z pt, Z window [60,120], tight (or loose) Z, low met, m3l>100, w1 mass < 30
            if analysis in ['WZ']:
                for z in ['Loose', 'Tight']:
                    fakeRegion = 'Z{0}Probe{1}{2}'.format(z,lepName[f],p)
                    denom = 'z1.Pt1>20. & z1.mass>60. & z1.mass<120. & z1.Pass{0}1 & z1.Pass{0}2 & finalstate.met<20. & finalstate.mass>100. & w1.mass<20. & w1.dR1_z1_1>0.1 & w1.dR1_z1_2>0.1 & w1Flv=="{1}"'.format(z,f)
                    numer = '{0} & w1.Pass{1}1'.format(denom,p)
                    fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
                    if p=='Tight':
                       fakeRegion += '_LooseProbe'
                       denom += ' & w1.PassLoose1'
                       numer += ' & w1.PassLoose1'
                       fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
            # select w lepton pt, z veto, met
            #'W' : 'w1.Pt1>20. & (z1.mass<60. | z1.mass>120.) & finalstate.met>30. & w1.mass>30.',
            if analysis in ['WZ_W']:
                for w in ['Loose','Tight']:
                    fakeRegion = 'W{0}Probe{1}{2}'.format(w,lepName[f],p)
                    denom = 'w1.Pt1>20. & w1.mass>30. & finalstate.met>30. & (z1.mass<60. | z1.mass>120.) & l1.Chg==l2.Chg & z1.dR>0.1 & w1.Pass{0}1 & w2Flv=="{1}"'.format(w,f)
                    numer = '{0} & w2.Pass{1}1'.format(denom,p)
                    fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
                    if p=='Tight':
                       fakeRegion += '_LooseProbe'
                       denom += ' & w2.PassLoose1'
                       numer += ' & w2.PassLoose1'
                       fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
            # ntuple cuts: zVeto 60-120, met vet 20, w veto 20, jet pt > 20, jet dr > 1.0
            if analysis in ['WZ_FakeRate']:
               fakeRegion = 'FakeRateProbe{0}{1}'.format(lepName[f],p)
               denom = 'l1Flv=="{0}"'.format(f)
               numer = '{0} & w1.Pass{1}1'.format(denom,p)
               fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
               if p=='Tight':
                  fakeRegion += '_LooseProbe'
                  denom += ' & w1.PassLoose1'
                  numer += ' & w1.PassLoose1'
                  fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}



    # setup selections
    ptBins = [10,20,40,100]
    etaBins = {
        'e': [0,1.479,2.5],
        'm': [0,1.2,2.4],
    }
    for fakeRegion in fakeRegions['WZ']:
        logger.info("%s:%s:%iTeV: Fake Region: %s" % (analysis,channel, runPeriod, fakeRegion))
        denom = fakeRegions['WZ'][fakeRegion]['denom']
        numer = fakeRegions['WZ'][fakeRegion]['numer']
        probe = fakeRegions['WZ'][fakeRegion]['probe']
        ptvar = fakeRegions['WZ'][fakeRegion]['ptVar']
        etavar = fakeRegions['WZ'][fakeRegion]['etaVar']

        if 'Muon' in fakeRegion: # prescale 23 workaround
            plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,rootName='{0}_fakeplots'.format(fakeRegion),scaleFactor='event.pu_weight*event.lep_scale*event.trig_scale*1./23.')
        elif 'Elec' in fakeRegion: # prescale 10 workaround
            plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,rootName='{0}_fakeplots'.format(fakeRegion),scaleFactor='event.pu_weight*event.lep_scale*event.trig_scale*1./23.')
        else:
            plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,mergeDict=mergeDict,rootName='{0}_fakeplots'.format(fakeRegion),scaleFactor='event.pu_weight*event.lep_scale*event.trig_scale')
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotMode = 'plotMCData'
        plotMethod = getattr(plotter,plotMode)
        logger.info("%s:%s:%iTeV: Plotting discriminating variables: All Probes" % (analysis,channel, runPeriod))
        plotDistributions(plotMethod,denom,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_all'.format(fakeRegion))
        logger.info("%s:%s:%iTeV: Plotting discriminating variables: Passing" % (analysis,channel, runPeriod))
        plotDistributions(plotMethod,numer,nl,isControl,analysis=analysis,savedir='fakeRate/{0}_prompts'.format(fakeRegion))
        logger.info("%s:%s:%iTeV: Plotting discriminating variables: Failing" % (analysis,channel, runPeriod))
        plotDistributions(plotMethod,'{0} & !({1})'.format(denom,numer),nl,isControl,analysis=analysis,savedir='fakeRate/{0}_fakes'.format(fakeRegion))

        # now plot the fake rates
        logger.info("%s:%s:%iTeV: Computing fake rates" % (analysis,channel, runPeriod))
        plotter = FakeRatePlotter(channel,ntupleDir=ntuples,saveDir=saves,period=runPeriod,rootName='{0}_fakerates'.format(fakeRegion),mergeDict=mergeDict,scaleFactor='event.pu_weight*event.lep_scale*event.trig_scale*1./23.')
        # this should be done in data... using mc until we get some!
        # subtract WZ, ZZ, ttbar contributions from data
        # only initialize Z... in MC?
        plotter.initializeBackgroundSamples([sigMap[runPeriod][x] for x in channelBackground[channel]])
        plotter.initializeDataSamples([sigMap[runPeriod]['data']])
        plotter.setIntLumi(intLumiMap[runPeriod])
        plotter.plotFakeRate(numer, denom, 'fakeRate/{0}_fakerate'.format(fakeRegion), ptBins=ptBins, etaBins=etaBins[probe], logx=1, ptVar=ptvar, etaVar=etavar)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['Z','WZ','WZ_W','WZ_FakeRate','Hpp2l','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['Z','WZ','W','FakeRate','TT','Hpp2l','Hpp3l','Hpp4l','FakeRate','LowMass'], help='Channel in analysis')
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
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

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
        logger.warning("7 TeV not implemented")
    elif args.doFakeRate:
        plotFakeRate(args.analysis,args.channel,args.period,mass=args.mass,loglevel=args.log)
    else:
        plotRegion(args.analysis,args.channel,args.period,plotFinalStates=args.plotFinalStates,runTau=args.runTau,blind=args.unblind,mass=args.mass,plotJetBins=args.plotJetBins,plotOverlay=args.plotOverlay,plotShapes=args.plotShapes,plotCutFlow=args.plotCutFlow,myCut=args.cut,finalStates=args.finalStates,nostack=args.nostack,normalize=args.normalize,scaleFactor=args.scaleFactor,loglevel=args.log)

    return 0


if __name__ == "__main__":
    main()
