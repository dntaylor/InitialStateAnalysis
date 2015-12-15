#!/usr/bin/env python

import argparse
import itertools
import sys
import os
import logging
import json
import pickle

from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS
from InitialStateAnalysis.Limits.limitUtils import getSignalStrength

loglevel = 'INFO'
logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

##########################
### Load the pkl files ###
##########################
with open('yields.pkl','rb') as f:
    yields = pickle.load(f)

with open('systematics.pkl','rb') as f:
    systematics = pickle.load(f)


fakesFile = 'InitialStateAnalysis/Analyzers/python/scale_factors/fakes_trigIso_dijet_13TeV.json'
#fakesFile = 'InitialStateAnalysis/Analyzers/python/scale_factors/fakes_trigIso_13TeV.json'
fakerates = {}
with open(fakesFile,'r') as f:
    faketemp = json.load(f)
pts = {}
for f in faketemp:
    key = 'm' if 'Muon' in f else 'e'
    fakerates[key] = {}
    pts = []
    for v in faketemp[f]:
        ptlow = v['pt_low']
        pthigh = v['pt_high']
        etalow = v['eta_low']
        etahigh = v['eta_high']
        fake = v['fakerate']
        err = v['error']
        if ptlow not in fakerates[key]:
            fakerates[key][ptlow] = {}
            pts += [ptlow]
        fakerates[key][ptlow]['endcap' if etalow else 'barrel'] = {'pthigh' : pthigh, 'fake': fake, 'err': err}


cut = 'select.passTight'
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 5 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 10 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 15 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 20 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 25 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 30 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS

sigStrengths = getSignalStrength(cut)
sigStrengthsExperimental = getSignalStrength(cut,mode='experimental')
sigStrengthsLumi = getSignalStrength(cut,mode='lumi')
sigStrengthsTheory = getSignalStrength(cut,mode='theory')
#sigStrengthsStatOnly = getSignalStrength(cut,statonly=True)
sigStrengthsStatOnly = getSignalStrength(cut,mode='stat')
print 'all         ', sigStrengths
print 'experimental', sigStrengthsExperimental
print 'lumi        ', sigStrengthsLumi
print 'theory      ', sigStrengthsTheory
print 'stat        ', sigStrengthsStatOnly
sigs = {}
for chan in sigStrengthsStatOnly:
    sigs[chan] = {}
    sigs[chan]['val'] = sigStrengths[chan][0]
    sigs[chan]['statdown'] = abs(sigStrengthsStatOnly[chan][1])
    sigs[chan]['statup'] = abs(sigStrengthsStatOnly[chan][2])
    sigs[chan]['stat'] = (sigs[chan]['statdown'] + sigs[chan]['statup'])/2
    sigs[chan]['systdown'] = abs(sigStrengthsExperimental[chan][1]**2 - sigStrengthsStatOnly[chan][1]**2)**0.5
    sigs[chan]['systup'] = abs(sigStrengthsExperimental[chan][2]**2 - sigStrengthsStatOnly[chan][2]**2)**0.5
    sigs[chan]['syst'] = (sigs[chan]['systdown'] + sigs[chan]['systup'])/2
    sigs[chan]['theodown'] = abs(sigStrengthsTheory[chan][1]**2 - sigStrengthsStatOnly[chan][1]**2)**0.5
    sigs[chan]['theoup'] = abs(sigStrengthsTheory[chan][2]**2 - sigStrengthsStatOnly[chan][2]**2)**0.5
    sigs[chan]['theo'] = (sigs[chan]['theodown'] + sigs[chan]['theoup'])/2
    sigs[chan]['lumidown'] = abs(sigStrengthsLumi[chan][1]**2 - sigStrengthsStatOnly[chan][1]**2)**0.5
    sigs[chan]['lumiup'] = abs(sigStrengthsLumi[chan][2]**2 - sigStrengthsStatOnly[chan][2]**2)**0.5
    sigs[chan]['lumi'] = (sigs[chan]['lumidown'] + sigs[chan]['lumiup'])/2

############################
### Templates for tables ###
############################
names = {
    'eee' : '\\eee',
    'eem' : '\\eem',
    'mme' : '\\mme',
    'mmm' : '\\mmm',
}

tables = {}

# fake rate values
tables['fakerate'] = {}
tables['fakerate']['header'] = '\\begin{table}[htbp]\n'\
                               '    \\centering\n'\
                               '    \\begin{tabular}{|c|cc|cc|}\n'\
                               '    \\hline\n'\
                               '                    & \\multicolumn{2}{c|}{Electron}           & \\multicolumn{2}{c|}{Muon}               \\\\\n'\
                               '    $\\pt$ bin       & Barrrel           & Endcap            & Barrrel           & Endcap            \\\\ \\hline\n'
tables['fakerate']['row'] = '$[{ptlo:4d},{pthi:4d}] \\GeV$  & ${eb:5.3f} \\pm {ebe:5.3f}$ & ${ee:5.3f} \\pm {eee:5.3f}$ & ${mb:5.3f} \\pm {mbe:5.3f}$ & ${me:5.3f} \\pm {mee:5.3f}$ \\\\'
tables['fakerate']['rows'] = []
tables['fakerate']['footer'] = '    \\hline\n'\
                               '    \\end{tabular}\n'\
                               '    \\caption{Fake rate values for electrons and muons in the barrel and endcap calculated in $\\cPZ$ + fake control region.}\n'\
                               '    \\label{tab:fakerates}\n'\
                               '\\end{table}\n'

# fake rate control yields
tables['controlyields'] = {}
tables['controlyields']['header'] = '\\begin{table}[htbp]\n'\
                                    '    \\centering\n'\
                                    '    \\begin{tabular}{|ll|cccc| }\n'\
                                    '\\hline %------------------------------------------------------------------------------------------ \n'\
                                    '      &     & \\multicolumn{4}{c|}{Yields in Control Regions}                                                                                                      \\\\ \n'\
                                    '\\multicolumn{2}{c|}{Control Region}   & $\\eee$    & $\\eem$   & $\\mme$  & $\\mmm$ \\\\ \n'\
                                    '\hline %------------------------------------------------------------------------------------------ \n'
tables['controlyields']['datarow'] = '{fr:4}    & Data & ${eeed:3d} \\pm {eeede:5.2f}$ & ${eemd:3d} \\pm {eemde:5.2f}$ & ${mmed:3d} \\pm {mmede:5.2f}$ & ${mmmd:3d} \\pm {mmmde:5.2f}$ \\\\'
tables['controlyields']['mcrow'] = '        & MC   & ${eeem:5.2f} \\pm {eeeme:5.2f}$ & ${eemm:5.2f} \\pm {eemme:5.2f}$ & ${mmem:5.2f} \\pm {mmeme:5.2f}$ & ${mmmm:5.2f} \\pm {mmmme:5.2f}$ \\\\ \\hline'
tables['controlyields']['rows'] = []
tables['controlyields']['footer'] = '\\hline %------------------------------------------------------------------------------------------ \n'\
                                    '     \\end{tabular}\n'\
                                    '    \\caption{ Control region yields in data and Monte Carlo. The Monte Carlo includes contributions from $\\WZ$, $\\ZZ$, $\\TTV$, and $\\VVV$ \n'\
                                    '              that must be subtracted from the data yields. }\n'\
                                    '    \\label{tab:control_regions3}\n'\
                                    '\\end{table}'


# fake rate control regions
tables['controlregions'] = {}
tables['controlregions']['header'] = '\\begin{table}[htbp]\n'\
                                     '    \\centering\n'\
                                     '    \\begin{tabular}{|l|cccc| }\n'\
                                     '\\hline %------------------------------------------------------------------------------------------ \n'\
                                     'Control &  \\multicolumn{4}{c|}{Estimated background in the signal region}                        \\\\\n'\
                                     'region  &  $\\eee$      & $\\eem$      & $\\mme$     & $\\mmm$            \\\\\n'\
                                     '\\hline %------------------------------------------------------------------------------------------ \n'
tables['controlregions']['row'] = '{fr:4}     & ${eee:5.2f} \\pm {eeee:4.2f}$ & ${eem:5.2f} \\pm {eeme:4.2f}$ & ${mme:5.2f} \\pm {mmee:4.2f}$ & ${mmm:5.2f} \\pm {mmme:4.2f}$ \\\\'
tables['controlregions']['rows'] = []
tables['controlregions']['footer'] = '\\hline %------------------------------------------------------------------------------------------ \n'\
                                     '     \\end{tabular}\n'\
                                     '    \\caption{ Estimated backgound contribution from different control regions.}\n'\
                                     '    \\label{tab:control_regions4}\n'\
                                     '\\end{table}\n'


# final cross section yields
tables['results'] = {}
tables['results']['header'] = '\\begin{table}[htbp]\n'\
                              '\\centering\n'\
                              '\\begin{tabular}{|l|c|cc|c|c|}\n'\
                              '\\hline\n'\
                              'Decay           & $N_{WZ}^{exp}$                &  Background                              & Total                        & Observed \\\\\n'\
                              'channel         &                               & Datadriven            & Monte Carlo      & expected                     &          \\\\\n'\
                              '\\hline\n'\
                              '\\hline\n'
#tables['results']['row'] = '${name:12} $ & $ {wz:5.2f} \\pm {wze:4.2f} \\pm {wzs:4.2f}$   & $ {bg:5.2f} \\pm {bge:4.2f} \\pm {bgs:4.2f}$ & $ {bgm:5.2f} \\pm {bgme:4.2f} \\pm {bgms:4.2f}$  & $ {ex:5.2f} \\pm {exe:4.2f} \\pm {exs:4.2f}$  & $ {obs:3d} $  \\\\'
# no syst
tables['results']['row'] = '${name:12} $ & $ {wz:5.2f} \\pm {wze:4.2f}$   & $ {bg:5.2f} \\pm {bge:4.2f}$ & $ {bgm:5.2f} \\pm {bgme:4.2f}$  & $ {ex:5.2f} \\pm {exe:4.2f}$  & $ {obs:3d} $  \\\\'
tables['results']['rows'] = []
tables['results']['footer'] = '\\hline\n'\
                              '\\end{tabular}\n'\
                              '\\vspace{0.5cm}\n'\
                              '\\caption{ The observed and expected yield of WZ events, and estimated yield of background events obtained \n'\
                              'from data are shown for each decay channel and are summed in the total expected yield.}\n'\
                              '\\label{table:results}\n'\
                              '\\end{table}\n'

# systematics TODO!!!

# signal strength
tables['crosssection'] = {}
tables['crosssection']['header'] = '\\begin{table}[htbp]\n'\
                                   '\\centering \n'\
                                   '\\begin{tabular}{|l|c|} \n'\
                                   '\\hline \n'\
                                   'Decay           & Signal Strength \\\\ \n'\
                                   'channel         &                 \\\\ \n'\
                                   '\\hline \n'
tables['crosssection']['row'] = '$ {name:12} $ & $ {val:4.2f} \\pm {stat:4.2f} \\mbox{{ (stat.) }} ^{{ {systup:+4.2f} }}_{{ {systdown:4.2f} }} \\mbox{{ (syst.) }}  ^{{ {theoup:+4.2f} }}_{{ {theodown:4.2f} }} \\mbox{{ (theory) }} \\pm {lumi:4.2f} \\mbox{{ (lumi.) }}$ \\\\'
tables['crosssection']['rows'] = []
tables['crosssection']['footer'] = '\\hline\n'\
                                   '\\end{tabular}\n'\
                                   '\\vspace{0.5cm}\n'\
                                   '\\caption{ The observed signal strength. }\n'\
                                   '\\label{table:signal_strength}\n'\
                                   '\\end{table}\n'

# table comparing TT and DY

tables['dyttbar'] = {}
tables['dyttbar']['header'] = '\\begin{table}[htbp]\n'\
                              '    \\centering\n'\
                              '    \\begin{tabular}{|l|l|cccc|c| }\n'\
                              '\\hline %------------------------------------------------------------------------------------------ \n'\
                              'Control &        & \\multicolumn{5}{c|}{Estimated background in the signal region}         \\\\\n'\
                              'region  & Sample &  $\\eee$      & $\\eem$      & $\\mme$     & $\\mmm$            & Total \\\\\n'\
                              '\\hline %------------------------------------------------------------------------------------------ \n'
tables['dyttbar']['row'] = '{name:8} & {sample:6} & ${eee:5.2f} \\pm {eeee:5.2f}$ & ${eem:5.2f} \\pm {eeme:5.2f}$ & ${mme:5.2f} \\pm {mmee:5.2f}$ & ${mmm:5.2f} \\pm {mmme:5.2f}$ & ${tot:5.2f} \\pm {tote:5.2f}$ \\\\'
tables['dyttbar']['rows'] = []
tables['dyttbar']['footer'] = '\\hline\n'\
                              '\\end{tabular}\n'\
                              '\\caption{ Drell-Yan and t\\bar{t} yields estimated from Monte Carlo. }\n'\
                              '\\label{table:dyttbar}\n'\
                              '\\end{table}\n'



# final cross section
equations = {}

fiducialTheory = 69.18+68.97+69.77+69.46

equations['fiducialcrosssection'] = {}
equations['fiducialcrosssection']['equation'] = '\\begin{{equation}}\n'\
                                                '\\sigma_{{\\rm{{fid}}}} (\\pp \\to \\PW\\PZ \\to \\ell\\nu\\ell\'\\ell\') = {xsec:4d} \\pm {stat:4d} \\stat ^{{ {systup:+4} }}_{{ {systdown:4d} }} \\syst \\pm {lumi:4d} \\lumi\\unit{{fb}}.\n'\
                                                '\\end{{equation}}\n'

totalTheory = 42.7

equations['totalcrosssection'] = {}
equations['totalcrosssection']['equation'] = '\\begin{{equation}}\n'\
                                             '\\sigma(\\pp \\to \\PW\\PZ) = {xsec:5.1f} \\pm {stat:5.1f} \\stat  ^{{ {systup:+5.1f} }}_{{ {systdown:5.1f} }} \\syst \\pm 0.6 \\theo \\pm {lumi:5.1f} \\lumi\\unit{{pb}}.\n'\
                                             '\\end{{equation}}\n'

########################
### Build the tables ###
########################

# fakerate
for ptlow in pts:
    if ptlow < 10: continue
    pthigh = fakerates['e'][ptlow]['barrel']['pthigh']
    eb  = fakerates['e'][ptlow]['barrel']['fake']
    ebe = fakerates['e'][ptlow]['barrel']['err']
    ee  = fakerates['e'][ptlow]['endcap']['fake']
    eee = fakerates['e'][ptlow]['endcap']['err']
    mb  = fakerates['m'][ptlow]['barrel']['fake']
    mbe = fakerates['m'][ptlow]['barrel']['err']
    me  = fakerates['m'][ptlow]['endcap']['fake']
    mee = fakerates['m'][ptlow]['endcap']['err']
    rowdict = {'ptlo' : ptlow, 'pthi' : pthigh, 'eb' : eb, 'ebe' : ebe, 'ee' : ee, 'eee' : eee, 'mb' : mb, 'mbe' : mbe, 'me' : me, 'mee' : mee}
    tables['fakerate']['rows'] += [tables['fakerate']['row'].format(**rowdict)]

# control regions
for r in ['PPF','PFP','FPP','PFF','FPF','FFP','FFF']:
    rowdict = {'fr':r}
    for c in ['eee','eem','mme','mmm']:
        rowdict['{0}d'.format(c)] = int(yields['yields'][c][r])
        rowdict['{0}de2'.format(c)] = yields['err2'][c][r]
        rowdict['{0}de'.format(c)] = rowdict['{0}de2'.format(c)]**0.5
        rowdict['{0}m'.format(c)] = 0.
        rowdict['{0}me2'.format(c)] = 0.
        for m in ['WZ','ZZ','TTV','VVV','ZG']:
             rowdict['{0}m'.format(c)] += yields['yieldsMC'][c][r][m]
             rowdict['{0}me2'.format(c)] += yields['err2MC'][c][r][m]
        rowdict['{0}me'.format(c)] = rowdict['{0}me2'.format(c)]**0.5
    tables['controlyields']['rows'] += [tables['controlyields']['datarow'].format(**rowdict)]
    tables['controlyields']['rows'] += [tables['controlyields']['mcrow'].format(**rowdict)]

sampMap = {'dy': {'name': 'Drell-Yan', 'samples': ['Z']}, 'tt': {'name': '$t\\bar{t} + t$', 'samples': ['T','TT']}}

for r in ['PPP','PPF','PFP','FPP','PFF','FPF','FFP','FFF']:
    for s in ['dy','tt']:
        rowdict = {}
        rowdict['name'] = r if s=='dy' else ''
        rowdict['sample'] = sampMap[s]['name']
        totval = 0.
        toterr2 = 0.
        for c in ['eee','eem','mme','mmm']:
            cval = 0.
            cerr2 = 0.
            for mc in sampMap[s]['samples']:
                cval += yields['yieldsMC'][c][r][mc]
                cerr2 += yields['err2MC'][c][r][mc]
            totval += cval
            toterr2 += cerr2
            cerr = cerr2**0.5
            rowdict[c] = cval
            rowdict[c+'e'] = cerr
        toterr = toterr2**0.5
        rowdict['tot'] = totval
        rowdict['tote'] = toterr
        tables['dyttbar']['rows'] += [tables['dyttbar']['row'].format(**rowdict)]


# final weighted control regions
for r in ['PPF','PFP','FPP','PFF','FPF','FFP','FFF']:
    eee = abs(yields['yieldsWeighted']['eee'][r])
    eem = abs(yields['yieldsWeighted']['eem'][r])
    mme = abs(yields['yieldsWeighted']['mme'][r])
    mmm = abs(yields['yieldsWeighted']['mmm'][r])
    eeee = yields['err2Weighted']['eee'][r]**0.5
    eeme = yields['err2Weighted']['eem'][r]**0.5
    mmee = yields['err2Weighted']['mme'][r]**0.5
    mmme = yields['err2Weighted']['mmm'][r]**0.5
    rowdict = {'fr' : r, 'eee' : eee, 'eeee': eeee, 'eem' : eem, 'eeme' : eeme, 'mme' : mme, 'mmee' : mmee, 'mmm' : mmm, 'mmme' : mmme}
    tables['controlregions']['rows'] += [tables['controlregions']['row'].format(**rowdict)]

# final cross sections
def getSystematic2(systs,sval,chan):
    su2 = 0.
    sd2 = 0.
    for syst in systs:
        su2 += (sval*systematics[syst][chan]['up'])**2
        sd2 += (sval*systematics[syst][chan]['down'])**2
    return [su2,sd2]

mcSyst = ['met','pileup','muon','electron']
#sigSyst = ['met','pileup','lepton','pdf','scale']
sigSyst = ['met','pileup','muon','electron']
lumiSyst = ['luminosity']
ddSyst = ['datadriven']
theoSyst = ['pdf','scale']
crosssectionSyst = {
    'TTV': 'ttv',
    'ZZ' : 'zz',
    'WW' : 'ww',
    'ZG' : 'zg',
    'VVV': 'vvv',
}


total = {'name' : 'Total', 'wz' : 0., 'wze2' : 0., 'wzs2' : 0., 'wzl2' : 0., 'bg' : 0., 'bge2' : 0., 'bgs2' : 0., 'bgm' : 0., 'bgme2' : 0., 'bgms2' : 0., 'bgml2' : 0., 'ex' : 0., 'exe2' : 0., 'exs2' : 0., 'exl2' : 0., 'obs' : 0., 'obse2' : 0.}
for c in ['eee','eem','mme','mmm']:
    name = names[c]
    wz = yields['yieldsMC'][c]['PPP']['WZ']
    wze2 = yields['err2MC'][c]['PPP']['WZ']
    wze = wze2**0.5
    wzs2 = max(getSystematic2(sigSyst,wz,c))
    wzs = wzs2**0.5
    wzl2 = max(getSystematic2(lumiSyst,wz,c))
    wzl = wzl2**0.5
    bgm = 0.
    bgme2 = 0.
    bgms2 = 0.
    bgml2 = 0.
    for mc in ['TTV','VVV','ZZ','ZG']:
        bgm += yields['yieldsMC'][c]['PPP'][mc]
        bgme2 += yields['err2MC'][c]['PPP'][mc]
        bgms2 += max(getSystematic2(mcSyst,yields['yieldsMC'][c]['PPP'][mc],c))
        bgml2 += max(getSystematic2(lumiSyst,yields['yieldsMC'][c]['PPP'][mc],c))
    bgme = bgme2**0.5
    bgms = bgms2**0.5
    bgml = bgml2**0.5
    bg = 0.
    bge2 = 0.
    bgs2 = 0.
    for fr in ['PPF','PFP','FPP','PFF','FPF','FFP','FFF']:
        bg += yields['yieldsWeighted'][c][fr]
        bge2 += yields['err2Weighted'][c][fr]
        #bgs2 += max(getSystematic2(ddSyst,yields['yieldsWeighted'][c][fr],c))
    bgs2 = max(getSystematic2(ddSyst,bg,c))
    bge = bge2**0.5
    bgs = bgs2**0.5
    ex = wz + bg + bgm
    exe2 = wze2 + bge2 + bgme2
    exs2 = wzs2 + bgs2 + bgms2
    exl2 = wzl2 + bgml2
    exe = exe2**0.5
    exs = exs2**0.5
    exl = exl2**0.5
    obs = yields['yields'][c]['PPP']
    obse2 = yields['err2'][c]['PPP']
    obse = obse2**0.5
    wzobs = obs - bg - bgm
    ssobs = wzobs/wz
    sse2 = ((obse2 + bge2 + bgme2)/wzobs**2 + wze2/wz**2) * ssobs**2
    sse = sse2**0.5
    #sigs[c]['stat'] = sse
    # sum them
    total['wz'] += wz
    total['wze2'] += wze2
    total['wzs2'] += wzs2
    total['wzl2'] += wzl2
    total['bg'] += bg
    total['bge2'] += bge2
    total['bgs2'] += bgs2
    total['bgm'] += bgm
    total['bgme2'] += bgme2
    total['bgms2'] += bgms2
    total['bgml2'] += bgml2
    total['ex'] += ex
    total['exe2'] += exe2
    total['exs2'] += exs2
    total['exl2'] += exl2
    total['obs'] += obs
    total['obse2'] += obse2
    rowdict = {'name':name, 'wz':wz, 'wze':wze, 'wzs':(wzs**2+wzl**2)**0.5, 'bg':bg, 'bge':bge, 'bgs':bgs, 'bgm':bgm, 'bgme':bgme, 'bgms':(bgms**2+bgml**2)**0.5, 'ex':ex, 'exe':exe, 'exs':(exs**2+exl**2)**0.5, 'obs':int(obs), 'obse':obse}
    tables['results']['rows'] += [tables['results']['row'].format(**rowdict)]
total['wze'] = total['wze2']**0.5
total['wzs'] = (total['wzs2']+total['wzl2'])**0.5
total['wzl'] = total['wzl2']**0.5
total['bge'] = total['bge2']**0.5
total['bgs'] = total['bgs2']**0.5
total['bgme'] = total['bgme2']**0.5
total['bgms'] = (total['bgms2']+total['bgml2'])**0.5
total['bgml'] = total['bgml2']**0.5
total['exe'] = total['exe2']**0.5
total['exs'] = (total['exs2']+ total['exl2'])**0.5
total['exl'] = total['exl2']**0.5
total['obse'] = total['obse2']**0.5
total['obs'] = int(total['obs'])
wzobs = total['obs'] - total['bg'] - total['bgm']
ssobs = wzobs/total['wz']
sse2 = ((total['obse2'] + total['bge2'] + total['bgme2'])/wzobs**2 + total['wze2']/total['wz']**2) * ssobs**2
sse = sse2**0.5
#sigs['wz']['stat'] = sse
tables['results']['rows'] += ['\\hline\n' + tables['results']['row'].format(**total)]

# generate signal strength table
for c in ['eee','eem','mme','mmm','wz']:
    # override lumi uncertainty
    sigs[c]['lumi'] = sigs[c]['val'] * systematics['luminosity']['eee']['up']
    if c=='wz':
        name = 'Total'
        tables['crosssection']['rows'] += ['\\hline']
    else: 
        name = names[c]
    sigs[c]['name'] = name
    tables['crosssection']['rows'] += [tables['crosssection']['row'].format(**sigs[c])]

# generate cross section equations
fidNums = {}
fidNums['xsec']     = int(fiducialTheory * sigs['wz']['val']) # sig strength * theory
fidNums['stat']     = int(fiducialTheory * sigs['wz']['stat'])
fidNums['systup']   = int(fiducialTheory * sigs['wz']['systup'])
fidNums['systdown'] = int(fiducialTheory * sigs['wz']['systdown'])
fidNums['lumi']     = int(fiducialTheory * sigs['wz']['lumi'])

totNums = {}
totNums['xsec']     = totalTheory * sigs['wz']['val'] # sig strength * theory
totNums['stat']     = totalTheory * sigs['wz']['stat']
totNums['systup']   = totalTheory * sigs['wz']['systup']
totNums['systdown'] = totalTheory * sigs['wz']['systdown']
totNums['lumi']     = totalTheory * sigs['wz']['lumi']



########################
### print everything ###
########################
for t in ['fakerate','controlyields','controlregions','results','crosssection','dyttbar']:
    print tables[t]['header']
    for row in tables[t]['rows']:
        print row
    print tables[t]['footer']
    print ''

print equations['fiducialcrosssection']['equation'].format(**fidNums)
print ''
print equations['totalcrosssection']['equation'].format(**totNums)
