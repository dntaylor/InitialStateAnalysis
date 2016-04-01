import itertools
import os
import sys
import errno
import argparse
from copy import deepcopy

_3L_MASSES = [170, 200, 250, 300, 350, 400, 450, 500, 600, 700]
_4L_MASSES = [130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]

ZMASS = 91.1876

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def getChannelSidebandSignalRegion(region,channel,**kwargs):
    mass = kwargs.pop('mass',500) # for higgs
    regionMap = {'Hpp3l':{}, 'Hpp4l':{}}
    regionMap['Hpp3l'][0] = {
        'sr' : 'hN.mass>0.9*%f && hN.mass<1.1*%f' %(mass,mass),
        'sb' : '((hN.mass<0.9*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
    }
    regionMap['Hpp3l'][1] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        'sb' : '((hN.mass<0.5*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
    }
    regionMap['Hpp3l'][2] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        'sb' : '((hN.mass<0.5*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
    }
    regionMap['Hpp4l'][0] = {
        'sr' : 'hN.mass>0.9*%f && hN.mass<1.1*%f' %(mass,mass),
        #'sb' : '((hN.mass<0.9*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
        'sb' : '(hN.mass>12. && hN.mass<800. && !(hN.mass>0.9*%f && hN.mass<1.1*%f))' %(mass,mass)
    }
    regionMap['Hpp4l'][1] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        #'sb' : '((hN.mass<0.5*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
        'sb' : '(hN.mass>12. && hN.mass<800. && !(hN.mass>0.5*%f && hN.mass<1.1*%f))' %(mass,mass)
    }
    regionMap['Hpp4l'][2] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        #'sb' : '((hN.mass<0.5*%f-20. && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
        'sb' : '(hN.mass>12. && hN.mass<800. && !(hN.mass>0.5*%f && hN.mass<1.1*%f))' %(mass,mass)
    }
    if region == 'Hpp3l':
        numTaus = channel[:2].count('t')
        theMap = regionMap['Hpp3l'][numTaus]
        cutMap = {'srcut': theMap['sr'].replace('hN','h1'), 'sbcut': theMap['sb'].replace('hN','h1')}
    if region == 'Hpp4l':
        h1Taus = channel[:2].count('t')
        h2Taus = channel[2:].count('t')
        h1Map = regionMap['Hpp4l'][h1Taus]
        h2Map = regionMap['Hpp4l'][h2Taus]
        cutMap = {
            'srcut': h1Map['sr'].replace('hN','h1') + ' && ' + h2Map['sr'].replace('hN','h2'),
            #'sbcut': h1Map['sb'].replace('hN','h1') + ' && ' + h2Map['sb'].replace('hN','h2'),
            'sbcut': '!(' + h1Map['sr'].replace('hN','h1') + ' && ' + h2Map['sr'].replace('hN','h2') + ')',
        }
    return cutMap


def getChannelCutFlowMap(region,channel,**kwargs):
    mass = kwargs.pop('mass',500) # for higgs
    regionMap = { 'Hpp3l' : {}, 'Hpp4l' : {}, 'WZ' : {}, 'Z' : {}, 'TT' : {}, }
    regionMap['Hpp3l'][0] = {
        #'st' : 'finalstate.sT>1.1*%f+60.' %mass,
        'st' : 'finalstate.sT>1.07*%f+45.' %mass,
        'zveto' : 'abs(z1.mass-%f)>80.' %ZMASS,
        'met' : None,
        'dphi' : 'abs(hN.dPhi)<%f/600.+1.95' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass,
        'dr' : 'hN.dR<%f/380.+2.06'%mass if mass<400 else 'hN.dR<%f/1200.+2.77'%mass,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp3l'][1] = {
        #'st' : 'finalstate.sT>0.85*%f+125.' %mass,
        'st' : 'finalstate.sT>0.72*%f+50.' %mass,
        'zveto' : 'abs(z1.mass-%f)>80.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'abs(hN.dPhi)<%f/200.+1.15' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.96'%mass if mass<400 else 'hN.dR<%f/1000.+2.6'%mass,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp3l'][2] = {
        #'st' : '(finalstate.sT>%f-10||finalstate.sT>200.)' %mass,
        'st' : 'finalstate.sT>0.44*%f+65' %mass,
        'zveto' : 'abs(z1.mass-%f)>50.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'abs(hN.dPhi)<2.1',
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.86'%mass if mass<400 else 'hN.dR<%f/750.+2.37'%mass,
        #'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass),
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp4l'][0] = {
        'st' : 'finalstate.sT>0.6*%f+130.' %mass,
        'zveto' : None,
        'dphi' : None,
        'dr' : None,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp4l'][1] = {
        'st' : '(finalstate.sT>%f+100.||finalstate.sT>400.)' %mass,
        'zveto' : 'abs(z1.mass-%f)>10.&&abs(z2.mass-%f)>10.' %(ZMASS,ZMASS),
        'dphi' : None,
        'dr' : None,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp4l'][2] = {
        'st' : 'finalstate.sT>120.',
        'zveto' : 'abs(z1.mass-%f)>50.&&abs(z2.mass-%f)>50.' %(ZMASS,ZMASS),
        'dphi' : 'abs(hN.dPhi)<2.5',
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        #'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass),
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }

    cutMap = {}
    if region == 'Hpp3l':
        numTaus = channel[:2].count('t')
        theMap = regionMap['Hpp3l'][numTaus]
        cuts = {}
        for cut in theMap:
            if theMap[cut]:
                cuts[cut] = theMap[cut].replace('hN','h1')
            else:
                cuts[cut] = '1'
        cuts['pre'] = '1'
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['met'], cuts['dphi'], cuts['mass']]
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
    if region == 'Hpp4l':
        h1Taus = channel[:2].count('t')
        h2Taus = channel[2:].count('t')
        h1Map = regionMap['Hpp4l'][h1Taus]
        h2Map = regionMap['Hpp4l'][h2Taus]
        cuts = {}
        for cut in h1Map:
            thisCut = []
            if h1Map[cut]:
                thisCut += [h1Map[cut].replace('hN','h1')]
            if h2Map[cut]:
                thisCut += [h2Map[cut].replace('hN','h2')]
            if thisCut:
                cuts[cut] = ' && '.join(thisCut)
            else:
                cuts[cut] = '1'
        cuts['pre'] = '1'
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
    return cutMap


def defineCutFlowMap(region,channels,mass):
    # define regions (based on number of taus in higgs candidate)
    regionMap = { 'Hpp3l' : {}, 'Hpp4l' : {}, 'WZ' : {}, 'Z' : {}, 'TT' : {}, }
    regionMap['Hpp3l'][0] = {
        #'st' : 'finalstate.sT>1.1*%f+60.' %mass,
        'st' : 'finalstate.sT>1.07*%f+45.' %mass,
        'zveto' : 'fabs(z1.mass-%f)>80.' %ZMASS,
        'met' : None,
        'dphi' : 'fabs(hN.dPhi)<%f/600.+1.95' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass,
        'dr' : 'hN.dR<%f/380.+2.06'%mass if mass<400 else 'hN.dR<%f/1200.+2.77'%mass,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp3l'][1] = {
        #'st' : 'finalstate.sT>0.85*%f+125.' %mass,
        'st' : 'finalstate.sT>0.72*%f+50.' %mass,
        'zveto' : 'fabs(z1.mass-%f)>80.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'fabs(hN.dPhi)<%f/200.+1.15' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.96'%mass if mass<400 else 'hN.dR<%f/1000.+2.6'%mass,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp3l'][2] = {
        #'st' : '(finalstate.sT>%f-10||finalstate.sT>200.)' %mass,
        'st' : 'finalstate.sT>0.44*%f+65' %mass,
        'zveto' : 'fabs(z1.mass-%f)>50.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'fabs(hN.dPhi)<2.1',
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.86'%mass if mass<400 else 'hN.dR<%f/750.+2.37'%mass,
        #'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass)
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][0] = {
        'st' : 'finalstate.sT>0.6*%f+130.' %mass,
        'zveto' : None,
        'dphi' : None,
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][1] = {
        'st' : '(finalstate.sT>%f+100.||finalstate.sT>400.)' %mass,
        'zveto' : 'fabs(z1.mass-%f)>10.&&fabs(z2.mass-%f)>10.' %(ZMASS,ZMASS),
        'dphi' : None,
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][2] = {
        'st' : 'finalstate.sT>120.',
        'zveto' : 'fabs(z1.mass-%f)>50.&&fabs(z2.mass-%f)>50.' %(ZMASS,ZMASS),
        'dphi' : 'fabs(hN.dPhi)<2.5',
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'mass' : None
    }
    regionMap['WZ'][0] = {
        'mass' : 'finalstate.mass>100.',
        'zpt' : '(z1.Pt1>20. && z1.Pt2>10.)',
        'zmass' : 'fabs(z1.mass-{0})<15.'.format(ZMASS),
        'bveto' : 'finalstate.bjetVeto20Tight==0',
        'wdr' : 'w1.dR1_z1_1>0.1 && w1.dR1_z1_2>0.1',
        'wmll' : 'w1.mll_z1_1>4. && w1.mll_z1_2>4.',
        'wpt' : 'w1.Pt1>20.',
        'met' : 'finalstate.met>30.',
    }
    regionMap['Z'][0] = {
        'zpt' : '(z1.Pt1>20.&z1.Pt2>10.)',
        'zmass' : 'fabs(z1.mass-%f)<20.' % ZMASS,
        'metveto' : 'met<30.',
        'bveto' : 'finalstate.bjetVeto30Medium==0',
    }
    regionMap['TT'][0] = {
        'zveto' : 'fabs(z1.mass-%f)>20.' % ZMASS,
        'met': 'met>30.',
        'jets' : 'finalstate.jetVeto30>1',
        'bjets' : 'finalstate.bjetVeto30Medium>0',
    }
    # define cutmap to be returned
    cutMap = { 'cuts' : [], 'labels': [], 'labels_simple': [], 'preselection': [] }
    if region == 'Hpp3l':
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','E_{T}^{miss}','#Delta#phi','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','MET','dPhi','Mass window']
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta#phi','Mass window']
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dPhi','Mass window']
        cutMap['labels_simple'] = ['Preselection','sT','ZVeto','dR','MassWindow']
        cutMap['preselection'] = ['All events', 'Three Lepton', 'Trigger', 'Fiducial',\
                                  'Trigger threshold', 'Lepton ID', 'Isolation', 'QCD Rejection',\
                                  'Lepton Charge', '4th Lepton Veto']
        cuts = { 'pre' : '',
                 'st' : '',
                 'zveto' : '',
                 'met' : '',
                 'dphi' : '',
                 'dr' : '',
                 'mass' : '' }
    elif region == 'Hpp4l':
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta#phi','Mass window']
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
        cutMap['labels'] = ['Preselection','s_{T}','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dPhi','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dR','Mass window']
        cutMap['labels_simple'] = ['Preselection','sT','MassWindow']
        cutMap['preselection'] = ['All events', 'Four Lepton', 'Trigger', 'Fiducial',\
                                  'Trigger threshold', 'Lepton ID', 'Isolation', 'QCD Rejection',\
                                  'Lepton Charge']
        cuts = { 'pre' : '',
                 'st' : '',
                 'zveto' : '',
                 'dphi' : '',
                 'dr' : '',
                 'mass' : '' }
    elif region == 'WZ':
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['labels'] = ['Preselection (ID)',
                            'Z lepton p_{T}',
                            'Z window',
                            'Mass 3l',
                            'b-jet Veto',
                            #'W #DeltaR to Z',
                            'M(W_\\ell,Z_\\ell)',
                            'W lepton p_{T}',
                            'E_{T}^{miss}',
                           ]
        cutMap['labels_simple'] = ['Preselection',
                                   'ZLepPt',
                                   'Zwindow',
                                   'mass3l',
                                   'bjetVeto',
                                   #'WDR',
                                   'WMll',
                                   'WLepPt',
                                   'MET',
                                  ]
        cutMap['cuts'] = ['1',
                          regionMap['WZ'][0]['zpt'],
                          regionMap['WZ'][0]['zmass'],
                          regionMap['WZ'][0]['mass'],
                          regionMap['WZ'][0]['bveto'],
                          #regionMap['WZ'][0]['wdr'],
                          regionMap['WZ'][0]['wmll'],
                          regionMap['WZ'][0]['wpt'],
                          regionMap['WZ'][0]['met'],
                         ]
        # this is 8 tev, add a period check later
        #cutMap['labels'] = ['Preselection (ID)', 'Z lepton p_{T}', 'Z window', 'Mass 3l', 'W #DeltaR to Z', 'W lepton p_{T}', 'E_{T}^{miss}']
        #cutMap['labels_simple'] = ['Presel (ID)', 'Z lep pt', 'Z window', 'mass3l', 'W dR', 'W lep pt', 'MET']
        #cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        #cutMap['cuts'] = ['1', regionMap['WZ'][0]['zpt'], regionMap['WZ'][0]['zmass'], regionMap['WZ'][0]['mass'],\
        #                  regionMap['WZ'][0]['wdr'], regionMap['WZ'][0]['wpt'], regionMap['WZ'][0]['met']]
    elif region=='Z':
        cutMap['labels'] = ['Preselection (ID)', 'Z lepton p_{T}', 'Z window', 'met veto', 'b veto']
        cutMap['labels_simple'] = ['Presel (ID)', 'Z lep pt', 'Z window', 'met veto', 'b veto']
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['cuts'] = ['1', regionMap['Z'][0]['zpt'], regionMap['Z'][0]['zmass'], regionMap['Z'][0]['metveto'], regionMap['Z'][0]['bveto']]
    elif region=='TT':
        cutMap['labels'] = ['Preselection (ID)', 'Z veto', 'met', 'jets', 'b jets']
        cutMap['labels_simple'] = ['Presel (ID)', 'Z veto', 'met', 'jets', 'b jets']
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['cuts'] = ['1', regionMap['TT'][0]['zveto'], regionMap['TT'][0]['met'], regionMap['TT'][0]['jets'], regionMap['TT'][0]['bjets']]
    else:
        cutMap['cuts'] = '1'
        cutMap['labels'] = ['%s Full Selection' %region]
        cutMap['labels_simple'] = [region]
    if region not in ['Hpp3l','Hpp4l']: return cutMap
    usedLepPairs = []
    for channel in channels:
        lepPairs = [channel[:2]]
        #if region=='Hpp4l': lepPairs += [channel[2:]]
        if lepPairs in usedLepPairs: continue
        for cut in cuts:
            if cuts[cut] and cuts[cut][-2:]!='||': cuts[cut] += '||'
        hNum = 0
        tempCut = {}
        for lepPair in lepPairs:
            hNum += 1
            numTau = lepPair.count('t')
            for cut in regionMap[region][numTau]:
                thisCut = regionMap[region][numTau][cut]
                if thisCut is not None:
                    if cut in tempCut:
                        tempCut[cut] += '&&'
                    else:
                        tempCut[cut] = ''
                    tempCut[cut] += 'h%iFlv=="%s"&&%s' % (hNum, lepPair, thisCut.replace('N',str(hNum)))
                else:
                    if cut in tempCut:
                        tempCut[cut] += '&&'
                    else:
                        tempCut[cut] = ''
                    tempCut[cut] += 'h%iFlv=="%s"' % (hNum, lepPair)
        for cut in cuts:
            if cut in tempCut:
                cuts[cut] += '(%s)' % tempCut[cut]
        usedLepPairs += [lepPairs]
    for cut in cuts:
        if not cuts[cut]:
            cuts[cut] = '('+'||'.join(['h1Flv=="%s"' %x[0] for x in usedLepPairs])+')'
        else:
            if cuts[cut][-2:]=='||': cuts[cut] = cuts[cut][:-2]
            cuts[cut] = '(%s)' % cuts[cut]

    if region == 'Hpp3l':
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['met'], cuts['dphi'], cuts['mass']]
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
    if region == 'Hpp4l':
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['mass']]
    return cutMap

def getChannels(numLeptons,**kwargs):
    '''Get channels for a given region.'''
    runTau = kwargs.pop('runTau',False)
    leptons = ['l%i' %(x+1) for x in range(numLeptons)]
    lepTypes = 'emt' if runTau else 'em'
    lepPairs = [x[0]+x[1] for x in itertools.product(lepTypes,repeat=2)]
    if numLeptons == 2:
        channels = [x[0]+x[1] for x in itertools.product(lepTypes,lepTypes)]
    elif numLeptons == 3:
        channels = [x[0]+x[1] for x in itertools.product(lepPairs,lepTypes)]
    else:
        channels = [x[0]+x[1] for x in itertools.product(lepPairs,lepPairs)]
    return channels,leptons

def getMergeDict(period):
    '''Return a dictionary of samples to merge in plotting.'''
    sampleMergeDict = {}
    # 13 tev
    if period==13:
        # phys14 samples
        #sampleMergeDict['SingleTop'] = {
        #     'TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola': '1',
        #     'TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola': '1',
        #     'TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'   : '1',
        #     'TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'   : '1',
        #     'T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola'          : '1',
        #     'Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola'       : '1',
        #}
        #sampleMergeDict['Diboson']   = {
        #    'WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola': '1',
        #    'ZZTo4L_Tune4C_13TeV-powheg-pythia8'       : '1',
        #}
        #sampleMergeDict['WZJets']    = {
        #    'WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola': '1',
        #}
        #sampleMergeDict['ZZJets']    = {
        #    'ZZTo4L_Tune4C_13TeV-powheg-pythia8': '1',
        #}
        #sampleMergeDict['TTJets']    = {
        #    'TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola': '1',
        #}
        #sampleMergeDict['ZJets']     = {
        #    'DYJetsToLL_M-50_13TeV-madgraph-pythia8': '1',
        #}
        #sampleMergeDict['WJets']     = {
        #    'WJetsToLNu_13TeV-madgraph-pythia8-tauola': '1',
        #}
        #sampleMergeDict['TTVJets']   = {
        #    'TTWJets_Tune4C_13TeV-madgraph-tauola': '1',
        #    'TTZJets_Tune4C_13TeV-madgraph-tauola': '1',
        #}
        # RunIISpring15Dr
        sampleMergeDict['SingleTop'] = {
            'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       : '1',
            #'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       : '1',
            #'ST_t-channel_5f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       : '1',
            'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : '1',
            'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1'     : '1',
            #'ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'  : '1',
            'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'     : '1',
            #'ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'      : '1',
            'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'         : '1',
        }
        sampleMergeDict['TTJets'] = {
            'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
            #'TTTo2L2Nu_13TeV-powheg'                         : '1',
            #'TT_TuneCUETP8M1_13TeV-powheg-pythia8'           : '1',
        }
        sampleMergeDict['ZJetsFiltered'] = {
            'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_filtered' : '1',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_filtered'     : '1',
        }
        sampleMergeDict['ZJets'] = {
            'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'     : '1',
        }
        sampleMergeDict['WJets'] = {
            'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
        }
        sampleMergeDict['WWJets'] = {
            'WWTo2L2Nu_13TeV-powheg'        : '1',
            #'WWTo4Q_13TeV-powheg'           : '1',
            #'WWToLNuQQ_13TeV-powheg'        : '1',
            #'WW_TuneCUETP8M1_13TeV-pythia8' : '1',
        }
        sampleMergeDict['WZJets'] = {
            #'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8' : '1',
            'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8'     : '1',
            #'WZ_TuneCUETP8M1_13TeV-pythia8'                  : '1',
        }
        sampleMergeDict['ZZJets'] = {
            #'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'   : '1',
            'ZZTo4L_13TeV_powheg_pythia8'                   : '1',
            #'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8'     : '1',
            #'ZZ_TuneCUETP8M1_13TeV-pythia8'                 : '1',
            'GluGluToZZTo2e2mu_BackgroundOnly_13TeV_MCFM'   : '1',
            'GluGluToZZTo2e2tau_BackgroundOnly_13TeV_MCFM'  : '1',
            'GluGluToZZTo2mu2tau_BackgroundOnly_13TeV_MCFM' : '1',
            'GluGluToZZTo4e_BackgroundOnly_13TeV_MCFM'      : '1',
            'GluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM'     : '1',
            'GluGluToZZTo4tau_BackgroundOnly_13TeV_MCFM'    : '1',
        }
        sampleMergeDict['TTZJets'] = {
            'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8'         : '1',
            #'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8'                  : '1',
        }
        sampleMergeDict['TTWJets'] = {
            'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8' : '1',
            #'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8'  : '1',
        }
        sampleMergeDict['TTVJets'] = {
            'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8'         : '1',
            #'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8'                  : '1',
            'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8' : '1',
            #'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8'  : '1',
        }
        sampleMergeDict['VVVJets']   = {
            'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8': '1',
            'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8': '1',
            'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8': '1',
            #'WWW_TuneCUETP8M1_13TeV-amcatnlo-pythia8': '1',
        }
        sampleMergeDict['ZGFiltered'] = {
            'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_filtered' : '1',
        }
        sampleMergeDict['ZG'] = {
            'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
        }
        sampleMergeDict['WG'] = {
            'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
        }
        sampleMergeDict['QCD'] = {
            #'QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8' : '1',
            #'QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8' : '1',

            #'QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8' : '1',

            #'QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8' : '1',
            #'QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8' : '1',
            #'QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8' : '1',

            'QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8'      : '1',
            'QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8'    : '1',
            'QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8'  : '1',
            'QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8'  : '1',
        }
        sampleMergeDict['data']      = {
            #'data_Run2015B': '1', # 50 ns
            'data_Run2015C': '1', # 25 ns and 50 ns
            'data_Run2015D': '1', # 25 ns
            #'data_DoubleEG_Run2015C_05Oct2015_25ns'       : '1',
            #'data_DoubleEG_Run2015D_05Oct2015_25ns'       : '1',
            #'data_DoubleEG_Run2015D_PromptReco-v4_25ns'   : '1',
            #'data_DoubleMuon_Run2015C_05Oct2015_25ns'     : '1',
            #'data_DoubleMuon_Run2015D_05Oct2015_25ns'     : '1',
            #'data_DoubleMuon_Run2015D_PromptReco-v4_25ns' : '1',
            #'data_MuonEG_Run2015C_05Oct2015_25ns'         : '1',
            #'data_MuonEG_Run2015D_05Oct2015_25ns'         : '1',
            #'data_MuonEG_Run2015D_PromptReco-v4_25ns'     : '1',
        }
    # 8 TeV sample aliases
    if period==8:
        sampleMergeDict['WWJets']    = {
            'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola': '1',
            #'WW_TuneZ2star_8TeV_pythia6_tauola': '1',
        }
        sampleMergeDict['WZJets']    = {
            'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola': '1',
            'WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola'    : '1',
            #'WZ_TuneZ2star_8TeV_pythia6_tauola': '1',
        }
        sampleMergeDict['ZZJets']    = {
            'ZZTo2e2mu_8TeV-powheg-pythia6'                : '1',
            'ZZTo2e2tau_8TeV-powheg-pythia6'               : '1',
            'ZZTo2mu2tau_8TeV-powheg-pythia6'               : '1',
            'ZZTo4e_8TeV-powheg-pythia6'                    : '1',
            'ZZTo4mu_8TeV-powheg-pythia6'                   : '1',
            'ZZTo4tau_8TeV-powheg-pythia6'                  : '1',
            #'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola'    : '1',
            'ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola'  : '1',
            'ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola' : '1',
            #'ZZ_TuneZ2star_8TeV_pythia6_tauola'             : '1',
            'GluGluToZZTo4L_8TeV-gg2zz-pythia6'             : '1',
            'GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6': '1',
        }
        sampleMergeDict['HZZ']       = {
            'GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6': '1',
        }
        sampleMergeDict['WGJets']    = {
            'WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola' : '1',
            'WGstarToLNu2Mu_TuneZ2star_8TeV-madgraph-tauola': '1',
        }
        sampleMergeDict['SingleTop'] = {
            'T_s-channel_TuneZ2star_8TeV-powheg-tauola'       : '1',
            'T_t-channel_TuneZ2star_8TeV-powheg-tauola'       : '1',
            'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'   : '1',
            'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'    : '1',
            'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'    : '1',
            'Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola': '1',
        }
        sampleMergeDict['TTJets']    = {
            'TTJetsFullLepMGDecays': '1',
            'TTJetsSemiLepMGDecays': '1',
        }
        sampleMergeDict['ZJets']     = {
            #'DYJetsToLL_M-10To50filter_8TeV-madgraph'         : '1',
            #'Z1jets_M50'                                      : '1',
            #'Z2jets_M50_S10'                                  : '1',
            #'Z3jets_M50'                                      : '1',
            #'Z4jets_M50'                                      : '1',
            #'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball': 'event.GenNUP==5',
            #'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball': '1',
            'DYJetsToLL_M-10To50filter_8TeV-madgraph_filtered'         : '1',
            'Z1jets_M50_filtered'                                      : '1',
            'Z2jets_M50_S10_filtered'                                  : '1',
            'Z3jets_M50_filtered'                                      : '1',
            'Z4jets_M50_filtered'                                      : '1',
            'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_filtered': 'event.GenNUP==5',
            #'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_filtered': '1',
        }
        sampleMergeDict['WJets']     = {
            'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v1': '1',
            'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v2': '1',
        }
        sampleMergeDict['ZG']        = {
            'ZGToLLG_8TeV-madgraph': '1',
        }
        sampleMergeDict['data']      = {
            'data_Run2012A': '1',
            'data_Run2012B': '1',
            'data_Run2012C': '1',
            'data_Run2012D': '1',
            #'data_DoubleElectron_Run2012A_22Jan2013_v1': '1',
            #'data_DoubleElectron_Run2012B_22Jan2013_v1': '1',
            #'data_DoubleElectron_Run2012C_22Jan2013_v1': '1',
            #'data_DoubleElectron_Run2012D_22Jan2013_v1': '1',
            #'data_DoubleMu_Run2012A_22Jan2013_v1': '1',
            #'data_DoubleMuParked_Run2012B_22Jan2013_v1': '1',
            #'data_DoubleMuParked_Run2012C_22Jan2013_v1': '1',
            #'data_DoubleMuParked_Run2012D_22Jan2013_v1': '1',
            #'data_MuEG_Run2012A_22Jan2013_v1': '1',
            #'data_MuEG_Run2012B_22Jan2013_v1': '1',
            #'data_MuEG_Run2012C_22Jan2013_v1': '1',
            #'data_MuEG_Run2012D_22Jan2013_v1': '1',
        }
        sampleMergeDict['TTVJets']   = {
            'TTZJets' : '1',
            'TTWJets' : '1',
            'TTWWJets': '1',
            'TTGJets' : '1',
        }
        sampleMergeDict['VVVJets']   = {
            'ZZZNoGstarJets': '1',
            'WZZNoGstarJets': '1',
            'WWZNoGstarJets': '1',
            'WWWJets'       : '1',
        }
    return sampleMergeDict


def getSigMap(numLeptons,mass=500):
    '''Return a signal map for a given running period'''
    sigMap = {
        8 : {
             'ZZ'  : 'ZZJets',
             'WZ'  : 'WZJets',
             'WW'  : 'WWJets',
             'WG'  : 'WGJets',
             'Z'   : 'ZJets',
             'Zlow': 'DYJetsToLL_M-10To50filter_8TeV-madgraph',
             'W'   : 'WJets',
             'ZG'  : 'ZG',
             'HZZ' : 'HZZ',
             'DB'  : 'Diboson',
             'TT'  : 'TTJets',
             'T'   : 'SingleTop',
             'TTV' : 'TTVJets',
             'TTZ' : 'TTZJets',
             'TTW' : 'TTWJets',
             'VVV' : 'VVVJets',
             'ZZZ' : 'ZZZNoGstarJets',
             'WWZ' : 'WWZNoGstarJets',
             'WWW' : 'WWWJets',
             'Sig' : 'HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % mass\
                      if numLeptons==3 else 'HPlusPlusHMinusMinusHTo4L_M-%i_8TeV-pythia6' % mass,
             'SigAP' : 'HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % mass,
             'SigPP' : 'HPlusPlusHMinusMinusHTo4L_M-%i_8TeV-pythia6' % mass,
             'data': 'data',
             'datadriven': 'datadriven',
        },
        13 : {
             'ZZ'  : 'ZZJets',
             'WZ'  : 'WZJets',
             'WW'  : 'WWJets',
             'Z'   : 'ZJets',
             'Zfiltered'   : 'ZJetsFiltered',
             'W'   : 'WJets',
             'TT'  : 'TTJets',
             'T'   : 'SingleTop',
             'TTV' : 'TTVJets',
             'VVV' : 'VVVJets',
             'QCD' : 'QCD',
             'WG'  : 'WG',
             'ZG'  : 'ZG',
             'ZGfiltered'  : 'ZGFiltered',
             'Sig' : 'DBLH_m500',
             'data': 'data',
             'datadriven': 'datadriven',
        }
    }
    if numLeptons==3:
        for m in _3L_MASSES:
            sigMap[8][m] = 'HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % m
    if numLeptons==4:
        for m in _4L_MASSES:
            sigMap[8][m] = 'HPlusPlusHMinusMinusHTo4L_M-%i_8TeV-pythia6' % m
    return sigMap

def getIntLumiMap():
    '''Get map of integrated luminosity to scale MC'''
    intLumiMap = {
        7 : 4900,
        8 : 19700,
        #13: 42, # 50ns
        #13: 1280, # 25ns, AN freeze
        # brilcalc lumi --normtag  ~lumipro/public/normtag_file/OfflineNormtagV2.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-258750_13TeV_PromptReco_Collisions15_25ns_JSON.txt
        #13: 1341, # 25ns, updated lumicalc 
        # brilcalc lumi --normtag  ~lumipro/public/normtag_file/OfflineNormtagV2.json -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt
        #13: 2263, # fule 25ns
        # moriond update
        # brilcalc lumi -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt --normtag /afs/cern.ch/user/l/lumipro/public/normtag_file/moriond16_normtag.json
        13: 2318, 
    }
    return intLumiMap

def getChannelStringsCuts(region,channels):
    channelCharMap = {'e':'e', 'm':'#mu', 't':'#tau'}
    channelStrings = []
    for channel in channels:
        channelString = ''
        for c in channel:
            channelString += channelCharMap[c]
        channelStrings += [channelString]
    channelCuts = ['channel=="%s"' % x for x in channels]
    channelsWZ = [['ee','e'],['ee','m'],['mm','e'],['mm','m']]
    channelStringsWZ = ['(ee)e','(ee)#mu','(#mu#mu)e','(#mu#mu)#mu']
    channelCutsWZ = ['z1Flv=="%s"&&w1Flv=="%s"' %(x[0],x[1]) for x in channelsWZ]
    channelsZ = [['ee'],['mm']]
    channelStringsZ = ['ee','#mu#mu']
    channelCutsZ = ['z1Flv=="%s"' %(x[0]) for x in channelsZ]
    plotChannelCuts = channelCuts
    if region in ['WZ']: plotChannelCuts = channelCutsWZ
    if region in ['Z']: plotChannelCuts = channelCutsZ
    plotChannelStrings = channelStrings
    if region in ['WZ']: plotChannelStrings = channelStringsWZ
    if region in ['Z']: plotChannelStrings = channelStringsZ
    return plotChannelStrings, plotChannelCuts

def getGenChannelStringsCuts(region,genChannels):
    channelCharMap = {'e':'e', 'm':'#mu', 't':'#tau'}
    channelStrings = []
    channelCuts = []
    for channel in genChannels:
        channelString = ''
        for c in channel:
            channelString += channelCharMap[c]
        channelStrings += [channelString]
        theGenChannels = [channel+c for c in ['e','m','t']]
        genCut = '(' + '||'.join(['genChannel=="%s"' %x for x in theGenChannels]) + ')'
        channelCuts += [genCut]
    return channelStrings, channelCuts
            
    

def getChannelBackgrounds(runPeriod):
    channelBackground = {
        'Hpp2l'   : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'TT'      : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Z'       : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Charge'  : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'W'       : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'FakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'TTFakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'HZZFakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Hpp3l'   : ['TT', 'TTV', 'Z', 'ZG', 'VVV', 'ZZ', 'WZ'],
        'WZ'      : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
        'NoVeto'  : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
        'LowMass' : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
        'Hpp4l'   : ['WZ', 'TTV', 'VVV', 'ZZ'],
        'ZZ'      : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
    }   
    if runPeriod==13:
        channelBackground = {
            'WZ'      : ['T', 'TT', 'TTV', 'Zfiltered', 'ZGfiltered', 'VVV', 'WW', 'ZZ', 'WZ'],
            #'WZMCClosure'  : ['TT','Z'],
            'WZMCClosure'  : ['QCD','T','TT','W','WW','Z'],
            'WZdatadriven' : ['TTV', 'ZG', 'VVV', 'ZZ', 'WZ'],
            'NoVeto'  : ['T', 'TT', 'TTV', 'Z', 'VVV', 'WW', 'ZZ', 'WZ'],
            'W'       : ['T', 'TT', 'TTV', 'W', 'Z', 'WW', 'ZZ', 'WZ'],
            'FakeRate': ['QCD','T', 'TT', 'TTV', 'W', 'Z', 'WW', 'VVV', 'ZZ', 'WZ'],
            'TTFakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'WW', 'VVV', 'ZZ', 'WZ'],
            'HZZFakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'WW', 'ZZ', 'WZ'],
            'Z'       : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'Charge'  : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'TT'      : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'Hpp2l'   : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'Hpp3l'   : ['TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'LowMass' : ['TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'Hpp4l'   : ['TT', 'Z', 'TTV', 'ZZ', 'WZ']
        }
    return channelBackground

def plotLepton(plotMethod,myCut,obj,**kwargs):
    savedir = kwargs.pop('savedir','')
    doDetailed = kwargs.pop('doDetailed',False)
    doMinimal = kwargs.pop('doMinimal',False)
    pre = kwargs.pop('pre','')
    post = kwargs.pop('post','')
    name = kwargs.pop('name','')
    pretty = kwargs.pop('pretty','')
    plotMethod('%s.%sPt%s' %(obj,pre,post), [20,0,200],           savedir+'%s/Pt' %name, yaxis='Events/10 GeV',xaxis='p_{T}^{%s} (GeV)' %pretty,       legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('%s.%sIso%s' %(obj,pre,post),[50,0,.5],            savedir+'%s/Iso' %name,yaxis='Events',        xaxis='Relative Isolation (%s)' %pretty,legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('%s.%sEta%s' %(obj,pre,post),[30,-3.0,3.0],        savedir+'%s/Eta' %name,yaxis='Events',        xaxis='\\eta^{%s}' %pretty,             legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
    plotMethod('%s.%sPhi%s' %(obj,pre,post),[30,-3.14159,3.14159],savedir+'%s/Phi' %name,yaxis='Events',        xaxis='\\phi^{%s}' %pretty,             legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
    if not doMinimal:
        #plotMethod('%s.ChargeConsistent' %name,[3,-1.5,1.5],savedir+'%sChargeId' %name,yaxis='Events',xaxis='Charge ID (\\ell%i)' %(l+1),legendpos=31,logy=0,cut=myCut,**kwargs)    
        if doDetailed:
            plotMethod('%s.%sDxy%s' %(obj,pre,post),                     [50,-0.1,0.1],savedir+'%s/Dxy' %name,             yaxis='Events',xaxis='d_{0}^{%s}' %pretty,                      legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sDz%s' %(obj,pre,post),                      [50,-0.1,0.1],savedir+'%s/Dz' %name,              yaxis='Events',xaxis='d_{Z}^{%s}' %pretty,                      legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sSigmaIEtaIEta%s' %(obj,pre,post),           [40,0.,0.08], savedir+'%s/SigmaIEtaIEta' %name,   yaxis='Events',xaxis='\\sigma_{i\\eta i\\eta}^{%s}' %pretty,    legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sDEtaIn%s' %(obj,pre,post),                  [50,-0.1,0.1],savedir+'%s/DEtaIn' %name,          yaxis='Events',xaxis='\\Delta\\eta_{In}^{%s}' %pretty,          legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sDPhiIn%s' %(obj,pre,post),                  [50,-0.1,0.1],savedir+'%s/DPhiIn' %name,          yaxis='Events',xaxis='\\Delta\\phi_{In}^{%s}' %pretty,          legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sHOverE%s' %(obj,pre,post),                  [50,0.,0.1],  savedir+'%s/HOverE' %name,          yaxis='Events',xaxis='\\frac{H}{E}^{%s}' %pretty,               legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sOoEmOoP%s' %(obj,pre,post),                 [50,0.,0.1],  savedir+'%s/OoEmOoP' %name,         yaxis='Events',xaxis='|\\frac{1}{E}-\\frac{1}{P}|^{%s}' %pretty,legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sTriggeringMVA%s' %(obj,pre,post),           [50,-1.,1.],  savedir+'%s/TriggeringMVA' %name,   yaxis='Events',xaxis='MVA_{Triggering} (%s)' %pretty,           legendpos=42,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sNonTriggeringMVA%s' %(obj,pre,post),        [50,-1.,1.],  savedir+'%s/NonTriggeringMVA' %name,yaxis='Events',xaxis='MVA_{NonTriggering} (%s)' %pretty,        legendpos=42,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sExpectedMissingInnerHits%s' %(obj,pre,post),[6,0,6],      savedir+'%s/MissingHits' %name,     yaxis='Events',xaxis='Missing Hits (%s)' %pretty,               legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sPassConversionVeto%s' %(obj,pre,post),      [2,0,2],      savedir+'%s/PassConversion' %name,  yaxis='Events',xaxis='Pass Conversion Veto (%s)' %pretty,       legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
            plotMethod('%s.%sNormalizedChi2%s' %(obj,pre,post),          [60,0.,15.],  savedir+'%s/NormalizedChi2' %name,  yaxis='Events',xaxis='\\chi_{Norm}^{2} (%s)' %pretty,           legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sIsGlobalMuon%s' %(obj,pre,post),            [2,0,2],      savedir+'%s/IsGlobalMuon' %name,    yaxis='Events',xaxis='Gloabl Muon (%s)' %pretty,                legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
            plotMethod('%s.%sIsPFMuon%s' %(obj,pre,post),                [2,0,2],      savedir+'%s/IsPFMuon' %name,        yaxis='Events',xaxis='PFMuon (%s)' %pretty,                     legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
            plotMethod('%s.%sIsTrackerMuon%s' %(obj,pre,post),           [2,0,2],      savedir+'%s/IsTrackerMuon' %name,   yaxis='Events',xaxis='Tracker Muon (%s)' %pretty,               legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
            plotMethod('%s.%sValidMuonHits%s' %(obj,pre,post),           [80,0,80],    savedir+'%s/ValidMuonHits' %name,   yaxis='Events',xaxis='Muon Hits (%s)' %pretty,                  legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sMatchedStations%s' %(obj,pre,post),         [6,0,6],      savedir+'%s/MatchedStations' %name, yaxis='Events',xaxis='Matched Stations (%s)' %pretty,           legendpos=31,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sValidPixelHits%s' %(obj,pre,post),          [8,0,8],      savedir+'%s/ValidPixelHits' %name,  yaxis='Events',xaxis='Valid Pixel Hits (%s)' %pretty,           legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sTrackerLayers%s' %(obj,pre,post),           [20,0,20],    savedir+'%s/TrackerLayers' %name,   yaxis='Events',xaxis='Tracker Layers (%s)' %pretty,             legendpos=31,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sGenPdgId%s' %(obj,pre,post),                [61,-30,31],  savedir+'%s/GenPdgId' %name,        yaxis='Events',xaxis='Gen PDG ID (%s)' %pretty,                 legendpos=31,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sMotherGenPdgId%s' %(obj,pre,post),          [61,-30,31],  savedir+'%s/MotherGenPdgId' %name,  yaxis='Events',xaxis='Mother Gen PDG ID (%s)' %pretty,          legendpos=31,logy=0,cut=myCut,**kwargs)
            plotMethod('%s.%sPassLoose%s' %(obj,pre,post),               [2,0,2],      savedir+'%s/PassLoose' %name,       yaxis='Events',xaxis='Pass Loose (%s)' %pretty,                 legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
            plotMethod('%s.%sPassTight%s' %(obj,pre,post),               [2,0,2],      savedir+'%s/PassTight' %name,       yaxis='Events',xaxis='Pass Tight (%s)' %pretty,                 legendpos=43,logy=0,cut=myCut,numcol=3,**kwargs)
            plotMethod('%s.%sJetPt%s' %(obj,pre,post),                   [40,0,200],   savedir+'%s/JetPt' %name,           yaxis='Events/5.0 GeV',xaxis='p_{T}^{%s jet} (GeV)' %pretty,    legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
            plotMethod('%s.%sJetBTag%s' %(obj,pre,post),                 [50,0.,1.],   savedir+'%s/JetBtag' %name,         yaxis='Events',xaxis='Jet b Tag (%s)' %pretty,                  legendpos=42,logy=0,cut=myCut,**kwargs)

def plotDistributions(plotMethod,myCut,nl,isControl,**kwargs):
    savedir = kwargs.pop('savedir','')
    analysis = kwargs.pop('analysis','')
    region = kwargs.pop('region','')
    mass = kwargs.pop('mass',500)
    doDetailed = kwargs.pop('doDetailed',False)
    doMinimal = kwargs.pop('doMinimal',False)
    doPAS = kwargs.pop('doPAS',False)
    if savedir: savedir += '/'

    # fast track pas plots
    if doPAS:
        if analysis in ['WZ']:
            mtext = 'm_{3\\ell} (GeV)'
            if analysis in ['Hpp4l']: mtext = 'm_{4\\ell} (GeV)'
            #if analysis in ['WZ']: mtext = 'm_{\\ell\\ell\\\'\\ell\\\'} (GeV)'
            plotMethod('finalstate.mass',[25,0,500],savedir+'mass',yaxis='Events/20 GeV',xaxis=mtext,logy=0,cut=myCut,overflow=True,**kwargs)
            plotMethod('z1.mass', [13,58.5,123.5],   savedir+'z1/Mass_newWidth_wide', yaxis='Events/5 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,**kwargs)
        if analysis in ['Hpp3l']:
            #plotMethod('h1.mass', [28,0,700],savedir+'hpp/Mass', yaxis='Events/25 GeV',xaxis='m_{\\ell^{\\pm}\\ell^{\\pm}} (GeV)', legendpos=43,yscale=3.5,logy=1,cut=myCut,overflow=True,numcol=2,ratiomin=0.,ratiomax=2.,**kwargs)
            plotMethod('h1.mass', [28,0,700],savedir+'hpp/Mass', yaxis='Events/25 GeV',xaxis='m_{l^{#pm}l^{#pm}} (GeV)', legendpos=43,yscale=3.5,logy=1,cut=myCut,overflow=True,numcol=2,ratiomin=0.,ratiomax=2.,**kwargs)
        if analysis in ['Hpp4l']:
            #plotMethod('h1.mass', [14,0,700],savedir+'hpp/Mass', yaxis='Events/50 GeV',xaxis='m_{\\ell^{+}\\ell^{+}} (GeV)', legendpos=43,yscale=3.5,logy=1,cut=myCut,overflow=True,numcol=2,ratiomin=0.,ratiomax=2.,**kwargs)
            plotMethod('h1.mass', [14,0,700],savedir+'hpp/Mass', yaxis='Events/50 GeV',xaxis='m_{l^{+}l^{+}} (GeV)', legendpos=43,yscale=3.5,logy=1,cut=myCut,overflow=True,numcol=2,ratiomin=0.,ratiomax=2.,**kwargs)
        return


    plotMethod('finalstate.sT',[40,0,1000],savedir+'sT',yaxis='Events/25 GeV',xaxis='S_{T} (GeV)',lumitext=33,logy=0,ymin=0.1,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.sT',[40,0,1000],savedir+'sT_log',yaxis='Events/25 GeV',xaxis='S_{T} (GeV)',lumitext=33,logy=1,ymin=0.1,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.elecVetoLoose',[8,0,8],savedir+'numElectronsLoose',yaxis='Events',xaxis='Number of Electrons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.muonVetoLoose',[8,0,8],savedir+'numMuonsLoose',yaxis='Events',xaxis='Number of Muons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.elecVetoLoose+finalstate.muonVetoLoose',[8,0,8],savedir+'numLeptonsLoose',yaxis='Events',xaxis='Number of Leptons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.jetVeto30',[8,0,8],savedir+'numJets30',yaxis='Events',xaxis='Number of Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.bjetVeto30Medium',[8,0,8],savedir+'numBJets30Medium',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    if not doMinimal:
        plotMethod('finalstate.elecVetoTight',[8,0,8],savedir+'numElectronsTight',yaxis='Events',xaxis='Number of Electrons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('finalstate.muonVetoTight',[8,0,8],savedir+'numMuonsTight',yaxis='Events',xaxis='Number of Muons (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.bjetVeto30Loose',[8,0,8],savedir+'numBJets30Loose',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.bjetVeto30Tight',[8,0,8],savedir+'numBJets30Tight',yaxis='Events',xaxis='Number of b Jets (p_{T}>30 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.jetVeto40',[8,0,8],savedir+'numJets40',yaxis='Events',xaxis='Number of Jets (p_{T}>40 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.muonVeto5',[8,0,8],savedir+'muonVeto5',yaxis='Events',xaxis='Muon Veto (p_{T}>5 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.muonVeto10Loose',[8,0,8],savedir+'muonVeto10',yaxis='Events',xaxis='Muon Veto (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.muonVeto15',[8,0,8],savedir+'muonVeto15',yaxis='Events',xaxis='Muon Veto (p_{T}>15 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
        #plotMethod('finalstate.elecVeto10',[8,0,8],savedir+'elecVeto10',yaxis='Events',xaxis='Electron Veto (p_{T}>10 GeV)',lumitext=33,logy=0,cut=myCut,**kwargs)
    plotMethod('finalstate.met',[20,0,200],savedir+'met',yaxis='Events/10 GeV',xaxis='E_{T}^{miss} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    mtext = 'm_{3\\ell} (GeV)'
    if analysis in ['Hpp4l']: mtext = 'm_{4\\ell} (GeV)'
    #if analysis in ['WZ']: mtext = 'm_{\\ell\\ell\\\'\\ell\\\'} (GeV)'
    plotMethod('finalstate.mass',[25,0,500],savedir+'mass',yaxis='Events/20 GeV',xaxis=mtext,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.mass',[100,0,500],savedir+'mass_zoom',yaxis='Events/5 GeV',xaxis=mtext,logy=0,cut=myCut,overflow=True,**kwargs)
    #plotMethod('finalstate.mass',[250,0,500],savedir+'mass_zoom',yaxis='Events/2.0 GeV',xaxis='M_{3\\ell} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    plotMethod('finalstate.mT',[50,0,1000],savedir+'mT',yaxis='Events/20 GeV',xaxis='m_T^{3\\ell+MET} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    #plotMethod('finalstate.mT',[150,0,300],savedir+'mT_zoom',yaxis='Events/2.0 GeV',xaxis='M_T^{3\\ell+MET} (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
    #plotMethod('finalstate.hT',[40,0,800],savedir+'hT',yaxis='Events/20.0 GeV',xaxis='H_{T} (GeV)',lumitext=33,logy=1,cut=myCut,overflow=True,**kwargs)
    plotMethod('event.nvtx',[50,0,50],savedir+'puVertices',yaxis='Events',xaxis='Number PU Vertices',legendpos=43,logy=0,cut=myCut,**kwargs)
    plotMethod('event.nvtx',[50,0,50],savedir+'puVertices_noreweight',yaxis='Events',xaxis='Number PU Vertices',legendpos=43,logy=0,cut=myCut, scalefactor='event.gen_weight*event.lep_scale*event.trig_scale',**kwargs)
    if analysis in ['WZ','WZ_Dijet']:
        plotMethod('finalstate.leadJetPt',[30,0,300],savedir+'JetPt',yaxis='Events/10 GeV',xaxis='p_{T}^{jet} (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('finalstate.leadJetEta',[50,-5.0,5.0],savedir+'JetEta',yaxis='Events',xaxis='\\eta^{jet}',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('finalstate.leadJetPhi',[30,-3.14159,3.14159],savedir+'JetPhi',yaxis='Events',xaxis='\\phi^{jet}',legendpos=43,logy=0,cut=myCut,**kwargs)
    # plot lepton kinematics
    if doDetailed:
        for l in range(nl):
            obj = 'l%i' % (l+1)
            plotLepton(plotMethod,myCut,obj,name=obj,pretty='\\ell%i'%(l+1),savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)

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
    #for l in ['e','m']:
    #    name = names[l]
    #    t = tex[l]
    #    cuts = ['%s & %s' %(myCut,'l%iFlv=="%s"' %((x+1),l)) for x in range(nl)]
    #    plotMethod(['l%i.Pt'  %(x+1) for x in range(nl)], [40,0,200],            savedir+'%s/Pt'  %name, yaxis='Events/5 GeV', xaxis='p_{T}^{%s} (GeV)' %t,        legendpos=43, logy=0, cut=cuts, overflow=True, **kwargs)
    #    plotMethod(['l%i.Iso' %(x+1) for x in range(nl)], [50,0,.5],             savedir+'%s/Iso' %name, yaxis='Events',         xaxis='Relative Isolation (%s)' %t, legendpos=43, logy=0, cut=cuts, overflow=True, **kwargs)
    #    plotMethod(['l%i.Eta' %(x+1) for x in range(nl)], [30,-3.0,3.0],         savedir+'%s/Eta' %name, yaxis='Events',         xaxis='\\eta^{%s}' %t,              legendpos=43, logy=0, cut=cuts, numcol=3, **kwargs)
    #    plotMethod(['l%i.Phi' %(x+1) for x in range(nl)], [30,-3.14159,3.14159], savedir+'%s/Phi' %name, yaxis='Events',         xaxis='\\phi^{%s}' %t,              legendpos=43, logy=0, cut=cuts, numcol=3, **kwargs)
    #    if not doMinimal:
    #        pass
    #        #plotMethod(['l%i.ChargeConsistent' %(x+1) for x in range(nl)], [3,-1.5,1.5], savedir+'%s/ChargeId' %name, yaxis='Events', xaxis='Charge ID (%s)' %t,         legendpos=43, logy=0, cut=cuts, numcol=3, **kwargs)

    # plot doubly charged higgs stuff
    if analysis in ['Hpp4l']:
        plotMethod('h1.mass', [24,0,600],savedir+'hpp/Mass',              yaxis='Events/25 GeV',xaxis='m_{\\ell^{+}\\ell^{+}} (GeV)',        legendpos=43,yscale=3.5,logy=1,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.mass', [32,0,800],savedir+'hpp/Mass_alpha',        yaxis='Events/25 GeV',xaxis='m_{\\ell^{+}\\ell^{+}} (GeV)', lumitext=33,logy=1,cut=myCut,boxes=[[12,0.9*mass,1],[1.1*mass,800,1],[0.9*mass,1.1*mass,2]],**kwargs)
        plotMethod('h1.dPhi', [32,0,3.2],savedir+'hpp/Dphi',              yaxis='Events/0.1 rad',       xaxis='\\Delta\\phi_{\\ell^{+}\\ell^{+}} (rad)',   legendpos=41,lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('h1.Pt',   [40,0,400],savedir+'hpp/Pt',                yaxis='Events/10 GeV',      xaxis='p_{T}^{\\Phi^{++}} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.dR',   [60,0,6], savedir+'hpp/dR',                 yaxis='Events',               xaxis='\\Delta R(\\ell^{+}\\ell^{+})',             legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('abs(h1.mass-h2.mass)', [30,0,300],savedir+'massdiff', yaxis='Events/10 GeV',xaxis='|m_{\\ell^{+}\\ell^{+}}-M_{\\ell^{-}\\ell^{-}}| (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h2.mass', [24,0,600],savedir+'hmm/Mass',              yaxis='Events/25 GeV',xaxis='m_{\\ell^{-}\\ell^{-}} (GeV)',        lumitext=33,logy=1,cut=myCut,overflow=True,**kwargs)
        plotMethod('h2.mass', [32,0,800],savedir+'hmm/Mass_alpha',        yaxis='Events/25 GeV',xaxis='m_{\\ell^{-}\\ell^{-}} (GeV)', lumitext=33,logy=1,cut=myCut,boxes=[[12,0.9*mass,1],[1.1*mass,800,1],[0.9*mass,1.1*mass,2]],**kwargs)
        plotMethod('h2.dPhi', [32,0,3.2],savedir+'hmm/Dphi',              yaxis='Events/0.1 rad',       xaxis='\\Delta\\phi_{\\ell^{-}\\ell^{-}} (rad)',   legendpos=41,lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('h2.Pt',   [40,0,400],savedir+'hmm/Pt',                yaxis='Events/10 GeV',      xaxis='p_{T}^{\\Phi^{--}} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h2.dR',   [60,0,6], savedir+'hmm/dR',                 yaxis='Events',               xaxis='\\Delta R(\\ell^{-}\\ell^{-})',             legendpos=43,logy=0,cut=myCut,**kwargs)
        plotLepton(plotMethod,myCut,'h1',post='1',name='hpp/Leading',pretty='\\Phi^{++} \\text{Leading Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'h1',post='2',name='hpp/SubLeading',pretty='\\Phi^{++} \\text{SubLeading Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'h2',post='1',name='hmm/Leading',pretty='\\Phi^{--} \\text{Leading Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'h2',post='2',name='hmm/SubLeading',pretty='\\Phi^{--} \\text{SubLeading Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
    if analysis in ['Hpp3l']:
        plotMethod('h1.mass', [24,0,600],savedir+'hpp/Mass',              yaxis='Events/25 GeV',xaxis='m_{\\ell^{\\pm}\\ell^{\\pm}} (GeV)',        legendpos=43,yscale=3.5,logy=1,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.mass', [32,0,800],savedir+'hpp/Mass_alpha',        yaxis='Events/25 GeV',xaxis='m_{\\ell^{\\pm}\\ell^{\\pm}} (GeV)', lumitext=33,logy=1,cut=myCut,boxes=[[12,0.9*mass,1],[1.1*mass,800,1],[0.9*mass,1.1*mass,2]],**kwargs)
        plotMethod('h1.dPhi', [32,0,3.2],savedir+'hpp/Dphi',              yaxis='Events/0.1 rad',       xaxis='\\Delta\\phi_{\\ell^{\\pm}\\ell^{\\pm}} (rad)',   legendpos=41,lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('h1.Pt',   [40,0,400],savedir+'hpp/Pt',                yaxis='Events/10 GeV',      xaxis='p_{T}^{\\Phi^{\\pm\\pm}} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h1.dR',   [60,0,6], savedir+'hpp/dR',                 yaxis='Events',               xaxis='\\Delta R(\\ell^{\\pm}\\ell^{\\pm})',             legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('abs(h1.mass-h2.mass)', [30,0,300],savedir+'massdiff',yaxis='Events/10 GeV',xaxis='|m_{\\ell^{\\pm}\\ell^{\\pm}}-M_{\\ell^{\\mp}}| (GeV)',lumitext=33,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h2.mass', [24,0,600],savedir+'hm/Mass',              yaxis='Events/25 GeV',xaxis='m_{\\ell^{\\mp},E_{T}^{miss}} (GeV)',        lumitext=33,logy=1,cut=myCut,overflow=True,**kwargs)
        plotMethod('h2.dPhi', [32,0,3.2],savedir+'hm/Dphi',              yaxis='Events/0.1 rad',       xaxis='\\Delta\\phi_{\\ell^{\\mp},E_{T}^{miss}} (rad)',   legendpos=41,lumitext=33,logy=0,cut=myCut,**kwargs)
        plotMethod('h2.Pt',   [40,0,400],savedir+'hm/Pt',                yaxis='Events/10 GeV',      xaxis='p_{T}^{\\Phi^{\\mp}} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('h2.dR',   [60,0,6], savedir+'hm/dR',                 yaxis='Events',               xaxis='\\Delta R(\\ell^{\\mp},E_{T}^{miss})',             legendpos=43,logy=0,cut=myCut,**kwargs)
        plotLepton(plotMethod,myCut,'h1',post='1',name='hpp/Leading',pretty='\\Phi^{\\pm\\pm} \\text{Leading Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'h1',post='2',name='hpp/SubLeading',pretty='\\Phi^{\\pm\\pm} \\text{SubLeading Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'h2',post='1',name='hm/Lepton',pretty='\\Phi^{\\mp} \\text{Lepton}',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
    # plot Z stuff
    if analysis in ['Z', 'Hpp3l', 'Hpp4l', 'WZ', 'WZ_W'] or region in ['Z', 'TT']:
        plotMethod('z1.mass', [60,60,120],   savedir+'z1/Mass_newWidth',      yaxis='Events/1.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',     logy=0,cut=myCut,**kwargs)
        plotMethod('abs(z1.mass-{0})'.format(ZMASS), [60,0,60],   savedir+'z1/mllMinusMZ',      yaxis='Events/1.0 GeV', xaxis='|M_{\\ell^{+}\\ell^{-}}-M_Z| (GeV)',     legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z1.mass', [13,58.5,123.5],   savedir+'z1/Mass_newWidth_wide', yaxis='Events/5 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',logy=0,cut=myCut,**kwargs)
        plotMethod('z1.Pt',   [40,0,400],    savedir+'z1/Pt',                 yaxis='Events/10 GeV',xaxis='p_{T}^{Z} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        if not doMinimal:
            plotMethod('z1.mass', [42,70,112],   savedir+'z1/Mass',               yaxis='Events/1.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('z1.mass', [7,80.5,101.5],savedir+'z1/Mass_wideBin',       yaxis='Events/3.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('z1.mass', [80,0,240],    savedir+'z1/Mass_fullWindow',    yaxis='Events/3.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
            plotMethod('z1.mass', [80,0,240],    savedir+'z1/Mass_fullWindow_log',yaxis='Events/3.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (GeV)',     legendpos=43,logy=1,cut=myCut,overflow=True,**kwargs)
            plotMethod('z1.dR',   [60,0,6],      savedir+'z1/dR',                 yaxis='Events',         xaxis='\\Delta R(\\ell^{+}\\ell^{-})',    legendpos=43,logy=0,cut=myCut,**kwargs)
        plotLepton(plotMethod,myCut,'z1',post='1',name='z1/Leading',pretty='Z Leading Lepton',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'z1',post='2',name='z1/SubLeading',pretty='Z SubLeading Lepton',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
    # plot second z stuff
    if analysis in ['Hpp4l']:
        plotMethod('z2.mass',[42,70,112],   savedir+'z2/Mass',               yaxis='Events/1.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=0,cut=myCut,**kwargs)
        #plotMethod('z2.mass',[7,80.5,101.5],savedir+'z2/Mass_wideBin',       yaxis='Events/3.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[80,0,240],    savedir+'z2/Mass_fullWindow',    yaxis='Events/3.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod('z2.mass',[80,0,240],    savedir+'z2/Mass_fullWindow_log',yaxis='Events/3.0 GeV', xaxis='m_{\\ell^{+}\\ell^{-}} (Z2) (GeV)', legendpos=43,logy=1,cut=myCut,**kwargs)
        plotMethod('z2.Pt',  [40,0,400],    savedir+'z2/Pt',                 yaxis='Events/10 GeV',xaxis='p_{T}^{Z2} (GeV)',                  legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotLepton(plotMethod,myCut,'z2',post='1',name='z2/Leading',pretty='Z2 Leading Lepton',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
        plotLepton(plotMethod,myCut,'z2',post='2',name='z2/SubLeading',pretty='Z2 SubLeading Lepton',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
    # plot W stuff
    if analysis in ['Hpp3l', 'WZ', 'WZ_W', 'WZ_Dijet']:
        plotMethod('w1.Pt',  [40,0,400],savedir+'w1/Pt',      yaxis='Events/10 GeV',xaxis='p_{T}^{W} (GeV)',                           legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.mass',[40,0,400],savedir+'w1/Mass',    yaxis='Events/10 GeV', xaxis='m_{T}^{W} (GeV)',                           legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w1.dPhi',[32,0,3.2],savedir+'w1/dPhi',    yaxis='Events/0.1 rad', xaxis='\\Delta\\phi(W lepton, E_{T}^{miss}) (rad)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotLepton(plotMethod,myCut,'w1',post='1',name='w1/Lepton',pretty='W Lepton',savedir=savedir,doDetailed=doDetailed,doMinimal=doMinimal,**kwargs)
    if analysis in ['Hpp3l', 'WZ'] and not doMinimal:
        if doDetailed:
            plotMethod('w1.mll_z1_1',[80,0,240],savedir+'w1/dilepton_mass_1_ss',yaxis='Events/3.0 GeV',xaxis='m(l^{#pm}l^{#pm}) (Z_{l1},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1==z1.Chg1 & %s' %myCut,overflow=True,**kwargs)
            plotMethod('w1.mll_z1_2',[80,0,240],savedir+'w1/dilepton_mass_2_ss',yaxis='Events/3.0 GeV',xaxis='m(l^{#pm}l^{#pm}) (Z_{l2},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1==z1.Chg2 & %s' %myCut,overflow=True,**kwargs)
            plotMethod('w1.mll_z1_1',[80,0,240],savedir+'w1/dilepton_mass_1_os',yaxis='Events/3.0 GeV',xaxis='m(l^{#pm}l^{#mp}) (Z_{l1},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1!=z1.Chg1 & %s' %myCut,overflow=True,**kwargs)
            plotMethod('w1.mll_z1_2',[80,0,240],savedir+'w1/dilepton_mass_2_os',yaxis='Events/3.0 GeV',xaxis='m(l^{#pm}l^{#mp}) (Z_{l2},W_{l}) (GeV)',legendpos=43,logy=0,cut='w1.Chg1!=z1.Chg2 & %s' %myCut,overflow=True,**kwargs)
        plotMethod(['w1.mll_z1_1','w1.mll_z1_2'],[80,0,240],savedir+'w1/dilepton_mass',yaxis='Events/3.0 GeV',xaxis='m(ll) (Z_{l},W_{l}) (GeV)',legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod(['w1.mll_z1_1','w1.mll_z1_2'],[80,0,240],savedir+'w1/dilepton_mass_ss',yaxis='Events/3.0 GeV',xaxis='m(ll) (Z_{l},W_{l}) (GeV)',legendpos=43,logy=0,cut=['w1.Chg1==z1.Chg1 & %s' %myCut,'w1.Chg1==z1.Chg2 & %s' %myCut],overflow=True,**kwargs)
        plotMethod(['w1.mll_z1_1','w1.mll_z1_2'],[80,0,240],savedir+'w1/dilepton_mass_os',yaxis='Events/3.0 GeV',xaxis='m(ll) (Z_{l},W_{l}) (GeV)',legendpos=43,logy=0,cut=['w1.Chg1!=z1.Chg1 & %s' %myCut,'w1.Chg1!=z1.Chg2 & %s' %myCut],overflow=True,**kwargs)
    if analysis in ['WZ_W']:
        plotMethod('w2.Pt',  [40,0,400],savedir+'w2/Pt',      yaxis='Events/10 GeV',xaxis='p_{T}^{W} (GeV)',                           legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w2.mass',[40,0,200],savedir+'w2/Mass',    yaxis='Events/5 GeV', xaxis='m_{T}^{W} (GeV)',                           legendpos=43,logy=0,cut=myCut,overflow=True,**kwargs)
        plotMethod('w2.dPhi',[32,0,3.2],savedir+'w2/dPhi',    yaxis='Events/0.1 rad', xaxis='\\Delta\\phi(W lepton, E_{T}^{miss}) (rad)',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotLepton(plotMethod,myCut,'w2',post='2',name='w2/Lepton',pretty='W Lepton',savedir=savedir,doDetailed=doDetailed,**kwargs)
    if analysis in ['Hpp3l', 'WZ'] and not doMinimal:
        if doDetailed:
            plotMethod('w1.dR1_z1_1',[60,0,6],savedir+'w1/dR_z1_1',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{leading lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)
            plotMethod('w1.dR1_z1_2',[60,0,6],savedir+'w1/dR_z1_2',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{subleading lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)
        plotMethod(['w1.dR1_z1_1','w1.dR1_z1_2'],[60,0,6],savedir+'w1/dR_z1_l',yaxis='Events',xaxis='#DeltaR(W^{lepton},Z^{lepton})',legendpos=43,logy=0,cut=myCut,**kwargs)

def getLeptonParams(obj,**kwargs):
    pre = kwargs.pop('pre','')
    post = kwargs.pop('post','')
    name = kwargs.pop('name','')
    pretty = kwargs.pop('pretty','')
    lepPlots = [
        #( variable,                  binning,              name,                   plotArgs
        ('%s.%sPt%s' %(obj,pre,post), [20,0,200],           '%s/Pt' %name, {'yaxis':'Events/10 GeV', 'xaxis':'p_{T}^{%s} (GeV)' %pretty,}),
        ('%s.%sIso%s' %(obj,pre,post),[50,0,.5],            '%s/Iso' %name,{'yaxis':'Events',        'xaxis':'Relative Isolation (%s)' %pretty,}),
        ('%s.%sEta%s' %(obj,pre,post),[30,-3.0,3.0],        '%s/Eta' %name,{'yaxis':'Events',        'xaxis':'\\eta^{%s}' %pretty,}),
        ('%s.%sPhi%s' %(obj,pre,post),[30,-3.14159,3.14159],'%s/Phi' %name,{'yaxis':'Events',        'xaxis':'\\phi^{%s}' %pretty,}),
    ]
    return lepPlots


def getPlotParams(plotMode,myCut,nl,plots,**kwargs):
    savedir = kwargs.pop('savedir','')
    analysis = kwargs.pop('analysis','')
    region = kwargs.pop('region','')
    mass = kwargs.pop('mass',500)
    doDetailed = kwargs.pop('doDetailed',False)
    doMinimal = kwargs.pop('doMinimal',False)
    doPAS = kwargs.pop('doPAS',False)
    if savedir: savedir += '/'

    basePlot = {
        'plotMode' : plotMode,
        'variable' : '1',
        'binning'  : [1,0,1],
        'savename' : savedir,
        'plotArgs' : {
            'cut' : myCut,
        },
    }

    # setup plot params
    allPlots = [
        #( variable,                                         binning,              name,                   plotArgs
        ('finalstate.sT',                                    [40,0,1000],          'sT',                   {'yaxis':'Events/25 GeV', 'xaxis':'S_{T} (GeV)',}),
        ('finalstate.elecVetoLoose',                         [4,0,4],              'numElectronsLoose',    {'yaxis':'Events',        'xaxis':'Number of Electrons (p_{T}>10 GeV)',}),
        ('finalstate.muonVetoLoose',                         [4,0,4],              'numMuonsLoose',        {'yaxis':'Events',        'xaxis':'Number of Muons (p_{T}>10 GeV)',}),
        ('finalstate.elecVetoLoose+finalstate.muonVetoLoose',[4,0,4],              'numLeptonsLoose',      {'yaxis':'Events',        'xaxis':'Number of Leptons (p_{T}>10 GeV)',}),
        ('finalstate.jetVeto30',                             [4,0,4],              'numJets30',            {'yaxis':'Events',        'xaxis':'Number of Jets (p_{T}>30 GeV)',}),
        ('finalstate.bjetVeto20Tight',                       [4,0,4],              'numBJets20Tight',      {'yaxis':'Events',        'xaxis':'Number of b Jets (p_{T}>20 GeV)',}),
        ('finalstate.met',                                   [20,0,200],           'met',                  {'yaxis':'Events/10 GeV', 'xaxis':'E_{T}^{miss} (GeV)',}),
        ('finalstate.mass',                                  [25,0,500],           'mass',                 {'yaxis':'Events/20 GeV', 'xaxis':'m_{3l}',}),
        ('finalstate.mT',                                    [15,0,600],           'mT',                   {'yaxis':'Events/40 GeV', 'xaxis':'m_{T}^{3l+MET} (GeV)',}),
        ('event.nvtx',                                       [50,0,50],            'puVertices',           {'yaxis':'Events',        'xaxis':'Number PU Vertices',}),
        ('event.nvtx',                                       [50,0,50],            'puVertices_noreweight',{'yaxis':'Events',        'xaxis':'Number PU Vertices', 'scalefactor':'event.gen_weight*event.lep_scale*event.trig_scale',}),
        ('finalstate.leadJetPt',                             [30,0,300],           'JetPt',                {'yaxis':'Events/10 GeV', 'xaxis':'p_{T}^{jet} (GeV)',}),
        ('finalstate.leadJetEta',                            [50,-5.0,5.0],        'JetEta',               {'yaxis':'Events',        'xaxis':'#eta^{jet}',}),
        ('finalstate.leadJetPhi',                            [30,-3.14159,3.14159],'JetPhi',               {'yaxis':'Events',        'xaxis':'#phi^{jet}',}),
        ('z1.mass',                                          [30,60,120],          'z1/Mass',              {'yaxis':'Events/2.0 GeV','xaxis':'m_{l^{+}l^{-}} (GeV)',}),
        ('abs(z1.mass-{0})'.format(ZMASS),                   [60,0,60],            'z1/mllMinusMZ',        {'yaxis':'Events/1.0 GeV','xaxis':'|m_{l^{+}l^{-}}-m_{Z}| (GeV)',}),
        ('z1.Pt',                                            [20,0,400],           'z1/Pt',                {'yaxis':'Events/20 GeV', 'xaxis':'p_{T}^{Z} (GeV)',}),
        ('w1.Pt',                                            [20,0,400],           'w1/Pt',                {'yaxis':'Events/20 GeV', 'xaxis':'p_{T}^{W} (GeV)',}),
        ('w1.mass',                                          [20,0,200],           'w1/Mass',              {'yaxis':'Events/10 GeV', 'xaxis':'m_{T}^{W} (GeV)',}),
        ('w1.dPhi',                                          [32,-3.2,3.2],        'w1/dPhi',              {'yaxis':'Events/0.2 rad','xaxis':'#Delta#phi(W lepton, E_{T}^{miss}) (rad)',}),
    ]
    # setup lepton params
    allPlots += getLeptonParams('z1',post='1',name='z1/Leading',pretty='Z Leading Lepton')
    allPlots += getLeptonParams('z1',post='2',name='z1/SubLeading',pretty='Z SubLeading Lepton')
    allPlots += getLeptonParams('w1',post='1',name='w1/Lepton',pretty='W Lepton')

    plotsMap = {}

    for plot in allPlots:
        variable, binning, name, plotArgs = plot
        if name not in plots and 'all' not in plots: continue
        plotsMap[name] = deepcopy(basePlot)
        plotsMap[name]['variable'] = variable
        plotsMap[name]['binning']  = binning
        plotsMap[name]['savename'] = savedir+name
        plotsMap[name]['plotArgs'].update(plotArgs)
  
    return plotsMap

def getCutflowParams(analysis,finalStates,cut,savedir,**kwargs):
    if savedir: savedir += '/'
    params = {}
    # channels on same plot
    plotChannelStrings, plotChannelCuts = getChannelStringsCuts(analysis,finalStates)
    params['individualChannels'] = {'args': [[cut] + ['{0} && {1}'.format(cut,c) for c in plotChannelCuts], '{0}individualChannels'.format(savedir)], 'kwargs': {'labels': ['Total'] + plotChannelStrings, 'nosum': True, 'numcol': 2, 'lumitext': 33}}

    return params


def getFakeParams(analysis,cut='1'):
    # define fake regions
    lepName = {'e': 'Elec', 'm': 'Muon', 't': 'Tau'}
    fakeRegions = {}
    fakeRegions[analysis] = {}
    for f in ['e', 'm']:
        #for p in ['Loose', 'Tight']:
        for p in ['Medium','Tight']:
            # select leading Z pt, Z window [60,120], tight (or loose) Z, low met, m3l>100, w1 mass < 30
            if analysis in ['WZ']:
                #for z in ['Loose', 'Tight']:
                for z in ['Tight']:
                    fakeRegion = 'Z{0}Probe{1}{2}'.format(z,lepName[f],p)
                    #denom = 'z1.Pass{0}1 && z1.Pass{0}2 && w1.dR1_z1_1>0.02 && w1.dR1_z1_2>0.02 && w1.mass<25. && w1Flv=="{1}"'.format(z,f)
                    basecut = 'z1.Pass{0}1 && z1.Pass{0}2 && w1.dR1_z1_1>0.02 && w1.dR1_z1_2>0.02 && finalstate.met<25. && w1.mass<30. && w1Flv=="{1}" && {2}'.format(z,f,cut)
                    numer = '{0} && w1.Pass{1}1==1'.format(basecut,p)
                    denom = '{0} && w1.Pass{1}1==0'.format(basecut,p)
                    fakeRegions[analysis][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
                    #if p=='Tight' and channel not in ['HZZFakeRate']:
                    #   fakeRegion += '_LooseProbe'
                    #   denom += ' && w1.PassLoose1'
                    #   numer += ' && w1.PassLoose1'
                    #   fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}
            # select w lepton pt, z veto, met
            #'W' : 'w1.Pt1>20. & (z1.mass<60. | z1.mass>120.) & finalstate.met>30. & w1.mass>30.',
            if analysis in ['WZ_W']:
                #for w in ['Loose','Tight']:
                for w in ['Tight']:
                    fakeRegion = 'W{0}Probe{1}{2}'.format(w,lepName[f],p)
                    basecut = 'w1.Pt1>20. && w1.mass>30. && finalstate.met>30. && (z1.mass<60. || z1.mass>120.) && l1.Chg==l2.Chg && z1.dR>0.1 && w1.Pass{0}1 && w2Flv=="{1}" && {2}'.format(w,f,cut)
                    numer = '{0} && w2.Pass{1}1==1'.format(basecut,p)
                    denom = '{0} && w2.Pass{1}1==0'.format(basecut,p)
                    fakeRegions[analysis][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
                    #if p=='Tight' and channel not in ['HZZFakeRate']:
                    #   fakeRegion += '_LooseProbe'
                    #   denom += ' && w2.PassLoose1'
                    #   numer += ' && w2.PassLoose1'
                    #   fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w2.Pt1', 'etaVar': 'w2.Eta1'}
            # ntuple cuts: zVeto 60-120, met vet 20, w veto 20, jet pt > 20, jet dr > 1.0
            if analysis in ['WZ_Dijet']:
               fakeRegion = 'FakeRateProbe{0}{1}'.format(lepName[f],p)
               basecut = 'l1Flv=="{0}" && finalstate.leadJetPt>{1} && {2}'.format(f,20 if f=='m' else 35,cut)
               numer = '{0} && w1.Pass{1}1==1'.format(basecut,p)
               #denom = '{0} && w1.Pass{1}1==0'.format(basecut,p)
               denom = basecut
               numerScale = 'event.gen_weight*event.pu_weight*w1.LepScale{0}1*(1./event.trig_prescale)'.format(p)
               denomScale = 'event.gen_weight*event.pu_weight*w1.LepScaleLoose1*(1./event.trig_prescale)'
               fakeRegions[analysis][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1', 'numerScale': numerScale, 'denomScale': denomScale}
               #if p=='Tight' and channel not in ['HZZFakeRate']:
               #   fakeRegion += '_LooseProbe'
               #   denom += ' && w1.PassLoose1'
               #   numer += ' && w1.PassLoose1'
               #   fakeRegions['WZ'][fakeRegion] = {'denom': denom, 'numer': numer, 'probe': f, 'ptVar': 'w1.Pt1', 'etaVar': 'w1.Eta1'}



    # setup selections
    ptBins = [10,1000]
    #if analysis in ['WZ_Dijet']: ptBins = [1,10,15,20,25,30,40,50,200]

    etaBins = {
        #'e': [0,1.479,2.5],
        #'m': [0,1.2,2.4],
        'e': [0,2.5],
        'm': [0,2.4],
    }

    return fakeRegions, ptBins, etaBins

def getWZGenYields():
    yields = {
        'fiducial' : {
            'eee' : 34092,
            'eem' : 33937,
            'mme' : 34480,
            'mmm' : 34276,
        },
        'zWindow' : {
            'eee' : 66171,
            'eem' : 66448,
            'mme' : 66954,
            'mmm' : 66398,
        },
    }

def getAcceptanceEfficiency(analysis):
    '''Return dictionary of acceptance times efficiency'''
    accEff = {
        'zWindow' : {
            'eee' : 0.1450,
            'eem' : 0.1637,
            'mme' : 0.1950,
            'mmm' : 0.2385,
        },
        'fiducial' : {
            'eee' : 0.2814,
            'eem' : 0.3205,
            'mme' : 0.3786,
            'mmm' : 0.4621,
        },
    }
    return accEff

def getNtupleDirectory(analysis,region,period):
    ntupleMap = {
        #('WZ','WZ',13) : '2015-11-14_WZ_WZ_13_isa_updatedNtuples_v1',
    }
    if (analysis,region,period) in ntupleMap:
        return '/hdfs/store/user/dntaylor/{0}'.format(ntupleMap[(analysis,region,period)])
    else:
        return 'ntuples/{0}_{1}TeV_{2}'.format(analysis,period,region)

def getPassTightDefinition(analysis,region,period):
    cutMap = {
        #('WZ','WZ',13) : 'finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && z1.mass>60. && z1.mass<120. && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.',
        ('WZ','WZ',13) : 'finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-91.1876)<15.   && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30. && finalstate.bjetVeto20Tight<1',
    }
    key = (analysis,region,period)
    if key in cutMap:
        return cutMap[key]
    else:
        return 'select.passTight'
